#!/usr/bin/env python3
"""LAGMX post-production analysis.

LAGMX takes a protein-ligand system from PDB to a finished production
trajectory. What it has never done is tell you anything about that trajectory:
every project ended with an .xtc and a fresh round of hand-typed `gmx rms`,
`gmx sasa`, `gmx covar` invocations, each with its own interactive group
prompts and its own chance of fitting on the wrong group.

This script closes that gap. It runs the standard analysis battery over every
complex_*/ directory that has a finished production run, writes CSV and PNG
per analysis, and reduces each complex to one row of summary numbers so that
several systems can be compared side by side.

    cd run_matrix && python3 ../analyze_md.py

Like LAGMX.py, every path is resolved against the directory you run it from,
not against the location of this file.

Analyses
    rmsd      complex stability, protein backbone and ligand separately
    rmsf      per-residue flexibility
    rg        radius of gyration, protein compactness
    sasa      solvent accessible surface, protein and complex
    hbond     protein-ligand hydrogen bonds over time
    contacts  per-residue protein-ligand contact occupancy
    pca       essential dynamics, PC1/PC2 projection
    fel       free energy landscape over PC1/PC2
    mmpbsa    binding free energy via gmx_MMPBSA (MM/GBSA and/or MM/PBSA)

Every one of these except fel reads the trajectory front to back on a single
core -- gmx analysis tools have no -nt -- so each has its own frame stride,
analysis_stride_<name>, set in frames. The defaults are not uniform because the
analyses are not: see the comment above STRIDABLE for what each costs and what
timescale it actually resolves.

Configuration comes from gmx_config.txt; see ANALYSIS_DEFAULTS below for the
keys this script adds. All of them are optional.
"""

import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Keys this script adds to gmx_config.txt. Everything has a working default so
# that an existing config file keeps working untouched.
ANALYSIS_DEFAULTS = {
    "analysis": "rmsd,rmsf,rg,sasa,hbond,contacts,pca,fel",
    "analysis_skip_ns": "0",          # drop this much from the start of the run
    "analysis_contact_cutoff": "0.4",  # nm, protein-ligand contact definition
    "analysis_gmx": "gmx",            # which gmx to analyse with
    # Frame stride per analysis: use one frame in every N, 1 meaning every
    # frame. The comment above STRIDABLE explains why each default is what it
    # is. fel has no entry because gmx sham bins the PCA projection rather than
    # reading the trajectory, so it inherits whatever stride pca used.
    "analysis_stride_rmsd": "1",
    "analysis_stride_rmsf": "1",
    "analysis_stride_rg": "1",
    "analysis_stride_sasa": "10",
    "analysis_stride_hbond": "5",
    "analysis_stride_contacts": "10",
    "analysis_stride_pca": "10",
    "mmpbsa_python": "",              # env holding gmx_MMPBSA; empty = disabled
    "mmpbsa_method": "gb",            # gb, pb, or both
    "mmpbsa_frames": "100",           # target frame count, if no interval given
    "mmpbsa_interval": "",            # explicit frame stride; overrides the target
    "mmpbsa_np": "0",                 # MPI ranks for gmx_MMPBSA; 0 = auto
    "mmpbsa_igb": "5",
    "mmpbsa_salt": "0.150",
}

ALL_ANALYSES = ["rmsd", "rmsf", "rg", "sasa", "hbond", "contacts", "pca", "fel", "mmpbsa"]

# Analyses whose stride is configurable, i.e. the ones that read the trajectory.
STRIDABLE = ["rmsd", "rmsf", "rg", "sasa", "hbond", "contacts", "pca"]

# Why the defaults differ. Every gmx analysis tool is single-threaded and reads
# the trajectory front to back, so cost is proportional to frame count and to
# nothing else -- there is no -nt to throw cores at. A 100 ns run written every
# 2 ps is 45000 frames, and on a real ThiM complex that cost 5.4 h for four
# analyses, of which SASA alone was 3.4 h.
#
# Striding is not a shortcut, it is a statement about which timescale each
# quantity lives on. Frames 2 ps apart are strongly correlated: they are close
# to the same sample counted many times, so they cost time without adding
# information. The defaults follow from that, plus from what each analysis
# actually costs:
#
#   rmsd, rmsf, rg   1   Cheap (10-19 min each). Striding would save little and
#                        these are the traces people read frame by frame.
#   sasa            10   By far the most expensive. Solvent-accessible surface
#                        is a slowly varying geometric quantity; 4500 samples
#                        over 90 ns resolve everything it can show.
#   hbond            5   Occupancy statistics only -- this script asks gmx hbond
#                        for -num and never for lifetimes or autocorrelation, so
#                        no result here needs 2 ps resolution. Kept finer than
#                        sasa because hydrogen bonds do break and reform on a
#                        ps-ns timescale and the count trace should keep its
#                        shape. Raising this is safe only while no lifetime
#                        analysis is added.
#   contacts        10   Per-residue occupancy is a fraction; 4500 frames put
#                        its standard error well under one percent.
#   pca             10   Covariance wants independent samples, which 2 ps frames
#                        are not. Striding improves the conditioning rather than
#                        harming it, and shrinks eigenvec.trr by the same factor.
#   fel                  Not listed: gmx sham bins the PCA projection and never
#                        opens the trajectory, so it inherits pca's stride.
#   mmpbsa               Has its own knob, mmpbsa_interval, for the same reason.

GMX = "gmx"


# --------------------------------------------------------------------------
# process helpers
# --------------------------------------------------------------------------

def gmx_run(args, stdin=None, cwd=None, quiet=True):
    """Call a gmx module with an argv list rather than a shell string.

    Group selection goes in on stdin as plain group numbers, which is what the
    interactive prompts read. Returns (ok, combined_output) instead of raising:
    one analysis failing must not take the other eight down with it.
    """
    cmd = [GMX] + [str(a) for a in args]
    env = dict(os.environ, GMX_MAXBACKUP="-1")
    try:
        proc = subprocess.run(
            cmd, input=stdin, cwd=cwd, capture_output=True, text=True,
            timeout=7200, env=env
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"{exc}"
    out = (proc.stdout or "") + (proc.stderr or "")
    if not quiet:
        print(out)
    return proc.returncode == 0, out


def say(msg, indent=0):
    print(f"{' ' * indent}{msg}", flush=True)


# --------------------------------------------------------------------------
# file format readers
# --------------------------------------------------------------------------

def read_xvg(path):
    """Return (data, legends) from a Grace .xvg file.

    Legends are pulled from the @ s0 legend lines so that plots and CSV
    headers carry the names GROMACS chose, instead of "column 1".
    """
    rows, legends, xlabel, ylabel = [], [], "", ""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("@"):
                m = re.match(r'@\s+s\d+\s+legend\s+"(.*)"', line)
                if m:
                    legends.append(m.group(1))
                m = re.match(r'@\s+xaxis\s+label\s+"(.*)"', line)
                if m:
                    xlabel = m.group(1)
                m = re.match(r'@\s+yaxis\s+label\s+"(.*)"', line)
                if m:
                    ylabel = m.group(1)
                continue
            if line.startswith(("#", "&")):
                continue
            parts = line.split()
            if not parts:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    if not rows:
        return np.empty((0, 0)), {"legends": legends, "xlabel": xlabel, "ylabel": ylabel}
    width = min(len(r) for r in rows)
    data = np.array([r[:width] for r in rows])
    return data, {"legends": legends, "xlabel": xlabel, "ylabel": ylabel}


def read_xpm(path):
    """Parse a GROMACS .xpm matrix into (values, x_axis, y_axis).

    gmx sham writes the free energy landscape as an XPM colour map. The kJ/mol
    value of each colour lives in the C comment after the colour definition --

        "A  c #000000 " /* "0" */,

    -- not inside the quoted colour string itself, so the file has to be
    decoded rather than read as an image. The header line gives the geometry
    and, crucially, how many characters encode one pixel.
    """
    header_re = re.compile(r'^"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*"')
    nx = ny = ncolours = nchar = None
    colours, rows, xaxis, yaxis = {}, [], [], []

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("/* x-axis:"):
                xaxis += [float(v) for v in re.findall(r"[-\d.eE+]+", line[10:])]
                continue
            if line.startswith("/* y-axis:"):
                yaxis += [float(v) for v in re.findall(r"[-\d.eE+]+", line[10:])]
                continue
            if not line.startswith('"'):
                continue

            if nx is None:
                m = header_re.match(line)
                if m:
                    nx, ny, ncolours, nchar = (int(g) for g in m.groups())
                continue

            if len(colours) < ncolours:
                m = re.match(r'^"(.{%d})\s+c\s+\S+\s*"\s*/\*\s*"([^"]*)"' % nchar, line)
                if m:
                    try:
                        colours[m.group(1)] = float(m.group(2))
                    except ValueError:
                        colours[m.group(1)] = np.nan
                continue

            body = line[1:line.rfind('"')]
            if len(body) >= nx * nchar:
                rows.append(body)

    if not rows or not colours:
        return None, None, None
    values = np.array([
        [colours.get(row[i:i + nchar], np.nan) for i in range(0, nx * nchar, nchar)]
        for row in rows
    ])
    # XPM rows run top to bottom, the y axis runs bottom to top.
    values = values[::-1]
    return values, np.array(xaxis), np.array(yaxis)


# --------------------------------------------------------------------------
# index groups
# --------------------------------------------------------------------------

def _index_group_names(ndx_path):
    """{group name: group number} read from an .ndx file.

    Group numbers are positions in the file, which is how gmx numbers them.
    """
    names = re.findall(r"^\s*\[\s*(.+?)\s*\]", open(ndx_path, errors="replace").read(), re.M)
    return {name: i for i, name in enumerate(names)}


def build_index(tpr, out_ndx, workdir, merged_group="Protein_LIG"):
    """Write an index file and return {group name: group number}.

    The default make_ndx groups already carry everything needed: "Other" is
    exactly the set of non-protein, non-water, non-ion residues, which for a
    LAGMX system is exactly the ligands, however many there are and whatever
    they are called. A merged Protein+ligand group is appended for centring
    and for MM/PBSA.
    """
    ok, out = gmx_run(["make_ndx", "-f", tpr, "-o", out_ndx], stdin="q\n", cwd=workdir)
    if not ok:
        return None
    groups = {name: int(num) for num, name in re.findall(r"^\s*(\d+)\s+(\S+)\s*:", out, re.M)}
    if "Other" not in groups:
        say("no 'Other' group: system has no ligand, ligand analyses will be skipped", 6)
    else:
        new_index = max(groups.values()) + 1
        ok, out2 = gmx_run(
            ["make_ndx", "-f", tpr, "-n", out_ndx, "-o", out_ndx],
            stdin=f'"Protein" | "Other"\nname {new_index} {merged_group}\nq\n', cwd=workdir,
        )
        if ok:
            groups[merged_group] = new_index

        # One group per ligand copy. A system with the same ligand in three
        # equivalent sites has three independent measurements in it; leaving
        # them fused in "Other" would report their sum as if it were one
        # binding event, and throw away the spread between the sites -- which
        # is the only internal check on how reproducible the number is.
        ok, _ = gmx_run(
            ["make_ndx", "-f", tpr, "-n", out_ndx, "-o", out_ndx],
            stdin=f'splitres {groups["Other"]}\nq\n', cwd=workdir)
        if ok:
            # Read the names back from the file, not from the command output.
            # make_ndx lists its groups when it starts and then quits after the
            # split without listing them again, so the new groups can never
            # appear in what the command printed -- parsing that output found
            # nothing, _ligan_terpisah stayed empty, and MM/PBSA silently fell
            # back to treating every ligand copy as one lumped "ligand". On the
            # ThiM trimer that meant a binding energy for three molecules at
            # once, reported as if it were one.
            for name, num in _index_group_names(out_ndx).items():
                if name not in groups:
                    groups[name] = num
                    groups.setdefault("_ligan_terpisah", []).append(name)
    return groups


# --------------------------------------------------------------------------
# trajectory preparation
# --------------------------------------------------------------------------

def frame_interval_ps(cdir):
    """Time between stored frames, read from the production .mdp.

    There is no cheap way to ask gmx: every tool that reports frame times reads
    the whole trajectory, which on a 12 GB file is exactly the cost striding
    exists to avoid. nstxout-compressed * dt is what produced the interval in
    the first place, so take it from there. Returns 0.0 if it cannot be read,
    which disables striding rather than guessing.
    """
    mdp = os.path.join(cdir, "md.mdp")
    if not os.path.exists(mdp):
        return 0.0
    nst, dt = None, None
    for line in open(mdp, errors="replace"):
        line = line.split(";", 1)[0]
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        key = key.lower().replace("_", "-")
        try:
            if key == "nstxout-compressed":
                nst = int(float(value.split()[0]))
            elif key == "dt":
                dt = float(value.split()[0])
        except (ValueError, IndexError):
            continue
    if not nst or not dt:
        return 0.0
    return nst * dt


def stride_args(ctx, key, tu="ps"):
    """gmx's -dt flag for one analysis, or nothing if it keeps every frame.

    gmx takes a time, not a frame count, so a stride of N frames is N times the
    interval between them. Configuration is in frames because that is the unit
    the user reasons in and the unit mmpbsa_interval already uses.

    `tu` must match the -tu the command is given, because -dt is interpreted in
    whatever unit -tu selects -- not in ps. Getting this wrong is silent and
    expensive: `gmx rms -tu ns -dt 10` was read as one frame per 10 ns rather
    than per 10 ps, and a 90 ns trajectory came back as ten data points instead
    of nine thousand. The file is valid, the plot draws, and nothing warns.
    """
    stride = ctx["strides"].get(key, 1)
    interval = ctx["frame_ps"]
    if stride <= 1 or interval <= 0:
        return []
    dt = stride * interval
    if tu == "ns":
        dt /= 1000.0
    return ["-dt", f"{dt:g}"]


# Identifies the PBC recipe in prepare_trajectory. Change it and every stored
# md_center.stamp stops matching, so prepared trajectories are rebuilt instead
# of silently reused with the treatment they were made under.
#   1 = whole -> nojump -> mol/compact/center   (broke multimers)
#   2 = whole -> cluster -> mol/compact/center
PREP_RECIPE = 2


def build_reference_tpr(cdir, adir, tpr):
    """A .tpr whose coordinates match the prepared trajectory.

    Every fitting analysis -- rms, covar, anaeig -- superimposes each frame on
    the structure held in the file given to -s. That file was md.tpr, whose
    coordinates are the raw ones grompp was handed: molecules split across the
    boundary and, for a multimer, protomers in different periodic images. The
    trajectory is corrected and the reference is not, and no rigid-body
    superposition can reconcile the two.

    Measured on the ThiM trimer: backbone RMSD came out at 3.56 nm against
    md.tpr and 0.222 nm against a reference built from the prepared first frame.
    The first number is stable, plausible-looking, and meaningless -- it is the
    size of a box vector, not a conformational change. This predates the
    nojump/cluster fix and affected every RMSD, RMSF and PCA produced here.

    A plain .pdb cannot serve as -s: gmx needs masses for the least-squares fit
    and a PDB carries none for ligand atoms. So grompp rebuilds a real .tpr from
    start.pdb, with the same topology and the same atom count.

    Returns the reference to use; falls back to the original tpr, saying so,
    rather than failing the whole analysis.
    """
    start = os.path.join(adir, "start.pdb")
    ref = os.path.join(adir, "ref.tpr")
    if not os.path.exists(start):
        return tpr
    if os.path.exists(ref) and os.path.getmtime(ref) >= os.path.getmtime(start):
        return ref

    mdp = os.path.join(cdir, "md.mdp")
    top = os.path.join(cdir, "topol.top")
    if not (os.path.exists(mdp) and os.path.exists(top)):
        say("md.mdp or topol.top missing; fitting against the raw tpr", 6)
        return tpr
    ok, out = gmx_run(["grompp", "-f", mdp, "-c", start, "-p", top,
                       "-n", os.path.join(adir, "analysis.ndx"),
                       "-o", ref, "-maxwarn", "10"], cwd=cdir)
    if not ok or not os.path.exists(ref):
        tail = out.strip().splitlines()[-1] if out.strip() else "?"
        say(f"reference tpr could not be built ({tail}); fitting against the raw tpr", 6)
        return tpr
    say("reference tpr rebuilt from the prepared first frame", 6)
    return ref


def prepare_trajectory(cdir, adir, tpr, xtc, groups, merged_group, skip_ps):
    """Undo periodic boundary artefacts before anything is measured.

    Three passes, in this order and no other: make molecules whole, gather the
    complex into one periodic image, then centre it in a compact box. Skipping
    this is the classic way to get an RMSD trace with a cliff in it that looks
    like unbinding and is really the ligand crossing the box edge.

    The middle pass is `-pbc cluster`, not `-pbc nojump`. nojump keeps each
    molecule continuous in time relative to the first frame, which is right for
    a single-chain protein and quietly wrong for a multimer: each protomer is
    its own moleculetype, so they drift into different images and the assembly
    comes apart. Measured on the ThiM trimer, whose three chains sit 32 A apart
    in the input: after nojump they were 74, 80 and 34 A apart -- one protomer
    displaced by roughly a box length -- in every one of 201 frames checked.
    Adding cluster before nojump does not help, because nojump undoes it.

    That break is invisible in the output and poisons everything computed on the
    assembly as a whole: the backbone fit behind every RMSD, the radius of
    gyration, the solvent-accessible surface (the A-B and A-C interfaces become
    exposed), the C-alpha covariance behind PCA and FEL, and the receptor
    MM/PBSA is handed. Per-site quantities -- contacts, hydrogen bonds -- survive,
    because each ligand travels with its own chain.

    With cluster in its place the same 201 frames hold at 32-33 A throughout,
    and each ligand stays 13-15 A from its own protomer.

    The result is reused if one is already here for the same source trajectory
    and the same analysis_skip_ns. Three passes over a 12 GB .xtc is an hour
    that should not be spent again just because a stride was changed -- and
    striding is per analysis, so it does not affect this file at all. The stamp
    records what the existing file was built from, so a changed skip_ns or a
    re-run production rebuilds it rather than being silently reused.
    """
    system = groups.get("System", 0)
    centre = groups.get(merged_group, groups.get("Protein", 1))

    whole = os.path.join(adir, "_whole.xtc")
    nojump = os.path.join(adir, "_clustered.xtc")
    final = os.path.join(adir, "md_center.xtc")
    start = os.path.join(adir, "start.pdb")
    stamp = os.path.join(adir, "md_center.stamp")

    # The stamp has to identify the pipeline, not just its inputs. It recorded
    # only skip_ps, the centring group and the source mtime, so changing the PBC
    # treatment left every existing stamp valid and the broken trajectories
    # would have been reused in silence. Bump PREP_RECIPE whenever the steps
    # below change.
    built_from = (f"recipe={PREP_RECIPE} skip_ps={skip_ps:g} centre={centre} "
                  f"src_mtime={os.path.getmtime(xtc):.0f}")
    if os.path.exists(final) and os.path.exists(start) and os.path.exists(stamp):
        if open(stamp, errors="replace").read().strip() == built_from:
            say("prepared trajectory reused", 6)
            return final

    steps = [
        (["trjconv", "-s", tpr, "-f", xtc, "-o", whole, "-pbc", "whole"], f"{system}\n"),
        # Cluster on the same group the box is centred on, so every chain and
        # every ligand copy is pulled into one image before centring.
        (["trjconv", "-s", tpr, "-f", whole, "-o", nojump, "-pbc", "cluster"],
         f"{centre}\n{system}\n"),
        (["trjconv", "-s", tpr, "-f", nojump, "-o", final,
          "-pbc", "mol", "-ur", "compact", "-center"], f"{centre}\n{system}\n"),
    ]
    if skip_ps > 0:
        steps[-1][0].extend(["-b", str(skip_ps)])

    for args, stdin in steps:
        ok, out = gmx_run(args + ["-n", os.path.join(adir, "analysis.ndx")],
                          stdin=stdin, cwd=cdir)
        if not ok:
            say(f"trjconv failed: {out.strip().splitlines()[-1] if out.strip() else '?'}", 6)
            return None

    for tmp in (whole, nojump):
        if os.path.exists(tmp):
            os.remove(tmp)

    # A single reference frame, used as the -s for analyses and for viewing.
    gmx_run(["trjconv", "-s", tpr, "-f", final, "-o", start,
             "-n", os.path.join(adir, "analysis.ndx"), "-dump", "0"],
            stdin=f"{system}\n", cwd=cdir)

    # Written last, so an interrupted preparation leaves no stamp to trust.
    if os.path.exists(start):
        with open(stamp, "w") as fh:
            fh.write(built_from + "\n")
    return final


# --------------------------------------------------------------------------
# plotting and export
# --------------------------------------------------------------------------

def export(data, meta, adir, name, title, xlabel=None, ylabel=None, columns=None):
    """Write one analysis to CSV and PNG."""
    if data.size == 0:
        return
    csv_path = os.path.join(adir, f"{name}.csv")
    header = columns or ([meta.get("xlabel") or "x"] +
                         (meta.get("legends") or
                          [f"y{i}" for i in range(1, data.shape[1])]))
    header = header[:data.shape[1]]
    while len(header) < data.shape[1]:
        header.append(f"y{len(header)}")
    np.savetxt(csv_path, data, delimiter=",", header=",".join(header), comments="", fmt="%.6g")

    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    for col in range(1, data.shape[1]):
        ax.plot(data[:, 0], data[:, col], lw=1.0, label=header[col])
    ax.set_xlabel(xlabel or meta.get("xlabel") or "")
    ax.set_ylabel(ylabel or meta.get("ylabel") or "")
    ax.set_title(title)
    if data.shape[1] > 2:
        ax.legend(fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(adir, f"{name}.png"), dpi=150)
    plt.close(fig)


def stat(data, col=1):
    """Mean and standard deviation of one column, ignoring empty input."""
    if data.size == 0 or data.shape[1] <= col:
        return None, None
    values = data[:, col]
    return float(np.mean(values)), float(np.std(values))


# --------------------------------------------------------------------------
# individual analyses
# --------------------------------------------------------------------------

def analyse_rmsd(ctx):
    """Backbone RMSD for stability, ligand RMSD for whether it stayed put.

    The ligand curve is deliberately fitted on the protein backbone, not on
    the ligand itself: fitting a ligand to its own reference measures internal
    conformational change and hides the thing you actually want to see, which
    is the ligand drifting out of the pocket.
    """
    results = {}
    backbone = ctx["groups"].get("Backbone")
    if backbone is None:
        return results

    # Both rms calls pass -tu ns, so -dt has to be in ns as well.
    step = stride_args(ctx, "rmsd", tu="ns")
    ok, _ = gmx_run(["rms", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-o", ctx["a"]("rmsd_protein.xvg"), "-tu", "ns"] + step,
                    stdin=f"{backbone}\n{backbone}\n", cwd=ctx["cdir"])
    if ok:
        data, meta = read_xvg(ctx["a"]("rmsd_protein.xvg"))
        export(data, meta, ctx["adir"], "rmsd_protein", "RMSD protein backbone",
               "waktu (ns)", "RMSD (nm)", ["waktu_ns", "rmsd_nm"])
        mean, sd = stat(data)
        results["rmsd_protein_nm"], results["rmsd_protein_sd"] = mean, sd

    ligand = ctx["groups"].get("Other")
    if ligand is not None:
        ok, _ = gmx_run(["rms", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                         "-o", ctx["a"]("rmsd_ligand.xvg"), "-tu", "ns"] + step,
                        stdin=f"{backbone}\n{ligand}\n", cwd=ctx["cdir"])
        if ok:
            data, meta = read_xvg(ctx["a"]("rmsd_ligand.xvg"))
            export(data, meta, ctx["adir"], "rmsd_ligand", "RMSD ligand (fit on protein)",
                   "waktu (ns)", "RMSD (nm)", ["waktu_ns", "rmsd_nm"])
            mean, sd = stat(data)
            results["rmsd_ligand_nm"], results["rmsd_ligand_sd"] = mean, sd
    return results


def analyse_rmsf(ctx):
    """Per-residue fluctuation, C-alpha only, averaged over each residue."""
    ca = ctx["groups"].get("C-alpha")
    if ca is None:
        return {}
    ok, _ = gmx_run(["rmsf", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-o", ctx["a"]("rmsf.xvg"), "-oq", ctx["a"]("rmsf_bfactor.pdb"),
                     "-res"] + stride_args(ctx, "rmsf"),
                    stdin=f"{ca}\n", cwd=ctx["cdir"])
    if not ok:
        return {}
    data, meta = read_xvg(ctx["a"]("rmsf.xvg"))
    export(data, meta, ctx["adir"], "rmsf", "RMSF per residu",
           "residu", "RMSF (nm)", ["residu", "rmsf_nm"])
    mean, _ = stat(data)
    out = {"rmsf_mean_nm": mean}
    if data.size:
        top = data[np.argsort(-data[:, 1])][:10]
        out["rmsf_residu_teratas"] = " ".join(f"{int(r[0])}:{r[1]:.2f}" for r in top)
    return out


def analyse_rg(ctx):
    """Radius of gyration. Falls back to the legacy tool: gmx gyrate was
    reimplemented in recent releases and older builds only ship one name."""
    protein = ctx["groups"].get("Protein")
    if protein is None:
        return {}
    for tool in ("gyrate", "gyrate-legacy"):
        ok, _ = gmx_run([tool, "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                         "-o", ctx["a"]("rg.xvg")] + stride_args(ctx, "rg"),
                        stdin=f"{protein}\n", cwd=ctx["cdir"])
        if ok:
            break
    else:
        return {}
    data, meta = read_xvg(ctx["a"]("rg.xvg"))
    if data.size:
        data = data[:, :2]
    export(data, meta, ctx["adir"], "rg", "Radius of gyration",
           "waktu (ps)", "Rg (nm)", ["waktu_ps", "rg_nm"])
    mean, sd = stat(data)
    return {"rg_nm": mean, "rg_sd": sd}


def analyse_sasa(ctx):
    """Solvent accessible surface of the protein and of the whole complex.

    The difference between the two, against the ligand's own free surface, is
    the buried area -- a cheap proxy for how deep the ligand sits.
    """
    out = {}
    jobs = [("sasa_protein", 'group "Protein"')]
    if "Other" in ctx["groups"]:
        jobs.append(("sasa_ligand", 'group "Other"'))
        jobs.append(("sasa_complex", 'group "Protein" or group "Other"'))
    for name, sel in jobs:
        ok, _ = gmx_run(["sasa", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                         "-surface", sel, "-o", ctx["a"](f"{name}.xvg")]
                        + stride_args(ctx, "sasa"), cwd=ctx["cdir"])
        if not ok:
            continue
        data, meta = read_xvg(ctx["a"](f"{name}.xvg"))
        if data.size:
            data = data[:, :2]
        export(data, meta, ctx["adir"], name, name.replace("_", " ").upper(),
               "waktu (ps)", "SASA (nm^2)", ["waktu_ps", "sasa_nm2"])
        mean, sd = stat(data)
        out[f"{name}_nm2"], out[f"{name}_sd"] = mean, sd

    if all(k in out for k in ("sasa_protein_nm2", "sasa_ligand_nm2", "sasa_complex_nm2")):
        out["sasa_terkubur_nm2"] = (out["sasa_protein_nm2"] + out["sasa_ligand_nm2"]
                                    - out["sasa_complex_nm2"])
    return out


def analyse_hbond(ctx):
    """Protein-ligand hydrogen bonds per frame.

    gmx hbond was rewritten with a selection interface in recent GROMACS, while
    older builds prompt for two index groups. Try the modern call first, fall
    back to the legacy one, so this works on both.
    """
    if "Other" not in ctx["groups"]:
        return {}
    target = ctx["a"]("hbond_num.xvg")

    step = stride_args(ctx, "hbond")
    ok, _ = gmx_run(["hbond", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-ref", 'group "Protein"', "-sel", 'group "Other"',
                     "-num", target] + step, cwd=ctx["cdir"])
    if not ok or not os.path.exists(target):
        protein, ligand = ctx["groups"]["Protein"], ctx["groups"]["Other"]
        for tool in ("hbond", "hbond-legacy"):
            ok, _ = gmx_run([tool, "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                             "-num", target] + step,
                            stdin=f"{protein}\n{ligand}\n", cwd=ctx["cdir"])
            if ok and os.path.exists(target):
                break
    if not os.path.exists(target):
        return {}

    data, meta = read_xvg(target)
    if data.size:
        data = data[:, :2]
    export(data, meta, ctx["adir"], "hbond", "Ikatan hidrogen protein-ligan",
           "waktu (ps)", "jumlah H-bond", ["waktu_ps", "n_hbond"])
    mean, sd = stat(data)
    out = {"hbond_rerata": mean, "hbond_sd": sd}
    if data.size:
        out["hbond_maks"] = float(np.max(data[:, 1]))
        out["hbond_frac_ada"] = float(np.mean(data[:, 1] > 0))
    return out


def _atom_residue_map(pdb_path):
    """Map 1-based atom serial to (chain, resid, resname) from a dumped frame."""
    mapping, serial = {}, 0
    with open(pdb_path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                serial += 1
                mapping[serial] = (line[21].strip() or "-", line[22:26].strip(),
                                   line[17:20].strip())
    return mapping


def analyse_contacts(ctx):
    """Which residues actually touch the ligand, and for what fraction of the run.

    A single contact map from the final frame says nothing about persistence.
    What matters for picking key residues is occupancy: the share of frames in
    which a residue sits within the cutoff of the ligand.
    """
    if "Other" not in ctx["groups"]:
        return {}
    cutoff = ctx["contact_cutoff"]
    sel = (f'group "Protein" and same residue as within {cutoff} of group "Other"')
    size_xvg, idx_dat = ctx["a"]("contacts_size.xvg"), ctx["a"]("contacts_index.dat")

    ok, _ = gmx_run(["select", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-select", sel, "-os", size_xvg, "-oi", idx_dat]
                    + stride_args(ctx, "contacts"), cwd=ctx["cdir"])
    if not ok or not os.path.exists(idx_dat):
        return {}

    data, meta = read_xvg(size_xvg)
    if data.size:
        data = data[:, :2]
    export(data, meta, ctx["adir"], "contacts_count",
           f"Atom protein dalam {cutoff} nm dari ligan",
           "waktu (ps)", "jumlah atom", ["waktu_ps", "n_atom"])

    amap = _atom_residue_map(ctx["a"]("start.pdb")) if os.path.exists(ctx["a"]("start.pdb")) else {}
    counts, frames = {}, 0
    with open(idx_dat, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            frames += 1
            seen = set()
            for token in parts[2:]:
                try:
                    key = amap.get(int(token))
                except ValueError:
                    continue
                if key and key not in seen:
                    seen.add(key)
                    counts[key] = counts.get(key, 0) + 1

    if not frames or not counts:
        return {}
    rows = sorted(counts.items(), key=lambda kv: -kv[1])
    with open(ctx["a"]("contacts_residue.csv"), "w") as fh:
        fh.write("rantai,nomor,residu,frame_kontak,okupansi\n")
        for (chain, resid, resname), n in rows:
            fh.write(f"{chain},{resid},{resname},{n},{n / frames:.4f}\n")

    top = rows[:15][::-1]
    fig, ax = plt.subplots(figsize=(7.5, max(3.0, 0.28 * len(top))))
    ax.barh([f"{r[0][2]}{r[0][1]}:{r[0][0]}" for r in top],
            [r[1] / frames for r in top], color="#4c78a8")
    ax.set_xlabel("okupansi kontak")
    ax.set_title(f"Residu kontak ligan (cutoff {cutoff} nm)")
    ax.grid(axis="x", alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(ctx["a"]("contacts_residue.png"), dpi=150)
    plt.close(fig)

    mean, _ = stat(data)
    return {
        "kontak_atom_rerata": mean,
        "kontak_residu_kunci": " ".join(
            f"{c[2]}{c[1]}:{c[0]}({n / frames:.2f})" for c, n in rows[:8]),
    }


def analyse_pca(ctx):
    """Essential dynamics: covariance of C-alpha motion, projected on PC1/PC2."""
    ca = ctx["groups"].get("C-alpha")
    if ca is None:
        return {}
    eigvec = ctx["a"]("eigenvec.trr")
    # covar and anaeig must see the same frames: the projection is only
    # meaningful against the covariance it came from.
    step = stride_args(ctx, "pca")
    ok, _ = gmx_run(["covar", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-o", ctx["a"]("eigenval.xvg"), "-v", eigvec,
                     "-av", ctx["a"]("average.pdb"), "-l", ctx["a"]("covar.log")] + step,
                    stdin=f"{ca}\n{ca}\n", cwd=ctx["cdir"])
    if not ok:
        return {}

    vals, meta = read_xvg(ctx["a"]("eigenval.xvg"))
    out = {}
    if vals.size:
        total = float(np.sum(vals[:, 1]))
        if total > 0:
            out["pca_pc1_persen"] = round(100 * vals[0, 1] / total, 2)
            out["pca_pc2_persen"] = round(100 * vals[1, 1] / total, 2) if len(vals) > 1 else None
            out["pca_pc1_pc2_persen"] = round(
                100 * float(np.sum(vals[:2, 1])) / total, 2)
        export(vals[:20], meta, ctx["adir"], "pca_eigenvalue",
               "Eigenvalue PCA (20 mode pertama)", "mode", "eigenvalue (nm^2)",
               ["mode", "eigenvalue_nm2"])

    proj = ctx["a"]("pca_proj_1_2.xvg")
    ok, _ = gmx_run(["anaeig", "-v", eigvec, "-s", ctx["tpr"], "-f", ctx["xtc"],
                     "-n", ctx["ndx"], "-first", "1", "-last", "2", "-2d", proj] + step,
                    stdin=f"{ca}\n{ca}\n", cwd=ctx["cdir"])
    if ok and os.path.exists(proj):
        data, _ = read_xvg(proj)
        if data.size >= 2:
            fig, ax = plt.subplots(figsize=(5.4, 5.0))
            sc = ax.scatter(data[:, 0], data[:, 1], c=np.arange(len(data)),
                            cmap="viridis", s=6)
            ax.set_xlabel("PC1 (nm)")
            ax.set_ylabel("PC2 (nm)")
            ax.set_title("Proyeksi PCA")
            fig.colorbar(sc, ax=ax, label="urutan frame")
            fig.tight_layout()
            fig.savefig(ctx["a"]("pca_proj.png"), dpi=150)
            plt.close(fig)
            np.savetxt(ctx["a"]("pca_proj_1_2.csv"), data[:, :2], delimiter=",",
                       header="pc1_nm,pc2_nm", comments="", fmt="%.6g")
    return out


def analyse_fel(ctx):
    """Free energy landscape over the PC1/PC2 projection.

    Runs only after PCA, because it is that projection that gets binned. The
    minimum of the surface is reported so that the lowest-energy conformer can
    be pulled out of the trajectory afterwards.
    """
    proj = ctx["a"]("pca_proj_1_2.xvg")
    if not os.path.exists(proj):
        return {}
    xpm = ctx["a"]("fel_gibbs.xpm")
    ok, _ = gmx_run(["sham", "-f", proj, "-ls", xpm, "-notime",
                     "-lsh", ctx["a"]("fel_enthalpy.xpm"),
                     "-lss", ctx["a"]("fel_entropy.xpm")], cwd=ctx["cdir"])
    if not ok or not os.path.exists(xpm):
        return {}

    values, xaxis, yaxis = read_xpm(xpm)
    if values is None:
        return {}
    np.savetxt(ctx["a"]("fel_gibbs.csv"), values, delimiter=",", fmt="%.6g")

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    extent = None
    if xaxis is not None and yaxis is not None and len(xaxis) and len(yaxis):
        extent = [xaxis.min(), xaxis.max(), yaxis.min(), yaxis.max()]
    im = ax.imshow(values, origin="lower", aspect="auto", extent=extent, cmap="jet")
    ax.set_xlabel("PC1 (nm)")
    ax.set_ylabel("PC2 (nm)")
    ax.set_title("Free energy landscape")
    fig.colorbar(im, ax=ax, label="G (kJ/mol)")
    fig.tight_layout()
    fig.savefig(ctx["a"]("fel_gibbs.png"), dpi=150)
    plt.close(fig)

    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {}
    return {"fel_min_kJmol": float(np.min(finite)),
            "fel_max_kJmol": float(np.max(finite))}


def _resolve_mmpbsa(setting):
    """Accept an env prefix, a python interpreter, or the executable itself."""
    if not setting:
        return shutil.which("gmx_MMPBSA")
    if os.path.isdir(setting):
        for candidate in (os.path.join(setting, "bin", "gmx_MMPBSA"),
                          os.path.join(setting, "gmx_MMPBSA")):
            if os.access(candidate, os.X_OK):
                return candidate
        return None
    if os.access(setting, os.X_OK):
        if os.path.basename(setting).startswith("python"):
            sibling = os.path.join(os.path.dirname(setting), "gmx_MMPBSA")
            return sibling if os.access(sibling, os.X_OK) else None
        return setting
    return None


def _auto_mmpbsa_np(n_frames):
    """Rank count to use when mmpbsa_np is 0.

    gmx_MMPBSA splits frames across ranks, so more ranks than frames is waste.
    Half the machine's cores, capped at eight, and never above the frame count --
    the same rule LADEEP's worker applies, so the two agree on any machine and
    not only on a large one.

    Capped rather than uncapped because the scaling flattens: the ranks share
    one filesystem and one trajectory reader, and past eight the coordination
    costs more than the extra parallelism returns.
    """
    return max(1, min(n_frames, max(1, (os.cpu_count() or 4) // 2), 8))


def _resolve_mpirun(mmpbsa_exe):
    """mpirun from the same environment as gmx_MMPBSA, or None.

    gmx_MMPBSA splits *frames* across MPI ranks, so this is the one knob that
    makes it finish sooner. It was never used: the command was built without
    mpirun at all, and one MM/GBSA pass over 181 frames took 98 minutes on a
    machine with 40 idle cores.

    The launcher has to come from gmx_MMPBSA's own environment. A system mpirun
    against a conda-built mpi4py is the classic way to get ranks that each think
    they are rank 0, which silently computes the same frames N times.
    """
    if not mmpbsa_exe:
        return None
    sibling = os.path.join(os.path.dirname(mmpbsa_exe), "mpirun")
    return sibling if os.access(sibling, os.X_OK) else None


def _frame_count(ctx):
    """Number of frames in the prepared trajectory, via the RMSD trace if it
    exists (cheap) and gmx check otherwise."""
    rmsd = ctx["a"]("rmsd_protein.xvg")
    if os.path.exists(rmsd):
        data, _ = read_xvg(rmsd)
        if data.size:
            return len(data)
    ok, out = gmx_run(["check", "-f", ctx["xtc"]], cwd=ctx["cdir"])
    m = re.search(r"Step\s+(\d+)", out or "")
    return int(m.group(1)) if m else 0


def analyse_mmpbsa(ctx):
    """Binding free energy with gmx_MMPBSA.

    Kept opt-in behind mmpbsa_python because it needs its own conda
    environment and, unlike everything else here, it costs real time: an
    end-state calculation over a few hundred frames is minutes to hours, not
    seconds.
    """
    exe = ctx["mmpbsa_exe"]
    if not exe:
        return {}
    if "Other" not in ctx["groups"] or "Protein" not in ctx["groups"]:
        return {}
    # Run once per ligand copy when the system holds several, so the answer is
    # a per-site binding energy with a spread, not one lumped number.
    copies = ctx["groups"].get("_ligan_terpisah") or []
    targets = [(name, ctx["groups"][name]) for name in copies] or [("Other", ctx["groups"]["Other"])]
    topol = os.path.join(ctx["cdir"], "topol.top")
    if not os.path.exists(topol):
        say("topol.top not found, MM/PBSA skipped", 6)
        return {}

    total = _frame_count(ctx)
    if total < 2:
        return {}
    # An explicit stride beats a target count, because the stride is the thing
    # with physical meaning: frames closer together than the system's
    # correlation time cost time without adding information, and a standard
    # error computed over them looks smaller than it is.
    if ctx["mmpbsa_interval"]:
        interval = max(1, ctx["mmpbsa_interval"])
    else:
        wanted = max(1, min(ctx["mmpbsa_frames"], total))
        interval = max(1, total // wanted)
    sampled = len(range(1, total + 1, interval))

    method = ctx["mmpbsa_method"].lower()
    say(f"gmx_MMPBSA: {total} frame tersedia, interval {interval} -> "
        f"{sampled} frame disampling; {len(targets)} salinan ligan "
        f"({', '.join(name for name, _ in targets)})", 6)
    per_copy, out = [], {}
    for label, group_num in targets:
        got = _run_one_mmpbsa(ctx, exe, label, group_num, total, interval, method)
        if got:
            per_copy.append((label, got))
    if not per_copy:
        return {}
    for key in sorted({k for _, g in per_copy for k in g}):
        vals = [g[key] for _, g in per_copy if key in g]
        if not vals or not key.endswith("_kcal"):
            continue
        out[key] = sum(vals) / len(vals)
        if len(vals) > 1:
            mean = out[key]
            out[key + "_sd_antar_situs"] = (
                sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
    for label, g in per_copy:
        for key, val in g.items():
            out[f"{key}_{label}"] = val
    return out


def _run_one_mmpbsa(ctx, exe, label, group_num, total, interval, method):
    """One gmx_MMPBSA run: the protein against a single ligand group."""
    blocks = [
        "&general",
        f'sys_name="{os.path.basename(ctx["cdir"])}_{label}",',
        f"startframe=1, endframe={total}, interval={interval},",
        "verbose=2,",
        "/",
    ]
    if method in ("gb", "both", "gbsa"):
        blocks += ["&gb", f"igb={ctx['mmpbsa_igb']}, saltcon={ctx['mmpbsa_salt']},", "/"]
    if method in ("pb", "both", "pbsa"):
        blocks += ["&pb", f"istrng={ctx['mmpbsa_salt']}, inp=2, radiopt=0,", "/"]
    infile = ctx["a"](f"mmpbsa_{label}.in")
    with open(infile, "w") as fh:
        fh.write("\n".join(blocks) + "\n")

    env = dict(os.environ)
    env["PATH"] = os.path.dirname(os.path.abspath(GMX)) + os.pathsep + env.get("PATH", "")
    cmd = [exe, "-O", "-i", infile,
           "-cs", ctx["tpr"], "-ci", ctx["ndx"],
           "-cg", str(ctx["groups"]["Protein"]), str(group_num),
           "-ct", ctx["xtc"], "-cp", "topol.top",
           "-o", ctx["a"](f"mmpbsa_{label}.dat"),
           "-eo", ctx["a"](f"mmpbsa_{label}.csv"), "-nogui"]
    # Frames are the unit of work, so more ranks than frames is waste.
    requested = int(ctx.get("mmpbsa_np", 0) or 0)
    nprocs = requested if requested > 0 else _auto_mmpbsa_np(total)
    nprocs = max(1, min(nprocs, total))
    if nprocs > 1:
        launcher = _resolve_mpirun(exe)
        if launcher:
            cmd = [launcher, "-np", str(nprocs)] + cmd
        else:
            say("mpirun not found next to gmx_MMPBSA; running serial", 8)
            nprocs = 1
    say(f"gmx_MMPBSA [{label}]: interval {interval}, metode {method}, "
        f"{nprocs} proses", 8)
    try:
        proc = subprocess.run(cmd, cwd=ctx["cdir"], capture_output=True,
                              text=True, env=env, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        say(f"gmx_MMPBSA [{label}] failed: {exc}", 6)
        return {}
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        say(f"gmx_MMPBSA [{label}] failed: " + " | ".join(tail), 6)
        return {}

    dat = ctx["a"](f"mmpbsa_{label}.dat")
    return parse_mmpbsa_dat(dat) if os.path.exists(dat) else {}


# Matches the delta row only. A .dat holds four TOTAL rows -- complex,
# receptor, ligand, then the difference -- and gmx_MMPBSA writes the last as
# "ΔTOTAL", one word, not "Delta TOTAL".
_DELTA_TOTAL = re.compile(r"^\s*(?:Δ|DELTA\s+|Delta\s+)TOTAL\s+([-\d.]+)\s+([-\d.]+)")


def parse_mmpbsa_dat(path):
    """Binding free energy and its standard deviation from a gmx_MMPBSA .dat.

    Returns {'mmpbsa_gb_dG_kcal', 'mmpbsa_gb_sd'} and/or the pb equivalents.

    Only the ΔTOTAL row counts. The previous pattern, `(?:Delta\\s+)?TOTAL`,
    matched the three plain TOTAL rows and missed the delta, and since each
    match overwrote the last, what was reported as the binding free energy was
    the ligand's own internal energy -- +70.81 kcal/mol where the answer was
    -35.19. Wrong quantity, wrong sign, and nothing looked broken: those
    numbers even ranked the candidates plausibly, because a ligand's internal
    energy grows with its size.
    """
    out, section = {}, None
    try:
        fh = open(path, errors="replace")
    except OSError:
        return out
    with fh:
        for line in fh:
            if "GENERALIZED BORN" in line:
                section = "gb"
            elif "POISSON BOLTZMANN" in line:
                section = "pb"
            m = _DELTA_TOTAL.match(line)
            if m and section:
                out[f"mmpbsa_{section}_dG_kcal"] = float(m.group(1))
                out[f"mmpbsa_{section}_sd"] = float(m.group(2))
    return out


ANALYSIS_FUNCS = {
    "rmsd": analyse_rmsd,
    "rmsf": analyse_rmsf,
    "rg": analyse_rg,
    "sasa": analyse_sasa,
    "hbond": analyse_hbond,
    "contacts": analyse_contacts,
    "pca": analyse_pca,
    "fel": analyse_fel,
    "mmpbsa": analyse_mmpbsa,
}


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def read_config(directory):
    """Same parser LAGMX.py uses, plus the analysis defaults."""
    values = dict(ANALYSIS_DEFAULTS)
    path = os.path.join(directory, "gmx_config.txt")
    if not os.path.exists(path):
        sys.exit(
            f"LAGMX analyze: no 'gmx_config.txt' in {directory}\n\n"
            "analyze_md.py reads its configuration and complex_*/ directories "
            "from the\ndirectory you run it from, exactly like LAGMX.py. Run it "
            "from the same\nplace you ran the simulation:\n\n"
            "    cd run_matrix && python3 ../analyze_md.py\n"
        )
    for line in open(path, errors="replace"):
        if line.strip().startswith("#"):
            continue
        try:
            name, value = line.strip().split(": ", 1)
        except ValueError:
            continue
        values[name.strip()] = value.strip()
    return values


def analyse_complex(cdir, cfg, requested):
    """Run the requested analyses on one complex directory."""
    name = os.path.basename(cdir)
    tpr, xtc = os.path.join(cdir, "md.tpr"), os.path.join(cdir, "md.xtc")
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        say(f"{name}: no finished production run (md.tpr/md.xtc missing), skipped")
        return None

    adir = os.path.join(cdir, "analysis")
    os.makedirs(adir, exist_ok=True)
    say(f"{name}")

    merged = cfg.get("merged_group", "Protein_LIG")
    groups = build_index(tpr, os.path.join(adir, "analysis.ndx"), cdir, merged)
    if not groups:
        say("make_ndx failed, skipped", 6)
        return None

    skip_ps = float(cfg["analysis_skip_ns"]) * 1000.0
    prepared = prepare_trajectory(cdir, adir, tpr, xtc, groups, merged, skip_ps)
    if not prepared:
        return None

    frame_ps = frame_interval_ps(cdir)
    strides = {}
    for key in STRIDABLE:
        try:
            strides[key] = max(1, int(float(cfg.get(f"analysis_stride_{key}", "1"))))
        except ValueError:
            strides[key] = 1
    thinned = {k: v for k, v in strides.items() if v > 1}
    if thinned and frame_ps <= 0:
        say("cannot read nstxout-compressed*dt from md.mdp; "
            "analysing every frame instead of striding", 6)
    elif thinned:
        say("stride: " + ", ".join(f"{k} 1/{v}" for k, v in thinned.items())
            + f" (frame {frame_ps:g} ps)", 6)

    # Fit against coordinates that match the trajectory, not the raw ones.
    ref_tpr = build_reference_tpr(cdir, adir, tpr)

    ctx = {
        "cdir": cdir, "adir": adir, "tpr": ref_tpr, "xtc": prepared,
        "ndx": os.path.join(adir, "analysis.ndx"), "groups": groups,
        "a": lambda f: os.path.join(adir, f),
        "strides": strides, "frame_ps": frame_ps,
        "contact_cutoff": float(cfg["analysis_contact_cutoff"]),
        "mmpbsa_exe": _resolve_mmpbsa(cfg.get("mmpbsa_python", "").strip()),
        "mmpbsa_method": cfg.get("mmpbsa_method", "gb"),
        "mmpbsa_frames": int(cfg.get("mmpbsa_frames", "100")),
        "mmpbsa_interval": int(cfg.get("mmpbsa_interval") or 0),
        "mmpbsa_np": int(cfg.get("mmpbsa_np", "1") or 1),
        "mmpbsa_igb": cfg.get("mmpbsa_igb", "5"),
        "mmpbsa_salt": cfg.get("mmpbsa_salt", "0.150"),
    }

    summary = {"complex": name}
    for key in requested:
        func = ANALYSIS_FUNCS.get(key)
        if func is None:
            continue
        say(f"{key} ...", 6)
        try:
            summary.update(func(ctx) or {})
        except Exception as exc:                                # noqa: BLE001
            say(f"{key} failed: {exc}", 8)

    with open(os.path.join(adir, "summary.csv"), "w") as fh:
        fh.write("besaran,nilai\n")
        for key, value in summary.items():
            fh.write(f"{key},{value}\n")
    return summary


def main():
    global GMX
    directory = os.getcwd()
    cfg = read_config(directory)

    GMX = cfg.get("analysis_gmx", "gmx").strip() or "gmx"
    resolved = shutil.which(GMX) or (GMX if os.access(GMX, os.X_OK) else None)
    if not resolved:
        sys.exit(f"LAGMX analyze: '{GMX}' not found on PATH. Set analysis_gmx in gmx_config.txt.")
    GMX = resolved

    requested = [a.strip().lower() for a in cfg["analysis"].split(",") if a.strip()]
    if "all" in requested:
        requested = list(ALL_ANALYSES)
    unknown = [a for a in requested if a not in ANALYSIS_FUNCS]
    if unknown:
        sys.exit(f"LAGMX analyze: unknown analysis {unknown}; choose from {ALL_ANALYSES}")

    ok, version = gmx_run(["--version"])
    banner = next((l.strip() for l in version.splitlines() if "GROMACS version" in l), "?")
    say(f"gmx      : {GMX}")
    say(f"           {banner}")
    say(f"analisis : {', '.join(requested)}")
    say(f"buang    : {cfg['analysis_skip_ns']} ns pertama")

    complex_dirs = sorted(d for d in glob.glob(os.path.join(directory, "complex*"))
                          if os.path.isdir(d))
    if not complex_dirs:
        sys.exit("LAGMX analyze: no complex*/ directories here.")
    say(f"kompleks : {', '.join(os.path.basename(d) for d in complex_dirs)}\n")

    rows = [r for r in (analyse_complex(d, cfg, requested) for d in complex_dirs) if r]
    if not rows:
        say("\nTidak ada kompleks yang bisa dianalisis.")
        return 1

    columns, seen = [], set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                columns.append(key)
    out_csv = os.path.join(directory, "analysis_summary.csv")
    with open(out_csv, "w") as fh:
        fh.write(",".join(columns) + "\n")
        for row in rows:
            fh.write(",".join(str(row.get(c, "")) for c in columns) + "\n")

    # Printed transposed: metrics down the side, complexes across. With four
    # or five systems a metric-per-row table is the one you can actually read.
    headline = [
        ("rmsd_protein_nm", "RMSD protein (nm)"),
        ("rmsd_ligand_nm", "RMSD ligan (nm)"),
        ("rmsf_mean_nm", "RMSF rerata (nm)"),
        ("rg_nm", "Rg (nm)"),
        ("sasa_terkubur_nm2", "SASA terkubur (nm2)"),
        ("hbond_rerata", "H-bond rerata"),
        ("kontak_atom_rerata", "kontak atom rerata"),
        ("pca_pc1_pc2_persen", "PC1+PC2 (%)"),
        ("fel_min_kJmol", "FEL min (kJ/mol)"),
        ("mmpbsa_gb_dG_kcal", "MM/GBSA dG (kcal/mol)"),
        ("mmpbsa_pb_dG_kcal", "MM/PBSA dG (kcal/mol)"),
    ]
    width = max(14, *(len(r["complex"]) + 2 for r in rows))
    say("\n=================== RINGKASAN ===================")
    say(f"{'besaran':<24}" + "".join(f"{r['complex']:>{width}}" for r in rows))
    say("-" * (24 + width * len(rows)))
    for key, label in headline:
        if key not in columns:
            continue
        cells = []
        for row in rows:
            v = row.get(key, "")
            cells.append(f"{v:>{width}.3f}" if isinstance(v, float) else f"{str(v):>{width}}")
        say(f"{label:<24}" + "".join(cells))

    for row in rows:
        key_res = row.get("kontak_residu_kunci")
        if key_res:
            say(f"\nresidu kunci {row['complex']}: {key_res}")
    say(f"\n-> {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
