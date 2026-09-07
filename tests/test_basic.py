"""Lightweight tests that don't need GROMACS/AmberTools installed.

These check repository structure, config/mdp sanity, and the pure-Python
logic in fix_ligand_valence.py. They do NOT run the actual MD pipeline --
that needs GROMACS/AmberTools and is verified manually against run_matrix/
(see CONTRIBUTING.md).
"""
import glob
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PY_FILES = [
    "LAGMX.py",
    "fix_structure.py",
    "analyze_md.py",
    "run_matrix/fix_ligand_valence.py",
]

# complex_6broken is intentionally malformed and expected to fail the real
# pipeline; it's excluded from checks that assume a working scenario.
KNOWN_BROKEN_SCENARIOS = {"complex_6broken"}


def test_python_files_exist_and_compile():
    for rel in PY_FILES:
        path = os.path.join(ROOT, rel)
        assert os.path.exists(path), f"missing {rel}"
        subprocess.run([sys.executable, "-m", "py_compile", path], check=True)


def test_only_one_copy_of_lagmx_exists():
    # LAGMX.py used to be duplicated into run_matrix/ and the two
    # copies drifted apart during development (one had checkpoint-resume,
    # the other had SEQRES curation, neither had both). There must be
    # exactly one canonical copy from here on -- run_matrix/run_lagmx.sh
    # invokes it via a relative path (../LAGMX.py) instead of a second
    # copy living next to the test scenarios.
    assert os.path.exists(os.path.join(ROOT, "LAGMX.py"))
    assert not os.path.exists(os.path.join(ROOT, "run_matrix", "LAGMX.py")), (
        "run_matrix/LAGMX.py should not exist -- there must be exactly "
        "one copy of the script (repo root); see run_lagmx.sh, which "
        "invokes ../LAGMX.py"
    )
    assert not os.path.exists(os.path.join(ROOT, "run_matrix", "fix_structure.py")), (
        "run_matrix/fix_structure.py should not exist -- fix_structure.py "
        "is resolved relative to LAGMX.py's own location, not cwd, so "
        "it only needs to exist once, at the repo root"
    )


def test_run_matrix_scenarios_have_inputs():
    matrix = os.path.join(ROOT, "run_matrix")
    complex_dirs = sorted(
        d for d in os.listdir(matrix)
        if os.path.isdir(os.path.join(matrix, d)) and d.startswith("complex_")
    )
    assert len(complex_dirs) >= 6, "expected at least 6 scenarios in run_matrix/"
    for d in complex_dirs:
        p = os.path.join(matrix, d)
        rec = glob.glob(os.path.join(p, "rec*.pdb"))
        lig = glob.glob(os.path.join(p, "lig*.pdb")) + glob.glob(os.path.join(p, "lig*.mol2"))
        assert rec, f"{d}: no receptor file (rec*.pdb)"
        assert lig, f"{d}: no ligand file (lig*.pdb / lig*.mol2)"


def test_gmx_config_has_analysis_keys():
    # analyze_md.py reads its settings from the same gmx_config.txt as the
    # simulation. Every one of them is optional at runtime, but the worked
    # example has to document them or nobody discovers the feature exists.
    required = ["analysis", "analysis_skip_ns", "analysis_contact_cutoff",
                "mmpbsa_method", "mmpbsa_frames"]
    path = os.path.join(ROOT, "run_matrix", "gmx_config.txt")
    with open(path) as f:
        text = f.read()
    for key in required:
        assert f"{key}:" in text, f"run_matrix/gmx_config.txt: missing '{key}'"


def test_gmx_config_has_required_keys():
    required = ["box_type", "distance", "solvent", "charge", "atom_type"]
    path = os.path.join(ROOT, "run_matrix", "gmx_config.txt")
    with open(path) as f:
        text = f.read()
    for key in required:
        assert f"{key}:" in text, f"run_matrix/gmx_config.txt: missing required key '{key}'"


def test_mdp_templates_have_required_params():
    for name in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp", "ions.mdp"]:
        path = os.path.join(ROOT, "run_matrix", name)
        assert os.path.exists(path), f"missing {path}"
        with open(path) as f:
            text = f.read()
        assert "integrator" in text, f"{path}: missing 'integrator ='"


def test_fix_ligand_valence_detects_missing_hydrogen(tmp_path):
    sys.path.insert(0, os.path.join(ROOT, "run_matrix"))
    import fix_ligand_valence as flv

    # Minimal synthetic MOL2: atom 1 is C.3 (expects 4 neighbours) but only
    # has 3 -- the same class of defect fix_ligand_valence.py was built to
    # repair (see README: complex_4mol2 needed exactly this fix).
    mol2 = """@<TRIPOS>MOLECULE
test
4 3 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1  C1        0.000    0.000    0.000 C.3     1  UNL1    0.0000
      2  C2        1.500    0.000    0.000 C.3     1  UNL1    0.0000
      3  N1       -0.750    1.299    0.000 N.3     1  UNL1    0.0000
      4  C3       -0.750   -1.299    0.000 C.3     1  UNL1    0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
     3     1     4    1
"""
    p = tmp_path / "synthetic.mol2"
    p.write_text(mol2)

    _, _, _, _, atom_lines, bond_lines = flv.parse_mol2(str(p))
    result = flv.find_missing_h_atom(atom_lines, bond_lines)

    assert result is not None, "expected atom 1 (C.3, only 3 neighbours) to be flagged"
    atom_id, deficit = result
    assert atom_id == 1
    assert deficit == 1


def _load_analyze_md():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "analyze_md", os.path.join(ROOT, "analyze_md.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_read_xvg_parses_data_and_legends(tmp_path):
    xvg = tmp_path / "rmsd.xvg"
    xvg.write_text(
        '# comment\n'
        '@    title "RMSD"\n'
        '@    xaxis  label "Time (ns)"\n'
        '@    yaxis  label "RMSD (nm)"\n'
        '@ s0 legend "Backbone"\n'
        '0.000    0.100\n'
        '1.000    0.150\n'
        '2.000    0.200\n'
    )
    am = _load_analyze_md()
    data, meta = am.read_xvg(str(xvg))
    assert data.shape == (3, 2)
    assert meta["legends"] == ["Backbone"]
    assert meta["xlabel"] == "Time (ns)"
    assert abs(float(data[2, 1]) - 0.2) < 1e-9


def test_read_xpm_decodes_values_from_the_comment(tmp_path):
    # gmx sham hides the kJ/mol value of each colour in the C comment after
    # the colour definition, not inside the quoted colour string. Reading the
    # quoted part instead returns a matrix of NaN that still plots, still
    # writes a CSV, and is entirely meaningless.
    xpm = tmp_path / "gibbs.xpm"
    xpm.write_text(
        '/* XPM */\n'
        '/* title:   "Gibbs Energy Landscape" */\n'
        'static char *gromacs_xpm[] = {\n'
        '"3 2   2 1",\n'
        '"A  c #000000 " /* "0" */,\n'
        '"B  c #FFFFFF " /* "12.5" */,\n'
        '/* x-axis: 1 2 3 */\n'
        '/* y-axis: 1 2 */\n'
        '"ABA",\n'
        '"BBB"\n'
        '};\n'
    )
    am = _load_analyze_md()
    values, xaxis, yaxis = am.read_xpm(str(xpm))
    assert values is not None, "XPM must decode, not come back as None"
    assert values.shape == (2, 3)
    # Rows are flipped: XPM runs top to bottom, the y axis bottom to top.
    assert list(values[0]) == [12.5, 12.5, 12.5]
    assert list(values[1]) == [0.0, 12.5, 0.0]
    assert list(xaxis) == [1.0, 2.0, 3.0]
    assert list(yaxis) == [1.0, 2.0]


def test_analysis_defaults_cover_every_named_analysis():
    am = _load_analyze_md()
    for name in am.ALL_ANALYSES:
        assert name in am.ANALYSIS_FUNCS, f"{name} has no implementation"
    # The default list must be runnable without extra setup, so mmpbsa -- which
    # needs its own conda environment -- stays opt-in.
    assert "mmpbsa" not in am.ANALYSIS_DEFAULTS["analysis"].split(",")


def test_mmpbsa_reads_the_delta_row_not_the_ligand(tmp_path):
    """A .dat holds four TOTAL rows: complex, receptor, ligand, then the
    difference. gmx_MMPBSA writes the last as "ΔTOTAL", one word. Matching
    `(?:Delta\\s+)?TOTAL` caught the first three and missed the fourth, and
    since each match overwrote the last, the ligand's own internal energy was
    reported as the binding free energy -- positive, and ranking by ligand
    size rather than by affinity."""
    dat = tmp_path / "mmpbsa_Other_K13_258.dat"
    dat.write_text(
        "GENERALIZED BORN:\n"
        "\n"
        "Complex:\n"
        "GGAS                  -5235.26        210.75\n"
        "TOTAL                -10461.67        260.21\n"
        "\n"
        "Receptor:\n"
        "TOTAL                -10497.29        259.70\n"
        "\n"
        "Ligand:\n"
        "TOTAL                    70.81          5.51\n"
        "\n"
        "Delta (Complex - Receptor - Ligand):\n"
        "ΔGGAS                   -50.71          1.17\n"
        "ΔTOTAL                  -35.19          1.20\n"
    )
    am = _load_analyze_md()
    out = am.parse_mmpbsa_dat(str(dat))
    assert out["mmpbsa_gb_dG_kcal"] == -35.19, "must be ΔTOTAL, not a plain TOTAL"
    assert out["mmpbsa_gb_sd"] == 1.20
    # complex - receptor - ligand reproduces the delta, so the row is the
    # binding energy and not one of the three absolute totals.
    assert round(-10461.67 - (-10497.29) - 70.81, 2) == out["mmpbsa_gb_dG_kcal"]


def test_mmpbsa_missing_file_is_empty_not_an_error():
    am = _load_analyze_md()
    assert am.parse_mmpbsa_dat("/nonexistent/mmpbsa.dat") == {}
