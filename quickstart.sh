#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# LAGMX quickstart -- "does this work on my machine?" in a few minutes.
#
#   ./quickstart.sh
#
# Runs the smallest scenario from the test matrix (one receptor chain, one
# ligand) end to end: preparation, minimization, NVT, NPT, and a short
# production run. Everything is deliberately shortened -- Gasteiger charges
# instead of AM1-BCC, 10 ps of equilibration instead of 100 ps, 10 ps of
# production instead of 1 ns -- so this is a smoke test, NOT a scientific
# result. For real runs use run_matrix/ and read the README.
#
# Output goes to quickstart_run/, which is gitignored and safe to delete.
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")"
ROOT="$PWD"
WORK="$ROOT/quickstart_run"

echo "=================== PREFLIGHT ==================="
fail=0
for tool in gmx antechamber parmchk2 tleap; do
    if command -v "$tool" >/dev/null 2>&1; then
        printf '  %-12s %s\n' "$tool" "$(command -v "$tool")"
    else
        printf '  %-12s MISSING\n' "$tool"; fail=1
    fi
done
if python3 -c "import parmed" 2>/dev/null; then
    printf '  %-12s OK\n' "parmed"
else
    printf '  %-12s MISSING\n' "parmed"; fail=1
fi
if [ "$fail" -ne 0 ]; then
    cat >&2 <<'MSG'

Some requirements are missing. The quickest way to get all of them:

    conda env create -f environment.yml
    conda activate lagmx

MSG
    exit 1
fi
echo "================================================="
echo

rm -rf "$WORK"
mkdir -p "$WORK/complex_smoke"

# Inputs: the single-receptor, single-ligand scenario from the test matrix.
cp run_matrix/complex_1single/rec_erbb.pdb "$WORK/complex_smoke/"
cp run_matrix/complex_1single/lig_erbb.pdb "$WORK/complex_smoke/"

# Templates, with equilibration shortened from 100 ps to 10 ps.
cp run_matrix/em.mdp run_matrix/ions.mdp run_matrix/md.mdp "$WORK/"
for f in nvt npt; do
    sed 's/^nsteps.*/nsteps                  = 5000      ; 2 * 5000 = 10 ps (quickstart)/' \
        "run_matrix/$f.mdp" > "$WORK/$f.mdp"
done
cp "$WORK"/*.mdp "$WORK/complex_smoke/"

cat > "$WORK/gmx_config.txt" <<'CFG'
# Quickstart configuration -- tuned for speed, not for accuracy.
box_type: Cubic
distance: 1.0
solvent: spc216.gro
charge: gas
atom_type: gaff2
forcefield: amber03
water: tip3p
merged_group: Protein_LIG
mdrun_options:
net_charge: auto
fix_structure: no
fixer_python:
seqres_reference:
production_ns: 0.01
CFG

echo "Running LAGMX in $WORK ..."
echo
cd "$WORK"
python3 -u "$ROOT/LAGMX.py" 2>&1 | tee quickstart.log
status=${PIPESTATUS[0]}

echo
echo "=================== RESULT ======================"
printf '%-14s %-7s %-7s %-7s %-7s\n' STAGE EM NVT NPT MD
printf '%-14s' "complex_smoke"
missing=0
for f in em.gro nvt.gro npt.gro md.gro; do
    if [ -f "complex_smoke/$f" ]; then printf ' %-7s' "OK"; else printf ' %-7s' "-"; missing=1; fi
done
echo
echo
if [ "$status" -eq 0 ] && [ "$missing" -eq 0 ]; then
    echo "PASS -- LAGMX ran the full pipeline on this machine."
    echo "Next: see README.md, and use run_matrix/ for real simulations."
else
    echo "FAIL -- exit code $status. See $WORK/quickstart.log for the error."
    echo "If you think this is a bug, please open an issue and attach that log:"
    echo "  https://github.com/laode-aman-ung/lagmx/issues"
fi
echo "================================================="
exit "$status"
