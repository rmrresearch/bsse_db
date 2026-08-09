#!/usr/bin/env bash
#
# make_tree_2_new_v2.sh
#
# Creates the directory structure for the Ne dimer BSSE study, version 2.
#
#   Grid  : 50 uniform points, R in [1.5, 3.5] A, dR = 0.040816 A
#   Bases : aug-cc-pvdz, aug-cc-pvtz, aug-cc-pvqz
#   Theory: CCSD(T)  (SCF, MP2, CCSD(T) all extracted from the same output)
#
# Run from:  <repo>/data/Ne/
# Usage   :  bash make_tree_2_new_v2.sh
#
# Idempotent: safe to re-run.  Never touches 2_new_v1/.

set -euo pipefail

ROOT="2_new_v2"
BASES=(aug-cc-pvdz aug-cc-pvtz aug-cc-pvqz)
LEVEL="CCSD_T"

if [[ ! -d "2_new_v1" ]]; then
    echo "ERROR: expected to find 2_new_v1/ in $(pwd)"
    echo "       Run this from the data/Ne/ directory."
    exit 1
fi

echo "Creating ${ROOT}/ under $(pwd)"

# --- top level -------------------------------------------------------------
mkdir -p "${ROOT}"/scripts
mkdir -p "${ROOT}"/R_config

# --- analysis folders ------------------------------------------------------
mkdir -p "${ROOT}"/2b_bsse_scf/fit_study/figures
mkdir -p "${ROOT}"/2b_bsse_scf/basis_study/figures
mkdir -p "${ROOT}"/2b_bsse_mp2c/figures
mkdir -p "${ROOT}"/2b_bsse_ccsd_tc/figures
mkdir -p "${ROOT}"/bsse_fits/figures

# --- geometry grid ---------------------------------------------------------
# 50 points: R_i = 1.5 + i * (2.0/49), i = 0..49
n=0
for i in $(seq 0 49); do
    n=$((n+1))
    R=$(python3 -c "print(f'{1.5 + $i * (2.0/49):.3f}')")
    LABEL=$(printf "%03d_R_%s" "$n" "${R/./p}")
    for B in "${BASES[@]}"; do
        mkdir -p "${ROOT}/R_config/${LABEL}/${LEVEL}/${B}"
    done
done

echo "  R_config/ : ${n} R-points x ${#BASES[@]} bases = $((n * ${#BASES[@]})) leaf directories"
echo
echo "Structure:"
echo "  ${ROOT}/"
echo "  |-- README.md            (write this next)"
echo "  |-- scripts/"
echo "  |-- R_config/NNN_R_XpXXX/${LEVEL}/<basis>/"
echo "  |-- 2b_bsse_scf/"
echo "  |     |-- fit_study/     (functional form from GPT/Boys)"
echo "  |     '-- basis_study/   (DZ/TZ/QZ ladder)"
echo "  |-- 2b_bsse_mp2c/"
echo "  |-- 2b_bsse_ccsd_tc/"
echo "  '-- bsse_fits/           (sum of components vs. single-form fits)"
echo
echo "Done."
