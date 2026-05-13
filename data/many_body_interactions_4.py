"""
many_body_interactions_4.py
----------------------------
Reads raw_data_{molecule}_4.pickle and computes the many-body interaction
energies for the Ne or HF tetramer using the VMFC (Valiron-Mayer) scheme.

Interaction terms computed
--------------------------
For each method (SCF, MP2, CCSD_T) and basis (aug-cc-pvdz):

  2b  — six dimer pairs:  AB, AC, AD, BC, BD, CD
  3b  — four trimer triples: ABC, ABD, ACD, BCD
  4b  — one tetramer:     ABCD

Each stored as {'no_vmfc': float, 'vmfc': float}

Output dictionary structure
---------------------------
result[mol][cluster][config][method][basis] = {
    '2b': {
        'AB': {'no_vmfc': float, 'vmfc': float},
        'AC': {'no_vmfc': float, 'vmfc': float},
        ...
    },
    '3b': {
        'ABC': {'no_vmfc': float, 'vmfc': float},
        ...
    },
    '4b': {'no_vmfc': float, 'vmfc': float},
}

Formulas (Valiron & Mayer 1997, PI task list Issue #46)
-------------------------------------------------------
2b_no_vmfc(IJ)  = E(IJ|IJ) - E(I|I) - E(J|J)
2b_vmfc(IJ)     = E(IJ|IJ) - E(I|IJ) - E(J|IJ)

3b_no_vmfc(IJK) = E(IJK|IJK) - E(IJ|IJ) - E(IK|IK) - E(JK|JK)
                + E(I|I) + E(J|J) + E(K|K)

3b_vmfc(IJK)    = E(IJK|IJK) - E(I|IJK) - E(J|IJK) - E(K|IJK)
                - eps2(IJ|IJK) - eps2(IK|IJK) - eps2(JK|IJK)
  where eps2(IJ|IJK) = E(IJ|IJK) - E(I|IJK) - E(J|IJK)

4b_no_vmfc      = E(ABCD|ABCD)
                - [E(ABC|ABC) + E(ABD|ABD) + E(ACD|ACD) + E(BCD|BCD)]
                + [E(AB|AB) + E(AC|AC) + E(AD|AD)
                   + E(BC|BC) + E(BD|BD) + E(CD|CD)]
                - [E(A|A) + E(B|B) + E(C|C) + E(D|D)]

4b_vmfc         = E(ABCD|ABCD)
                - [E(A|ABCD) + E(B|ABCD) + E(C|ABCD) + E(D|ABCD)]
                - sum of eps2(IJ|ABCD) for all 6 pairs
                - sum of eps3(IJK|ABCD) for all 4 trimers
  where eps2(IJ|ABCD)  = E(IJ|ABCD) - E(I|ABCD) - E(J|ABCD)
        eps3(IJK|ABCD) = E(IJK|ABCD) - E(I|ABCD) - E(J|ABCD) - E(K|ABCD)
                       - eps2(IJ|ABCD) - eps2(IK|ABCD) - eps2(JK|ABCD)

Notation: E(X|Y) means fragment X computed in basis Y, i.e. file key X_Y.

Molecule-specific notes
-----------------------
    Ne  : B_B = C_C = D_D = A_A was injected at the build_raw_data stage.
          No special handling needed here — all 65 keys are present.

    HF  : All four monomers are independent. No special handling needed.

Usage
-----
    python3 many_body_interactions_4.py --mol Ne
    python3 many_body_interactions_4.py --mol HF
    # produces data/many_body_interactions_{mol}_4.pickle
"""

import pickle
import argparse
from pathlib import Path
from itertools import combinations


# ---------------------------------------------------------------------------
# Topology — same for Ne and HF
# ---------------------------------------------------------------------------

MONOMERS = ['A', 'B', 'C', 'D']
PAIRS    = [''.join(p) for p in combinations(MONOMERS, 2)]   # AB AC AD BC BD CD
TRIMERS  = [''.join(t) for t in combinations(MONOMERS, 3)]   # ABC ABD ACD BCD
TETRA    = 'ABCD'


# ---------------------------------------------------------------------------
# Energy helper
# ---------------------------------------------------------------------------

def total_energy(d, method, basis, key):
    """
    Reconstruct the total energy at the requested method level.

    raw_data stores:
        SCF    -> total SCF energy
        MP2    -> MP2 correlation (E_mp2 - E_scf)
        CCSD_T -> CCSD(T) component (E_ccsdt - E_mp2)

    Parameters
    ----------
    d      : dict  — raw_data[mol][cluster][config]
    method : str   — 'SCF', 'MP2', or 'CCSD_T'
    basis  : str
    key    : str   — file key, e.g. 'A_AB'
    """
    e = d['SCF'][basis][key]
    if method in ('MP2', 'CCSD_T'):
        e += d['MP2'][basis][key]
    if method == 'CCSD_T':
        e += d['CCSD_T'][basis][key]
    return e


# ---------------------------------------------------------------------------
# Many-body interaction builders
# ---------------------------------------------------------------------------

def compute_2b(e, pairs):
    """Compute 2b no_vmfc and vmfc for all pairs."""
    result = {}
    for ij in pairs:
        i, j = ij[0], ij[1]
        no_vmfc = e(f'{ij}_{ij}') - e(f'{i}_{i}') - e(f'{j}_{j}')
        vmfc    = e(f'{ij}_{ij}') - e(f'{i}_{ij}') - e(f'{j}_{ij}')
        result[ij] = {'no_vmfc': no_vmfc, 'vmfc': vmfc}
    return result


def compute_3b(e, trimers):
    """Compute 3b no_vmfc and vmfc for all trimer triples."""
    result = {}
    for ijk in trimers:
        i, j, k = ijk[0], ijk[1], ijk[2]
        ij, ik, jk = i+j, i+k, j+k

        # no_vmfc — own basis recursion
        no_vmfc = (
            e(f'{ijk}_{ijk}')
            - e(f'{ij}_{ij}') - e(f'{ik}_{ik}') - e(f'{jk}_{jk}')
            + e(f'{i}_{i}')   + e(f'{j}_{j}')   + e(f'{k}_{k}')
        )

        # vmfc — trimer basis CP corrections
        eps2_ij = e(f'{ij}_{ijk}') - e(f'{i}_{ijk}') - e(f'{j}_{ijk}')
        eps2_ik = e(f'{ik}_{ijk}') - e(f'{i}_{ijk}') - e(f'{k}_{ijk}')
        eps2_jk = e(f'{jk}_{ijk}') - e(f'{j}_{ijk}') - e(f'{k}_{ijk}')

        vmfc = (
            e(f'{ijk}_{ijk}')
            - e(f'{i}_{ijk}') - e(f'{j}_{ijk}') - e(f'{k}_{ijk}')
            - eps2_ij - eps2_ik - eps2_jk
        )

        result[ijk] = {'no_vmfc': no_vmfc, 'vmfc': vmfc}
    return result


def compute_4b(e, pairs, trimers):
    """Compute 4b no_vmfc and vmfc for the full tetramer."""
    # no_vmfc — inclusion-exclusion over own-basis energies
    no_vmfc = e('ABCD_ABCD')
    for ijk in trimers:
        no_vmfc -= e(f'{ijk}_{ijk}')
    for ij in pairs:
        no_vmfc += e(f'{ij}_{ij}')
    for i in MONOMERS:
        no_vmfc -= e(f'{i}_{i}')

    # vmfc — all sub-clusters evaluated in full tetramer basis
    # eps2(IJ|ABCD)
    eps2 = {}
    for ij in pairs:
        i, j = ij[0], ij[1]
        eps2[ij] = (
            e(f'{ij}_ABCD') - e(f'{i}_ABCD') - e(f'{j}_ABCD')
        )

    # eps3(IJK|ABCD)
    eps3 = {}
    for ijk in trimers:
        i, j, k = ijk[0], ijk[1], ijk[2]
        ij, ik, jk = i+j, i+k, j+k
        eps3[ijk] = (
            e(f'{ijk}_ABCD')
            - e(f'{i}_ABCD') - e(f'{j}_ABCD') - e(f'{k}_ABCD')
            - eps2[ij] - eps2[ik] - eps2[jk]
        )

    vmfc = (
        e('ABCD_ABCD')
        - sum(e(f'{i}_ABCD') for i in MONOMERS)
        - sum(eps2.values())
        - sum(eps3.values())
    )

    return {'no_vmfc': no_vmfc, 'vmfc': vmfc}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_interactions(molecule, cluster='4', config='0',
                       raw_data_path=None):
    here = Path(__file__).resolve().parent

    if raw_data_path is None:
        raw_data_path = here / f'raw_data_{molecule}_{cluster}.pickle'

    with open(raw_data_path, 'rb') as f:
        bundle = pickle.load(f)

    raw_data = bundle['raw_data']
    d_cfg    = raw_data[molecule][cluster][config]

    result = {}

    for method in d_cfg.keys():                        # SCF, MP2, CCSD_T
        for basis in d_cfg[method].keys():             # aug-cc-pvdz

            def e(key, _d=d_cfg, _m=method, _b=basis):
                return total_energy(_d, _m, _b, key)

            two_b   = compute_2b(e, PAIRS)
            three_b = compute_3b(e, TRIMERS)
            four_b  = compute_4b(e, PAIRS, TRIMERS)

            result \
                .setdefault(molecule, {}) \
                .setdefault(cluster, {}) \
                .setdefault(config, {}) \
                .setdefault(method, {}) \
                [basis] = {
                    '2b': two_b,
                    '3b': three_b,
                    '4b': four_b,
                }

    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute many-body interactions for Ne or HF tetramer'
    )
    parser.add_argument(
        '--mol', required=True, choices=['Ne', 'HF'],
        help='Molecule: Ne or HF'
    )
    args = parser.parse_args()

    molecule = args.mol
    cluster  = '4'
    config   = '0'
    here     = Path(__file__).resolve().parent

    print(f'Computing {molecule} tetramer many-body interactions...')
    result = build_interactions(molecule, cluster, config)

    # --- Sample checks ---
    basis = 'aug-cc-pvdz'
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d = result[molecule][cluster][config][method][basis]
        print(f'\n{method}:')
        for pair in PAIRS:
            print(f"  2b_{pair}  no_vmfc={d['2b'][pair]['no_vmfc']:+.6e}"
                  f"  vmfc={d['2b'][pair]['vmfc']:+.6e}")
        for tri in TRIMERS:
            print(f"  3b_{tri} no_vmfc={d['3b'][tri]['no_vmfc']:+.6e}"
                  f"  vmfc={d['3b'][tri]['vmfc']:+.6e}")
        print(f"  4b       no_vmfc={d['4b']['no_vmfc']:+.6e}"
              f"  vmfc={d['4b']['vmfc']:+.6e}")

    # --- Sanity check for Ne: all 2b pairs should be equal (identical monomers) ---
    if molecule == 'Ne':
        print('\nNe symmetry check — all 2b pairs should be equal:')
        for method in ['SCF', 'MP2', 'CCSD_T']:
            d = result[molecule][cluster][config][method][basis]
            vals = [d['2b'][p]['vmfc'] for p in PAIRS]
            spread = max(vals) - min(vals)
            status = '✅' if spread < 1e-10 else f'⚠️  spread={spread:.2e}'
            print(f'  {method}: max spread across pairs = {spread:.2e}  {status}')

    # --- Pickle ---
    out_path = here / f'many_body_interactions_{molecule}_{cluster}.pickle'
    with open(out_path, 'wb') as f:
        pickle.dump(result, f)
    print(f'\nPickled to: {out_path}')