"""
four_body_bsse.py
------------------
Reads many_body_interactions_{molecule}_4.pickle and raw_data_{molecule}_4.pickle
and computes the BSSE decomposition for the Ne or HF tetramer.

FORM(I) — grouped by interaction order (PI's convention)
---------------------------------------------------------
    2b_bsse(IJ)  = 2b_no_vmfc(IJ)  - 2b_vmfc(IJ)    for 6 pairs
    3b_bsse(IJK) = 3b_no_vmfc(IJK) - 3b_vmfc(IJK)   for 4 triples
    4b_bsse      = 4b_no_vmfc      - 4b_vmfc
    total_bsse   = Σ 2b_bsse + Σ 3b_bsse + 4b_bsse

FORM(II) — 50 individual delta terms
---------------------------------------------------------
Each term: delta[X_in_Y] = E(X|Y) - E(X|X)
Sign rule: (-1)^(|Y| - |X| + 1)

    Group          | Count | Sign | Keys
    ---------------|-------|------|----------------------------
    1b-in-2b       |  12   |  +   | I_in_IJ  (each pair)
    1b-in-3b       |  12   |  -   | I_in_IJK (each triple)
    2b-in-3b       |  12   |  +   | IJ_in_IJK (each pair in triple)
    1b-in-4b       |   4   |  +   | I_in_ABCD
    2b-in-4b       |   6   |  -   | IJ_in_ABCD
    3b-in-4b       |   4   |  +   | IJK_in_ABCD
    Total          |  50   |

Consistency check: FORM(I) total_bsse == FORM(II) total_bsse_form2

Output structure
----------------
result[mol][cluster][config][method][basis] = {
    '2b'              : {'AB': {'no_vmfc', 'vmfc', 'bsse'}, ...},
    '3b'              : {'ABC': {'no_vmfc', 'vmfc', 'bsse'}, ...},
    '4b'              : {'no_vmfc', 'vmfc', 'bsse'},
    'total_bsse'      : float,          # FORM(I)
    'delta_terms'     : {50 keys},      # FORM(II) raw deltas (unsigned)
    'total_bsse_form2': float,          # FORM(II) signed sum
}

Molecule-specific notes
-----------------------
    Ne  : B_B = C_C = D_D = A_A was injected at the build_raw_data stage.
          No special handling needed here — all 65 keys are present.

    HF  : All four monomers are independent. No special handling needed.

Usage
-----
    python3 four_body_bsse.py --mol Ne
    python3 four_body_bsse.py --mol HF
    # produces data/four_body_bsse_{mol}_4.pickle
"""

import pickle
import argparse
from pathlib import Path
from itertools import combinations


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------

MONOMERS = ['A', 'B', 'C', 'D']
PAIRS    = [''.join(p) for p in combinations(MONOMERS, 2)]
TRIMERS  = [''.join(t) for t in combinations(MONOMERS, 3)]
TETRA    = 'ABCD'


# ---------------------------------------------------------------------------
# Energy helper (same as in many_body_interactions_4.py)
# ---------------------------------------------------------------------------

def total_energy(d, method, basis, key):
    e = d['SCF'][basis][key]
    if method in ('MP2', 'CCSD_T'):
        e += d['MP2'][basis][key]
    if method == 'CCSD_T':
        e += d['CCSD_T'][basis][key]
    return e


# ---------------------------------------------------------------------------
# FORM(II) delta terms
# ---------------------------------------------------------------------------

def compute_delta_terms(e):
    """
    Compute all 50 FORM(II) delta terms.
    Each value is the raw difference: E(X|Y) - E(X|X)  (unsigned, sign applied later)

    Parameters
    ----------
    e : callable — e(key) returns total energy for file key at current method/basis

    Returns
    -------
    dict : {term_key: float}
    """
    deltas = {}

    # --- 1b-in-2b (12 terms, sign +) ---
    for ij in PAIRS:
        for i in ij:
            deltas[f'{i}_in_{ij}'] = e(f'{i}_{ij}') - e(f'{i}_{i}')

    # --- 1b-in-3b (12 terms, sign -) ---
    for ijk in TRIMERS:
        for i in ijk:
            deltas[f'{i}_in_{ijk}'] = e(f'{i}_{ijk}') - e(f'{i}_{i}')

    # --- 2b-in-3b (12 terms, sign +) ---
    for ijk in TRIMERS:
        for r in range(len(ijk)):
            for s in range(r + 1, len(ijk)):
                ij = ijk[r] + ijk[s]
                deltas[f'{ij}_in_{ijk}'] = e(f'{ij}_{ijk}') - e(f'{ij}_{ij}')

    # --- 1b-in-4b (4 terms, sign +) ---
    for i in MONOMERS:
        deltas[f'{i}_in_{TETRA}'] = e(f'{i}_{TETRA}') - e(f'{i}_{i}')

    # --- 2b-in-4b (6 terms, sign -) ---
    for ij in PAIRS:
        deltas[f'{ij}_in_{TETRA}'] = e(f'{ij}_{TETRA}') - e(f'{ij}_{ij}')

    # --- 3b-in-4b (4 terms, sign +) ---
    for ijk in TRIMERS:
        deltas[f'{ijk}_in_{TETRA}'] = e(f'{ijk}_{TETRA}') - e(f'{ijk}_{ijk}')

    return deltas


def bsse_form2(deltas, pairs, trimers):
    """
    Apply signs to delta terms and sum to get FORM(II) total BSSE.

    Sign rule: (-1)^(|Y| - |X| + 1)
        1b-in-2b  |Y|-|X|=1  sign = +
        1b-in-3b  |Y|-|X|=2  sign = -
        2b-in-3b  |Y|-|X|=1  sign = +
        1b-in-4b  |Y|-|X|=3  sign = +
        2b-in-4b  |Y|-|X|=2  sign = -
        3b-in-4b  |Y|-|X|=1  sign = +
    """
    total = 0.0

    # 1b-in-2b (+)
    for ij in pairs:
        for i in ij:
            total += deltas[f'{i}_in_{ij}']

    # 1b-in-3b (-)
    for ijk in trimers:
        for i in ijk:
            total -= deltas[f'{i}_in_{ijk}']

    # 2b-in-3b (+)
    for ijk in trimers:
        for r in range(len(ijk)):
            for s in range(r + 1, len(ijk)):
                ij = ijk[r] + ijk[s]
                total += deltas[f'{ij}_in_{ijk}']

    # 1b-in-4b (+)
    for i in MONOMERS:
        total += deltas[f'{i}_in_{TETRA}']

    # 2b-in-4b (-)
    for ij in pairs:
        total -= deltas[f'{ij}_in_{TETRA}']

    # 3b-in-4b (+)
    for ijk in trimers:
        total += deltas[f'{ijk}_in_{TETRA}']

    return total


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_bsse(molecule, cluster='4', config='0',
               interactions_path=None, raw_data_path=None):

    here = Path(__file__).resolve().parent

    if interactions_path is None:
        interactions_path = here / f'many_body_interactions_{molecule}_{cluster}.pickle'
    if raw_data_path is None:
        raw_data_path = here / f'raw_data_{molecule}_{cluster}.pickle'

    with open(interactions_path, 'rb') as f:
        interactions = pickle.load(f)
    with open(raw_data_path, 'rb') as f:
        bundle = pickle.load(f)

    raw_data = bundle['raw_data']
    d_int    = interactions[molecule][cluster][config]
    d_raw    = raw_data[molecule][cluster][config]

    result = {}

    for method in d_int.keys():
        for basis in d_int[method].keys():

            mb = d_int[method][basis]

            def e(key, _d=d_raw, _m=method, _b=basis):
                return total_energy(_d, _m, _b, key)

            # --- FORM(I) ---
            two_b      = {}
            total_bsse = 0.0

            for ij in PAIRS:
                bsse_ij = mb['2b'][ij]['no_vmfc'] - mb['2b'][ij]['vmfc']
                two_b[ij] = {
                    'no_vmfc': mb['2b'][ij]['no_vmfc'],
                    'vmfc'   : mb['2b'][ij]['vmfc'],
                    'bsse'   : bsse_ij,
                }
                total_bsse += bsse_ij

            three_b = {}
            for ijk in TRIMERS:
                bsse_ijk = mb['3b'][ijk]['no_vmfc'] - mb['3b'][ijk]['vmfc']
                three_b[ijk] = {
                    'no_vmfc': mb['3b'][ijk]['no_vmfc'],
                    'vmfc'   : mb['3b'][ijk]['vmfc'],
                    'bsse'   : bsse_ijk,
                }
                total_bsse += bsse_ijk

            bsse_4b = mb['4b']['no_vmfc'] - mb['4b']['vmfc']
            four_b  = {
                'no_vmfc': mb['4b']['no_vmfc'],
                'vmfc'   : mb['4b']['vmfc'],
                'bsse'   : bsse_4b,
            }
            total_bsse += bsse_4b

            # --- FORM(II) ---
            deltas      = compute_delta_terms(e)
            total_form2 = bsse_form2(deltas, PAIRS, TRIMERS)

            # --- Consistency check ---
            diff = abs(total_bsse - total_form2)
            status = '✅' if diff < 1e-10 else f'⚠️  MISMATCH'
            if diff > 1e-10:
                print(f'WARNING: FORM(I) vs FORM(II) mismatch '
                      f'({method}/{basis}): |diff| = {diff:.2e} Ha')

            result \
                .setdefault(molecule, {}) \
                .setdefault(cluster, {}) \
                .setdefault(config, {}) \
                .setdefault(method, {}) \
                [basis] = {
                    '2b'              : two_b,
                    '3b'              : three_b,
                    '4b'              : four_b,
                    'total_bsse'      : total_bsse,
                    'delta_terms'     : deltas,
                    'total_bsse_form2': total_form2,
                }

    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute tetramer BSSE decomposition for Ne or HF'
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

    print(f'Computing {molecule} tetramer BSSE decomposition...\n')
    result = build_bsse(molecule, cluster, config)

    basis = 'aug-cc-pvdz'
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d = result[molecule][cluster][config][method][basis]
        print(f'--- {method} ---')

        # FORM(I) — per-pair, per-triple, per-4b
        for ij in PAIRS:
            print(f"  2b_bsse({ij})  = {d['2b'][ij]['bsse']:+.6e} Ha")
        for ijk in TRIMERS:
            print(f"  3b_bsse({ijk}) = {d['3b'][ijk]['bsse']:+.6e} Ha")
        print(f"  4b_bsse       = {d['4b']['bsse']:+.6e} Ha")

        # FORM(I) total vs FORM(II) total
        print(f"  total_bsse    = {d['total_bsse']:+.6e} Ha  [FORM I]")
        print(f"  total_form2   = {d['total_bsse_form2']:+.6e} Ha  [FORM II]")
        diff = abs(d['total_bsse'] - d['total_bsse_form2'])
        print(f"  |FORM I - FORM II| = {diff:.2e} Ha  "
              f"{'✅' if diff < 1e-10 else '⚠️  MISMATCH'}\n")

    # --- Body-order scaling summary ---
    print('--- Body-order scaling summary ---')
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d = result[molecule][cluster][config][method][basis]
        sum_2b = sum(d['2b'][ij]['bsse'] for ij in PAIRS)
        sum_3b = sum(d['3b'][ijk]['bsse'] for ijk in TRIMERS)
        bsse_4b = d['4b']['bsse']
        total   = d['total_bsse']
        print(f'  {method}:  Σ2b={sum_2b:+.4e}  Σ3b={sum_3b:+.4e}'
              f'  4b={bsse_4b:+.4e}  total={total:+.4e} Ha')

    # --- Pickle ---
    out_path = here / f'four_body_bsse_{molecule}_{cluster}.pickle'
    with open(out_path, 'wb') as f:
        pickle.dump(result, f)
    print(f'\nPickled to: {out_path}')