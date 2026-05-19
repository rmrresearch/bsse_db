"""
four_body_bsse_10.py
---------------------
Reads many_body_interactions_HF_10.pickle and raw_data_HF_10.pickle
and computes the BSSE decomposition for the HF decamer (truncated at 4b).

FORM(I) -- grouped by body order
---------------------------------
    2b_bsse(IJ)   = 2b_no_vmfc(IJ)   - 2b_vmfc(IJ)    for 45 pairs
    3b_bsse(IJK)  = 3b_no_vmfc(IJK)  - 3b_vmfc(IJK)   for 120 triples
    4b_bsse(IJKL) = 4b_no_vmfc(IJKL) - 4b_vmfc(IJKL)  for 210 quadruples

    a_0 = sum of 2b_bsse  (45 terms)
    a_1 = sum of 3b_bsse  (120 terms)
    a_2 = sum of 4b_bsse  (210 terms)
    total_bsse = a_0 + a_1 + a_2

FORM(II) -- 3750 individual delta terms
-----------------------------------------
Each term: delta[X_in_Y] = E(X|Y) - E(X|X)
Sign rule: (-1)^(|Y| - |X| + 1)

    Group       | Count | Sign
    ------------|-------|------
    1b-in-2b    |    90 |  +
    1b-in-3b    |   360 |  -
    2b-in-3b    |   360 |  +
    1b-in-4b    |   840 |  +
    2b-in-4b    |  1260 |  -
    3b-in-4b    |   840 |  +
    Total       |  3750 |

Consistency check: FORM(I) total_bsse == FORM(II) total_bsse_form2

Output structure
----------------
result['HF']['10']['0'][method]['aug-cc-pvdz'] = {
    '2b'              : {'AB': {'no_vmfc', 'vmfc', 'bsse'}, ...},  # 45
    '3b'              : {'ABC': {'no_vmfc', 'vmfc', 'bsse'}, ...}, # 120
    '4b'              : {'ABCD': {'no_vmfc', 'vmfc', 'bsse'}, ...},# 210
    'a_0'             : float,   # sum of 2b_bsse over 45 pairs
    'a_1'             : float,   # sum of 3b_bsse over 120 triples
    'a_2'             : float,   # sum of 4b_bsse over 210 quadruples
    'total_bsse'      : float,   # FORM(I): a_0 + a_1 + a_2
    'delta_terms'     : {3750 keys},
    'total_bsse_form2': float,   # FORM(II) signed sum
}

Usage
-----
    python3 four_body_bsse_10.py
    # produces  data/four_body_bsse_HF_10.pickle
"""

import pickle
from pathlib import Path
from itertools import combinations


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------

MONOMERS   = list('ABCDEFGHIJ')
PAIRS      = [''.join(p) for p in combinations(MONOMERS, 2)]   # 45
TRIMERS    = [''.join(t) for t in combinations(MONOMERS, 3)]   # 120
QUADRUPLES = [''.join(q) for q in combinations(MONOMERS, 4)]   # 210

MOLECULE = 'HF'
CLUSTER  = '10'
CONFIG   = '0'
BASIS    = 'aug-cc-pvdz'

EXPECTED_DELTA_TERMS = 3750


# ---------------------------------------------------------------------------
# Energy helper  (same as many_body_interactions_10.py)
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
    Compute all 3750 FORM(II) delta terms.
    Each value: E(X|Y) - E(X|X)  (unsigned; sign applied in bsse_form2)

    Key format: '{real}_in_{basis}'  e.g. 'A_in_AB', 'AB_in_ABCD'

    Returns
    -------
    dict : {term_key: float}  -- 3750 entries
    """
    deltas = {}

    # --- 1b-in-2b (90 terms) ---
    for ij in PAIRS:
        for i in ij:
            deltas[f'{i}_in_{ij}'] = e(f'{i}_{ij}') - e(f'{i}_{i}')

    # --- 1b-in-3b (360 terms) ---
    for ijk in TRIMERS:
        for i in ijk:
            deltas[f'{i}_in_{ijk}'] = e(f'{i}_{ijk}') - e(f'{i}_{i}')

    # --- 2b-in-3b (360 terms) ---
    for ijk in TRIMERS:
        for r in range(len(ijk)):
            for s in range(r + 1, len(ijk)):
                ij = ijk[r] + ijk[s]
                deltas[f'{ij}_in_{ijk}'] = e(f'{ij}_{ijk}') - e(f'{ij}_{ij}')

    # --- 1b-in-4b (840 terms) ---
    for ijkl in QUADRUPLES:
        for i in ijkl:
            deltas[f'{i}_in_{ijkl}'] = e(f'{i}_{ijkl}') - e(f'{i}_{i}')

    # --- 2b-in-4b (1260 terms) ---
    for ijkl in QUADRUPLES:
        for r in range(len(ijkl)):
            for s in range(r + 1, len(ijkl)):
                ij = ijkl[r] + ijkl[s]
                deltas[f'{ij}_in_{ijkl}'] = e(f'{ij}_{ijkl}') - e(f'{ij}_{ij}')

    # --- 3b-in-4b (840 terms) ---
    for ijkl in QUADRUPLES:
        for r in range(len(ijkl)):
            for s in range(r + 1, len(ijkl)):
                for t in range(s + 1, len(ijkl)):
                    ijk = ijkl[r] + ijkl[s] + ijkl[t]
                    deltas[f'{ijk}_in_{ijkl}'] = (
                        e(f'{ijk}_{ijkl}') - e(f'{ijk}_{ijk}')
                    )

    return deltas


def bsse_form2(deltas):
    """
    Apply signs to delta terms and sum to get FORM(II) total BSSE.

    Sign rule: (-1)^(|Y| - |X| + 1)
        1b-in-2b  (+) : |Y|-|X| = 1
        1b-in-3b  (-) : |Y|-|X| = 2
        2b-in-3b  (+) : |Y|-|X| = 1
        1b-in-4b  (+) : |Y|-|X| = 3
        2b-in-4b  (-) : |Y|-|X| = 2
        3b-in-4b  (+) : |Y|-|X| = 1
    """
    total = 0.0

    # 1b-in-2b (+)
    for ij in PAIRS:
        for i in ij:
            total += deltas[f'{i}_in_{ij}']

    # 1b-in-3b (-)
    for ijk in TRIMERS:
        for i in ijk:
            total -= deltas[f'{i}_in_{ijk}']

    # 2b-in-3b (+)
    for ijk in TRIMERS:
        for r in range(len(ijk)):
            for s in range(r + 1, len(ijk)):
                ij = ijk[r] + ijk[s]
                total += deltas[f'{ij}_in_{ijk}']

    # 1b-in-4b (+)
    for ijkl in QUADRUPLES:
        for i in ijkl:
            total += deltas[f'{i}_in_{ijkl}']

    # 2b-in-4b (-)
    for ijkl in QUADRUPLES:
        for r in range(len(ijkl)):
            for s in range(r + 1, len(ijkl)):
                ij = ijkl[r] + ijkl[s]
                total -= deltas[f'{ij}_in_{ijkl}']

    # 3b-in-4b (+)
    for ijkl in QUADRUPLES:
        for r in range(len(ijkl)):
            for s in range(r + 1, len(ijkl)):
                for t in range(s + 1, len(ijkl)):
                    ijk = ijkl[r] + ijkl[s] + ijkl[t]
                    total += deltas[f'{ijk}_in_{ijkl}']

    return total


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_bsse(molecule=MOLECULE, cluster=CLUSTER, config=CONFIG,
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
            two_b = {}
            for ij in PAIRS:
                bsse_ij = mb['2b'][ij]['no_vmfc'] - mb['2b'][ij]['vmfc']
                two_b[ij] = {
                    'no_vmfc': mb['2b'][ij]['no_vmfc'],
                    'vmfc'   : mb['2b'][ij]['vmfc'],
                    'bsse'   : bsse_ij,
                }

            three_b = {}
            for ijk in TRIMERS:
                bsse_ijk = mb['3b'][ijk]['no_vmfc'] - mb['3b'][ijk]['vmfc']
                three_b[ijk] = {
                    'no_vmfc': mb['3b'][ijk]['no_vmfc'],
                    'vmfc'   : mb['3b'][ijk]['vmfc'],
                    'bsse'   : bsse_ijk,
                }

            four_b = {}
            for ijkl in QUADRUPLES:
                bsse_ijkl = mb['4b'][ijkl]['no_vmfc'] - mb['4b'][ijkl]['vmfc']
                four_b[ijkl] = {
                    'no_vmfc': mb['4b'][ijkl]['no_vmfc'],
                    'vmfc'   : mb['4b'][ijkl]['vmfc'],
                    'bsse'   : bsse_ijkl,
                }

            a_0        = sum(v['bsse'] for v in two_b.values())
            a_1        = sum(v['bsse'] for v in three_b.values())
            a_2        = sum(v['bsse'] for v in four_b.values())
            total_bsse = a_0 + a_1 + a_2

            # --- FORM(II) ---
            deltas      = compute_delta_terms(e)
            total_form2 = bsse_form2(deltas)

            # --- Consistency check ---
            diff = abs(total_bsse - total_form2)
            if diff > 1e-10:
                print(f'WARNING: FORM(I) vs FORM(II) mismatch '
                      f'({method}/{basis}): |diff| = {diff:.2e} Ha')

            (result
             .setdefault(molecule, {})
             .setdefault(cluster, {})
             .setdefault(config, {})
             .setdefault(method, {}))[basis] = {
                '2b'              : two_b,
                '3b'              : three_b,
                '4b'              : four_b,
                'a_0'             : a_0,
                'a_1'             : a_1,
                'a_2'             : a_2,
                'total_bsse'      : total_bsse,
                'delta_terms'     : deltas,
                'total_bsse_form2': total_form2,
            }

    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    here = Path(__file__).resolve().parent

    print('Computing HF decamer BSSE decomposition...\n')
    result = build_bsse()

    basis = BASIS

    # --- Delta term count ---
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d      = result[MOLECULE][CLUSTER][CONFIG][method][basis]
        n_deltas = len(d['delta_terms'])
        ok = '✅' if n_deltas == EXPECTED_DELTA_TERMS else f'⚠️  expected {EXPECTED_DELTA_TERMS}'
        print(f'{method:8s}  delta terms: {n_deltas}  {ok}')

    print()

    # --- FORM(I) vs FORM(II) consistency ---
    print('FORM(I) vs FORM(II) consistency:')
    all_pass = True
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d    = result[MOLECULE][CLUSTER][CONFIG][method][basis]
        diff = abs(d['total_bsse'] - d['total_bsse_form2'])
        ok   = diff < 1e-10
        flag = '✅' if ok else '⚠️  MISMATCH'
        print(f'  {method:8s}  |FORM I - FORM II| = {diff:.2e} Ha  {flag}')
        if not ok:
            all_pass = False

    print()

    # --- Body-order decomposition ---
    print('Body-order decomposition (Ha):')
    KCAL = 627.509474
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d = result[MOLECULE][CLUSTER][CONFIG][method][basis]
        print(f'  {method}:')
        print(f'    a_0 (Σ 2b, 45 pairs)    = {d["a_0"]:+.6e} Ha'
              f'  ({d["a_0"]*KCAL:+.4f} kcal/mol)')
        print(f'    a_1 (Σ 3b, 120 triples) = {d["a_1"]:+.6e} Ha'
              f'  ({d["a_1"]*KCAL:+.4f} kcal/mol)')
        print(f'    a_2 (Σ 4b, 210 quads)   = {d["a_2"]:+.6e} Ha'
              f'  ({d["a_2"]*KCAL:+.4f} kcal/mol)')
        print(f'    total                   = {d["total_bsse"]:+.6e} Ha'
              f'  ({d["total_bsse"]*KCAL:+.4f} kcal/mol)')
        print()

    # --- Pickle ---
    out_path = here / f'four_body_bsse_{MOLECULE}_{CLUSTER}.pickle'
    with open(out_path, 'wb') as f:
        pickle.dump(result, f)
    print(f'Pickled to: {out_path}')