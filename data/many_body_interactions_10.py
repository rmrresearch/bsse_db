"""
many_body_interactions_10.py
-----------------------------
Reads raw_data_HF_10.pickle and computes the many-body interaction
energies for the HF decamer using the VMFC (Valiron-Mayer) scheme.
Truncated at 4-body: only 2b, 3b, and 4b terms are computed.

Interaction terms computed
--------------------------
For each method (SCF, MP2, CCSD_T) and basis (aug-cc-pvdz):

  2b  -- C(10,2) = 45  dimer pairs
  3b  -- C(10,3) = 120 trimer triples
  4b  -- C(10,4) = 210 quadruple subsystems

Each stored as {'no_vmfc': float, 'vmfc': float}

Output dictionary structure
---------------------------
result['HF']['10']['0'][method]['aug-cc-pvdz'] = {
    '2b': {
        'AB':   {'no_vmfc': float, 'vmfc': float},
        'AC':   {'no_vmfc': float, 'vmfc': float},
        ...                                          # 45 entries
    },
    '3b': {
        'ABC':  {'no_vmfc': float, 'vmfc': float},
        ...                                          # 120 entries
    },
    '4b': {
        'ABCD': {'no_vmfc': float, 'vmfc': float},
        ...                                          # 210 entries
    },
}

NOTE on '4b' structure
----------------------
Unlike the tetramer (where '4b' is a single dict), here '4b' is keyed
by each of the 210 quadruple labels -- consistent with how '2b' and '3b'
are already keyed by subsystem label.

Formulas (Valiron & Mayer 1997, PI task list Issue #46)
-------------------------------------------------------
2b_no_vmfc(IJ)   = E(IJ|IJ) - E(I|I) - E(J|J)
2b_vmfc(IJ)      = E(IJ|IJ) - E(I|IJ) - E(J|IJ)

3b_no_vmfc(IJK)  = E(IJK|IJK) - E(IJ|IJ) - E(IK|IK) - E(JK|JK)
                 + E(I|I) + E(J|J) + E(K|K)
3b_vmfc(IJK)     = E(IJK|IJK) - E(I|IJK) - E(J|IJK) - E(K|IJK)
                 - eps2(IJ|IJK) - eps2(IK|IJK) - eps2(JK|IJK)
  where eps2(IJ|IJK) = E(IJ|IJK) - E(I|IJK) - E(J|IJK)

4b_no_vmfc(IJKL) = E(IJKL|IJKL)
                 - [E(IJK|IJK) + E(IJL|IJL) + E(IKL|IKL) + E(JKL|JKL)]
                 + [E(IJ|IJ) + E(IK|IK) + E(IL|IL)
                    + E(JK|JK) + E(JL|JL) + E(KL|KL)]
                 - [E(I|I) + E(J|J) + E(K|K) + E(L|L)]

4b_vmfc(IJKL)    = E(IJKL|IJKL)
                 - [E(I|IJKL) + E(J|IJKL) + E(K|IJKL) + E(L|IJKL)]
                 - sum of eps2(XY|IJKL)  for all 6 pairs XY in IJKL
                 - sum of eps3(XYZ|IJKL) for all 4 triples XYZ in IJKL
  where eps2(XY|IJKL)  = E(XY|IJKL) - E(X|IJKL) - E(Y|IJKL)
        eps3(XYZ|IJKL) = E(XYZ|IJKL) - E(X|IJKL) - E(Y|IJKL) - E(Z|IJKL)
                       - eps2(XY|IJKL) - eps2(XZ|IJKL) - eps2(YZ|IJKL)

Notation: E(X|Y) = fragment X computed in basis Y = file key X_Y.

Usage
-----
    python3 many_body_interactions_10.py
    # produces  data/many_body_interactions_HF_10.pickle
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


# ---------------------------------------------------------------------------
# Energy helper  (same as many_body_interactions_4.py)
# ---------------------------------------------------------------------------

def total_energy(d, method, basis, key):
    """
    Reconstruct the total energy at the requested method level.

    raw_data stores incremental components:
        SCF    -> total SCF energy
        MP2    -> MP2 correlation (E_mp2 - E_scf)
        CCSD_T -> CCSD(T) component (E_ccsdt - E_mp2)

    Parameters
    ----------
    d      : dict  -- raw_data[mol][cluster][config]
    method : str   -- 'SCF', 'MP2', or 'CCSD_T'
    basis  : str
    key    : str   -- file key, e.g. 'A_ABCD'
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
    """Compute 2b no_vmfc and vmfc for all 45 pairs."""
    result = {}
    for ij in pairs:
        i, j = ij[0], ij[1]
        no_vmfc = e(f'{ij}_{ij}') - e(f'{i}_{i}') - e(f'{j}_{j}')
        vmfc    = e(f'{ij}_{ij}') - e(f'{i}_{ij}') - e(f'{j}_{ij}')
        result[ij] = {'no_vmfc': no_vmfc, 'vmfc': vmfc}
    return result


def compute_3b(e, trimers):
    """Compute 3b no_vmfc and vmfc for all 120 triples."""
    result = {}
    for ijk in trimers:
        i, j, k = ijk[0], ijk[1], ijk[2]
        ij, ik, jk = i+j, i+k, j+k

        no_vmfc = (
            e(f'{ijk}_{ijk}')
            - e(f'{ij}_{ij}') - e(f'{ik}_{ik}') - e(f'{jk}_{jk}')
            + e(f'{i}_{i}')   + e(f'{j}_{j}')   + e(f'{k}_{k}')
        )

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


def compute_4b_for_quad(e, quad):
    """
    Compute 4b no_vmfc and vmfc for a single quadruple subsystem.

    All sub-cluster energies for vmfc are evaluated in the quadruple
    basis (IJKL), not in the full decamer basis -- consistent with the
    truncated VMFC scheme.

    Parameters
    ----------
    e    : callable  -- e(key) returns total energy for file key
    quad : str       -- 4-letter label, e.g. 'ABCD'

    Returns
    -------
    {'no_vmfc': float, 'vmfc': float}
    """
    mons  = list(quad)                                        # 4 monomers
    pairs = [''.join(p) for p in combinations(mons, 2)]      # 6 pairs
    tris  = [''.join(t) for t in combinations(mons, 3)]      # 4 triples

    # no_vmfc -- inclusion-exclusion over own-basis energies
    no_vmfc = e(f'{quad}_{quad}')
    for ijk in tris:
        no_vmfc -= e(f'{ijk}_{ijk}')
    for ij in pairs:
        no_vmfc += e(f'{ij}_{ij}')
    for i in mons:
        no_vmfc -= e(f'{i}_{i}')

    # vmfc -- all sub-clusters evaluated in the quadruple basis
    eps2 = {}
    for ij in pairs:
        i, j     = ij[0], ij[1]
        eps2[ij] = e(f'{ij}_{quad}') - e(f'{i}_{quad}') - e(f'{j}_{quad}')

    eps3 = {}
    for ijk in tris:
        i, j, k = ijk[0], ijk[1], ijk[2]
        ij, ik, jk = i+j, i+k, j+k
        eps3[ijk] = (
            e(f'{ijk}_{quad}')
            - e(f'{i}_{quad}') - e(f'{j}_{quad}') - e(f'{k}_{quad}')
            - eps2[ij] - eps2[ik] - eps2[jk]
        )

    vmfc = (
        e(f'{quad}_{quad}')
        - sum(e(f'{i}_{quad}') for i in mons)
        - sum(eps2.values())
        - sum(eps3.values())
    )

    return {'no_vmfc': no_vmfc, 'vmfc': vmfc}


def compute_4b(e, quadruples):
    """Compute 4b no_vmfc and vmfc for all 210 quadruples."""
    return {quad: compute_4b_for_quad(e, quad) for quad in quadruples}


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_interactions(molecule=MOLECULE, cluster=CLUSTER, config=CONFIG,
                       raw_data_path=None):
    here = Path(__file__).resolve().parent

    if raw_data_path is None:
        raw_data_path = here / f'raw_data_{molecule}_{cluster}.pickle'

    with open(raw_data_path, 'rb') as f:
        bundle = pickle.load(f)

    raw_data = bundle['raw_data']
    d_cfg    = raw_data[molecule][cluster][config]

    result = {}

    for method in d_cfg.keys():                    # SCF, MP2, CCSD_T
        for basis in d_cfg[method].keys():         # aug-cc-pvdz

            def e(key, _d=d_cfg, _m=method, _b=basis):
                return total_energy(_d, _m, _b, key)

            two_b   = compute_2b(e, PAIRS)
            three_b = compute_3b(e, TRIMERS)
            four_b  = compute_4b(e, QUADRUPLES)

            (result
             .setdefault(molecule, {})
             .setdefault(cluster, {})
             .setdefault(config, {})
             .setdefault(method, {}))[basis] = {
                '2b': two_b,
                '3b': three_b,
                '4b': four_b,
            }

    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    here = Path(__file__).resolve().parent

    print('Computing HF decamer many-body interactions...')
    print(f'  Pairs: {len(PAIRS)}  |  Triples: {len(TRIMERS)}'
          f'  |  Quadruples: {len(QUADRUPLES)}\n')

    result = build_interactions()

    basis = BASIS

    # --- Entry count verification ---
    print('Entry counts:')
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d = result[MOLECULE][CLUSTER][CONFIG][method][basis]
        n2 = len(d['2b'])
        n3 = len(d['3b'])
        n4 = len(d['4b'])
        ok2 = '✅' if n2 == 45  else f'⚠️  expected 45'
        ok3 = '✅' if n3 == 120 else f'⚠️  expected 120'
        ok4 = '✅' if n4 == 210 else f'⚠️  expected 210'
        print(f'  {method:8s}  2b:{n2} {ok2}  3b:{n3} {ok3}  4b:{n4} {ok4}')

    # --- Sample 2b values (first 5 pairs) ---
    print('\nSample 2b interactions (CCSD_T, Ha):')
    d = result[MOLECULE][CLUSTER][CONFIG]['CCSD_T'][basis]
    for pair in PAIRS[:5]:
        v = d['2b'][pair]
        print(f"  {pair}  no_vmfc={v['no_vmfc']:+.6e}  "
              f"vmfc={v['vmfc']:+.6e}  "
              f"bsse={v['no_vmfc']-v['vmfc']:+.6e}")

    # --- Sample 3b values (first 5 triples) ---
    print('\nSample 3b interactions (CCSD_T, Ha):')
    for tri in TRIMERS[:5]:
        v = d['3b'][tri]
        print(f"  {tri}  no_vmfc={v['no_vmfc']:+.6e}  "
              f"vmfc={v['vmfc']:+.6e}  "
              f"bsse={v['no_vmfc']-v['vmfc']:+.6e}")

    # --- Sample 4b values (first 3 quadruples) ---
    print('\nSample 4b interactions (CCSD_T, Ha):')
    for quad in QUADRUPLES[:3]:
        v = d['4b'][quad]
        print(f"  {quad}  no_vmfc={v['no_vmfc']:+.6e}  "
              f"vmfc={v['vmfc']:+.6e}  "
              f"bsse={v['no_vmfc']-v['vmfc']:+.6e}")

    # --- Body-order sums (all methods) ---
    print('\nBody-order sums (Ha):')
    for method in ['SCF', 'MP2', 'CCSD_T']:
        d = result[MOLECULE][CLUSTER][CONFIG][method][basis]
        a0 = sum(v['no_vmfc'] - v['vmfc'] for v in d['2b'].values())
        a1 = sum(v['no_vmfc'] - v['vmfc'] for v in d['3b'].values())
        a2 = sum(v['no_vmfc'] - v['vmfc'] for v in d['4b'].values())
        total = a0 + a1 + a2
        print(f'  {method:8s}  a0={a0:+.6e}  a1={a1:+.6e}  '
              f'a2={a2:+.6e}  total={total:+.6e}')

    # --- Pickle ---
    out_path = here / f'many_body_interactions_{MOLECULE}_{CLUSTER}.pickle'
    with open(out_path, 'wb') as f:
        pickle.dump(result, f)
    print(f'\nPickled to: {out_path}')