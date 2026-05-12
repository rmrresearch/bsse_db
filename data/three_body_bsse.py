"""
three_body_bsse.py  (generalized)
==================================
Trimer analogue of two_body_bsse.py. Handles H2O, Ne, and HF trimers
via --mol CLI argument.

Reads:   data/many_body_interactions_{molecule}_3.pickle
Writes:  data/three_body_bsse_{molecule}.pickle

Key differences from two_body_bsse.py
--------------------------------------
1. No `motion` level — trimer configurations are numbered (0, 1, 2, ...)
2. Two-body terms have three pairs (AB, AC, BC) instead of one scalar.
3. Three-body term is a single scalar per config/method/basis.
4. FORM(II) delta terms — 12 individual basis-sharing contributions.

Output dictionary structure
----------------------------
result[mol][cluster][config][method][basis] = {

    # FORM(I) — two-body terms per pair
    '2b': {
        'AB': {'no_vmfc': float, 'vmfc': float, 'bsse': float},
        'AC': {'no_vmfc': float, 'vmfc': float, 'bsse': float},
        'BC': {'no_vmfc': float, 'vmfc': float, 'bsse': float},
    },

    # FORM(I) — three-body term
    '3b': {'no_vmfc': float, 'vmfc': float, 'bsse': float},

    # FORM(I) — total BSSE = sum of 2b_bsse(AB,AC,BC) + 3b_bsse
    'total_bsse': float,

    # FORM(II) — 12 individual delta terms (raw positive differences)
    'delta_terms': {
        # 1-body in 2-body basis (positive in BSSE sum)
        'A_in_AB': float, 'B_in_AB': float,
        'A_in_AC': float, 'C_in_AC': float,
        'B_in_BC': float, 'C_in_BC': float,
        # 1-body in 3-body basis (negative in BSSE sum)
        'A_in_ABC': float, 'B_in_ABC': float, 'C_in_ABC': float,
        # 2-body in 3-body basis (positive in BSSE sum)
        'AB_in_ABC': float, 'AC_in_ABC': float, 'BC_in_ABC': float,
    },

    # FORM(II) — total BSSE recovered from delta terms (must equal total_bsse)
    'total_bsse_form2': float,
}

Usage
-----
    python three_body_bsse.py --mol Ne
    python three_body_bsse.py --mol HF
    python three_body_bsse.py --mol H2O
"""

import pickle
import argparse
from pathlib import Path
from open_pickle_data import open_pickle


# Signs for each delta term in the BSSE sum (FORM II)
DELTA_SIGNS = {
    'A_in_AB'  : +1, 'B_in_AB'  : +1,
    'A_in_AC'  : +1, 'C_in_AC'  : +1,
    'B_in_BC'  : +1, 'C_in_BC'  : +1,
    'A_in_ABC' : -1, 'B_in_ABC' : -1, 'C_in_ABC' : -1,
    'AB_in_ABC': +1, 'AC_in_ABC': +1, 'BC_in_ABC': +1,
}


def threebody_bsse(pickle_file, molecule):
    """
    Compute BSSE quantities for all trimer configurations.

    Parameters
    ----------
    pickle_file : str or Path
        Path to many_body_interactions_{molecule}_3.pickle
    molecule : str
        e.g. 'Ne', 'HF', 'H2O'

    Returns
    -------
    result : nested dict
        Structure described in module docstring.
    """
    raw_data = open_pickle(pickle_file)

    two_b_no   = raw_data['2b_no_vmfc']
    two_b_vf   = raw_data['2b_vmfc']
    three_b_no = raw_data['3b_no_vmfc']
    three_b_vf = raw_data['3b_vmfc']
    deltas     = raw_data['bsse_delta_terms']

    result = {molecule: {}}

    for cluster in two_b_no[molecule]:
        result[molecule][cluster] = {}

        for config in sorted(two_b_no[molecule][cluster], key=int):
            result[molecule][cluster][config] = {}

            for method in two_b_no[molecule][cluster][config]:
                result[molecule][cluster][config][method] = {}

                for basis in two_b_no[molecule][cluster][config][method]:

                    # FORM(I) two-body terms per pair
                    two_b = {}
                    for pair in ('AB', 'AC', 'BC'):
                        no_v = two_b_no[molecule][cluster][config][method][basis][pair]
                        vf   = two_b_vf [molecule][cluster][config][method][basis][pair]
                        two_b[pair] = {
                            'no_vmfc': no_v,
                            'vmfc'   : vf,
                            'bsse'   : no_v - vf,
                        }

                    # FORM(I) three-body term
                    no_v3 = three_b_no[molecule][cluster][config][method][basis]
                    vf3   = three_b_vf [molecule][cluster][config][method][basis]
                    three_b = {
                        'no_vmfc': no_v3,
                        'vmfc'   : vf3,
                        'bsse'   : no_v3 - vf3,
                    }

                    # FORM(I) total BSSE
                    total_bsse_f1 = (
                          two_b['AB']['bsse']
                        + two_b['AC']['bsse']
                        + two_b['BC']['bsse']
                        + three_b['bsse']
                    )

                    # FORM(II) delta terms
                    delta_entry = deltas[molecule][cluster][config][method][basis]

                    total_bsse_f2 = sum(
                        DELTA_SIGNS[key] * val
                        for key, val in delta_entry.items()
                    )

                    # Store everything
                    result[molecule][cluster][config][method][basis] = {
                        '2b'              : two_b,
                        '3b'              : three_b,
                        'total_bsse'      : total_bsse_f1,
                        'delta_terms'     : delta_entry,
                        'total_bsse_form2': total_bsse_f2,
                    }

    return result


def print_summary(result, molecule, params_pickle=None):
    """
    Print a compact summary table of total BSSE per config at CCSD_T level.

    Parameters
    ----------
    result       : the nested BSSE result dict
    molecule     : e.g. 'Ne', 'HF', 'H2O'
    params_pickle: path to raw_data_{molecule}_3.pickle (to get real distances).
                   If None, falls back to showing config index.
    """
    # Load distance lookup from parameters pickle if provided
    dist_map = {}
    if params_pickle is not None:
        stored = open_pickle(params_pickle)
        params = stored.get('parameters', {})
        for cfg_key, cfg_data in params.get(molecule, {}).get('3', {}).items():
            dist_map[cfg_key] = cfg_data['params']['distance']

    print()
    print(f"{'='*65}")
    print(f"  {molecule} Trimer BSSE Summary")
    print(f"{'='*65}")
    print(f"  {'Config':<8}  {'R (Å)':<10}  {'BSSE CCSD(T) (Ha)':>20}  {'BSSE (µHa)':>12}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*20}  {'-'*12}")

    cluster = '3'
    configs = sorted(result[molecule][cluster], key=int)
    for cfg in configs:
        methods = list(result[molecule][cluster][cfg].keys())
        method  = 'CCSD_T' if 'CCSD_T' in methods else methods[0]
        bases   = list(result[molecule][cluster][cfg][method].keys())
        basis   = bases[0]

        entry    = result[molecule][cluster][cfg][method][basis]
        bsse     = entry['total_bsse']
        dist_str = f"{dist_map[cfg]:.4f}" if cfg in dist_map else cfg
        print(f"  {cfg:<8}  {dist_str:<10}  {bsse:>20.10f}  {bsse*1e6:>12.4f}")

    print(f"{'='*65}")
    print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compute three-body BSSE for trimer configurations."
    )
    parser.add_argument(
        "--mol",
        type=str,
        required=True,
        choices=["H2O", "Ne", "HF"],
        help="Molecule name: H2O, Ne, or HF",
    )
    args = parser.parse_args()
    molecule = args.mol

    here        = Path(__file__).resolve().parent
    pickle_file = here / f'many_body_interactions_{molecule}_3.pickle'
    out_file    = here / f'three_body_bsse_{molecule}.pickle'

    print(f'Reading : {pickle_file}')
    result = threebody_bsse(pickle_file, molecule)

    print(f'Writing : {out_file}')
    with open(out_file, 'wb') as fh:
        pickle.dump(result, fh)
    print('Done.')

    params_pickle = here / f'raw_data_{molecule}_3.pickle'
    print_summary(result, molecule, params_pickle=params_pickle)