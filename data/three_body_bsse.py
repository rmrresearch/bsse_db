"""
three_body_bsse.py
==================
Trimer analogue of two_body_bsse.py.

Reads many_body_interactions_H2O_3.pickle and computes BSSE quantities
for both FORM(I) and FORM(II).

Lives in data/ alongside two_body_bsse.py, as per the PI's convention.

Key differences from two_body_bsse.py
--------------------------------------
1. No `motion` level — trimer configurations are just numbered (0, 1, 2, ...)
   Dimer nesting : [mol][cluster][motion][config][method][basis]
   Trimer nesting: [mol][cluster][config][method][basis]

2. Two-body terms have three pairs (AB, AC, BC) instead of one scalar.
   Each pair stores: no_vmfc, vmfc, bsse

3. Three-body term is a single scalar per config/method/basis.
   Stores: no_vmfc, vmfc, bsse

4. FORM(II) delta terms — 12 individual basis-sharing contributions
   organized in 3 physical groups:
       1b_in_2b  : A_in_AB, B_in_AB, A_in_AC, C_in_AC, B_in_BC, C_in_BC
       1b_in_3b  : A_in_ABC, B_in_ABC, C_in_ABC
       2b_in_3b  : AB_in_ABC, AC_in_ABC, BC_in_ABC

   Signs in total BSSE:
       1b_in_2b  → positive  (ghost functions inflate interaction)
       1b_in_3b  → negative  (correction for overcounting)
       2b_in_3b  → positive  (dimer gains trimer ghost functions)

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
        # 1-body in 2-body basis (positive contribution to BSSE)
        'A_in_AB': float, 'B_in_AB': float,
        'A_in_AC': float, 'C_in_AC': float,
        'B_in_BC': float, 'C_in_BC': float,
        # 1-body in 3-body basis (negative contribution to BSSE)
        'A_in_ABC': float, 'B_in_ABC': float, 'C_in_ABC': float,
        # 2-body in 3-body basis (positive contribution to BSSE)
        'AB_in_ABC': float, 'AC_in_ABC': float, 'BC_in_ABC': float,
    },

    # FORM(II) — total BSSE recovered from delta terms (must equal total_bsse)
    'total_bsse_form2': float,
}

Usage
-----
From data/ or any script that adds data/ to sys.path:

    from three_body_bsse import threebody_bsse
    result = threebody_bsse(pickle_file, 'H2O')

    # Access example
    entry = result['H2O']['3']['0']['CCSD_T']['aug-cc-pvdz']
    print(entry['total_bsse'])
    print(entry['delta_terms']['AB_in_ABC'])
"""

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
        Path to many_body_interactions_H2O_3.pickle
    molecule : str
        e.g. 'H2O'

    Returns
    -------
    result : nested dict
        Structure described in module docstring.
    """
    raw_data = open_pickle(pickle_file)

    # Pull the four FORM(I) quantities and the 12 delta terms
    two_b_no   = raw_data['2b_no_vmfc']
    two_b_vf   = raw_data['2b_vmfc']
    three_b_no = raw_data['3b_no_vmfc']
    three_b_vf = raw_data['3b_vmfc']
    deltas     = raw_data['bsse_delta_terms']

    three_body_bsse = {molecule: {}}

    for cluster in two_b_no[molecule]:
        three_body_bsse[molecule][cluster] = {}

        # No motion loop — trimer has no rot/trnl distinction
        for config in two_b_no[molecule][cluster]:
            three_body_bsse[molecule][cluster][config] = {}

            for method in two_b_no[molecule][cluster][config]:
                three_body_bsse[molecule][cluster][config][method] = {}

                for basis in two_b_no[molecule][cluster][config][method]:

                    # ---- FORM(I) two-body terms per pair ----------------
                    two_b = {}
                    for pair in ('AB', 'AC', 'BC'):
                        no_v = two_b_no[molecule][cluster][config][method][basis][pair]
                        vf   = two_b_vf [molecule][cluster][config][method][basis][pair]
                        two_b[pair] = {
                            'no_vmfc': no_v,
                            'vmfc'   : vf,
                            'bsse'   : no_v - vf,
                        }

                    # ---- FORM(I) three-body term ------------------------
                    no_v3 = three_b_no[molecule][cluster][config][method][basis]
                    vf3   = three_b_vf [molecule][cluster][config][method][basis]
                    three_b = {
                        'no_vmfc': no_v3,
                        'vmfc'   : vf3,
                        'bsse'   : no_v3 - vf3,
                    }

                    # ---- FORM(I) total BSSE -----------------------------
                    total_bsse_f1 = (
                        two_b['AB']['bsse']
                      + two_b['AC']['bsse']
                      + two_b['BC']['bsse']
                      + three_b['bsse']
                    )

                    # ---- FORM(II) delta terms ---------------------------
                    delta_entry = deltas[molecule][cluster][config][method][basis]

                    # Recover total BSSE from signed sum of delta terms
                    total_bsse_f2 = sum(
                        DELTA_SIGNS[key] * val
                        for key, val in delta_entry.items()
                    )

                    # ---- Store everything -------------------------------
                    three_body_bsse[molecule][cluster][config][method][basis] = {
                        '2b'              : two_b,
                        '3b'              : three_b,
                        'total_bsse'      : total_bsse_f1,
                        'delta_terms'     : delta_entry,
                        'total_bsse_form2': total_bsse_f2,
                    }

    return three_body_bsse