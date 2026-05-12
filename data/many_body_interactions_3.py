"""
many_body_interactions_3.py  (generalized)
==========================================
Step 4 of Issue #46 — Compute many-body interactions for trimers.
Handles H2O, Ne, and HF trimers via --mol CLI argument.

Reads:   data/raw_data_{molecule}_3.pickle
Writes:  data/many_body_interactions_{molecule}_3.pickle

=============================================================================
WHAT THIS SCRIPT COMPUTES
=============================================================================

PART A — PI's Step 4 quantities (FORM I building blocks)
---------------------------------------------------------
  2b_no_vmfc[pair] : Two-body interaction energy without CP correction
      E^(2)_AB      = EAB(AB) − EA(A) − EB(B)

  2b_vmfc[pair]    : Two-body interaction energy WITH CP correction (VMFC)
      E^(2)_AB(AB)  = EAB(AB) − EA(AB) − EB(AB)

  3b_no_vmfc       : Three-body interaction energy without CP correction
      E^(3)_ABC     = EABC(ABC) − EA(A) − EB(B) − EC(C)
                    − E^(2)_AB − E^(2)_AC − E^(2)_BC

  3b_vmfc          : Three-body interaction energy WITH CP correction (VMFC)
      E^(3)_ABC(ABC)= EABC(ABC) − EA(ABC) − EB(ABC) − EC(ABC)
                    − E^(2)_AB(ABC) − E^(2)_AC(ABC) − E^(2)_BC(ABC)

PART B — FORM(II) extension (12 individual δ terms)
----------------------------------------------------
Each δ[X,Y] = E(fragment X in basis Y) − E(fragment X in own basis).
Signs in the total BSSE sum:
    POSITIVE : δ[monomer, dimer basis]   — 6 terms  (1b_in_2b)
    NEGATIVE : δ[monomer, trimer basis]  — 3 terms  (1b_in_3b)
    POSITIVE : δ[dimer,   trimer basis]  — 3 terms  (2b_in_3b)

CONSISTENCY CHECK
-----------------
At the end of every config/method/basis triple we verify:
    BSSE_FORM1 == BSSE_FORM2  to within 1e-10 Hartree.

NOTES ON Ne AND HF TRIMERS
---------------------------
All 19 file keys (A_A, B_B, C_C, AB_AB, ..., BC_ABC) are present as
separate .out files for both Ne and HF trimers. No aliasing is required.
The formulas are identical to H2O — full generality confirmed.

=============================================================================
ENERGY CONVENTIONS (inherited from build_raw_data_3.py)
=============================================================================
raw_data stores COMPONENTS, not total energies:
    raw_data[mol][cluster][config]['SCF'][basis][key]    = E_scf (total SCF)
    raw_data[mol][cluster][config]['MP2'][basis][key]    = E_mp2_corr  (= E_mp2 − E_scf)
    raw_data[mol][cluster][config]['CCSD_T'][basis][key] = E_ccsdt_comp(= E_ccsdt − E_mp2)

=============================================================================
OUTPUT — many_body_interactions_{molecule}_3.pickle
=============================================================================
A dict with keys:
    'molecule'        : str, e.g. 'Ne'
    'cluster'         : str, '3'
    '2b_no_vmfc'      : nested dict [mol][cluster][config][method][basis][pair]
    '2b_vmfc'         : same structure
    '3b_no_vmfc'      : nested dict [mol][cluster][config][method][basis] → scalar
    '3b_vmfc'         : same structure
    'bsse_delta_terms': nested dict [mol][cluster][config][method][basis][term_key]

Usage
-----
    python many_body_interactions_3.py --mol Ne
    python many_body_interactions_3.py --mol HF
    python many_body_interactions_3.py --mol H2O
"""

import pickle
import argparse
from pathlib import Path


# =============================================================================
# SECTION 1 — Energy reconstruction helper
# =============================================================================

def get_total_energy(raw_data, molecule, cluster, config, method, basis, file_key):
    """
    Reconstruct the total energy for a given fragment/basis combination.

    build_raw_data_3.py stores CORRELATION COMPONENTS, so we add them up:
        SCF    → E_scf
        MP2    → E_scf + mp2_corr
        CCSD_T → E_scf + mp2_corr + ccsdt_component
    """
    cfg_data = raw_data[molecule][cluster][config]

    e_scf = cfg_data['SCF'][basis][file_key]
    if method == 'SCF':
        return e_scf

    e_mp2_total = e_scf + cfg_data['MP2'][basis][file_key]
    if method == 'MP2':
        return e_mp2_total

    # method == 'CCSD_T'
    return e_mp2_total + cfg_data['CCSD_T'][basis][file_key]


# =============================================================================
# SECTION 2 — FORM(I) quantities
# =============================================================================

def compute_2b_no_vmfc(e):
    """
    Two-body interaction energies WITHOUT counterpoise correction.
    E^(2)_IJ = E_IJ(IJ) − E_I(I) − E_J(J)
    Returns dict with keys 'AB', 'AC', 'BC'.
    """
    return {
        'AB': e('AB_AB') - e('A_A') - e('B_B'),
        'AC': e('AC_AC') - e('A_A') - e('C_C'),
        'BC': e('BC_BC') - e('B_B') - e('C_C'),
    }


def compute_2b_vmfc(e):
    """
    Two-body interaction energies WITH CP correction (VMFC).
    E^(2)_IJ(IJ) = E_IJ(IJ) − E_I(IJ) − E_J(IJ)
    Returns dict with keys 'AB', 'AC', 'BC'.
    """
    return {
        'AB': e('AB_AB') - e('A_AB') - e('B_AB'),
        'AC': e('AC_AC') - e('A_AC') - e('C_AC'),
        'BC': e('BC_BC') - e('B_BC') - e('C_BC'),
    }


def compute_3b_no_vmfc(e, two_b_no):
    """
    Three-body interaction energy WITHOUT counterpoise correction.
    E^(3)_ABC = EABC(ABC) − EA(A) − EB(B) − EC(C) − E^(2)_AB − E^(2)_AC − E^(2)_BC
    Returns scalar.
    """
    return (
        e('ABC_ABC')
        - e('A_A') - e('B_B') - e('C_C')
        - two_b_no['AB'] - two_b_no['AC'] - two_b_no['BC']
    )


def compute_3b_vmfc(e, two_b_vmfc):
    """
    Three-body interaction energy WITH CP correction (VMFC).
    All fragments computed in the full trimer basis.
    Returns scalar.
    """
    e2_AB_ABC = e('AB_ABC') - e('A_ABC') - e('B_ABC')
    e2_AC_ABC = e('AC_ABC') - e('A_ABC') - e('C_ABC')
    e2_BC_ABC = e('BC_ABC') - e('B_ABC') - e('C_ABC')

    return (
        e('ABC_ABC')
        - e('A_ABC') - e('B_ABC') - e('C_ABC')
        - e2_AB_ABC - e2_AC_ABC - e2_BC_ABC
    )


# =============================================================================
# SECTION 3 — FORM(II) delta terms
# =============================================================================

def compute_bsse_delta_terms(e):
    """
    Compute the 12 individual basis-sharing (delta) terms for FORM(II).
    Each δ[X,Y] = E(fragment X in basis Y) − E(fragment X in own basis).
    Values stored as raw positive differences; signs applied at BSSE sum time.
    """
    return {
        # 1-body in 2-body basis (positive in BSSE sum)
        'A_in_AB' : e('A_AB')  - e('A_A'),
        'B_in_AB' : e('B_AB')  - e('B_B'),
        'A_in_AC' : e('A_AC')  - e('A_A'),
        'C_in_AC' : e('C_AC')  - e('C_C'),
        'B_in_BC' : e('B_BC')  - e('B_B'),
        'C_in_BC' : e('C_BC')  - e('C_C'),
        # 1-body in 3-body basis (negative in BSSE sum)
        'A_in_ABC': e('A_ABC') - e('A_A'),
        'B_in_ABC': e('B_ABC') - e('B_B'),
        'C_in_ABC': e('C_ABC') - e('C_C'),
        # 2-body in 3-body basis (positive in BSSE sum)
        'AB_in_ABC': e('AB_ABC') - e('AB_AB'),
        'AC_in_ABC': e('AC_ABC') - e('AC_AC'),
        'BC_in_ABC': e('BC_ABC') - e('BC_BC'),
    }


def bsse_from_delta_terms(delta):
    """
    Recover total BSSE from the 12 delta terms using FORM(II) signs.
    Used in the consistency check.
    """
    return (
        + delta['A_in_AB']   + delta['B_in_AB']
        + delta['A_in_AC']   + delta['C_in_AC']
        + delta['B_in_BC']   + delta['C_in_BC']
        - delta['A_in_ABC']  - delta['B_in_ABC']  - delta['C_in_ABC']
        + delta['AB_in_ABC'] + delta['AC_in_ABC'] + delta['BC_in_ABC']
    )


# =============================================================================
# SECTION 4 — Consistency check
# =============================================================================

def check_consistency(two_b_no, two_b_vf, three_b_no, three_b_vf,
                      delta, config, method, basis, tol=1e-10):
    """
    Verify FORM(I) and FORM(II) give the same total BSSE.
    Raises ValueError if they disagree beyond tol (default 1e-10 Hartree).
    Returns the BSSE value (verified consistent).
    """
    bsse_form1 = (
          (two_b_no['AB'] - two_b_vf['AB'])
        + (two_b_no['AC'] - two_b_vf['AC'])
        + (two_b_no['BC'] - two_b_vf['BC'])
        + (three_b_no - three_b_vf)
    )
    bsse_form2 = bsse_from_delta_terms(delta)

    diff = abs(bsse_form1 - bsse_form2)
    if diff > tol:
        raise ValueError(
            f"CONSISTENCY FAIL — config={config}, method={method}, basis={basis}\n"
            f"  BSSE FORM(I)  = {bsse_form1:.12f}\n"
            f"  BSSE FORM(II) = {bsse_form2:.12f}\n"
            f"  |difference|  = {diff:.2e}  (tolerance = {tol:.2e})"
        )
    return bsse_form1


# =============================================================================
# SECTION 5 — Main driver
# =============================================================================

def compute_many_body(molecule, cluster='3', data_dir=None):
    """
    Main function. Loops over all configs, methods, and bases in the pickle
    and computes FORM(I) and FORM(II) quantities for every combination.

    Parameters
    ----------
    molecule : str   e.g. 'Ne', 'HF', 'H2O'
    cluster  : str   '3' (trimer — this script is trimer-specific)
    data_dir : Path  directory containing input pickle and where output is
                     written. Defaults to the directory where this script lives.

    Returns
    -------
    bundle : dict  the full output dictionary (also pickled to disk)
    """
    here = Path(__file__).resolve().parent
    if data_dir is None:
        data_dir = here

    in_pickle  = data_dir / f'raw_data_{molecule}_{cluster}.pickle'
    out_pickle = data_dir / f'many_body_interactions_{molecule}_{cluster}.pickle'

    print(f'Reading : {in_pickle}')
    with open(in_pickle, 'rb') as fh:
        stored = pickle.load(fh)
    raw_data = stored['raw_data']

    # Dynamically discover methods and bases from first config
    mol_data  = raw_data[molecule][cluster]
    first_cfg = next(iter(mol_data))
    methods   = list(mol_data[first_cfg].keys())
    bases     = list(mol_data[first_cfg][methods[0]].keys())
    configs   = list(mol_data.keys())

    print(f'Molecule : {molecule}   Cluster : {cluster}')
    print(f'Configs  : {sorted(configs, key=int)}')
    print(f'Methods  : {methods}')
    print(f'Bases    : {bases}')
    print()

    # Initialise output containers
    two_b_no  = {molecule: {cluster: {}}}
    two_b_vf  = {molecule: {cluster: {}}}
    three_b_no= {molecule: {cluster: {}}}
    three_b_vf= {molecule: {cluster: {}}}
    delta_out = {molecule: {cluster: {}}}

    # Main loop
    for cfg in sorted(configs, key=int):
        two_b_no  [molecule][cluster][cfg] = {}
        two_b_vf  [molecule][cluster][cfg] = {}
        three_b_no[molecule][cluster][cfg] = {}
        three_b_vf[molecule][cluster][cfg] = {}
        delta_out [molecule][cluster][cfg] = {}

        for method in methods:
            two_b_no  [molecule][cluster][cfg][method] = {}
            two_b_vf  [molecule][cluster][cfg][method] = {}
            three_b_no[molecule][cluster][cfg][method] = {}
            three_b_vf[molecule][cluster][cfg][method] = {}
            delta_out [molecule][cluster][cfg][method] = {}

            for basis in bases:

                def e(file_key, _cfg=cfg, _method=method, _basis=basis):
                    return get_total_energy(
                        raw_data, molecule, cluster, _cfg, _method, _basis, file_key
                    )

                # FORM(I)
                two_b_no_result   = compute_2b_no_vmfc(e)
                two_b_vf_result   = compute_2b_vmfc(e)
                three_b_no_result = compute_3b_no_vmfc(e, two_b_no_result)
                three_b_vf_result = compute_3b_vmfc(e, two_b_vf_result)

                # FORM(II)
                delta_result = compute_bsse_delta_terms(e)

                # Consistency check
                bsse_val = check_consistency(
                    two_b_no_result, two_b_vf_result,
                    three_b_no_result, three_b_vf_result,
                    delta_result, cfg, method, basis
                )

                # Store results
                two_b_no  [molecule][cluster][cfg][method][basis] = two_b_no_result
                two_b_vf  [molecule][cluster][cfg][method][basis] = two_b_vf_result
                three_b_no[molecule][cluster][cfg][method][basis] = three_b_no_result
                three_b_vf[molecule][cluster][cfg][method][basis] = three_b_vf_result
                delta_out [molecule][cluster][cfg][method][basis] = delta_result

                print(f'  ✅ cfg={cfg}  method={method}  basis={basis}'
                      f'  BSSE={bsse_val:.8f} Ha')

    print()
    print(f'All {len(configs)} configurations processed. FORM(I) == FORM(II) for all entries.')

    # Assemble and pickle output bundle
    bundle = {
        'molecule'        : molecule,
        'cluster'         : cluster,
        '2b_no_vmfc'      : two_b_no,
        '2b_vmfc'         : two_b_vf,
        '3b_no_vmfc'      : three_b_no,
        '3b_vmfc'         : three_b_vf,
        'bsse_delta_terms': delta_out,
    }

    print(f'Writing : {out_pickle}')
    with open(out_pickle, 'wb') as fh:
        pickle.dump(bundle, fh)
    print('Done.')

    return bundle


# =============================================================================
# SECTION 6 — Entry point
# =============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compute many-body interactions for trimers."
    )
    parser.add_argument(
        "--mol",
        type=str,
        required=True,
        choices=["H2O", "Ne", "HF"],
        help="Molecule name: H2O, Ne, or HF",
    )
    args = parser.parse_args()

    bundle = compute_many_body(molecule=args.mol, cluster='3')

    # Quick sanity peek
    mol, cl = bundle['molecule'], bundle['cluster']
    first_cfg    = next(iter(bundle['2b_no_vmfc'][mol][cl]))
    first_method = next(iter(bundle['2b_no_vmfc'][mol][cl][first_cfg]))
    first_basis  = next(iter(bundle['2b_no_vmfc'][mol][cl][first_cfg][first_method]))

    print()
    print('--- Quick peek at first entry ---')
    print(f'Config : {first_cfg}   Method : {first_method}   Basis : {first_basis}')
    print()

    two_b_no = bundle['2b_no_vmfc'][mol][cl][first_cfg][first_method][first_basis]
    two_b_vf = bundle['2b_vmfc']   [mol][cl][first_cfg][first_method][first_basis]
    three_no = bundle['3b_no_vmfc'][mol][cl][first_cfg][first_method][first_basis]
    three_vf = bundle['3b_vmfc']   [mol][cl][first_cfg][first_method][first_basis]
    deltas   = bundle['bsse_delta_terms'][mol][cl][first_cfg][first_method][first_basis]

    print('FORM(I) quantities:')
    for pair in ['AB', 'AC', 'BC']:
        print(f'  2b_no_vmfc({pair}) = {two_b_no[pair]:.10f} Ha')
        print(f'  2b_vmfc({pair})    = {two_b_vf[pair]:.10f} Ha')
        print(f'  2b_bsse({pair})    = {two_b_no[pair]-two_b_vf[pair]:.10f} Ha')
    print(f'  3b_no_vmfc        = {three_no:.10f} Ha')
    print(f'  3b_vmfc           = {three_vf:.10f} Ha')
    print(f'  3b_bsse           = {three_no-three_vf:.10f} Ha')
    total_bsse_f1 = (
          (two_b_no['AB'] - two_b_vf['AB'])
        + (two_b_no['AC'] - two_b_vf['AC'])
        + (two_b_no['BC'] - two_b_vf['BC'])
        + (three_no - three_vf)
    )
    print(f'  Total BSSE FORM(I)= {total_bsse_f1:.10f} Ha')

    print()
    print('FORM(II) delta terms:')
    signs = {
        'A_in_AB': +1, 'B_in_AB': +1,
        'A_in_AC': +1, 'C_in_AC': +1,
        'B_in_BC': +1, 'C_in_BC': +1,
        'A_in_ABC': -1, 'B_in_ABC': -1, 'C_in_ABC': -1,
        'AB_in_ABC': +1, 'AC_in_ABC': +1, 'BC_in_ABC': +1,
    }
    total_bsse_f2 = 0.0
    for key, val in deltas.items():
        sign = signs[key]
        contrib = sign * val
        total_bsse_f2 += contrib
        print(f'  {key:12s}  raw={val:+.10f}  sign={sign:+d}  contribution={contrib:+.10f} Ha')
    print(f'  Total BSSE FORM(II)= {total_bsse_f2:.10f} Ha')
    print(f'  |FORM(I) − FORM(II)| = {abs(total_bsse_f1 - total_bsse_f2):.2e} Ha')