"""
many_body_interactions_3.py
============================
Step 4 (extended) of Issue #46 — Compute many-body interactions for trimers.

Naming convention follows the codebase (_3 = cluster size 3, not molecule-specific):
    build_raw_data_3.py          → produces raw_data_H2O_3.pickle
    tar_output_discovery_3.py    → streams .out files from 3.tar.gz
    many_body_interactions_3.py  → produces many_body_interactions_3.pickle (THIS FILE)

Currently wired for H2O/3. Will serve HF/3 and Ne/3 by changing the `molecule`
argument once those pickles exist.

Reads:   data/raw_data_{molecule}_3.pickle
Writes:  data/many_body_interactions_3.pickle

=============================================================================
WHAT THIS SCRIPT COMPUTES
=============================================================================

PART A — PI's Step 4 quantities (FORM I building blocks)
---------------------------------------------------------
These are the interaction energies needed to compute BSSE in FORM(I).
FORM(I) groups the BSSE into interaction-order terms:

    BSSE = 2b_bsse(AB) + 2b_bsse(AC) + 2b_bsse(BC) + 3b_bsse(ABC)

To get there, we first compute the interaction energies WITHOUT and WITH
counterpoise correction:

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

      where E^(2)_AB(ABC) = EAB(ABC) − EA(ABC) − EB(ABC)  [dimer in trimer basis]

PART B — Your FORM(II) extension (12 individual δ terms)
---------------------------------------------------------
FORM(II) fully expands the BSSE into 12 individual basis-sharing contributions.
Each δ[X,Y] = E(fragment X computed in Y basis) − E(fragment X in own basis).

The 12 terms with their signs in the BSSE sum are:

  POSITIVE (monomer gains ghost functions from dimer basis):
    δ[A,AB]  = EA(AB)  − EA(A)     file keys: A_AB  − A_A
    δ[B,AB]  = EB(AB)  − EB(B)     file keys: B_AB  − B_B
    δ[A,AC]  = EA(AC)  − EA(A)     file keys: A_AC  − A_A
    δ[C,AC]  = EC(AC)  − EC(C)     file keys: C_AC  − C_C
    δ[B,BC]  = EB(BC)  − EB(B)     file keys: B_BC  − B_B
    δ[C,BC]  = EC(BC)  − EC(C)     file keys: C_BC  − C_C

  NEGATIVE (correction: monomer already has trimer basis — over-counted above):
    δ[A,ABC] = EA(ABC) − EA(A)     file keys: A_ABC − A_A
    δ[B,ABC] = EB(ABC) − EB(B)     file keys: B_ABC − B_B
    δ[C,ABC] = EC(ABC) − EC(C)     file keys: C_ABC − C_C

  POSITIVE (dimer gains ghost functions from trimer basis):
    δ[AB,ABC]= EAB(ABC)− EAB(AB)   file keys: AB_ABC − AB_AB
    δ[AC,ABC]= EAC(ABC)− EAC(AC)   file keys: AC_ABC − AC_AC
    δ[BC,ABC]= EBC(ABC)− EBC(BC)   file keys: BC_ABC − BC_BC

Full FORM(II) BSSE expression (Image 3 in your notes):
    BSSE = +δ[A,AB]  + δ[B,AB]
           +δ[A,AC]  + δ[C,AC]
           +δ[B,BC]  + δ[C,BC]
           −δ[A,ABC] − δ[B,ABC] − δ[C,ABC]
           +δ[AB,ABC]+ δ[AC,ABC]+ δ[BC,ABC]

The δ values are stored as RAW POSITIVE DIFFERENCES (always larger_basis − own_basis).
The signs above are applied when summing to recover total BSSE.
This way each δ value on its own tells you the magnitude of that basis-sharing effect,
which is what you need for the threshold/approximation analysis in Step 10.

CONSISTENCY CHECK
-----------------
At the end of the script we verify:
    BSSE_FORM1 = (2b_no_vmfc − 2b_vmfc) summed over pairs + (3b_no_vmfc − 3b_vmfc)
    BSSE_FORM2 = signed sum of the 12 δ terms
    |BSSE_FORM1 − BSSE_FORM2| < 1e-10   for every config/method/basis

=============================================================================
ENERGY CONVENTIONS (inherited from build_raw_data_3.py)
=============================================================================
raw_data stores COMPONENTS, not total energies:
    raw_data[mol][cluster][config]['SCF'][basis][key]    = E_scf (total SCF)
    raw_data[mol][cluster][config]['MP2'][basis][key]    = E_mp2_corr  (= E_mp2 − E_scf)
    raw_data[mol][cluster][config]['CCSD_T'][basis][key] = E_ccsdt_comp(= E_ccsdt − E_mp2)

The helper function get_total_energy() reconstructs total energies from these.

=============================================================================
OUTPUT — many_body_interactions_3.pickle
=============================================================================
A dict with keys:
    'molecule'      : str, e.g. 'H2O'
    'cluster'       : str, e.g. '3'
    '2b_no_vmfc'    : nested dict [mol][cluster][config][method][basis][pair]
    '2b_vmfc'       : same structure
    '3b_no_vmfc'    : nested dict [mol][cluster][config][method][basis] → scalar
    '3b_vmfc'       : same structure
    'bsse_delta_terms': nested dict [mol][cluster][config][method][basis][term_key]

                       The 12 term_keys are organized in 3 physical groups:

                       1-body in 2-body basis (monomer gains dimer ghost functions):
                           'A_in_AB'   'B_in_AB'          ← from dimer AB
                           'A_in_AC'   'C_in_AC'          ← from dimer AC
                           'B_in_BC'   'C_in_BC'          ← from dimer BC

                       1-body in 3-body basis (monomer gains trimer ghost functions):
                           'A_in_ABC'  'B_in_ABC'  'C_in_ABC'

                       2-body in 3-body basis (dimer gains trimer ghost functions):
                           'AB_in_ABC' 'AC_in_ABC' 'BC_in_ABC'

                       Physical meaning of each group:
                           1b_in_2b  → cost of BSSE at the two-body level       (+)
                           1b_in_3b  → correction for overcounting in 1b_in_2b  (-)
                           2b_in_3b  → additional BSSE from three-body basis     (+)
"""

import pickle
from pathlib import Path


# =============================================================================
# SECTION 1 — Energy reconstruction helper
# =============================================================================

def get_total_energy(raw_data, molecule, cluster, config, method, basis, file_key):
    """
    Reconstruct the total energy for a given fragment/basis combination.

    Because build_raw_data_3.py stores CORRELATION COMPONENTS (not totals),
    we need to add them up to recover the total energy at a given level of theory.

    Parameters
    ----------
    raw_data  : the nested dictionary from raw_data_H2O_3.pickle
    molecule  : e.g. 'H2O'
    cluster   : e.g. '3'
    config    : e.g. '0'
    method    : one of 'SCF', 'MP2', 'CCSD_T'
    basis     : e.g. 'aug-cc-pvdz'
    file_key  : e.g. 'A_AB', 'AB_ABC', 'ABC_ABC'

    Returns
    -------
    float : total energy in Hartree

    Examples
    --------
    get_total_energy(..., 'SCF',   ..., 'A_A')   → E_scf of monomer A in basis A
    get_total_energy(..., 'MP2',   ..., 'A_AB')  → E_mp2_total of monomer A in basis AB
    get_total_energy(..., 'CCSD_T',..., 'AB_AB') → E_ccsdt_total of dimer AB in basis AB
    """
    cfg_data = raw_data[molecule][cluster][config]

    e_scf = cfg_data['SCF'][basis][file_key]
    if method == 'SCF':
        return e_scf

    e_mp2_total = e_scf + cfg_data['MP2'][basis][file_key]
    if method == 'MP2':
        return e_mp2_total

    # method == 'CCSD_T'
    e_ccsdt_total = e_mp2_total + cfg_data['CCSD_T'][basis][file_key]
    return e_ccsdt_total


# =============================================================================
# SECTION 2 — FORM(I) quantities
# =============================================================================

def compute_2b_no_vmfc(e):
    """
    Two-body interaction energies WITHOUT counterpoise correction.

    Formula:  E^(2)_IJ = E_IJ(IJ) − E_I(I) − E_J(J)
    This is the standard (uncorrected) two-body interaction energy.
    Each monomer is computed in its OWN basis, no ghost functions.

    Returns a dict with keys 'AB', 'AC', 'BC'.
    `e` is a shorthand function: e(file_key) → total energy for that file.
    """
    return {
        'AB': e('AB_AB') - e('A_A') - e('B_B'),   # EAB(AB) − EA(A) − EB(B)
        'AC': e('AC_AC') - e('A_A') - e('C_C'),   # EAC(AC) − EA(A) − EC(C)
        'BC': e('BC_BC') - e('B_B') - e('C_C'),   # EBC(BC) − EB(B) − EC(C)
    }


def compute_2b_vmfc(e):
    """
    Two-body interaction energies WITH counterpoise correction (VMFC).

    Formula:  E^(2)_IJ(IJ) = E_IJ(IJ) − E_I(IJ) − E_J(IJ)
    Each monomer is computed in the DIMER basis (ghost functions included).
    This removes the BSSE at the two-body level.

    Returns a dict with keys 'AB', 'AC', 'BC'.
    """
    return {
        'AB': e('AB_AB') - e('A_AB') - e('B_AB'),  # EAB(AB) − EA(AB) − EB(AB)
        'AC': e('AC_AC') - e('A_AC') - e('C_AC'),  # EAC(AC) − EA(AC) − EC(AC)
        'BC': e('BC_BC') - e('B_BC') - e('C_BC'),  # EBC(BC) − EB(BC) − EC(BC)
    }


def compute_3b_no_vmfc(e, two_b_no):
    """
    Three-body interaction energy WITHOUT counterpoise correction.

    Formula:
        E^(3)_ABC = EABC(ABC) − EA(A) − EB(B) − EC(C)
                  − E^(2)_AB − E^(2)_AC − E^(2)_BC

    In words: take the full trimer energy, subtract all monomer energies
    (in their own basis), then subtract all two-body interactions (also
    without CP). What remains is the pure three-body contribution.

    `two_b_no` is the dict already computed by compute_2b_no_vmfc().
    Returns a scalar.
    """
    return (
        e('ABC_ABC')
        - e('A_A') - e('B_B') - e('C_C')
        - two_b_no['AB'] - two_b_no['AC'] - two_b_no['BC']
    )


def compute_3b_vmfc(e, two_b_vmfc):
    """
    Three-body interaction energy WITH counterpoise correction (VMFC).

    Formula:
        E^(3)_ABC(ABC) = EABC(ABC) − EA(ABC) − EB(ABC) − EC(ABC)
                       − E^(2)_AB(ABC) − E^(2)_AC(ABC) − E^(2)_BC(ABC)

    where the dimer-in-trimer-basis interactions are:
        E^(2)_AB(ABC) = EAB(ABC) − EA(ABC) − EB(ABC)
        E^(2)_AC(ABC) = EAC(ABC) − EA(ABC) − EC(ABC)
        E^(2)_BC(ABC) = EBC(ABC) − EB(ABC) − EC(ABC)

    All fragments are computed in the FULL TRIMER BASIS (ghost functions from
    all three monomers included). This is the fully CP-corrected three-body term.

    Returns a scalar.
    """
    # Dimer interactions evaluated in the trimer basis
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

    Each δ[X,Y] = E(fragment X in basis Y) − E(fragment X in its own basis).
    This measures how much a fragment's energy changes when ghost functions
    from a larger system are added. This is the raw cost of BSSE for that
    specific fragment-basis combination.

    The δ values are stored as RAW POSITIVE DIFFERENCES. Their signs in the
    total BSSE expression are:
        POSITIVE : δ[monomer, dimer basis]   — 6 terms
        NEGATIVE : δ[monomer, trimer basis]  — 3 terms
        POSITIVE : δ[dimer,   trimer basis]  — 3 terms

    Knowing the sign of each contribution tells you whether that term
    inflates or deflates the BSSE, and by how much.

    Returns a dict with 12 keys as described below.
    """
    return {
        # --- Monomer in dimer basis (6 terms) ---
        # These are POSITIVE in BSSE: ghost functions from the dimer lower
        # the monomer energy, making the non-CP interaction look too favorable.
        'A_in_AB' : e('A_AB')  - e('A_A'),   # EA(AB)  − EA(A)
        'B_in_AB' : e('B_AB')  - e('B_B'),   # EB(AB)  − EB(B)
        'A_in_AC' : e('A_AC')  - e('A_A'),   # EA(AC)  − EA(A)
        'C_in_AC' : e('C_AC')  - e('C_C'),   # EC(AC)  − EC(C)
        'B_in_BC' : e('B_BC')  - e('B_B'),   # EB(BC)  − EB(B)
        'C_in_BC' : e('C_BC')  - e('C_C'),   # EC(BC)  − EC(C)

        # --- Monomer in trimer basis (3 terms) ---
        # These are NEGATIVE in BSSE: they correct for the fact that in the
        # dimer-basis terms above, we didn't yet account for the full trimer
        # basis. Subtracting these avoids double-counting.
        'A_in_ABC': e('A_ABC') - e('A_A'),   # EA(ABC) − EA(A)
        'B_in_ABC': e('B_ABC') - e('B_B'),   # EB(ABC) − EB(B)
        'C_in_ABC': e('C_ABC') - e('C_C'),   # EC(ABC) − EC(C)

        # --- Dimer in trimer basis (3 terms) ---
        # These are POSITIVE in BSSE: ghost functions from the third monomer
        # lower the dimer energy beyond what the dimer-only calculation sees.
        'AB_in_ABC': e('AB_ABC') - e('AB_AB'),  # EAB(ABC) − EAB(AB)
        'AC_in_ABC': e('AC_ABC') - e('AC_AC'),  # EAC(ABC) − EAC(AC)
        'BC_in_ABC': e('BC_ABC') - e('BC_BC'),  # EBC(ABC) − EBC(BC)
    }


def bsse_from_delta_terms(delta):
    """
    Recover total BSSE from the 12 delta terms using FORM(II) signs.

    BSSE = +δ[A,AB]  + δ[B,AB]  + δ[A,AC]  + δ[C,AC]  + δ[B,BC]  + δ[C,BC]
           −δ[A,ABC] − δ[B,ABC] − δ[C,ABC]
           +δ[AB,ABC]+ δ[AC,ABC]+ δ[BC,ABC]

    Used in the consistency check.
    """
    return (
        + delta['A_in_AB']  + delta['B_in_AB']
        + delta['A_in_AC']  + delta['C_in_AC']
        + delta['B_in_BC']  + delta['C_in_BC']
        - delta['A_in_ABC'] - delta['B_in_ABC'] - delta['C_in_ABC']
        + delta['AB_in_ABC']+ delta['AC_in_ABC']+ delta['BC_in_ABC']
    )


# =============================================================================
# SECTION 4 — Consistency check
# =============================================================================

def check_consistency(two_b_no, two_b_vf, three_b_no, three_b_vf,
                      delta, config, method, basis, tol=1e-10):
    """
    Verify that FORM(I) and FORM(II) give the same total BSSE.

    BSSE_FORM1 = (2b_no_vmfc - 2b_vmfc) summed over AB, AC, BC
               + (3b_no_vmfc - 3b_vmfc)

    BSSE_FORM2 = signed sum of the 12 delta terms

    Raises ValueError if they disagree beyond `tol` (default 1e-10 Hartree).
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
    return bsse_form1  # return the value since we've verified both agree


# =============================================================================
# SECTION 5 — Main driver
# =============================================================================

def compute_many_body(molecule='H2O', cluster='3', data_dir=None):
    """
    Main function. Loops over all configs, methods, and bases found in the
    pickle and computes FORM(I) and FORM(II) quantities for every combination.

    Parameters
    ----------
    molecule : str   e.g. 'H2O'  (will extend to 'HF', 'Ne' later)
    cluster  : str   '3'          (trimer — this script is trimer-specific)
    data_dir : Path  directory containing the input pickle and where the
                     output pickle will be written. Defaults to the directory
                     where this script lives (i.e. bsse_db/data/).

    Returns
    -------
    bundle : dict  the full output dictionary (also pickled to disk)
    """
    # ---- resolve paths -------------------------------------------------------
    here = Path(__file__).resolve().parent
    if data_dir is None:
        data_dir = here

    in_pickle  = data_dir / f'raw_data_{molecule}_{cluster}.pickle'
    out_pickle = data_dir / f'many_body_interactions_{molecule}_{cluster}.pickle'

    print(f'Reading : {in_pickle}')
    with open(in_pickle, 'rb') as fh:
        stored = pickle.load(fh)
    raw_data = stored['raw_data']

    # ---- discover methods and bases from first config ------------------------
    # We don't hard-code 'aug-cc-pvdz' or ['SCF','MP2','CCSD_T'] — we discover
    # them dynamically so the script works if new bases or methods are added.
    mol_data  = raw_data[molecule][cluster]
    first_cfg = next(iter(mol_data))
    methods   = list(mol_data[first_cfg].keys())
    bases     = list(mol_data[first_cfg][methods[0]].keys())
    configs   = list(mol_data.keys())

    print(f'Molecule : {molecule}   Cluster : {cluster}')
    print(f'Configs  : {configs}')
    print(f'Methods  : {methods}')
    print(f'Bases    : {bases}')
    print()

    # ---- initialise output containers ----------------------------------------
    # Each container mirrors the nesting: [mol][cluster][config][method][basis]
    two_b_no  = {molecule: {cluster: {}}}
    two_b_vf  = {molecule: {cluster: {}}}
    three_b_no= {molecule: {cluster: {}}}
    three_b_vf= {molecule: {cluster: {}}}
    delta_out = {molecule: {cluster: {}}}

    # ---- main loop -----------------------------------------------------------
    for cfg in configs:
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

                # Shorthand: e(file_key) returns total energy for this
                # (cfg, method, basis, file_key) combination.
                # This keeps all the formulas below clean and readable.
                def e(file_key, _cfg=cfg, _method=method, _basis=basis):
                    return get_total_energy(
                        raw_data, molecule, cluster, _cfg, _method, _basis, file_key
                    )

                # -- FORM(I) ---------------------------------------------------
                two_b_no_result = compute_2b_no_vmfc(e)
                two_b_vf_result = compute_2b_vmfc(e)
                three_b_no_result = compute_3b_no_vmfc(e, two_b_no_result)
                three_b_vf_result = compute_3b_vmfc(e, two_b_vf_result)

                # -- FORM(II) --------------------------------------------------
                delta_result = compute_bsse_delta_terms(e)

                # -- Consistency check -----------------------------------------
                bsse_val = check_consistency(
                    two_b_no_result, two_b_vf_result,
                    three_b_no_result, three_b_vf_result,
                    delta_result, cfg, method, basis
                )

                # -- Store results ---------------------------------------------
                two_b_no  [molecule][cluster][cfg][method][basis] = two_b_no_result
                two_b_vf  [molecule][cluster][cfg][method][basis] = two_b_vf_result
                three_b_no[molecule][cluster][cfg][method][basis] = three_b_no_result
                three_b_vf[molecule][cluster][cfg][method][basis] = three_b_vf_result
                delta_out [molecule][cluster][cfg][method][basis] = delta_result

                print(f'  ✅ cfg={cfg}  method={method}  basis={basis}'
                      f'  BSSE={bsse_val:.8f} Ha')

    print()
    print('All configurations processed. FORM(I) == FORM(II) for all entries.')

    # ---- assemble and pickle the output bundle -------------------------------
    bundle = {
        'molecule' : molecule,
        'cluster'  : cluster,

        # FORM(I) — PI's step 4 quantities
        # Use these for step 5: 2b_bsse = 2b_no_vmfc − 2b_vmfc (per pair)
        #                        3b_bsse = 3b_no_vmfc − 3b_vmfc
        '2b_no_vmfc' : two_b_no,
        '2b_vmfc'    : two_b_vf,
        '3b_no_vmfc' : three_b_no,
        '3b_vmfc'    : three_b_vf,

        # FORM(II) — your extension: 12 individual delta terms
        # Each value is a RAW POSITIVE DIFFERENCE (larger_basis − own_basis).
        # Apply signs from bsse_from_delta_terms() to recover total BSSE.
        # Use these for Step 10 threshold/approximation analysis.
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
    # Default: H2O trimer. Change molecule='HF' or molecule='Ne' when ready.
    bundle = compute_many_body(molecule='H2O', cluster='3')

    # Quick sanity peek at the output
    mol, cl = bundle['molecule'], bundle['cluster']
    first_cfg = next(iter(bundle['2b_no_vmfc'][mol][cl]))
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
        (two_b_no['AB']-two_b_vf['AB'])
      + (two_b_no['AC']-two_b_vf['AC'])
      + (two_b_no['BC']-two_b_vf['BC'])
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
    print(f'  |FORM(I) − FORM(II)| = {abs(total_bsse_f1-total_bsse_f2):.2e} Ha')