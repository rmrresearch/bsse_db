import pickle
from pathlib import Path
from open_pickle_data import open_pickle

"""
Construction of two-body interaction energy dictionaries for BSSE analysis.

This module builds two-body (dimer) interaction energies from raw quantum
chemistry output that has already been parsed and serialized into a pickle
file. Two interaction definitions are constructed:

1. Standard (no-VMFC) interaction energies:
   E_int = E_AB(AB) - E_A(A) - E_B(B)

2. Counterpoise / VMFC interaction energies:
   E_int^CP = E_AB(AB) - E_A(AB) - E_B(AB)

The resulting data structures are organized by:
    molecule → cluster size → motion → configuration → method → basis

and serialized into a single pickle file suitable for subsequent BSSE
analysis (Step 5).

This design is intentionally modular and supports extension to higher
many-body terms (3-body, 4-body, etc.) and multiple molecules.
"""


def _dimmer_no_vmfc(pickle_file,molecule):
    """
    Construct two-body (dimer) standard interaction energies (no VMFC).

    This function computes the standard dimer interaction energy for each
    configuration using monomer energies evaluated in their own basis sets:

        E_int = E_AB(AB) - E_A(A) - E_B(B)

    For homodimers (e.g., H2O–H2O), E_B(B) is not explicitly present in the
    raw data and is assumed equal to E_A(A).

    Parameters
    ----------
    pickle_file : str or pathlib.Path
        Path to the pickle file containing parsed raw total energies.
        The pickle is expected to have the structure:
            raw_data[molecule][cluster][motion][config][method][basis][component]

    molecule : str
        Molecular system label (e.g., "H2O").

    Returns
    -------
    dict
        Nested dictionary containing standard two-body interaction energies:
            result[molecule][cluster][motion][config][method][basis]["energy"]

    Notes
    -----
    - Metadata keys (e.g., '_params', '_config_label') are ignored.
    - Only configurations that converged and were included in the raw pickle
      are processed.
    - This function does not perform any file I/O.
    """

    raw_data = open_pickle(pickle_file)['raw_data']
    two_body_no_vmfc={}
    two_body_no_vmfc[molecule] = {}
    for cluster in raw_data[molecule]:
        two_body_no_vmfc[molecule][cluster] = {}
        #--
        for motion in raw_data[molecule][cluster]:
            two_body_no_vmfc[molecule][cluster][motion] = {}
            #-- 
            for config in raw_data[molecule][cluster][motion]:
                two_body_no_vmfc[molecule][cluster][motion][config] = {}
                #-- 
                for method in raw_data[molecule][cluster][motion][config]:
                    two_body_no_vmfc[molecule][cluster][motion][config][method] = {}
                    #--
                    for basis in raw_data[molecule][cluster][motion][config][method]:
                        two_body_no_vmfc[molecule][cluster][motion][config][method][basis] = {}
                        #--
                        E_A_A  = raw_data[molecule][cluster][motion][config][method][basis]['A_A']
                        E_AB_AB = raw_data[molecule][cluster][motion][config][method][basis]['AB_AB']
                        if 'B_B' not in raw_data[molecule][cluster][motion][config][method][basis]:
                            E_B_B = E_A_A
                        else:
                            E_B_B = raw_data[molecule][cluster][motion][config][method][basis]['B_B'] 
                        energy = E_AB_AB -(E_A_A + E_B_B)
                        #--
                        two_body_no_vmfc[molecule][cluster][motion][config][method][basis]['energy']= energy
    
    return two_body_no_vmfc

def _dimmer_vmfc(pickle_file,molecule):
    """
    Construct two-body (dimer) counterpoise / VMFC interaction energies.

    This function computes the counterpoise-corrected interaction energy
    using monomer energies evaluated in the full dimer basis:

        E_int^CP = E_AB(AB) - E_A(AB) - E_B(AB)

    Parameters
    ----------
    pickle_file : str or pathlib.Path
        Path to the pickle file containing parsed raw total energies.

    molecule : str
        Molecular system label (e.g., "H2O").

    Returns
    -------
    dict
        Nested dictionary containing VMFC (counterpoise) two-body interaction
        energies:
            result[molecule][cluster][motion][config][method][basis]["energy"]

    Notes
    -----
    - Metadata keys (e.g., '_params', '_config_label') are ignored.
    - This function assumes that A_AB and B_AB energies are present in the
      raw data.
    - This function does not perform any file I/O.
    """
    raw_data = open_pickle(pickle_file)['raw_data']
    two_body_vmfc={}
    two_body_vmfc[molecule] = {}
    for cluster in raw_data[molecule]:
        two_body_vmfc[molecule][cluster] = {}
        #--
        for motion in raw_data[molecule][cluster]:
            two_body_vmfc[molecule][cluster][motion] = {}
            #-- 
            for config in raw_data[molecule][cluster][motion]:
                two_body_vmfc[molecule][cluster][motion][config] = {}
                #-- 
                for method in raw_data[molecule][cluster][motion][config]:
                    two_body_vmfc[molecule][cluster][motion][config][method] = {}
                    #--
                    for basis in raw_data[molecule][cluster][motion][config][method]:
                        two_body_vmfc[molecule][cluster][motion][config][method][basis] = {}
                        #--
                        E_A_AB  = raw_data[molecule][cluster][motion][config][method][basis]['A_AB']
                        E_B_AB = raw_data[molecule][cluster][motion][config][method][basis]['B_AB'] 
                        E_AB_AB = raw_data[molecule][cluster][motion][config][method][basis]['AB_AB']
                        energy = E_AB_AB -(E_A_AB + E_B_AB)
                        #--
                        two_body_vmfc[molecule][cluster][motion][config][method][basis]['energy']= energy
    
    return two_body_vmfc      
    
def many_body_int_pickle(pickle_file, molecule, outdir="."):
    """
    Generate and serialize many-body interaction energy data.

    This function builds the two-body interaction energy dictionaries
    (standard and VMFC) for a given molecule and writes them to a single
    pickle file named:

        many_body_interactions.pickle

    The output file is intended to serve as a central data source for
    subsequent BSSE analysis and can be extended to include higher-order
    many-body terms and additional molecules.

    Parameters
    ----------
    pickle_file : str or pathlib.Path
        Path to the raw-data pickle file for the given molecule.

    molecule : str
        Molecular system label (e.g., "H2O").

    outdir : str or pathlib.Path, optional
        Output directory where the pickle file will be written.
        Defaults to the current directory.

    Output
    ------
    Writes a pickle file containing:
        {
            "2b_no_vmfc": {molecule: ...},
            "2b_vmfc":    {molecule: ...}
        }

    Notes
    -----
    - The output file overwrites any existing file with the same name.
    - The file format is designed to accommodate multiple molecules
      in future extensions.
    """
    two_body_no_vmfc = _dimmer_no_vmfc(pickle_file,molecule)
    two_body_vmfc = _dimmer_vmfc(pickle_file,molecule)
    body_interactions = {
    "2b_no_vmfc": two_body_no_vmfc,
    "2b_vmfc": two_body_vmfc,
}
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    raw_file = outdir / f"many_body_interactions.pickle"

    with raw_file.open("wb") as f:
        pickle.dump(body_interactions, f)
    



