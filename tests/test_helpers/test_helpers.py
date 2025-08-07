import os


def get_assets_dir():
    '''Works out the path to the NwChem files used for testing'''

    test_bsse_calc_dir = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.dirname(test_bsse_calc_dir)
    return os.path.join(tests_dir, 'test_assets', 'nwchem')


def get_output_files():
    '''Returns a list of output files that reside in the assets dir'''

    return ['output_A_AB.txt', 'output_B_AB.txt']


# ******************************************************************************
# Correct values that should be extracted from the output files
# ******************************************************************************


def corr_input_geometries():
    return {
        'output_A_AB.txt': [('H', '0.00000000', '0.00000000', '-0.90000000'),
                            ('F', '0.00000000', '0.00000000', '0.10000000'),
                            ('bqF', '0.00000000', '0.00000000', '4.10000000'),
                            ('bqH', '0.00000000', '0.00000000', '5.10000000')],
        'output_B_AB.txt': [('bqH', '0.00000000', '0.00000000', '-5.10000000'),
                            ('bqF', '0.00000000', '0.00000000', '-4.10000000'),
                            ('F', '0.00000000', '0.00000000', '-0.10000000'),
                            ('H', '0.00000000', '0.00000000', '0.90000000')]
    }


def corr_ao_basis_sets():
    return {'output_A_AB.txt': 'aug-cc-pvdz', 'output_B_AB.txt': 'aug-cc-pvdz'}


def corr_scf_energies():
    return {
        'output_A_AB.txt': '-100.023680528337',
        'output_B_AB.txt': '-100.023680528357'
    }


def corr_mp2_energies():
    return {
        'output_A_AB.txt': '-100.252984659592130',
        'output_B_AB.txt': '-100.252984659611371'
    }


def corr_ccsd_energies():
    return {
        'output_A_AB.txt': '-100.256011486273863',
        'output_B_AB.txt': '-100.256011486309021'
    }


def corr_ccsd_t_energies():
    return {
        'output_A_AB.txt': '-100.260535167198398',
        'output_B_AB.txt': '-100.260535167257615'
    }


def corr_nwchem_results():
    bases = corr_ao_basis_sets()
    geoms = corr_input_geometries()
    scf_egys = corr_scf_energies()
    mp2_egys = corr_mp2_energies()
    ccsd_egys = corr_ccsd_energies()
    ccsd_t_egys = corr_ccsd_t_energies()

    rv = {}

    for output_file in get_output_files():
        rv[output_file] = {
            'AO Basis Set': bases[output_file],
            'Input Geometry (angstroms)': geoms[output_file],
            'Total SCF Energy (a.u.)': scf_egys[output_file],
            'Total MP2 Energy (a.u.)': mp2_egys[output_file],
            'Total CCSD Energy (a.u.)': ccsd_egys[output_file],
            'Total CCSD(T) Energy (a.u.)': ccsd_t_egys[output_file]
        }

    return rv
