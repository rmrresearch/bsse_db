from nwchem_helpers.extract_energy import *
from test_helpers import *
import os

# ******************************************************************************
# Tests of finding the start of energy blocks
# ******************************************************************************


def test_is_total_scf_energy_start():

    assert is_total_scf_energy_start('Total SCF energy')
    assert not is_total_scf_energy_start('not the start')


def test_is_total_mp2_energy_start():

    assert is_total_mp2_energy_start('Total MP2 energy:')
    assert not is_total_mp2_energy_start('not the start')


def test_is_total_ccsd_energy_start():

    assert is_total_ccsd_energy_start('Total CCSD energy:')
    assert not is_total_ccsd_energy_start('not the start')


def test_is_total_ccsd_t_energy_start():

    assert is_total_ccsd_t_energy_start('Total CCSD(T) energy:')
    assert not is_total_ccsd_t_energy_start('not the start')


def test_is_total_energy_start():

    assert is_total_energy_start('Total SCF energy') == 'SCF'
    assert is_total_energy_start('Total MP2 energy:') == 'MP2'
    assert is_total_energy_start('Total CCSD energy:') == 'CCSD'
    assert is_total_energy_start('Total CCSD(T) energy:') == 'CCSD(T)'
    assert not is_total_energy_start('not the start')


# ******************************************************************************
# Tests of extracting the energies
# ******************************************************************************


def test_extract_energies():
    corr_energies = {
        'SCF': corr_scf_energies(),
        'MP2': corr_mp2_energies(),
        'CCSD': corr_ccsd_energies(),
        'CCSD(T)': corr_ccsd_t_energies()
    }

    for output_file in get_output_files():
        with open(os.path.join(get_assets_dir(), output_file)) as f:
            for line in f:
                if is_total_scf_energy_start(line):
                    egy = extract_total_scf_energy(line, f)
                    assert egy == extract_total_energy('SCF', line, f)
                    assert egy == corr_energies['SCF'][output_file]

                elif is_total_mp2_energy_start(line):
                    egy = extract_total_mp2_energy(line, f)
                    assert egy == extract_total_energy('MP2', line, f)
                    assert egy == corr_energies['MP2'][output_file]

                elif is_total_ccsd_energy_start(line):
                    egy = extract_total_ccsd_energy(line, f)
                    assert egy == extract_total_energy('CCSD', line, f)
                    assert egy == corr_energies['CCSD'][output_file]

                elif is_total_ccsd_t_energy_start(line):
                    egy = extract_total_ccsd_t_energy(line, f)
                    assert egy == extract_total_energy('CCSD(T)', line, f)
                    assert egy == corr_energies['CCSD(T)'][output_file]
