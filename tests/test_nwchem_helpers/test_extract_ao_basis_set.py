from nwchem_helpers.extract_ao_basis_set import *
from test_helpers import *
import os


def test_is_ao_basis_set_start():

    assert is_ao_basis_set_start('Summary of "ao basis" -> ""')
    assert not is_ao_basis_set_start('not the start')


def test_extract_ao_basis_set():
    corr_bases = corr_ao_basis_sets()

    for output_file in get_output_files():
        with open(os.path.join(get_assets_dir(), output_file)) as f:
            for line in f:
                if is_ao_basis_set_start(line):
                    basis = extract_ao_basis_set(line, f)

                    assert basis == corr_bases[output_file]
