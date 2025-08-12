from nwchem_helpers.parse_nwchem_output import *
from test_helpers import *
import pytest


class TestCheckValue:

    def parsed_values(self):
        return {'existing key': 'Hello World'}

    def test_add_new_value_dne(self):
        values = self.parsed_values()
        check_value('a key', values, 1.23)

        corr_values = self.parsed_values()
        corr_values['a key'] = 1.23

        assert values == corr_values

    def test_add_same_existing_value(self):
        values = self.parsed_values()
        check_value('existing key', values, 'Hello World')

        corr_values = self.parsed_values()
        assert values == corr_values

    def test_add_different_existing_value(self):
        values = self.parsed_values()
        with pytest.raises(Exception):
            check_value('existing key', values, 1.23)


def test_parse_nwchem_output():

    corr = corr_nwchem_results()

    for output_file in get_output_files():
        with open(os.path.join(get_assets_dir(), output_file)) as f:
            assert parse_nwchem_output(f) == corr[output_file]
