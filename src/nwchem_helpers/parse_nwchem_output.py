'''
Driver for extracting useful information out of an NWChem output file
'''

from nwchem_helpers.extract_energy import *
from nwchem_helpers.extract_geometry import is_geometry_start, extract_geometry


def check_value(key, parsed_values, new_value):
    '''
    Wraps the process of adding a parsed result to a dictionary.

    This function will check if `parsed_values` contains a "key, value" pair
    with the key `key`. If the pair does not exist, one will be created using
    `key` as the key and `new_value` as the value. If the pair already exists,
    this function will verify that `new_value` is the same as the value of the
    existing pair, raising an exception if this is not the case.

    :param key: The label for the result
    :type key: str
    :param parsed_values: The dictionary which may or may not already contain a
        result labeled `key`.
    :type parsed_values: dict
    :param new_value: The value to associate with `key`
    :type new_value: any
    '''

    if key not in parsed_values:
        parsed_values[key] = new_value

    elif parsed_values[key] != new_value:
        msg = 'Value for {} does not match previous value.'.format(key)
        raise Exception(msg)


def parse_nwchem_output(file):
    '''
    Parses the NWChem output that the interable `file` refers to. Generally
    speaking `file` is created by doing `with open(path_to_file) as file:`.

    N.b., all results come back as strings to preserve the representation found
    in the output file. Conversions to 
    '''
    parsed_values = {}

    for line in file:

        # These are "signals" indicating we found a result of interest
        egy_type = is_total_energy_start(line)
        is_geom = is_geometry_start(line)

        # These are the "signal handlers" that extract the result of interest
        if egy_type:
            egy_value = extract_total_energy(egy_type, line, file)
            egy_key = 'Total {} Energy (a.u.)'.format(egy_type)
            check_value(egy_key, parsed_values, egy_value)

        elif is_geom:
            geom = extract_geometry(line, file)
            geom_key = 'Input Geometry (angstroms)'
            check_value(geom_key, parsed_values, geom)

    return parsed_values
