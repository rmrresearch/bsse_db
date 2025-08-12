'''
Functions to facilitate extracting computed energies from NWChem outputs.

If you want to add the ability to extract a new energy type, please add a 
function for detecting the line containing the energy (the "energy block") and a
function for extracting the energy. Please also add both functions to their
respective convenience wrappers.

N.b., even though the individual functions presently share much of the same 
logic, please do not try to combine them. In practice, many of the energies
appear in multiple places in the output file and more fleshed out versions would
detect all such places and extract all such values (for sanity checks). This
would cause the implementations to diverge justifying the current design.
'''

# ******************************************************************************
# Functions for detecting the start of an energy block
# ******************************************************************************


def is_total_scf_energy_start(line):
    '''Wraps the logic for detecting the start of the "SCF Energy Block"'''

    return 'Total SCF energy' in line


def is_total_mp2_energy_start(line):
    '''Wraps the logic for detecting the start of the "MP2 Energy Block"'''

    return 'Total MP2 energy:' in line


def is_total_ccsd_energy_start(line):
    '''Wraps the logic for detecting the start of the "CCSD Energy Block"'''

    return 'Total CCSD energy:' in line


def is_total_ccsd_t_energy_start(line):
    '''Wraps the logic for detecting the start of the "CCSD(T) Energy Block"'''

    return 'Total CCSD(T) energy:' in line


def is_total_energy_start(line):
    '''
    Convenience driver for determining if `line` is the start of an energy
    block. 
    
    If `line` is the start of an energy block, this function returns a
    non-empty string containing the type of energy. If the line is NOT the start
    of such a block this function returns False.
    '''

    if is_total_scf_energy_start(line):
        return 'SCF'

    elif is_total_mp2_energy_start(line):
        return 'MP2'

    elif is_total_ccsd_energy_start(line):
        return 'CCSD'

    elif is_total_ccsd_t_energy_start(line):
        return 'CCSD(T)'

    else:
        return False


# ******************************************************************************
# Functions for extracting the energies
# ******************************************************************************


def extract_total_scf_energy(line, file):
    '''
    
    This function assumes that the iterator `file` is pointing at `line` for
    which `is_total_scf_energy_start(line)` returned true.
    '''

    return line.split()[-1]


def extract_total_mp2_energy(line, file):
    '''
    
    This function assumes that the iterator `file` is pointing at `line` for
    which `is_total_mp2_energy_start(line)` returned true.
    '''

    return line.split()[-1]


def extract_total_ccsd_energy(line, file):
    '''
    
    This function assumes that the iterator `file` is pointing at `line` for
    which `is_total_ccsd_energy_start(line)` returned true.
    '''

    return line.split()[-1]


def extract_total_ccsd_t_energy(line, file):
    '''
    
    This function assumes that the iterator `file` is pointing at `line` for
    which `is_total_ccsd_t_energy_start(line)` returned true.
    '''

    return line.split()[-1]


def extract_total_energy(egy_name, line, file):
    '''
    Convenience driver for dispatching to the appropriate energy extraction
    function based on a string.

    
    '''

    if egy_name == 'SCF':
        return extract_total_scf_energy(line, file)

    elif egy_name == 'MP2':
        return extract_total_mp2_energy(line, file)

    elif egy_name == 'CCSD':
        return extract_total_ccsd_energy(line, file)

    elif egy_name == 'CCSD(T)':
        return extract_total_ccsd_t_energy(line, file)
