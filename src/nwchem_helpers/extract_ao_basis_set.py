'''
Functions for extracting the AO basis set from an NWChem output file.
'''


def is_ao_basis_set_start(line):
    '''Wraps the logic for detecting the start of the "AO basis set block"'''

    return 'Summary of "ao basis" -> ""' in line


def extract_ao_basis_set(line, file):
    '''
    Extracts tha AO basis set used for a calculation.

    This function assumes that ``line`` is pointing at a line of ``file`` for
    which ``is_ao_basis_set_start`` returns true.
    '''

    bases = set()

    for _ in range(4):
        line = next(file)

    while line.split():
        bases.add(line.split()[1])
        line = next(file)

    # For now expecting the same basis set on all atoms
    assert len(bases) == 1

    return list(bases)[0]
