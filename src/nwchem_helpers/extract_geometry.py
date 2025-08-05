'''
Functions to facilitate extracting input geometries from NWchem outputs.
'''


def is_geometry_start(line):
    '''
    Wraps the logic for determining if the current line starts the geometry
    block.
    '''

    return 'Geometry "geometry" -> ""' in line

def extract_geometry(file):
    """
    Extracts the geometry from an NWChem file.

    Given a file iterator `file` currently pointing at the line:

    'Geometry "geometry" -> ""

    This function will advance the iterator 7 times (skips over dashes, blank
    lines, conversion factor, and column headings) to get to the first atom.

    """

    geom = []

    line = None
    for _ in range(7):
        line = next(file)
    
    # line points to the first atom, iterate until we hit a blank line
    while line.split():
        (idx, sym, q, x, y, z) = line.split()
        geom.append((sym, x, y, z))
        line = next(file)

    return geom


    
