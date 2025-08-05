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
        print(line)
    

    while not line.strip():
        print(line)
        line = next(file)
        geom.append(line)

    return geom


    
