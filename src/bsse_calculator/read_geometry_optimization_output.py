import re

table_entry_pattern = re.compile(
    r"^\s*(\d+)\s+(\w+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)"
)


def find_last_line_number(filepath):
    current_latest_line = None
    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            if "Output coordinates in angstroms" in line:
                current_latest_line = i
    return current_latest_line


def read_geometry_optimization_output(filename):
    """Reads the geometry optimization output file and returns the geometry as a list of strings"""
    start = find_last_line_number(filename)
    with open(filename, "r") as f:
        lines = f.readlines()[start + 1 :]
    line = lines.pop(0)
    geometry = []
    while "Atomic Mass" not in line:
        match = table_entry_pattern.match(line)
        if match:
            no, tag, charge, x, y, z = match.groups()
            geometry.append(f"{tag} {x} {y} {z}")
        line = lines.pop(0)
    return geometry
