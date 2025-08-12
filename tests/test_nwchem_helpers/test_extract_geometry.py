from nwchem_helpers.extract_geometry import *
from test_helpers import *
import os


def test_is_geometry_start():

    assert is_geometry_start('Geometry "geometry" -> ""')
    assert not is_geometry_start('not the start')


def test_extract_geometry():
    corr_geoms = corr_input_geometries()

    for output_file in get_output_files():
        with open(os.path.join(get_assets_dir(), output_file)) as f:
            for line in f:
                if is_geometry_start(line):
                    geom = extract_geometry(line, f)

                    assert geom == corr_geoms[output_file]
