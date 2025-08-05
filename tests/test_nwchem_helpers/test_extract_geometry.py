from nwchem_helpers.extract_geometry import *
from test_bsse_calculator.test_bsse_helpers import get_assets_dir
import os

corr_geoms = {
    'output_A_AB.txt' : 
    [('H', '0.00000000', '0.00000000', '-0.90000000'), 
     ('F', '0.00000000', '0.00000000', '0.10000000'), 
     ('bqF', '0.00000000', '0.00000000', '4.10000000'), 
     ('bqH', '0.00000000', '0.00000000', '5.10000000')],
    'output_B_AB.txt' :
    [('bqH', '0.00000000', '0.00000000', '-5.10000000'),
    ('bqF','0.00000000', '0.00000000', '-4.10000000'),
    ('F', '0.00000000', '0.00000000', '-0.10000000'),
    ('H', '0.00000000', '0.00000000', '0.90000000')]
}

def test_is_geometry_start():

    assert is_geometry_start('Geometry "geometry" -> ""')
    assert not is_geometry_start('not the start')

def test_extract_geometry():

    assets_dir = get_assets_dir()
    output_files = ['output_A_AB.txt', 'output_B_AB.txt']
    
    for output_file in output_files:
        with open(os.path.join(assets_dir, output_file)) as f:
            for line in f:
                if is_geometry_start(line):
                    geom = extract_geometry(f)
                 
                    assert geom == corr_geoms[output_file]
