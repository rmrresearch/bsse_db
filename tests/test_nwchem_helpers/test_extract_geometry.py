# from nwchem_helpers.extract_geometry import extract_geometry
# import os

# def test_extract_geometry():

#     tests_dir = os.path.dirname(os.path.abspath(__file__))
#     assets_dir = os.path.join(tests_dir, 'test_assets', 'nwchem')
#     output_files = ['output_A_AB.txt']
    
#     for output_file in output_files:
#         with open(os.path.join(assets_dir, output_file)) as f:
#             for line in f:
#                 if 'Geometry "geometry" -> ""' in line:
#                     geom = extract_geometry(f)
