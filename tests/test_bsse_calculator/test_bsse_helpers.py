import os

def get_assets_dir():
    '''Works out the path to the NwChem files used for testing'''
    
    test_bsse_calc_dir = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.dirname(test_bsse_calc_dir)
    return os.path.join(tests_dir, 'test_assets', 'nwchem')