"""
check_raw_data_keys.py
======================
Quick sanity check before running many_body_interactions_3.py.

Verifies that all 19 file keys required by many_body_interactions_3.py
are present in raw_data_H2O_3.pickle for every config, method, and basis.

Run from data/:
    python check_raw_data_keys.py
"""

import pickle
from pathlib import Path

# ---- The 19 required file keys, grouped by physical meaning ----
REQUIRED_KEYS = {
    'Monomers in own basis (3)':    ['A_A',   'B_B',   'C_C'],
    'Monomers in dimer basis (6)':  ['A_AB',  'B_AB',  'A_AC',
                                     'C_AC',  'B_BC',  'C_BC'],
    'Monomers in trimer basis (3)': ['A_ABC', 'B_ABC', 'C_ABC'],
    'Dimers in own basis (3)':      ['AB_AB', 'AC_AC', 'BC_BC'],
    'Dimers in trimer basis (3)':   ['AB_ABC','AC_ABC','BC_ABC'],
    'Full trimer (1)':              ['ABC_ABC'],
}

ALL_REQUIRED = [k for keys in REQUIRED_KEYS.values() for k in keys]

# ---- Load pickle ----
here = Path(__file__).resolve().parent
pickle_path = here / 'raw_data_H2O_3.pickle'

print(f'Reading: {pickle_path}\n')
with open(pickle_path, 'rb') as f:
    stored = pickle.load(f)

raw_data = stored['raw_data']

# ---- Discover structure ----
mol      = 'H2O'
cluster  = '3'
mol_data = raw_data[mol][cluster]
configs  = sorted(mol_data.keys(), key=int)
methods  = list(mol_data[configs[0]].keys())
bases    = list(mol_data[configs[0]][methods[0]].keys())

print(f'Configs  : {configs}')
print(f'Methods  : {methods}')
print(f'Bases    : {bases}')
print(f'Total required file keys: {len(ALL_REQUIRED)}')
print()

# ---- Check each config / method / basis ----
all_good = True

for cfg in configs:
    for method in methods:
        for basis in bases:
            keys_found = set(mol_data[cfg][method][basis].keys())
            missing    = [k for k in ALL_REQUIRED if k not in keys_found]
            extra      = sorted(keys_found - set(ALL_REQUIRED))

            status = '✅' if not missing else '❌'
            print(f'{status} cfg={cfg}  method={method}  basis={basis}'
                  f'  →  {len(keys_found)} keys found')

            if missing:
                all_good = False
                print(f'   ⚠️  MISSING : {missing}')
            if extra:
                print(f'   ℹ️  EXTRA   : {extra}  (not required but present)')

# ---- Summary ----
print()
if all_good:
    print('✅ All 19 required keys present for every config/method/basis.')
    print('   many_body_interactions_3.py is safe to run.')
else:
    print('❌ Some keys are missing. Check the output above.')
    print('   Do not run many_body_interactions_3.py until resolved.')

# ---- Bonus: print all keys found for first config as reference ----
print()
print(f'--- All keys in cfg={configs[0]}, method={methods[0]}, basis={bases[0]} ---')
all_keys = sorted(mol_data[configs[0]][methods[0]][bases[0]].keys())
for group_name, group_keys in REQUIRED_KEYS.items():
    print(f'  {group_name}:')
    for k in group_keys:
        val = mol_data[configs[0]][methods[0]][bases[0]].get(k, 'MISSING')
        if val == 'MISSING':
            print(f'    {k:12s} = MISSING ❌')
        else:
            print(f'    {k:12s} = {val:.10f}')