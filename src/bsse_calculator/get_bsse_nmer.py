from itertools import chain, combinations
from bsse_calculator.add_ghosts import add_ghosts
from bsse_calculator.prepare_input import prepare_input
from bsse_calculator.get_bsse import run_software
from bsse_calculator.read_total_energy import read_energies
import os
import glob
from tqdm import tqdm

letter_mapping = {
    "A": 0,
    "B": 1,
    "C": 2,
    "D": 3,
    "E": 4,
    "F": 5,
    "G": 6,
    "H": 7,
    "I": 8,
    "J": 9,
}


def get_powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def get_indices(letters, letter_mapping=letter_mapping):
    return [letter_mapping[atom] for atom in letters]


def get_bsse_nmer(
    n, nmer_input, scripts_dir, basis_set, force=True, max_nbody_interaction=5
):
    output = {}
    letter_mapping = {
        "A": 0,
        "B": 1,
        "C": 2,
        "D": 3,
        "E": 4,
        "F": 5,
        "G": 6,
        "H": 7,
        "I": 8,
        "J": 9,
        "K": 10,
    }
    powerset = list(get_powerset(list(letter_mapping.keys())[:n]))
    for nonghost_atoms in tqdm(powerset):
        if len(nonghost_atoms) == 0:
            continue
        if len(nonghost_atoms) > max_nbody_interaction:
            continue
        possible_basis_sets = [
            all_atoms
            for all_atoms in powerset
            if set(nonghost_atoms).issubset(set(all_atoms))
        ]
        for basis in tqdm(possible_basis_sets):
            if len(basis) > max_nbody_interaction:
                continue
            calculation_name = f"E({''.join(basis)})_{''.join(nonghost_atoms)}"
            file_name = f"E_{''.join(basis)}_{''.join(nonghost_atoms)}"
            print(f"Starting calculation: {calculation_name}")
            prepare_input(
                add_ghosts(
                    nmer_input,
                    get_indices(nonghost_atoms, letter_mapping),
                    get_indices(basis, letter_mapping),
                    basis_set,
                ),
                os.path.join(scripts_dir, f"input_{file_name}.txt"),
                "energy",
            )
            if force:
                run_software([os.path.join(scripts_dir, f"input_{file_name}.txt")])
            else:
                if not os.path.exists(
                    os.path.join(scripts_dir, f"output_{file_name}.txt")
                ):
                    run_software([os.path.join(scripts_dir, f"input_{file_name}.txt")])
            energies = read_energies(
                os.path.join(scripts_dir, f"output_{file_name}.txt")
            )
            output[calculation_name] = energies
            for file_path in glob.glob("input_E_*"):
                os.remove(file_path)
    return output
