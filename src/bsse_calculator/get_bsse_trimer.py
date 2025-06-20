from itertools import chain, combinations
from bsse_calculator.add_ghosts import add_ghosts
from bsse_calculator.prepare_input import prepare_input
from bsse_calculator.get_bsse import run_software
from bsse_calculator.read_total_energy import read_energies
import os


def get_powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def get_indices(atoms):
    letter_mapping = {
        "A": 0,
        "B": 1,
        "C": 2,
    }
    return [letter_mapping[atom] for atom in atoms]


def get_bsse_trimer(trimer_input, scripts_dir):
    output = {}
    powerset = list(get_powerset(["A", "B", "C"]))
    for nonghost_atoms in powerset:
        if len(nonghost_atoms) == 0:
            continue
        possible_basis_sets = [
            all_atoms
            for all_atoms in powerset
            if set(nonghost_atoms).issubset(set(all_atoms))
        ]
        for basis in possible_basis_sets:
            calculation_name = f"E({'_'.join(nonghost_atoms)})_{'_'.join(basis)}"
            print(f"Starting calculation: {calculation_name}")
            prepare_input(
                add_ghosts(
                    trimer_input, get_indices(nonghost_atoms), get_indices(basis)
                ),
                os.path.join(scripts_dir, f"input_{calculation_name}.txt"),
                "energy",
            )
            run_software([os.path.join(scripts_dir, f"input_{calculation_name}.txt")])
            energies = read_energies(
                os.path.join(scripts_dir, f"output_{calculation_name}.txt")
            )
            output[calculation_name] = energies
    return output
