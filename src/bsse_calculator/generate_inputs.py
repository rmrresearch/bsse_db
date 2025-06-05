import os
class GenerateInputs:
    def __init__(self):
        pass

    def get_geometry(self, HF_bond_length, **kwargs):
        if "FF_distance" in kwargs:
            FF_distance = kwargs["FF_distance"]
            return [
                [
                    f"H 0 0 -{HF_bond_length}",
                    "F 0 0 0",
                ],
                [f"F 0 0 {FF_distance}", f"H 0 0 {FF_distance + HF_bond_length}"],
            ]
        elif "F2_coords" in kwargs:
            F2_coords = kwargs["F2_coords"]
            return [
                [f"H 0 0 -{HF_bond_length}", "F 0 0 0"],
                [
                    f"F {F2_coords[0]} {F2_coords[1]} {F2_coords[2]}",
                    f"H {F2_coords[0]} {F2_coords[1]} {(F2_coords[2] + HF_bond_length)}",
                ],
            ]
        else:
            raise ValueError("Please specify either FF_distance or F2_coords")

    def prepare_input(self, geometry_str, input_file_path):
        raise NotImplementedError

    def dimer_geometry_str(self, geometry):
        raise NotImplementedError

    def bsse_corrected_monomer_geometry_str(self, geometry, ghost_monomer):
        raise NotImplementedError

    def dimer_input(self, geometry, input_file_path):
        """
        Creates psi4 input file for geometry to calculate E(AB, AB)
        """
        self.prepare_input(self.dimer_geometry_str(geometry), input_file_path)

    def bsse_corrected_monomer_input(self, geometry, input_A, input_B):
        """
        Creates psi4 input files for geometry to calculate E(A, AB) and E(B, AB)
        """
        self.prepare_input(
            self.bsse_corrected_monomer_geometry_str(geometry, 0), input_A
        )
        self.prepare_input(
            self.bsse_corrected_monomer_geometry_str(geometry, 1), input_B
        )

    def monomer_geometry_str(self, geometry):
        raise NotImplementedError

    def monomer_input(self, geometry, input_file_path):
        """Creates psi4 input file for geometry to calculate E(A, A)"""
        self.prepare_input(self.monomer_geometry_str(geometry), input_file_path)

    def run_software(self, input_files, force_run=True):
        raise NotImplementedError

    def read_energies(self, file_path):
        raise NotImplementedError

    def get_bsse(self, geometry, scripts_dir="scripts", force_rerun=True) -> dict:
        """
        Calculates BSSE for geometry (of dimer)
        Returns a dictonary of the form
        {
        "Delta_E_AB_AB": float,
        "BSSE_A": float,
        "BSSE_B": float
        }
        where Delta_E_AB_AB = ΔE(AB, AB) and BSSE_A = ε(A, AB) and BSSE_B = ε(B, AB)
        """
        os.makedirs(scripts_dir, exist_ok=True)
        self.dimer_input(geometry, os.path.join(scripts_dir, "input_dimer.txt"))
        self.bsse_corrected_monomer_input(
            geometry,
            os.path.join(scripts_dir, "input_A_AB.txt"),
            os.path.join(scripts_dir, "input_B_AB.txt"),
        )
        self.monomer_input(geometry, os.path.join(scripts_dir, "input_monomer.txt"))
        if force_rerun:
            self.run_software(
                [
                    os.path.join(scripts_dir, "input_dimer.txt"),
                    os.path.join(scripts_dir, "input_A_AB.txt"),
                    os.path.join(scripts_dir, "input_B_AB.txt"),
                    os.path.join(scripts_dir, "input_monomer.txt"),
                ]
            )
        else:
            needed_files = [
                os.path.join(scripts_dir, "input_dimer.txt"),
                os.path.join(scripts_dir, "input_A_AB.txt"),
                os.path.join(scripts_dir, "input_B_AB.txt"),
                os.path.join(scripts_dir, "input_monomer.txt"),
            ]
            if os.path.isdir(scripts_dir) and set(needed_files).issubset(
                [os.path.join(scripts_dir, file) for file in os.listdir(scripts_dir)]
            ):
                pass
            else:
                self.run_software(
                    [
                        os.path.join(scripts_dir, "input_dimer.txt"),
                        os.path.join(scripts_dir, "input_A_AB.txt"),
                        os.path.join(scripts_dir, "input_B_AB.txt"),
                        os.path.join(scripts_dir, "input_monomer.txt"),
                    ]
                )
        total_energy = self.read_energies(
            os.path.join(scripts_dir, "output_dimer.txt")
        )["Total_energy"]
        E_A_AB = self.read_energies(os.path.join(scripts_dir, "output_A_AB.txt"))[
            "Total_energy"
        ]
        E_B_AB = self.read_energies(os.path.join(scripts_dir, "output_B_AB.txt"))[
            "Total_energy"
        ]
        E_monomer = self.read_energies(os.path.join(scripts_dir, "output_monomer.txt"))[
            "Total_energy"
        ]
        return {
            "Delta_E_AB_AB": total_energy - E_A_AB - E_B_AB,
            "BSSE_A": E_monomer - E_A_AB,
            "BSSE_B": E_monomer - E_B_AB,
        }
