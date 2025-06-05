from bsse_calculator.generate_inputs import GenerateInputs
import subprocess, os, shutil


class GenerateInputsPsi4(GenerateInputs):
    def prepare_input(self, geometry_str, input_file_path):
        with open(input_file_path, "w") as f:
            f.write("memory 600 mb\n")
            f.write("\n")
            f.write("molecule HF {\n")
            f.write(geometry_str)
            f.write("}\n")
            f.write("\n")
            f.write("set basis aug-cc-pVDZ\n")
            f.write("energy('ccsd(t)')\n")

    def dimer_geometry_str(self, geometry):
        ret = ""
        for i, molecule in enumerate(geometry):
            for line in molecule:
                ret += f"  {line}\n"
            if i != len(geometry) - 1:
                ret += "  --\n"
        return ret

    def bsse_corrected_monomer_geometry_str(self, geometry, ghost_monomer):
        ret = ""
        for i, molecule in enumerate(geometry):
            for line in molecule:
                ret += f"  {line}\n" if i == ghost_monomer else f"  @{line}\n"
            if i != len(geometry) - 1:
                ret += f"  --\n"
        return ret

    def monomer_geometry_str(self, geometry):
        monomer = geometry[0]
        ret = ""
        for line in monomer:
            ret += f"  {line}\n"
        return ret

    def read_energies(self, file_path):
        """
        Returns a dictionary of the form
        {
        "Total_energy": float,
        "MP2_correlation_energy": float,
        "CCSD(T)_correlation_energy": float
        }
        output
        """
        with open(file_path, "r") as f:
            for line in f:
                if "CCSD(T) total energy  " in line:
                    total_energy = float(line.split()[-1])
                if "MP2 correlation energy" in line:
                    mp2_correlation_energy = float(line.split()[-1])
                if "CCSD(T) total energy" in line:
                    ccsd_total_energy = float(line.split()[-1])
                if "SCF energy" in line:
                    scf_energy = float(line.split()[-1])
            return {
                "Total_energy": total_energy,
                "MP2_correlation_energy": mp2_correlation_energy,
                "CCSD(T)_correlation_energy": ccsd_total_energy
                - scf_energy
                - mp2_correlation_energy,
                "SCF_energy": scf_energy,
            }

    def run_software(self, input_files, force_run=True):
        """
        Runs psi4 on the input files.

        Args:
            input_files (list): list of input files
            force_run (bool): whether to run psi4 or not if output files already exist
        Returns the directory of the output path"""
        if not force_run:
            needed_files = [
                "output_dimer.txt",
                "output_A_AB.txt",
                "output_B_AB.txt",
                "output_monomer.txt",
            ]
            if set(needed_files).issubset(os.listdir(os.getcwd())):
                return
            else:
                for input_file in input_files:
                    subprocess.run(
                        ["psi4", input_file, input_file.replace("input", "output")]
                    )
        else:
            for input_file in input_files:
                subprocess.run(
                    ["psi4", input_file, input_file.replace("input", "output")]
                )
