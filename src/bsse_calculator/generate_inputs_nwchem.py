from bsse_calculator.generate_inputs import GenerateInputs
import subprocess, os


class GenerateInputsNwchem(GenerateInputs):
    def prepare_input(self, geometry_str, input_file_path):
        with open(input_file_path, "w") as f:
            f.write("geometry\n")
            f.write(geometry_str)
            f.write("end\n")
            f.write("basis spherical\n")
            f.write("  F library aug-cc-pvdz\n")
            f.write("  H library aug-cc-pvdz\n")
            f.write("end\n")
            f.write("task ccsd(t) energy\n")

    def dimer_geometry_str(self, geometry):
        ret = ""
        for molecule in geometry:
            for line in molecule:
                ret += f"  {line}\n"
        return ret

    def bsse_corrected_monomer_geometry_str(self, geometry, ghost_monomer):
        ret = ""
        for i, molecule in enumerate(geometry):
            for line in molecule:
                ret += f"  {line}\n" if i == ghost_monomer else f"  bq{line}\n"
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
                if "Total CCSD(T) energy" in line:
                    total_energy = float(line.split()[-1])
                if "(T) corr. energy" in line:
                    paren_t_correlation_energy = float(line.split()[-1])
                if "CCSD corr. energy" in line:
                    CSSD_corr_energy = float(line.split()[-1])
                if "MP2 Corr. energy:" in line:
                    mp2_correlation_energy = float(line.split()[-1])
                if "Total SCF energy" in line:
                    scf_energy = float(line.split()[-1])
            return {
                "Total_energy": total_energy,
                "MP2_correlation_energy": mp2_correlation_energy,
                "(T)_correlation_energy": paren_t_correlation_energy,
                "CCSD corr.energy": CSSD_corr_energy,
                "SCF_energy": scf_energy,
            }

    def monomer_geometry_str(self, geometry):
        monomer = geometry[0]
        ret = ""
        for line in monomer:
            ret += f"  {line}\n"
        return ret

    def run_software(self, input_files, force_run=True):
        """
        Runs nwchem on the input files.

        Args:
            input_files (list): list of input files
            force_run (bool): whether to run nwchem or not if output files already exist
        Returns the directory of the output path"""
        if force_run:
            for input_file in input_files:
                output_file = input_file.replace("input", "output")
                with open(output_file, "w") as out:
                    subprocess.run(["nwchem", input_file], stdout=out)
        else:
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
                    output_file = input_file.replace("input", "output")
                    with open(output_file, "w") as out:
                        subprocess.run(["nwchem", input_file], stdout=out)
