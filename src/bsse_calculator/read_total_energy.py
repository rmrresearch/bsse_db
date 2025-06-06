def read_energies(file_path):
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
