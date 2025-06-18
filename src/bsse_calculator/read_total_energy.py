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
        total_energy = None
        paren_t_correlation_energy = None
        CSSD_corr_energy = None
        mp2_correlation_energy = None
        scf_energy = None
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
        for energy in [
            total_energy,
            paren_t_correlation_energy,
            CSSD_corr_energy,
            mp2_correlation_energy,
            scf_energy,
        ]:
            if energy is None:
                raise ValueError(
                    f"The output ({file_path}) does not contain all energy values"
                )
        if total_energy > 0:
            raise ValueError("The total energy is positive, something went wrong")
        return {
            "Total_energy": total_energy,
            "MP2_correlation_energy": mp2_correlation_energy,
            "(T)_correlation_energy": paren_t_correlation_energy,
            "CCSD_correlation_energy": CSSD_corr_energy - mp2_correlation_energy,
            "SCF_energy": scf_energy,
        }
