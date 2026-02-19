# build_raw_data_simple.py
import re
from pathlib import PurePosixPath
from tar_output_discovery import iter_out_texts

# --- Regex patterns (same idea as your sanity check) ---
SCF_PATTERN = re.compile(r"Total\s+SCF\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")
MP2_PATTERN = re.compile(r"Total\s+MP2\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")
CCSDT_PATTERN = re.compile(r"Total\s+CCSD\(T\)\s+energy\s*[:=]\s*([\-0-9.+Ee]+)")

def extract_last(pattern, text):
    """Return the last float matched by pattern in text, or None."""
    value = None
    for line in text.splitlines():
        m = pattern.search(line)
        if m:
            try:
                value = float(m.group(1))
            except ValueError:
                pass
    return value

def ensure_path(d, keys):
    """Ensure nested dict structure exists and return the deepest dict."""
    cur = d
    for k in keys:
        cur = cur.setdefault(k, {})
    return cur

def motion_label(component):
    """
    Map the second path component to a short label.
    - translational -> trnl
    - rotational*   -> rot
    If it's something else, keep it as-is (still dynamic).
    """
    c = component.lower()
    if c.startswith("trans"):
        return "trnl"
    if c.startswith("rot"):
        return "rot"
    return component

def parse_params(motion, config_label):
    """
    Parse configuration parameters from a config_label string, dispatching by motion.

    motion: "trnl" or "rot"
    config_label examples:
      translational:      "0_0.25_0.25_1.5"         -> x,y,z
      rotational_radial:  "3_0.0_0.0_90.0_2.76"     -> alpha,beta,gamma,R
    """
    parts = config_label.split("_")


    if motion == "trnl":
        if len(parts) != 4:
            return {}
        _, x, y, z = parts
        return {"x": float(x), "y": float(y), "z": float(z)}

    if motion == "rot":
        if len(parts) != 5:
            return {}
        _, a, b, g, r = parts
        return {"alpha": float(a), "beta": float(b), "gamma": float(g), "R": float(r)}

    return {}

def missing_energies(E_scf, E_mp2, E_ccsdt):
    """
    Identify which energy components are missing for a given output file.

    Parameters
    ----------
    E_scf : float or None
        Total SCF energy extracted from the output file.
    E_mp2 : float or None
        Total MP2 energy extracted from the output file.
    E_ccsdt : float or None
        Total CCSD(T) energy extracted from the output file.

    Returns
    -------
    list of str
        A list of energy labels that were not found in the output file.
        Possible entries are "SCF", "MP2", and "CCSD(T)".

    Notes
    -----
    This function is used to flag non-converged or incomplete calculations.
    If the returned list is non-empty, the corresponding configuration
    is excluded from the raw_data dictionary and recorded in raw_data_noconv.
    """
    
    missing = []
    if E_scf is None:
        missing.append("SCF")
    if E_mp2 is None:
        missing.append("MP2")
    if E_ccsdt is None:
        missing.append("CCSD(T)")
    return missing


def build_raw_data(molecule, verbose=True):
    """
    Build raw_data dictionary for a molecule (energies only) plus a separate
    parameters mapping.

    raw_data is keyed by:
      raw_data[molecule][cluster][motion][config_id][method_key][basis][file_key] = energy

    parameters is keyed by:
      parameters[molecule][cluster][motion][config_id] = {
          "config_label": <full config folder label>,
          "params": <dict from parse_params(...)>
      }

    Exclusion rule:
    If ANY .out file in a given (cluster/motion/config_id/method_dir/basis) folder is
    missing required energies (SCF/MP2/CCSD(T)), then that entire folder is skipped.
    """
    raw_data = {}
    raw_data_noconv = {}
    parameters = {}

    bucket = {}
    bad = set()

    # NEW: store labels seen, but don't commit params yet
    config_label_by_cfg = {}  # (cluster, motion, config_key) -> config_label

    for tar_path, internal_path, text in iter_out_texts(molecule):
        internal_path = PurePosixPath(str(internal_path))

        if len(internal_path.parts) < 2:
            continue

        cluster = internal_path.parts[0]
        motion = motion_label(internal_path.parts[1])

        config_label = internal_path.parent.parent.parent.name
        config_key = config_label.split("_", 1)[0]

        # record latest label for this config (should be consistent)
        config_label_by_cfg[(cluster, motion, config_key)] = config_label

        file_key = internal_path.stem
        basis = internal_path.parent.name
        method_dir = internal_path.parent.parent.name

        basis_folder_id = (cluster, motion, config_key, method_dir, basis)

        E_scf = extract_last(SCF_PATTERN, text)
        E_mp2 = extract_last(MP2_PATTERN, text)
        E_ccsdt = extract_last(CCSDT_PATTERN, text)

        if (E_scf is None) or (E_mp2 is None) or (E_ccsdt is None):
            bad.add(basis_folder_id)

            miss = missing_energies(E_scf, E_mp2, E_ccsdt)
            cfg_bad = ensure_path(raw_data_noconv, (molecule, cluster, motion, config_key))

            failures = cfg_bad.setdefault("failures", [])
            failures.append({
                "missing": miss,
                "method_dir": method_dir,
                "basis": basis,
                "file_key": file_key,
                "tar_path": str(tar_path),
                "internal_path": str(internal_path),
            })
            continue

        mp2_corr = E_mp2 - E_scf
        ccsdt_component = E_ccsdt - E_mp2

        bucket.setdefault(basis_folder_id, []).append(
            (file_key, E_scf, mp2_corr, ccsdt_component)
        )

    # --- Commit only good folders ---
    converged_cfgs = set()  # NEW: which (cluster,motion,config_key) actually make it in

    for (cluster, motion, config_key, method_dir, basis), records in bucket.items():
        if (cluster, motion, config_key, method_dir, basis) in bad:
            continue

        # mark config as converged (at least one good basis/method_dir survived)
        converged_cfgs.add((cluster, motion, config_key))

        for file_key, E_scf, mp2_corr, ccsdt_component in records:
            cfg_node = ensure_path(raw_data, (molecule, cluster, motion, config_key))
            ensure_path(cfg_node, ("SCF", basis))[file_key] = E_scf
            ensure_path(cfg_node, ("MP2", basis))[file_key] = mp2_corr
            ensure_path(cfg_node, ("CCSD_T", basis))[file_key] = ccsdt_component

    # NEW: build parameters only for converged configs
    for (cluster, motion, config_key) in converged_cfgs:
        config_label = config_label_by_cfg.get((cluster, motion, config_key))
        if config_label is None:
            # shouldn't happen, but safe guard
            continue
        node = ensure_path(parameters, (molecule, cluster, motion))
        node[config_key] = {
            "config_label": config_label,
            "params": parse_params(motion, config_label),
        }

    if verbose:
        if raw_data_noconv:
            print("\n⚠️  WARNING: Some configurations did NOT converge.")
            print("    Please inspect `raw_data_noconv` for details.\n")
            for mol, mol_data in raw_data_noconv.items():
                for cluster, cluster_data in mol_data.items():
                    for motion, motion_data in cluster_data.items():
                        bad_cfgs = sorted(motion_data.keys(), key=int)
                        print(f"    {mol} | cluster {cluster} | {motion} | failed configs: {bad_cfgs}")
        else:
            print("\n✅ All configurations converged successfully.\n")

    return raw_data, raw_data_noconv, parameters





if __name__ == "__main__":
    raw_data, raw_data_noconv, parameters = build_raw_data("H2O")


   

