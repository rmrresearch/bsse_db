from pathlib import Path
import pickle

def save_raw_bundle(molecule, raw_data, parameters, outdir=".", cluster=None):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    # If cluster is provided, name the file raw_data_{molecule}_{cluster}.pickle
    # Otherwise fall back to raw_data_{molecule}.pickle (backwards compatible)
    if cluster is not None:
        bundle_file = outdir / f"raw_data_{molecule}_{cluster}.pickle"
    else:
        bundle_file = outdir / f"raw_data_{molecule}.pickle"

    bundle = {
        "raw_data":   raw_data,
        "parameters": parameters,
    }

    with bundle_file.open("wb") as f:
        pickle.dump(bundle, f)


    
