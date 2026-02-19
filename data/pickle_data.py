from pathlib import Path
import pickle

def save_raw_bundle(molecule, raw_data, parameters, outdir="."):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    bundle_file = outdir / f"raw_data_{molecule}.pickle"
    bundle = {
        "raw_data": raw_data,
        "parameters": parameters,
    }

    with bundle_file.open("wb") as f:
        pickle.dump(bundle, f)


    
