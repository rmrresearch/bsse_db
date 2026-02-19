# build_trnl_dataframe.py

# ---- path setup ----
import sys
from pathlib import Path

here = Path(__file__).resolve()
data_dir = here.parent.parent.parent  # bsse_db/data
sys.path.insert(0, str(data_dir))

# ---- imports that depend on that path ----
from open_pickle_data import open_pickle
from two_body_bsse import twobody_bsse

# ---- standard imports ----
import pandas as pd


def main():
    # pickle files
    raw_pickle = data_dir / "raw_data_H2O.pickle"
    many_body_pickle = data_dir / "many_body_interactions.pickle"

    # Data dictionaries from pickle files
    two_body_bsse = twobody_bsse(many_body_pickle, "H2O")
    parameters = open_pickle(raw_pickle)["parameters"]

    # Rotational information
    rot_bsse = two_body_bsse["H2O"]["2"]["rot"]
    rot_pars = parameters["H2O"]["2"]["rot"]

    methods = ["SCF", "MP2", "CCSD_T"]
    basis = "aug-cc-pvdz"

    for method in methods:
        rows = []

        for cfg in sorted(rot_bsse.keys(), key=int):
            p = rot_pars[cfg]["params"]
            e = rot_bsse[cfg][method][basis]

            rows.append(
                {
                    "config": cfg,
                    "alpha": p["alpha"],
                    "beta": p["beta"],
                    "gamma": p["gamma"],                    
                    "R":p["R"],
                    "method": method,
                    "basis": basis,
                    "no_vmfc": e["no_vmfc"],
                    "vmfc": e["vmfc"],
                    "bsse": e["bsse"],
                }
            )

        df = pd.DataFrame(rows)

        out_file = here.parent / f"rot_{method}_{basis}.csv"
        df.to_csv(out_file, index=False)

        print(method, df.shape, "saved to:", out_file)


if __name__ == "__main__":
    main()
