import tarfile
import glob
import os

parent_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

paths = [
    os.path.join(parent_folder, "examples", "data", "nwchem"),
    os.path.join(parent_folder, "examples", "data", "psi4"),
]
for path in paths:
    files = glob.glob("output_*.txt", root_dir=path)
    with tarfile.open(os.path.join(path, "outputs.tar"), "w") as tar:
        for file in [os.path.join(path, file) for file in files]:
            tar.add(file)
            os.remove(file)
