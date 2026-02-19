import tarfile
from pathlib import Path, PurePosixPath

def iter_out_texts(molecule):
    """
    Parameters
    ----------
    molecule : str
        Name of the molecule directory (e.g., H2O, Ne, HF).

    Yields
    ------
    (Path, PurePosixPath, str)
        The tar archive path, internal path inside the archive,
        and the decoded text of the .out file.
    """
    molecule_dir = Path(molecule)

    for tar_path in molecule_dir.rglob("*.tar.gz"):
        with tarfile.open(tar_path, "r:gz") as tar:
            for tarinfo in tar:
                if tarinfo.isfile() and tarinfo.name.endswith(".out"):
                    f = tar.extractfile(tarinfo)
                    if f is None:
                        continue
                    text = f.read().decode("utf-8", errors="replace")
                    yield (tar_path, PurePosixPath(tarinfo.name), text)
