# tar_output_discovery_3.py
import tarfile
from pathlib import Path, PurePosixPath

def iter_out_texts(data_dir):
    """
    Parameters
    ----------
    data_dir : Path
        Path to the molecule data directory (e.g., .../data/H2O/).
        Only searches for *.tar.gz files directly in this directory
        (non-recursive) to avoid picking up nested backups like
        H20_trimer_original.tar.gz inside 3/.

    Yields
    ------
    (Path, PurePosixPath, str)
        The tar archive path, internal path inside the archive,
        and the decoded text of the .out file.
    """
    for tar_path in data_dir.glob("*.tar.gz"):
        with tarfile.open(tar_path, "r:gz") as tar:
            for tarinfo in tar:
                if tarinfo.isfile() and tarinfo.name.endswith(".out"):
                    f = tar.extractfile(tarinfo)
                    if f is None:
                        continue
                    text = f.read().decode("utf-8", errors="replace")
                    yield (tar_path, PurePosixPath(tarinfo.name), text)