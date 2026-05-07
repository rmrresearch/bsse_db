"""
tar_output_discovery_4.py
--------------------------
Streams the text content of every .out file found inside *.tar.gz archives
located at the molecule directory level (e.g. data/H2O/).

Uses glob (non-recursive) to avoid accidentally reading nested backup archives
such as original_data/ or any .tar.gz files inside the cluster folders.

Yields
------
(tar_path, internal_path, text)
    tar_path      : Path  — full path to the tarball on disk
    internal_path : PurePosixPath — path of the file inside the archive
    text          : str   — decoded text content of the .out file

Usage
-----
    from tar_output_discovery_4 import iter_out_texts
    from pathlib import Path

    data_dir = Path(__file__).resolve().parent / "H2O"
    for tar_path, internal_path, text in iter_out_texts(data_dir):
        print(internal_path)
"""

import tarfile
from pathlib import Path, PurePosixPath


def iter_out_texts(data_dir):
    """
    Iterate over all .out files inside *.tar.gz archives found in data_dir.

    Parameters
    ----------
    data_dir : Path
        Directory containing the molecule-level tarballs (e.g. data/H2O/).
        Uses glob (non-recursive) — only searches one level deep.
    """
    for tar_path in sorted(data_dir.glob("*.tar.gz")):
        with tarfile.open(tar_path, "r:gz") as tar:
            for tarinfo in tar:
                if tarinfo.isfile() and tarinfo.name.endswith(".out"):
                    f = tar.extractfile(tarinfo)
                    if f is None:
                        continue
                    text = f.read().decode("utf-8", errors="replace")
                    yield (tar_path, PurePosixPath(tarinfo.name), text)


if __name__ == "__main__":
    here     = Path(__file__).resolve().parent   # data/
    data_dir = here / "H2O"

    count = 0
    for tar_path, internal_path, text in iter_out_texts(data_dir):
        if internal_path.parts[0] == "4":
            print(f"{internal_path}  ({len(text)} chars)")
            count += 1

    print(f"\nTotal .out files found in tetramer tarball: {count}")
    print(f"Expected: 65")