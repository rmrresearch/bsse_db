from pathlib import Path
import pickle

def open_pickle(pickle_file):

    path = Path(pickle_file)
    with path.open('rb') as f:
        raw_data = pickle.load(f)
        return raw_data