from utils import get_bsse, get_geometry
import itertools
import os

if __name__ == "__main__":
    BSSE_As = []
    BSSE_Bs = []
    coordinate_list = []
    for coords in itertools.product([1, 2, 3, 4], repeat=3):
        geometry = get_geometry(0.924, F2_coords=coords)
        output = get_bsse(geometry, f"FF_distance_{'_'.join(map(str, coords))}")
        BSSE_As.append(output["BSSE_A"])
        BSSE_Bs.append(output["BSSE_B"])
        coordinate_list.append(coords)
    os.makedirs("data", exist_ok=True)
    with open(os.path.join("data", "BSSE_by_FF_distance.csv"), "w") as f:
        f.write("coordinates,BSSE_A,BSSE_B\n")
        for coord, BSSE_A, BSSE_B in zip(coordinate_list, BSSE_As, BSSE_Bs):
            f.write(f"{coord},{BSSE_A},{BSSE_B}\n")
