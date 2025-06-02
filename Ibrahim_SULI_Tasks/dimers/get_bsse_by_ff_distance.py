from utils import get_bsse, get_geometry
import os

if __name__ == "__main__":
    BSSE_As = []
    BSSE_Bs = []
    FF_distances = [3, 4, 5, 6]
    for FF_distance in FF_distances:
        geometry = get_geometry(0.924, FF_distance)
        output = get_bsse(geometry, f"FF_distance_{FF_distance}")
        BSSE_As.append(output["BSSE_A"])
        BSSE_Bs.append(output["BSSE_B"])
    os.makedirs("data", exist_ok=True)
    with open(os.path.join("data", "BSSE_by_FF_distance.csv"), "w") as f:
        f.write("FF_distance,BSSE_A,BSSE_B\n")
        for BSSE_A, BSSE_B, FF_distance in zip(BSSE_As, BSSE_Bs, FF_distances):
            f.write(f"{FF_distance},{BSSE_A},{BSSE_B}\n")
