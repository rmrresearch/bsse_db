from bsse_calculator.make_diagram import make_diagram
import subprocess, os


def test_make_diagram(tmp_path):
    geometry = [["H 0.0 0.0 -0.92", "F 0.0 0.0 0.0"], ["F 0.0 0.0 2", "H 0.0 0.0 2.92"]]
    make_diagram(geometry, os.path.join(tmp_path, "diagram.html"))
    if os.environ.get("OPEN_IMAGES"):
        subprocess.run(["open", os.path.join(tmp_path, "diagram.html")])
