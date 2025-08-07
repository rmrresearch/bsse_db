from nwchem_helpers.compare_results import *
from test_helpers import *


class TestAreSimilarGeometries:

    def test_is_same(self):
        geoms = corr_input_geometries()
        for k, v in geoms.items():
            assert are_similar_geometries(v, v, 1E-6)

    def test_are_different(self):
        geoms = corr_input_geometries()
        keys = list(geoms.keys())
        assert not are_similar_geometries(geoms[keys[0]], geoms[keys[1]], 1E-6)

    def test_similar_within_tol(self):
        geoms = corr_input_geometries()
        diff_geoms = corr_input_geometries()
        for k, v in diff_geoms.items():
            v[0] = (v[0][0], v[0][1], v[0][2], str(float(v[0][3]) + 1E-7))
            assert are_similar_geometries(geoms[k], v, 1E-6)

    def test_differ_more_than_tol(self):
        geoms = corr_input_geometries()
        diff_geoms = corr_input_geometries()
        for k, v in diff_geoms.items():
            v[0] = (v[0][0], v[0][1], v[0][2], str(float(v[0][3]) + 1E-5))
            assert not are_similar_geometries(geoms[k], v, 1E-6)


class TestAreSimilarEnergies:

    def test_is_same(self):
        egys = corr_scf_energies()
        for k, v in egys.items():
            assert are_similar_energies(v, v, 1E-6)

    def test_are_different(self):
        egys = corr_scf_energies()
        for k, v in egys.items():
            # They are all non-zero
            assert not are_similar_energies(v, 0.0, 1E-6)

    def test_similar_within_tol(self):
        egys = corr_scf_energies()
        diff_egys = corr_scf_energies()
        for k, v in diff_egys.items():
            v = str(float(v) + 1E-7)
            assert are_similar_energies(egys[k], v, 1E-6)

    def test_differ_more_than_tol(self):
        egys = corr_scf_energies()
        diff_egys = corr_scf_energies()
        for k, v in diff_egys.items():
            v = str(float(v) + 1E-5)
            assert not are_similar_energies(egys[k], v, 1E-6)


class TestSimilarNWChemRuns:

    def test_is_same(self):
        for k, v in corr_nwchem_results().items():
            assert similar_nwchem_runs(v, v)

    def test_different_sizes(self):
        for k, v in corr_nwchem_results().items():
            assert not similar_nwchem_runs(v, {})

    def test_different_keys(self):
        results0 = {'Total SCF Energy (a.u.)': -0.5}
        results1 = {'Total MP2 Energy (a.u.)': -0.5}

        assert not similar_nwchem_runs(results0, results1)

    def test_different_bases(self):
        results0 = {'AO Basis Set': 'aug-cc-pvdz'}
        results1 = {'AO Basis Set': 'aug-cc-pvtz'}

        assert not similar_nwchem_runs(results0, results1)

    def test_slightly_different_geometries(self):
        keys = get_output_files()
        r0 = corr_nwchem_results()[keys[0]]
        r1 = corr_nwchem_results()[keys[0]]
        geom = r1['Input Geometry (angstroms)']
        newz = str(float(geom[0][3]) + 1.0e-7)
        geom[0] = (geom[0][0], geom[0][1], geom[0][2], newz)

        assert similar_nwchem_runs(r0, r1)  # Passes with default
        assert not similar_nwchem_runs(r0, r1, geometry_tolerance=1E-8)

    def test_slightly_different_energies(self):
        keys = get_output_files()
        r0 = corr_nwchem_results()[keys[0]]
        r1 = corr_nwchem_results()[keys[0]]
        scf_key = 'Total SCF Energy (a.u.)'
        r1[scf_key] = str(float(r1[scf_key]) + 1E-6)

        assert similar_nwchem_runs(r0, r1)  # Passes with default
        assert not similar_nwchem_runs(r0, r1, energy_tolerance=1E-7)
