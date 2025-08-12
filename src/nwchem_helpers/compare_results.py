'''
Methods for comparing the parsed results of two NWChem output files.
'''

import math


def are_similar_geometries(geom0, geom1, geom_tol):
    '''
    Compares two geometries for similarity.

    This function assumes that it has been given two geometries, each obtained 
    by calling ``parse_nwchem_output`` on an NWChem output file. It is assumed
    that the geometries should have the same Cartesian coordinates (either
    because they have been generated the same way or because they have been
    reoriented according to the same convention).

    :param geom0: The first geometry
    :param geom1: The second geometry
    :param geom_tol: How 
    '''
    if len(geom0) != len(geom1):
        return False

    for a0, a1 in zip(geom0, geom1):
        if a0[0] != a1[0]:  # Different atomic symbols
            return False

        for i in range(1, 4):
            if not math.isclose(float(a0[i]), float(a1[i]), abs_tol=geom_tol):
                return False

    return True


def are_similar_energies(egy0, egy1, egy_tol):
    '''
    Determines if two energy values are similar.

    :param egy0: The first energy (in a.u.)
    :type egy0: float
    :param egy1: The second energy (in a.u.)
    :type egy1: float
    :param egy_tol: How much ``egy0`` and ``egy1`` can differ by and still be
        similar (in a.u.).
    :type egy_tol: float
    '''
    return math.isclose(float(egy0), float(egy1), abs_tol=egy_tol)


def similar_nwchem_runs(results0, results1, **kwargs):
    '''
    Determines if two runs of NWChem are "similar."

    This function takes the parsed results returned by two (potentially) 
    different calls of ``parse_nwchem_output`` and determines if they are 
    similar. Here similar is defined as:

    - Has the same list of properties.
    - Geometries are similar up to ``geometry_tolerance``.
    - Components of the energy are similar up to ``energy_tolerance``.

    This function relies on ``are_similar_geometries`` and 
    ``are_similar_energies`` for respectively comparing the geometries and
    energies in ``results0`` and ``results1``.

    :param results0: The first set of results.
    :param results1: The second set of results.
    :param **kwargs: Arbitrary keyword arguments.
    :key geometry_tolerance: The absolute value by which any geometric 
        coordinate (e.g., the x-coordinate of the third atom) may differ between
        outputs and still be considered similar (in angstroms). Default value is
        1e-7 because NWChem only prints 8 decimal places.
    :type geometry_tolerance: float
    :key energy_tolerance: The absolute value by which a component of the total
        energy may differ between outputs and still be considered similar (in
        a.u.). Default value is 1e-6 which is the default SCF convergence.
    :type energy_tolerance: float
    '''

    egy_tol = kwargs.get('energy_tolerance', 1e-6)
    geom_tol = kwargs.get('geometry_tolerance', 1e-7)

    if len(results0) != len(results1):
        return False

    keys = results0.keys()
    for key in keys:
        if not key in results1.keys():
            return False

    for k, v in results0.items():
        if k == 'Input Geometry (angstroms)':
            if not are_similar_geometries(v, results1[k], geom_tol):
                return False

        elif k == 'AO Basis Set':
            if v != results1[k]:
                return False

        else:
            if not are_similar_energies(v, results1[k], egy_tol):
                return False

    return True
