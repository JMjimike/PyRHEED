import numpy as np
from pymatgen.core import structure as pgStructure
from pymatgen.core.lattice import Lattice
from pymatgen.io.cif import CifParser

_AXIS_INDEX = {"a": 0, "b": 1, "c": 2}
_ALL_AXES = ("a", "b", "c")


def _dominant_species(site):
    if site.is_ordered:
        return site.specie
    return max(site.species.items(), key=lambda x: x[1])[0]


def load_ordered_structure(cif_path, structure=None):
    if structure is None:
        structure = CifParser(cif_path).parse_structures(primitive=False)[0]
    species = []
    coords = []
    for site in structure.sites:
        species.append(_dominant_species(site))
        coords.append(site.frac_coords)
    return pgStructure.Structure(structure.lattice, species, coords)


def structure_to_molecule(structure):
    species = [_dominant_species(site) for site in structure.sites]
    return pgStructure.Molecule(species, structure.cart_coords)


def permute_crystal_axes(structure, z_axis, x_axis):
    """Remap CIF axes so z_axis becomes out-of-plane (new c) and x_axis becomes new a."""
    if z_axis not in _AXIS_INDEX or x_axis not in _AXIS_INDEX:
        raise ValueError("z_axis and x_axis must each be one of 'a', 'b', 'c'")
    if z_axis == x_axis:
        raise ValueError("z_axis and x_axis must be different")

    z_i = _AXIS_INDEX[z_axis]
    x_i = _AXIS_INDEX[x_axis]
    y_i = 3 - z_i - x_i
    perm = [x_i, y_i, z_i]

    old_matrix = structure.lattice.matrix
    new_matrix = old_matrix[perm]
    frac_coords = np.asarray(structure.frac_coords)
    cart = frac_coords @ old_matrix
    new_frac = cart @ np.linalg.inv(new_matrix)

    species = [_dominant_species(site) for site in structure.sites]
    return pgStructure.Structure(
        Lattice(new_matrix), species, new_frac, site_properties=structure.site_properties
    )


def _lab_rotation_matrix(lattice_matrix):
    """Orthonormal rotation mapping lattice c -> lab +Z and lattice a -> lab +X."""
    matrix = np.asarray(lattice_matrix, dtype=float)
    a_vec = matrix[0]
    c_vec = matrix[2]
    ez = c_vec / np.linalg.norm(c_vec)
    ax = a_vec - np.dot(a_vec, ez) * ez
    norm_ax = np.linalg.norm(ax)
    if norm_ax < 1e-10:
        b_vec = matrix[1]
        ax = b_vec - np.dot(b_vec, ez) * ez
        norm_ax = np.linalg.norm(ax)
    if norm_ax < 1e-10:
        return np.eye(3)
    ex = ax / norm_ax
    ey = np.cross(ez, ex)
    norm_ey = np.linalg.norm(ey)
    if norm_ey < 1e-10:
        return np.eye(3)
    ey = ey / norm_ey
    return np.vstack([ex, ey, ez])


def orient_structure_to_lab(structure):
    """Rigid-body rotation: new c along lab +Z, new a along lab +X."""
    matrix = np.asarray(structure.lattice.matrix, dtype=float)
    rotation = _lab_rotation_matrix(matrix)
    cart = np.asarray(structure.cart_coords, dtype=float)
    cart_lab = (rotation @ cart.T).T
    new_matrix = matrix @ rotation.T
    new_frac = cart_lab @ np.linalg.inv(new_matrix)
    species = [_dominant_species(site) for site in structure.sites]
    return pgStructure.Structure(
        Lattice(new_matrix), species, new_frac, site_properties=structure.site_properties
    )


def get_mount_matrix(structure):
    """Return 3x3 rotation mapping oriented crystal axes to the laboratory frame."""
    return _lab_rotation_matrix(np.asarray(structure.lattice.matrix, dtype=float))


def orient_crystal_axes(structure, z_axis, x_axis):
    """Permute CIF axes then rotate the crystal into the fixed laboratory frame."""
    return orient_structure_to_lab(permute_crystal_axes(structure, z_axis, x_axis))


def in_plane_axis_options(z_axis):
    return [axis for axis in _ALL_AXES if axis != z_axis]
