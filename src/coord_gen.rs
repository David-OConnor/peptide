//! This module contains code for calculating atom coordinates.

use crate::{
    chem_definitions::{AminoAcidType, BackboneRole},
    lin_alg::{Quaternion, Vec3},
};

// Double bond len of C' to N.
// todo: 1.089000
const LEN_CP_N: f64 = 1.32; // angstrom
const LEN_N_CALPHA: f64 = 1.; // angstrom // todo
const LEN_CALPHA_CP: f64 = 1.; // angstrom // todo

// todo: These may change under diff circumstances, and be diff for the diff atoms
// and diff things on the atom.
const BACKBONE_BOND_ANGLES: f64 = 1.911; // radians

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub aas: Vec<Residue>,
}

#[derive(Debug)]
/// Describes the sequence of atoms that make up a protein backbone, with worldspace coordinates.
/// this is what is needed for the render, and spacial manipulations.
pub struct ProteinCoords {
    pub atoms_backbone: Vec<AtomCoords>,
}

impl ProteinCoords {
    pub fn from_descrip(descrip: &ProteinDescription) -> Self {
        let mut atoms_backbone = Vec::new();

        let mut id = 0;

        // N-terminus nitrogen, at the *start* of our chain.
        let starting_n = AtomCoords {
            aa_id: id,
            role: BackboneRole::N,
            position: Vec3::new(0., 0., 0.),
            orientation: Quaternion::new_identity(),
        };

        let mut prev_n = starting_n.clone();

        atoms_backbone.push(starting_n);
        id += 1;

        for aa in &descrip.aas {
            let aa_coords = aa.backbone_cart_coords(prev_n.position, prev_n.orientation);

            atoms_backbone.push(AtomCoords {
                aa_id: id,
                role: BackboneRole::Cα,
                position: aa_coords.cα,
                orientation: aa_coords.cα_orientation,
            });
            id += 1;

            atoms_backbone.push(AtomCoords {
                aa_id: id,
                role: BackboneRole::Cp,
                position: aa_coords.cp,
                orientation: aa_coords.cp_orientation,
            });
            id += 1;

            let n = AtomCoords {
                aa_id: id,
                role: BackboneRole::N,
                // atom: AtomType::N,
                position: aa_coords.n_next,
                orientation: aa_coords.n_next_orientation,
            };

            prev_n = n.clone();

            atoms_backbone.push(n);
            id += 1;
        }

        Self { atoms_backbone }
    }
}

/// Location of an atom, in the worldspace coordinate system.
#[derive(Clone, Debug)]
pub struct AtomCoords {
    /// id of the Amino Acid this atom is part of
    pub aa_id: usize, // todo: Do we want this id, or use an index?
    pub role: BackboneRole,
    pub position: Vec3,
    pub orientation: Quaternion,
}

/// Holds backbone atom coordinates and orientations, relative to the alpha carbon, for
/// a single amino acid. Generated from an AA's bond angles.
///
/// We will use the convention of the chain starting at NH2 (amine) end (N terminus), and ending in
/// COOH (carboxyl) (C terminus).
/// Flow is N -> C_alpha -> C'.
/// Coordinates are in relation to the Nitrogen molecule at the start of the segement, using the
/// directional convention described above. (todo: How do we orient axes?)
/// todo: To start, convention is the previous C' to this AA's starting A is on the positive X axis.
#[derive(Debug)]
pub struct BackboneCoordsAa {
    // ///Nitrogen atom bonded to C_alpha // defined as 0.
    // pub n: f64,
    /// Cα
    pub cα: Vec3,
    /// Carbon' atom bound to C_alpha
    pub cp: Vec3,
    /// Nitrogen atom of the next module.
    pub n_next: Vec3,
    pub cα_orientation: Quaternion,
    pub cp_orientation: Quaternion,
    pub n_next_orientation: Quaternion,
}

/// An amino acid in a protein structure, including position information.
#[derive(Debug)]
pub struct Residue {
    pub aa: AminoAcidType,
    /// Dihedral angle between C' and N
    /// Tor (Cα, C, N, Cα) is the ω torsion angle
    /// Assumed to be TAU/2 for most cases
    pub ω: f64,
    /// Dihedral angle between Cα and N.
    /// Tor (C, N, Cα, C) is the φ torsion angle
    pub φ: f64,
    /// Dihedral angle, between Cα and C'
    ///  Tor (N, Cα, C, N) is the ψ torsion angle
    pub ψ: f64,
    // todo: Include bond lengths here if they're not constant.
}

impl Residue {
    /// Generate cartesian coordinates of points from diahedral angles and bond lengths. Starts with
    /// N, and ends with C'.
    /// Accepts position, and orientation of the N atom that starts this segment.
    /// Also returns orientation of the N atom, for use when calculating coordinates for the next
    /// AA in the chain.
    pub fn backbone_cart_coords(&self, pos_n: Vec3, n_orientation: Quaternion) -> BackboneCoordsAa {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.
        let bond_prev = Vec3::new(-1., 1., 1.).to_normalized();
        let bond_next = Vec3::new(1., 1., -1.).to_normalized();

        // The other bonds, ie 2xH for C' and N, and H + side chain for Cα.
        let _bond_a = Vec3::new(1., -1., 1.).to_normalized();
        let _bond_b = Vec3::new(-1., -1., -1.).to_normalized();

        // todo: You must anchor the dihedral angle in terms of the plane formed by the atom pairs
        // todo on each side of the angle; that's how it's defined.
        let cα_orientation = find_backbone_atom_orientation(n_orientation, bond_prev, self.ψ);
        let cp_orientation = find_backbone_atom_orientation(cα_orientation, bond_prev, self.φ);
        let n_next_orientation = find_backbone_atom_orientation(cp_orientation, bond_prev, self.ω);

        // Orientations for the various atoms, Z up. todo ?

        // Calculate the orientation of each atom from the orientation, the angle of the
        // bond from the previous atom, and the dihedral angle between the two.

        // Set up atom-oriented vectors encoding the angle and len of each bond.
        let bond_n_cα = bond_next * LEN_N_CALPHA;
        let bond_cα_cp = bond_next * LEN_CALPHA_CP;
        let bond_cp_n = bond_next * LEN_CP_N;

        // Calculate the position of each atom from the position, orientation and bond angle
        // from the previous atom.
        // todo: Start by describing each term in terms of its dependencies.
        let cα = pos_n + n_orientation.rotate_vec(bond_n_cα);
        let cp = cα + cα_orientation.rotate_vec(bond_cα_cp);
        let n_next = cp + cp_orientation.rotate_vec(bond_cp_n);

        BackboneCoordsAa {
            cα,
            cp,
            n_next,
            cα_orientation,
            cp_orientation,
            n_next_orientation,
        }
    }
}

/// Calculate the orientation, as a quaternion, of a backbone atom, given the orientation of the
/// previous atom, and the bond angle. `bond_angle` is the vector representing the bond to
/// this atom from the previous atom's orientation.
pub fn find_backbone_atom_orientation(
    q_prev: Quaternion,
    bond_prev: Vec3,
    dihedral_angle: f64,
) -> Quaternion {
    // We split up determining the orientation of each backbone atom into 3 steps:

    // #1: Find the bond axis. This is the world-space vector of the backbone bond by rotating the
    // bond angle with the previous atom's orientation quaternion.
    let bond_axis = q_prev.rotate_vec(bond_prev);

    // todo: Rotate from prev bond to next bond.

    // #2: Calculate the orientation quaternion of the current atom so its bond to the previous atom
    // is aligned with the bond axis found above, and rotated per the dihedral angle.
    let current_atom_orientation = Quaternion::from_axis_angle(bond_axis, dihedral_angle);

    // #3: Rotate the current atom's orientation by the dihedral angle between the two atoms. (eg φ, or ψ).

    // todo: Instead of actual dihedral angles, we are using a somewhat arbitrary rotation. todo: Change
    // todo this later to be actual dihedral angle. To do so, you probably need to pass the plane formed.
    // todo: Consider changing this variable name to `bond_rotation` in the meanwhile.

    // by the prev 2 angles, and subtract it from the plane formed by this and the next.

    // current_atom_orientation * dihedral_angle

    current_atom_orientation
}

/// Used to represent one atom in a system built of atoms.
pub struct _ZmatrixItem {
    atomic_number: u8,
    bond_len: f64,       // Angstrom
    bond_angle: f64,     // radians
    dihedral_angle: f64, // radians
}
