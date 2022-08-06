//! This module contains code for calculating atom coordinates.

// todo: QC quaternion error creep, and re-normalize A/R

use crate::{
    chem_definitions::{AminoAcidType, BackboneRole},
    lin_alg::{Quaternion, Vec3},
};

// Double bond len of C' to N.
// todo: 1.089000
pub const LEN_CP_N: f64 = 1.32; // angstrom
pub const LEN_N_CALPHA: f64 = 1.; // angstrom // todo
pub const LEN_CALPHA_CP: f64 = 1.; // angstrom // todo

// todo: These may change under diff circumstances, and be diff for the diff atoms
// and diff things on the atom.
const BACKBONE_BOND_ANGLES: f64 = 1.911; // radians

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub residues: Vec<Residue>,
}

#[derive(Debug)]
/// Describes the sequence of atoms that make up a protein backbone, with worldspace coordinates.
/// this is what is needed for the render, and spacial manipulations.
pub struct ProteinCoords {
    pub atoms_backbone: Vec<AtomCoords>,
}

impl ProteinCoords {
    pub fn from_descrip(descrip: &ProteinDescription) -> Self {
        let mut backbone = Vec::new();

        let mut id = 0;

        // N-terminus nitrogen, at the *start* of our chain.
        let starting_n = AtomCoords {
            residue_id: id,
            role: BackboneRole::N,
            position: Vec3::new(0., 0., 0.),
            orientation: Quaternion::new_identity(),
        };

        let mut prev_position = starting_n.position;
        // todo: This may need adjustment. to match physical reality.
        let mut prev_2_position = Vec3::new(1., 0., 0.);
        let mut prev_orientation = starting_n.orientation;

        backbone.push(starting_n);
        id += 1;

        for aa in &descrip.residues {
            let aa_coords =
                aa.backbone_cart_coords(prev_position, prev_orientation, prev_2_position);

            backbone.push(AtomCoords {
                residue_id: id,
                role: BackboneRole::Cα,
                position: aa_coords.cα,
                orientation: aa_coords.cα_orientation,
            });
            id += 1;

            backbone.push(AtomCoords {
                residue_id: id,
                role: BackboneRole::Cp,
                position: aa_coords.cp,
                orientation: aa_coords.cp_orientation,
            });
            id += 1;

            let n = AtomCoords {
                residue_id: id,
                role: BackboneRole::N,
                position: aa_coords.n_next,
                orientation: aa_coords.n_next_orientation,
            };

            prev_position = n.position;
            prev_orientation = n.orientation;

            backbone.push(n);
            id += 1;
        }

        Self {
            atoms_backbone: backbone,
        }
    }
}

/// Location of an atom, in the worldspace coordinate system.
#[derive(Debug)]
pub struct AtomCoords {
    /// id of the Amino Acid this atom is part of
    pub residue_id: usize, // todo: Do we want this id, or use an index?
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

/// Calculate the orientation, as a quaternion, of a backbone atom, given the orientation of the
/// previous atom, and the bond angle. `bond_angle` is the vector representing the bond to
/// this atom from the previous atom's orientation. `bond_prev` and `bond_next` are in the atom's
/// coordinates; not worldspace.
pub fn find_backbone_atom_orientation(
    q_prev: Quaternion,
    bond_to_prev: Vec3,         // Local space
    bond_to_next: Vec3,         // Local space
    prev_bond_worldspace: Vec3, // World space
    dihedral_angle: f64,
) -> Quaternion {
    // #1: Align the prev atom's bond vector to world space based on the prev atom's orientation.
    let bond_to_this_worldspace = q_prev.rotate_vec(bond_to_next);

    // #2: Find the rotation quaternion that aligns the (inverse of) the local-space bond to
    // the prev atom with the world-space "to" bond of the previous atom. This is also the
    // orientation of our atom, without applying the dihedral angle.
    let bond_alignment_rotation =
        Quaternion::from_unit_vecs(bond_to_prev * -1., bond_to_this_worldspace);

    // #3: Rotate the orientation along the dihedral angle.
    // let dihedral_rotation = Quaternion::from_axis_angle(bond_to_this_worldspace, dihedral_angle);
    // let orientation = dihedral_rotation * bond_alignment_rotation;
    let orientation = bond_alignment_rotation;

    // #4: Adjust the dihedral angle to be in reference to the previous 2 atoms, per the convention.
    let next_bond_worldspace = orientation.rotate_vec(bond_to_next);

    // Change the basis so that one axis (we choose z) is aligned with the bond we're rotating.
    // Then remove the z component of the two flanking vectors. Their angle is the dihedral angle.
    // let basis_change = Quaternion::from_unit_vecs(Vec3::new(0., 0., 1.), bond_to_this_worldspace); // todo: Negate?
    //
    // let next_bond_2d = basis_change.rotate_vec(bond_to_next_worldspace);
    // let prev_bond_2d = basis_change.rotate_vec(vec_prev_2_atoms);
    //

    // Project the next and previous bonds onto the plane that has this
    // bond as its normal.
    let prev_bond_on_plane = prev_bond_worldspace.project_to_plane(bond_to_this_worldspace);
    // todo: Do we want to/need to reverse the plane axis to make the angles work correctly? Consider an
    // todo experiment in python.
    let next_bond_on_plane = next_bond_worldspace.project_to_plane(bond_to_this_worldspace * -1.);

    // // todo: There's a lot of room for sign, order etc errors.
    let dihedral_angle_current = (next_bond_on_plane.dot(prev_bond_on_plane)).acos();

    let dihedral_rotation = Quaternion::from_axis_angle(
        bond_to_this_worldspace,
        dihedral_angle - dihedral_angle_current,
    );

    dihedral_rotation * orientation
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
    /// `prev_cp_pos` is usd to anchor the dihedral angle properly, since it's defined by planes
    ///  of 3 atoms.
    pub fn backbone_cart_coords(
        &self,
        n_pos: Vec3,
        n_orientation: Quaternion,
        prev_cp_pos: Vec3,
    ) -> BackboneCoordsAa {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.
        let bond_to_prev = Vec3::new(-1., 1., 1.).to_normalized();
        let bond_to_next = Vec3::new(1., 1., -1.).to_normalized();

        // The other bonds, ie 2xH for C' and N, and H + side chain for Cα.
        let _bond_a = Vec3::new(1., -1., 1.).to_normalized();
        let _bond_b = Vec3::new(-1., -1., -1.).to_normalized();

        // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
        // (world space)
        let bond_prev_worldspace = n_pos - prev_cp_pos; // todo: Direction?

        // todo: You must anchor the dihedral angle in terms of the plane formed by the atom pairs
        // todo on each side of the angle; that's how it's defined.
        let cα_orientation = find_backbone_atom_orientation(
            n_orientation,
            bond_to_prev,
            bond_to_next,
            bond_prev_worldspace,
            self.φ,
        );
        let cp_orientation = find_backbone_atom_orientation(
            cα_orientation,
            bond_to_prev,
            bond_to_next,
            bond_prev_worldspace,
            self.ψ,
        );
        let n_next_orientation = find_backbone_atom_orientation(
            cp_orientation,
            bond_to_prev,
            bond_to_next,
            bond_prev_worldspace,
            self.ω,
        );

        // Calculate the orientation of each atom from the orientation, the angle of the
        // bond from the previous atom, and the dihedral angle between the two.

        // Set up atom-oriented vectors encoding the angle and len of each bond.
        let bond_n_cα = bond_to_next * LEN_N_CALPHA;
        let bond_cα_cp = bond_to_next * LEN_CALPHA_CP;
        let bond_cp_n = bond_to_next * LEN_CP_N;

        // Calculate the position of each atom from the position, orientation and bond angle
        // from the previous atom.
        // todo: Start by describing each term in terms of its dependencies.
        let cα = n_pos + n_orientation.rotate_vec(bond_n_cα);
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

/// Used to represent one atom in a system built of atoms.
pub struct _ZmatrixItem {
    atomic_number: u8,
    bond_len: f64,       // Angstrom
    bond_angle: f64,     // radians
    dihedral_angle: f64, // radians
}
