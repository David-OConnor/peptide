//! This module contains code for calculating atom coordinates from dihedral angles.
//! It solves the *forward kinematics problem*.
//!
//! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2810841/

use std::{f64::consts::TAU, fmt};

use crate::{
    bond_vecs::*,
    chem_definitions::{AminoAcidType, AtomType, BackboneRole},
    sidechain::{self, Sidechain},
};

use lin_alg2::f64::{self, Quaternion, Vec3};

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub name: String,
    pub pdb_ident: String,
    pub residues: Vec<Residue>,
}

impl fmt::Display for ProteinDescription {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Name: {} PDB ident: {}, Residues:\n\n", self.name, self.pdb_ident);

        for (i, residue) in self.residues.iter().enumerate() {
            write!(f, "Id {} - {}\n", i + 1, residue);
        }

        Ok(())
    }
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
pub struct BackboneCoords {
    /// Cα
    pub cα: Vec3,
    /// Carbon' atom bound to C_alpha
    pub cp: Vec3,
    /// Nitrogen atom of the next module.
    pub n_next: Vec3,
    /// Oxygen atom bonded to c'
    pub o: Vec3,
    pub cα_orientation: Quaternion,
    pub cp_orientation: Quaternion,
    pub n_next_orientation: Quaternion,
    pub o_orientation: Quaternion,
}

/// Calculate the dihedral angle between 4 atoms.
fn calc_dihedral_angle(bond_middle: Vec3, bond_adjacent1: Vec3, bond_adjacent2: Vec3) -> f64 {
    // Project the next and previous bonds onto the plane that has this bond as its normal.
    // Re-normalize after projecting.
    let bond1_on_plane = bond_adjacent1.project_to_plane(bond_middle).to_normalized();
    let bond2_on_plane = bond_adjacent2.project_to_plane(bond_middle).to_normalized();

    // Not sure why we need to offset by 𝜏/2 here, but it seems to be the case
    let result = (bond1_on_plane.dot(bond2_on_plane)).acos() + TAU / 2.;

    // The dot product approach to angles between vectors only covers half of possible
    // rotations; use a determinant of the 3 vectors as matrix columns to determine if we
    // need to modify to be on the second half.
    let det = lin_alg2::f64::det_from_cols(bond1_on_plane, bond2_on_plane, bond_middle);

    // todo: Exception if vecs are the same??
    if det < 0. {
        result
    } else {
        TAU - result
    }
}

/// Calculate the orientation, as a quaternion, and position, as a vector, of an atom, given the orientation of a
/// previous atom, and the bond angle. `bond_angle` is the vector representing the bond to
/// this atom from the previous atom's orientation. `bond_prev` and `bond_next` are in the atom's
/// coordinates; not worldspace. This is solving an iteration of the *forward kinematics problem*.
pub fn find_atom_placement(
    o_prev: Quaternion,
    bond_to_prev_local: Vec3, // Local space
    bond_to_next_local: Vec3, // Local space
    dihedral_angle: f64,
    posit_prev: Vec3,
    posit_2_back: Vec3,
    bond_to_this_local: Vec3, // direction-only unit vec
    bond_to_this_len: f64,
) -> (Vec3, Quaternion) {
    let prev_bond_world = posit_prev - posit_2_back;

    // Find the position:
    let position = posit_prev + o_prev.rotate_vec(bond_to_this_local) * bond_to_this_len;

    // #1: Align the prev atom's bond vector to world space based on the prev atom's orientation.
    let bond_to_this_worldspace = o_prev.rotate_vec(bond_to_next_local);

    // #2: Find the rotation quaternion that aligns the (inverse of) the local(world?)-space bond to
    // the prev atom with the world-space "to" bond of the previous atom. This is also the
    // orientation of our atom, without applying the dihedral angle.
    let bond_alignment_rotation =
        Quaternion::from_unit_vecs(bond_to_prev_local * -1., bond_to_this_worldspace);

    // #3: Rotate the orientation around the dihedral angle. We must do this so that our
    // dihedral angle is in relation to the previous and next bonds.
    // Adjust the dihedral angle to be in reference to the previous 2 atoms, per the convention.
    let next_bond_worldspace = bond_alignment_rotation.rotate_vec(bond_to_next_local);

    let dihedral_angle_current = calc_dihedral_angle(
        bond_to_this_worldspace,
        prev_bond_world,
        next_bond_worldspace,
    );

    let mut angle_dif = dihedral_angle - dihedral_angle_current;

    let rotate_axis = bond_to_this_worldspace;

    let mut dihedral_rotation = Quaternion::from_axis_angle(rotate_axis, angle_dif);

    let dihedral_angle_current2 = calc_dihedral_angle(
        bond_to_this_worldspace,
        prev_bond_world,
        (dihedral_rotation * bond_alignment_rotation).rotate_vec(bond_to_next_local),
    );

    // todo: Quick and dirty approach here. You can perhaps come up with something
    // todo more efficient, ie in one shot.
    if (dihedral_angle - dihedral_angle_current2).abs() > 0.0001 {
        dihedral_rotation = Quaternion::from_axis_angle(rotate_axis, -angle_dif + TAU);
    }

    (position, dihedral_rotation * bond_alignment_rotation)
}

/// An amino acid in a protein structure, including all dihedral angles required to determine
/// the conformation. Includes backbone and side chain dihedral angles. Doesn't store coordinates,
/// but coordinates can be generated using forward kinematics from the angles.
#[derive(Debug)]
pub struct Residue {
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
    /// Contains the χ angles that define t
    pub sidechain: Sidechain,
    pub dipole: Vec3,
}

impl fmt::Display for Residue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\nω: {:.2}τ, φ: {:.2}τ ψ: {:.2}τ\n",
            self.sidechain,
            self.ω / TAU,
            self.φ / TAU,
            self.ψ / TAU
        )
    }
}

impl Residue {
    pub fn new(ω: f64, φ: f64, ψ: f64, sidechain: Sidechain) -> Self {
        Self {
            ω,
            φ,
            ψ,
            sidechain,
            dipole: Vec3::new_zero(),
        }
    }

    /// Generate cartesian coordinates of points from diahedral angles and bond lengths. Starts with
    /// N, and ends with C'. This is solving the *forward kinematics problem*.
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
    ) -> BackboneCoords {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.

        let (cα, cα_orientation) = find_atom_placement(
            n_orientation,
            unsafe { CALPHA_N_BOND },
            CALPHA_CP_BOND,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.φ,
            n_pos,
            prev_cp_pos,
            N_CALPHA_BOND,
            LEN_N_CALPHA,
        );

        let (cp, cp_orientation) = find_atom_placement(
            cα_orientation,
            unsafe { CP_CALPHA_BOND },
            CP_N_BOND,
            self.ψ,
            cα,
            n_pos,
            CALPHA_CP_BOND,
            LEN_CALPHA_CP,
        );

        let (n_next, n_next_orientation) = find_atom_placement(
            cp_orientation,
            unsafe { N_CP_BOND },
            N_CALPHA_BOND,
            self.ω,
            cp,
            cα,
            CP_N_BOND,
            LEN_CP_N,
        );

        // Oxygen isn't part of the backbone chain; it's attached to C' with a double-bond.
        // The vectors we use are fudged; we just need to keep the orientation conistent relative
        // tod c'.
        let (o, o_orientation) = find_atom_placement(
            cp_orientation,
            O_CP_BOND,
            Vec3::new(0., 1., 0.), // arbitrary, since O doesn't continue the chain
            0.,                    // arbitrary
            cp,
            cα,
            unsafe { CP_O_BOND },
            LEN_CP_O,
        );

        BackboneCoords {
            cα,
            cp,
            n_next,
            o,
            cα_orientation,
            cp_orientation,
            n_next_orientation,
            o_orientation,
        }
    }
}

/// Describes a water molecule. These aren't directly part of a protein, but may play a role in its
/// folding, among other potential roles.
#[derive(Debug)]
// todo: Consider if you want this to be a struct, a const of some other struct etc.
pub struct WaterMolecule {
    /// Worldspace coordinates of the O atom.
    position_o_world: Vec3,
    /// Using the same orientation ref as protein atoms.
    orientation: Quaternion,
}
