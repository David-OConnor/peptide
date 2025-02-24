//! This module contains code for calculating atom coordinates from dihedral angles.
//! It solves the *forward kinematics problem*.

use std::f64::consts::TAU;

use lin_alg::f64::{self, Quaternion, Vec3};

use crate::{
    bond_vecs::*,
    types::{BackboneCoords, Residue},
};

/// Calculate the dihedral angle between 4 atoms.
fn calc_dihedral_angle(bond_middle: Vec3, bond_adjacent1: Vec3, bond_adjacent2: Vec3) -> f64 {
    // Project the next and previous bonds onto the plane that has this bond as its normal.
    // Re-normalize after projecting.
    let bond1_on_plane = bond_adjacent1.project_to_plane(bond_middle).to_normalized();
    let bond2_on_plane = bond_adjacent2.project_to_plane(bond_middle).to_normalized();

    // Not sure why we need to offset by 𝜏/2 here, but it seems to be the case
    let result = bond1_on_plane.dot(bond2_on_plane).acos() + TAU / 2.;

    // The dot product approach to angles between vectors only covers half of possible
    // rotations; use a determinant of the 3 vectors as matrix columns to determine if we
    // need to modify to be on the second half.
    let det = lin_alg::f64::det_from_cols(bond1_on_plane, bond2_on_plane, bond_middle);

    // todo: Exception if vecs are the same??
    if det < 0. {
        result
    } else {
        TAU - result
    }
}

/// Calculate the orientation, as a quaternion, and position, as a vector, of an atom, given the orientation of a
/// previous atom, and the bond angle. `bond_to_prev_local` is the vector representing the bond to
/// this atom from the previous atom's orientation. `bond_to_prev_local`, `bond_to_next_local`, and
/// `bond_to_this_local` are in the atom's local coordinate space; not worldspace. They are all unit vectors.
///
/// This is solving an iteration of the *forward kinematics problem*.
pub fn find_atom_placement(
    or_prev: Quaternion,
    bond_to_prev_local: Vec3,
    bond_to_next_local: Vec3,
    dihedral_angle: f64,
    posit_prev: Vec3,
    posit_2_back: Vec3,
    bond_to_this_local: Vec3,
    bond_to_this_len: f64,
) -> (Vec3, Quaternion) {
    let prev_bond_world = posit_prev - posit_2_back;

    // Find the position; this is passed directly to the output, and isn't used for further
    // calcualtions within this function.
    let position = posit_prev + or_prev.rotate_vec(bond_to_this_local) * bond_to_this_len;

    // #1: Align the prev atom's bond vector to world space based on the prev atom's orientation.
    let bond_to_this_worldspace = or_prev.rotate_vec(bond_to_this_local);

    // #2: Find the rotation quaternion that aligns the (inverse of) the local-space bond to
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

    let angle_dif = dihedral_angle - dihedral_angle_current;

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

impl Residue {
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

        let (h_cα, h_cα_orientation) = find_atom_placement(
            cα_orientation,
            H_CALPHA_BOND,
            Vec3::new(0., 1., 0.), // arbitrary, since H doesn't continue the chain
            0.,                    // arbitrary
            cα,
            n_pos,
            unsafe { CALPHA_H_BOND },
            LEN_CALPHA_H,
        );

        let (h_n, h_n_orientation) = find_atom_placement(
            n_orientation,
            H_N_BOND,
            Vec3::new(0., 1., 0.), // arbitrary, since H doesn't continue the chain
            0.,                    // arbitrary
            n_pos,
            // todo: n_next, or n_pos?
            // cp,
            prev_cp_pos, // todo: prev cp, or cp?
            unsafe { N_H_BOND },
            LEN_N_H,
        );

        BackboneCoords {
            cα,
            cp,
            n_next,
            o,
            h_cα,
            h_n,
            cα_orientation,
            cp_orientation,
            n_next_orientation,
            o_orientation,
            h_cα_orientation,
            h_n_orientation,
        }
    }
}
