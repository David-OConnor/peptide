//! This module contains code for calculating atom coordinates from dihedral angles.
//! It solves the *forward kinematics problem*.
//!
//! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2810841/

use core::f64::consts::TAU;

use crate::{
    chem_definitions::{AminoAcidType, AtomType, BackboneRole},
    lin_alg::{self, Quaternion, Vec3},
    sidechain::{self, Sidechain},
};

// Double bond len of C' to N.
pub const LEN_CP_N: f64 = 1.33; // angstrom
pub const LEN_N_CALPHA: f64 = 1.46; // angstrom
pub const LEN_CALPHA_CP: f64 = 1.53; // angstrom

pub const LEN_CP_O: f64 = 1.2; // angstrom // todo!

// Ideal bond angles. There are an approximation; from averages. Consider replacing with something
// more robust later. All angles are in radians. We use degrees with math to match common sources.
// R indicates the side chain.
const BOND_ANGLE_N_CALPHA_CP: f64 = 121.7 * TAU / 360.; // This is to Calpha and C'.
                                                        // todo: n hydrogen
                                                        // Bond from the Calpha atom
const BOND_ANGLE_CALPHA_N_R: f64 = 110.6 * TAU / 360.;
const BOND_ANGLE_CALPHA_R_CP: f64 = 110.6 * TAU / 360.;
const BOND_ANGLE_CALPHA_N_CP: f64 = 111.0 * TAU / 360.;
// todo: calpha hydrogen

// Bonds from the C' atom
// Note that these bonds add up to exactly 360, so these must be along
// a circle (in the same plane).
const BOND_ANGLE_CP_CALPHA_O: f64 = 120.1 * TAU / 360.;
const BOND_ANGLE_CP_CALPHA_N: f64 = 117.2 * TAU / 360.;
const BOND_ANGLE_CP_O_N: f64 = 122.7 * TAU / 360.;
// todo: c' hydrogen

// An arbitrary vector that anchors the others
const INIT_BOND_VEC: Vec3 = Vec3 {
    x: 1.,
    y: 0.,
    z: 0.,
};

// These bonds are unit vecs populated by init_local_bond_vecs.

// We use `static mut` here instead of constant, since we need non-const fns (like sin and cos, and
// the linear algebra operations that operate on them) in their construction.
const CALPHA_CP_BOND: Vec3 = INIT_BOND_VEC;

pub static mut CALPHA_N_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub static mut CALPHA_R_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

const CP_N_BOND: Vec3 = INIT_BOND_VEC;

static mut CP_CALPHA_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

// The O bond on CP turns out to be [-0.5402403204776551, 0, 0.841510781945306], given our calculated
// anchors for the N and Calpha bonds on it.
static mut CP_O_BOND: Vec3 = Vec3 {
    // todo: These are arbitrary values
    x: 0.,
    y: 0.,
    z: 0.,
};

const N_CALPHA_BOND: Vec3 = INIT_BOND_VEC;

static mut N_CP_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

const O_CP_BOND: Vec3 = INIT_BOND_VEC;

/// Calculate local bond vectors based on relative angles, and arbitrary constants.
/// The absolute bonds used are arbitrary; their positions relative to each other are
/// defined by the bond angles.
/// As an arbitrary convention, we'll make the first vector the one to the next atom
/// in the chain, and the second to the previous. The third is for C'oxygen, or Cα side chain.
pub fn init_local_bond_vecs() {
    // Calculate (arbitrary) vectors normal to the anchor vectors for each atom.
    // Find the second bond vector by rotating the first around this by the angle
    // between the two.
    // todo: Given we're anchoring the initial vecs to a specific vector, we can
    // todo skip this and use a known orthonormal vec to it like 0, 1, 0.
    // let normal_cα = Vec3::new(0., 1., 0.).cross(CALPHA_CP_BOND);
    let normal_cα = Vec3::new(0., 1., 0.);
    let normal_cp = Vec3::new(0., 1., 0.);
    let normal_n = Vec3::new(0., 1., 0.);

    //
    let rotation_cα = Quaternion::from_axis_angle(normal_cα, BOND_ANGLE_CALPHA_N_CP);
    let rotation_cp = Quaternion::from_axis_angle(normal_cp, BOND_ANGLE_CP_CALPHA_N);
    let rotation_n = Quaternion::from_axis_angle(normal_n, BOND_ANGLE_N_CALPHA_CP);

    // todo: Experimenting with bonds to R.
    // let rotation_n = Quaternion::from_axis_angle(normal_cα, BOND_ANGLE_CALPHA_R_CP);

    unsafe {
        CALPHA_N_BOND = rotation_cα.rotate_vec(CALPHA_CP_BOND);
        CP_CALPHA_BOND = rotation_cp.rotate_vec(CP_N_BOND);
        N_CP_BOND = rotation_n.rotate_vec(N_CALPHA_BOND);

        // CALPHA_R_BOND = rotation_cα.rotate_vec(CALPHA_R_BOND); // todo: QC etc.

        // println!("cp_n: [{}, {}, {}]", CP_N_BOND.x, CP_N_BOND.y, CP_N_BOND.z);
        // println!("cp_calpha: [{}, {}, {}]", CP_CALPHA_BOND.x, CP_CALPHA_BOND.y, CP_CALPHA_BOND.z);

        // todo: To find the 3rd (and later 4th, ie hydrogen-to-atom) bonds, we we
        // todo taking an iterative approach based on that above. Long-temr, you should
        // todo be able to calculate one of 2 valid choices (for carbon).

        // The CP ON angle must be within this (radians) to stop the process.
        // The other 2 angles should be ~exactly.
        const EPS: f64 = 0.01;

        // This approach performs a number of calculations at runtime, at init only. Correct
        // solution for now, since axis=0 turns out to be correct, given the bond from the CP atom
        // are in a circle.
        for axis in 0..10 {
            let r1 = Quaternion::from_axis_angle(CP_CALPHA_BOND, axis as f64 / 20.);
            let normal = r1.rotate_vec(normal_cp);

            // Set exactly to one vector; find the other through iterating the normal angles.
            let r2 = Quaternion::from_axis_angle(normal, BOND_ANGLE_CP_CALPHA_O);

            CP_O_BOND = r2.rotate_vec(CP_CALPHA_BOND);

            // todo: Quick/sloppy. Also, exits on satisfying CP_O_bond; not this. For now,
            // todo appears to produce the correct result.
            CALPHA_R_BOND = r2.rotate_vec(CALPHA_N_BOND);

            // let angle_calpha_o = (CP_CALPHA_BOND.dot(CP_O_BOND)).acos();
            // let angle_calpha_n = (CP_CALPHA_BOND.dot(CP_N_BOND)).acos();
            let angle_n_o = (CP_N_BOND.dot(CP_O_BOND)).acos();

            // Of note, our first value (axis = 0) seems to be the answer here?!
            if (angle_n_o - BOND_ANGLE_CP_O_N).abs() < EPS {
                break;
            }
            //
            // println!("cp_o: [{}, {}, {}]", CP_O_BOND.x, CP_O_BOND.y, CP_O_BOND.z);
            // // println!("angle calpha O: {angle_calpha_o}");
            // // println!("angle calpha N: {angle_calpha_n}");
            // println!("angle N O: {angle_n_o}");
        }
    }

    // Find vectors from C' to O, and Cα to R, given the previous 2 bonds for each.
}

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub residues: Vec<Residue>,
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
    let det = lin_alg::det_from_cols(bond1_on_plane, bond2_on_plane, bond_middle);

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
}

impl Residue {
    pub fn new(ω: f64, φ: f64, ψ: f64, sidechain: Sidechain) -> Self {
        Self {
            ω, φ, ψ, sidechain
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
