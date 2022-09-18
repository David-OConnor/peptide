//! Code for generating bond vectors for an atom.

use std::f64::consts::TAU;

use lin_alg2::f64::{Quaternion, Vec3};

// Double bond len of C' to N.
pub const LEN_CP_N: f64 = 1.33; // angstrom
pub const LEN_N_CALPHA: f64 = 1.46; // angstrom
pub const LEN_CALPHA_CP: f64 = 1.53; // angstrom

pub const LEN_CP_O: f64 = 1.2; // angstrom // todo placeholder!
pub const LEN_N_H: f64 = 1.0; // angstrom // todo placeholder!

// Ideal bond angles. There are an approximation; from averages. Consider replacing with something
// more robust later. All angles are in radians. We use degrees with math to match common sources.
// R indicates the side chain.
pub const BOND_ANGLE_N_CALPHA_CP: f64 = 121.7 * TAU / 360.; // This is to Calpha and C'.
pub const BOND_ANGLE_N_CALPHA_H: f64 = 120. * TAU / 360.; // todo: Placeholder; not real data

// todo: n hydrogen
// Bond from the Calpha atom

pub const BOND_ANGLE_CALPHA_N_R: f64 = 110.6 * TAU / 360.;
pub const BOND_ANGLE_CALPHA_R_CP: f64 = 110.6 * TAU / 360.;
pub const BOND_ANGLE_CALPHA_N_CP: f64 = 111.0 * TAU / 360.;
// todo: calpha hydrogen

// Bonds from the C' atom
// Note that these bonds add up to exactly 360, so these must be along
// a circle (in the same plane).
pub const BOND_ANGLE_CP_CALPHA_O: f64 = 120.1 * TAU / 360.;
pub const BOND_ANGLE_CP_CALPHA_N: f64 = 117.2 * TAU / 360.;
pub const BOND_ANGLE_CP_O_N: f64 = 122.7 * TAU / 360.;
// todo: c' hydrogen

// An arbitrary vector that anchors the others
pub const ANCHOR_BOND_VEC: Vec3 = Vec3 {
    x: 1.,
    y: 0.,
    z: 0.,
};

// These bonds are unit vecs populated by init_local_bond_vecs.

// We use `static mut` here instead of constant, since we need non-const fns (like sin and cos, and
// the linear algebra operations that operate on them) in their construction.
pub const CALPHA_CP_BOND: Vec3 = ANCHOR_BOND_VEC;

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

pub const CP_N_BOND: Vec3 = ANCHOR_BOND_VEC;

pub static mut CP_CALPHA_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

// The O bond on CP turns out to be [-0.5402403204776551, 0, 0.841510781945306], given our calculated
// anchors for the N and Calpha bonds on it.
pub static mut CP_O_BOND: Vec3 = Vec3 {
    // todo: These are arbitrary values
    x: 0.,
    y: 0.,
    z: 0.,
};

pub const N_CALPHA_BOND: Vec3 = ANCHOR_BOND_VEC;

pub static mut N_CP_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub static mut N_H_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub const O_CP_BOND: Vec3 = ANCHOR_BOND_VEC;
pub const H_N_BOND: Vec3 = ANCHOR_BOND_VEC;

// todo: temp
// These are updated in `init_local_bond_vecs`.
pub static mut SIDECHAIN_BOND_TO_PREV: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub static mut SIDECHAIN_BOND_OUT1: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub static mut SIDECHAIN_BOND_OUT2: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

/// Find the next vector in a plane, given a starting vector, and rotation angle.
fn find_vec_in_plane(start: Vec3, angle: f64, plane_norm: Vec3) -> Vec3 {
    let rotation = Quaternion::from_axis_angle(plane_norm, angle);
    rotation.rotate_vec(start)
}

/// Calculate local bond vectors based on relative angles, and arbitrary constants.
/// The absolute bonds used are arbitrary; their positions relative to each other are
/// defined by the bond angles.
/// As an arbitrary convention, we'll make the first vector the one to the next atom
/// in the chain, and the second to the previous. The third is for C'oxygen, or Cα side chain.
pub fn init_local_bond_vecs() {
    // Calculate (arbitrary) vectors normal to the anchor vectors for each atom.
    // Find the second bond vector by rotating the first around this by the angle
    // between the two.
    // Given we're anchoring the initial vecs to a specific vector
    // ANCHOR_BOND_VEC = (1, 0, 0), we can
    // skip this and use a known orthonormal vec to it like 0, 1, 0.

    // let normal_cα = Vec3::new(0., 1., 0.).cross(CALPHA_CP_BOND);
    let normal_plane = Vec3::new(0., 1., 0.);

    unsafe {
        CALPHA_N_BOND = find_vec_in_plane(CALPHA_CP_BOND, BOND_ANGLE_CALPHA_N_CP, normal_plane);
        CP_CALPHA_BOND = find_vec_in_plane(CP_N_BOND, BOND_ANGLE_CP_CALPHA_N, normal_plane);
        N_CP_BOND = find_vec_in_plane(N_CALPHA_BOND, BOND_ANGLE_N_CALPHA_CP, normal_plane);

        // todo: We don't have actual H measurements; using a dummy angle for now.
        N_H_BOND = find_vec_in_plane(N_CALPHA_BOND, TAU * 2. / 3., normal_plane);

        SIDECHAIN_BOND_OUT1 = ANCHOR_BOND_VEC;
        SIDECHAIN_BOND_OUT2 = find_vec_in_plane(ANCHOR_BOND_VEC, TAU / 3., normal_plane);
        SIDECHAIN_BOND_TO_PREV = find_vec_in_plane(ANCHOR_BOND_VEC, TAU * 2. / 3., normal_plane);

        // todo: To find the 3rd (and later 4th, ie hydrogen-to-atom) bonds, we we
        // todo taking an iterative approach based on that above. Long-term, you should
        // todo be able to calculate one of 2 valid choices (for carbon).

        // The CP ON angle must be within this (radians) to stop the process.
        // The other 2 angles should be ~exactly.
        const EPS: f64 = 0.01;

        // This approach performs a number of calculations at runtime, at init only. Correct
        // solution for now, since axis=0 turns out to be correct, given the bonds from the CP atom
        // are in a circle.
        for axis in 0..10 {
            let r1 = Quaternion::from_axis_angle(CP_CALPHA_BOND, axis as f64 / 20.);
            let normal = r1.rotate_vec(normal_plane);

            // Set exactly to one vector; find the other through iterating the normal angles.
            let r2 = Quaternion::from_axis_angle(normal, BOND_ANGLE_CP_CALPHA_O);

            CP_O_BOND = r2.rotate_vec(CP_CALPHA_BOND);

            // todo: Quick/sloppy. Also, exits on satisfying CP_O_bond; not this. For now,
            // todo appears to produce the correct result.
            CALPHA_R_BOND = r2.rotate_vec(CALPHA_N_BOND);

            // let angle_calpha_o = (CP_CALPHA_BOND.dot(CP_O_BOND)).acos();
            // let angle_calpha_n = (CP_CALPHA_BOND.dot(CP_N_BOND)).acos();
            let angle_n_o = (CP_N_BOND.dot(CP_O_BOND)).acos();

            // SIDECHAIN_BOND_OUT1 = CALPHA_R_BOND;

            // Of note, our first value (axis = 0) seems to be the answer here?!
            if (angle_n_o - BOND_ANGLE_CP_O_N).abs() < EPS {
                break;
            }
        }
    }

    // Find vectors from C' to O, and Cα to R, given the previous 2 bonds for each.
}
