//! Code for generating bond vectors for an atom.

use std::f64::consts::TAU;

use lin_alg2::f64::{Quaternion, Vec3};

// Double bond len of C' to N.
pub const LEN_CP_N: f64 = 1.33; // angstrom
pub const LEN_N_CALPHA: f64 = 1.46; // angstrom
pub const LEN_CALPHA_CP: f64 = 1.53; // angstrom

pub const LEN_CP_O: f64 = 1.2; // angstrom // todo!

// Ideal bond angles. There are an approximation; from averages. Consider replacing with something
// more robust later. All angles are in radians. We use degrees with math to match common sources.
// R indicates the side chain.
pub const BOND_ANGLE_N_CALPHA_CP: f64 = 121.7 * TAU / 360.; // This is to Calpha and C'.
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
pub const INIT_BOND_VEC: Vec3 = Vec3 {
    x: 1.,
    y: 0.,
    z: 0.,
};

// These bonds are unit vecs populated by init_local_bond_vecs.

// We use `static mut` here instead of constant, since we need non-const fns (like sin and cos, and
// the linear algebra operations that operate on them) in their construction.
pub const CALPHA_CP_BOND: Vec3 = INIT_BOND_VEC;

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

pub const CP_N_BOND: Vec3 = INIT_BOND_VEC;

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

pub const N_CALPHA_BOND: Vec3 = INIT_BOND_VEC;

pub static mut N_CP_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub const O_CP_BOND: Vec3 = INIT_BOND_VEC;

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

        SIDECHAIN_BOND_OUT2 = CALPHA_CP_BOND;
        SIDECHAIN_BOND_TO_PREV = CALPHA_N_BOND;
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

            SIDECHAIN_BOND_OUT1 = CALPHA_R_BOND;

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
