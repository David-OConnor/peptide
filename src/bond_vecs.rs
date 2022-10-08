//! Code for generating bond vectors for an atom.

// [Includes some common bond angles](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2810841/)

use std::f64::consts::TAU;

use lin_alg2::f64::{Quaternion, Vec3};

// Double bond len of C' to N.
pub const LEN_CP_N: f64 = 1.33; // angstrom
pub const LEN_N_CALPHA: f64 = 1.46; // angstrom
pub const LEN_CALPHA_CP: f64 = 1.53; // angstrom

pub const LEN_CP_O: f64 = 1.2; // angstrom // todo placeholder!
pub const LEN_CALPHA_H: f64 = 1.0; // angstrom // todo placeholder!
pub const LEN_N_H: f64 = 1.0; // angstrom // todo placeholder!

// Ideal bond angles. There are an approximation; from averages. Consider replacing with something
// more robust later. All angles are in radians. We use degrees with math to match common sources.
// R indicates the side chain.
// Our convention is to base each angle off the anchor bond.

// todo: N is listed as 121.7 or 122.7 depending on which side??
pub const BOND_ANGLE_N_CALPHA_CP: f64 = 122.7 * TAU / 360.; // This is to Calpha and C'.
pub const BOND_ANGLE_N_CALPHA_H: f64 = 120. * TAU / 360.; // todo: Placeholder; not real data

// Bond from the Calpha atom
pub const BOND_ANGLE_CALPHA_CP_R: f64 = 110.6 * TAU / 360.;
pub const BOND_ANGLE_CALPHA_CP_N: f64 = 111.0 * TAU / 360.;
pub const BOND_ANGLE_CALPHA_CP_H: f64 = 109.5 * TAU / 360.; // todo: Placeholder; Find this.
pub const BOND_ANGLE_CALPHA_N_R: f64 = 110.6 * TAU / 360.; // todo: Placeholder; Find this.

// todo: It would be more consistent to base all bonds off one (eg N), but the value we
// todo found for the sidechain is from CP. Ideally, get CALPHA_N_R.

// Bonds from the C' atom
// Note that these bonds add up to exactly 360, so these must be along
// a circle (in the same plane).
pub const BOND_ANGLE_CP_N_CALPHA: f64 = 117.2 * TAU / 360.;
pub const BOND_ANGLE_CP_N_O: f64 = 122.7 * TAU / 360.;
pub const BOND_ANGLE_CP_CALPHA_O: f64 = 120.1 * TAU / 360.;

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

pub static mut CALPHA_H_BOND: Vec3 = Vec3 {
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
pub const H_CALPHA_BOND: Vec3 = ANCHOR_BOND_VEC;
pub const H_N_BOND: Vec3 = ANCHOR_BOND_VEC;

// todo: temp
// These are updated in `init_local_bond_vecs`.
pub static mut SIDECHAIN_BOND_TO_PREV: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

pub const SIDECHAIN_BOND_OUT1: Vec3 = ANCHOR_BOND_VEC;

pub static mut SIDECHAIN_BOND_OUT2: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

/// Rotate a vector in 3d a given angle, given a rotation plane.
fn rotate_vec(start: Vec3, angle: f64, rot_plane_norm: Vec3) -> Vec3 {
    let rotation = Quaternion::from_axis_angle(rot_plane_norm, angle);
    rotation.rotate_vec(start)
}

/// Find the third bond by iterating over rotation axes, and using the one that provides
/// the closest match. For each rotation axis, we apply one bond angle as a constraint, and
/// attempt to minimize the second one.
/// The initial rotation plane norm is the plane used to find bond2 by rotating bond1 around it.
fn find_third_bond_vec(
    bond1: Vec3,
    bond2: Vec3,
    angle_1_3: f64,
    angle_2_3: f64,
    bond2_rot_plane_norm: Vec3,
) -> Vec3 {
    // The angle between bond 2 and 3 must be within this (radians) to find our result.
    const EPS: f64 = 0.01;
    // The number of rotation plane to try.
    const N_PLANES: usize = 120;

    let rot_plane_incr = TAU / N_PLANES as f64;

    let mut result = Vec3::new_zero();

    // Rotate vec1 along different rotation planes, starting with the one we used to find bond2.
    // This approximates the set of vectors orthonormal to the anchor vec (bond1)
    for i in 0..N_PLANES {
        let iterated_rot_plane_norm =
            rotate_vec(bond2_rot_plane_norm, rot_plane_incr * i as f64, bond1);

        // Rotate bond 1 around the angle between bond 1 and 3, along our iterated rotation angle.
        result = rotate_vec(bond1, angle_1_3, iterated_rot_plane_norm);

        // Measure the other constraint; the angle between bond 1 and our bond3 candidate.
        let angle_2_3_meas = (bond2.dot(result)).acos();

        // Compare the other constraint against the target; the bond 1-3 angle.
        if (angle_2_3_meas - angle_2_3).abs() < EPS {
            break;
        }
    }

    result
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

    // We use this normal plane for all rotations if the bond are in plane. We use it for
    // the first rotation if not.
    let rot_plane_norm = Vec3::new(0., 1., 0.);

    // The first bond vectors are defined as the anchor vec. These are Calpha's CP bond, Cp's N bond,
    // and N's Calpha bond. This is also the first sidechain bond out.

    unsafe {
        // Find the second bond vectors. The initial anchor bonds (eg `CALPHA_CP_BOND`) are rotated
        // along the (underconstrained) normal plane above.
        CALPHA_N_BOND = rotate_vec(CALPHA_CP_BOND, BOND_ANGLE_CALPHA_CP_N, rot_plane_norm);
        CP_CALPHA_BOND = rotate_vec(CP_N_BOND, BOND_ANGLE_CP_N_CALPHA, rot_plane_norm);
        N_CP_BOND = rotate_vec(N_CALPHA_BOND, BOND_ANGLE_N_CALPHA_CP, rot_plane_norm);

        // todo: We don't have actual H measurements; using a dummy angle for now.
        // todo: Rather than assume it's in plane, we may need to apply the `find_third_bond_vec`
        // todo logic below to find H atom locations.
        N_H_BOND = rotate_vec(N_CALPHA_BOND, TAU * 2. / 3., rot_plane_norm);

        SIDECHAIN_BOND_OUT2 = rotate_vec(ANCHOR_BOND_VEC, TAU / 3., rot_plane_norm);
        SIDECHAIN_BOND_TO_PREV = rotate_vec(ANCHOR_BOND_VEC, TAU * 2. / 3., rot_plane_norm);

        // todo: To find the 3rd (and later 4th, ie hydrogen-to-atom) bonds, we we
        // todo taking an iterative approach based on that above. Long-term, you should
        // todo be able to calculate one of 2 valid choices (for carbon).

        CALPHA_R_BOND = find_third_bond_vec(
            CALPHA_CP_BOND,
            CALPHA_N_BOND,
            BOND_ANGLE_CALPHA_CP_R,
            BOND_ANGLE_CALPHA_N_R,
            rot_plane_norm,
        );

        // println!("CALPHA_R_BOND {:?}", CALPHA_R_BOND);
        // println!("CALPHA_N_BOND {:?}", CALPHA_N_BOND);
        // println!("CALPHA_CP_BOND {:?}", CALPHA_CP_BOND);

        // todo: This is a 4th bond, and you don't have angles to it... this is a temp hack.
        CALPHA_H_BOND = find_third_bond_vec(
            CALPHA_CP_BOND,
            CALPHA_N_BOND,
            // todo: These value are placeholders; they produce rough, but OK results.
            -120. * TAU / 360.,
            120. * TAU / 360.,
            rot_plane_norm,
        );

        // todo temp
        // CALPHA_H_BOND = Vec3::new(0.7, -1., 1.).to_normalized();

        CP_O_BOND = find_third_bond_vec(
            CP_N_BOND,
            CP_CALPHA_BOND,
            BOND_ANGLE_CP_N_O,
            BOND_ANGLE_CP_CALPHA_O,
            rot_plane_norm,
        );

        // todo: The other ones.
    }

    // Find vectors from C' to O, and Cα to R, given the previous 2 bonds for each.
}
