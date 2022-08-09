//! This file contains render-related code that's engine agnostic.

use crate::lin_alg::{Quaternion, Vec3};

pub const BACKGROUND_COLOR: (f32, f32, f32) = (0.9, 0.9, 0.9);

pub const BOND_COLOR: (f32, f32, f32) = (0.2, 0.2, 0.2);

// Atom colors
pub const CALPHA_COLOR: (f32, f32, f32) = (0.7, 0.86, 0.66);
pub const CP_COLOR: (f32, f32, f32) = (0.86, 0.82, 0.68);
pub const N_COLOR: (f32, f32, f32) = (0., 0., 1.);
pub const O_COLOR: (f32, f32, f32) = (1., 0., 0.);
pub const C_OTHER_COLOR: (f32, f32, f32) = (0.3, 0.3, 0.3);

pub const LIGHT_INTENSITY: f32 = 1500.0;

pub const BOND_RADIUS: f32 = 0.03;
pub const BOND_N_SIDES: usize = 10;

pub const CAM_MOVE_SENS: f64 = 1.1;
pub const CAM_ROTATE_SENS: f64 = 0.3;
pub const CAM_ROTATE_KEY_SENS: f64 = 0.1;

pub const DT: f64 = 1. / 60.;

pub const UP_VEC: Vec3 = Vec3 {
    x: 0.,
    y: 1.,
    z: 0.,
};
pub const RIGHT_VEC: Vec3 = Vec3 {
    x: 1.,
    y: 0.,
    z: 0.,
};
pub const FWD_VEC: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 1.,
};

pub struct Camera {
    pub position: Vec3,
    pub orientation: Quaternion,
}

// impl Camera {
//     /// Find the camera's forward direction.
//     pub fn fwd(&self) -> Vec3 {
//         // todo: Add this pattern to your article.
//         self.orientation.rotate_vec(FWD_VEC)
//     }
//     /// Find the camera's upforward direction.
//     pub fn up(&self) -> Vec3 {
//         self.orientation.rotate_vec(UP_VEC)
//     }
//
// }
