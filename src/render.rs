//! This file contains render-related code that's engine agnostic.

use crate::lin_alg::{Quaternion, Vec3};

pub const BACKGROUND_COLOR: (f32, f32, f32) = (0.9, 0.9, 0.9);
pub const BOND_COLOR: (f32, f32, f32) = (0.2, 0.2, 0.2);

pub const LIGHT_INTENSITY: f32 = 1500.0;

pub const BOND_RADIUS: f32 = 0.03;
pub const BOND_N_SIDES: usize = 10;

pub const CAM_MOVE_SENS: f64 = 0.6;
pub const CAM_ROTATE_SENS: f64 = 0.4;
pub const CAM_ROTATE_KEY_SENS: f64 = 1.;

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
