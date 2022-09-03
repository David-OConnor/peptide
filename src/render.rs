//! This file contains render-related code that's engine agnostic.

use lin_alg2::f64::{Quaternion, Vec3};

pub const BACKGROUND_COLOR: (f32, f32, f32) = (0.9, 0.9, 0.9);

pub const BOND_COLOR_BACKBONE: (f32, f32, f32) = (0.2, 0.2, 0.2);
pub const BOND_COLOR_SIDECHAIN: (f32, f32, f32) = (0.6, 0.6, 0.8);

// Atom colors
pub const CALPHA_COLOR: (f32, f32, f32) = (0.7, 0.86, 0.66);
pub const CP_COLOR: (f32, f32, f32) = (0.86, 0.82, 0.68);
pub const N_COLOR: (f32, f32, f32) = (0., 0., 1.);
pub const O_COLOR: (f32, f32, f32) = (1., 0., 0.);
pub const C_SIDECHAIN_COLOR: (f32, f32, f32) = (0.5, 0.5, 0.5);

pub const LIGHT_INTENSITY: f32 = 5_000.0;

pub const BOND_RADIUS_BACKBONE: f32 = 0.07;
pub const BOND_RADIUS_SIDECHAIN: f32 = 0.03;
pub const BOND_N_SIDES: usize = 10;

// These sensitivities are in units (position), or radians (orientation) per second.
pub const CAM_MOVE_SENS: f64 = 1.1;
pub const CAM_ROTATE_SENS: f64 = 0.3;
pub const CAM_ROTATE_KEY_SENS: f64 = 0.5;
// Move speed multiplier when the run modifier key is held.
pub const RUN_FACTOR: f64 = 5.;

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
