//! Misc utility functionality

use rand;

use lin_alg2::{
    self,
    f64::{Quaternion, Vec3},
};

/// Return a random value in the range of -1 to +1: This is what we need for a random quaternion.
fn rand_for_orientation() -> f64 {
    (rand::random::<f64>() - 0.5) * 2.
}

/// Return a random value in the range of -1 to +1: This is what we need for a random quaternion.
pub fn rand_orientation() -> Quaternion {
    Quaternion::new(
        rand_for_orientation(),
        rand_for_orientation(),
        rand_for_orientation(),
        rand_for_orientation(),
    )
    .to_normalized()
}

pub fn vec3_to_f32(v: Vec3) -> lin_alg2::f32::Vec3 {
    lin_alg2::f32::Vec3::new(v.x as f32, v.y as f32, v.z as f32)
}

pub fn quat_to_f32(q: Quaternion) -> lin_alg2::f32::Quaternion {
    lin_alg2::f32::Quaternion::new(q.w as f32, q.x as f32, q.y as f32, q.z as f32)
}
