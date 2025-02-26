//! Misc utility functionality

use lin_alg::{
    self,
    complex_nums::Cplx,
    f64::{Quaternion, Vec3},
};
use rand;

pub type Arr2dReal = Vec<Vec<f64>>;
pub type Arr3dReal = Vec<Vec<Vec<f64>>>;

pub type Arr2d = Vec<Vec<Cplx>>;
pub type Arr3d = Vec<Vec<Vec<Cplx>>>;
pub type Arr4d = Vec<Vec<Vec<Vec<Cplx>>>>;

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

pub fn vec3_to_f32(v: Vec3) -> lin_alg::f32::Vec3 {
    lin_alg::f32::Vec3::new(v.x as f32, v.y as f32, v.z as f32)
}

pub fn quat_to_f32(q: Quaternion) -> lin_alg::f32::Quaternion {
    lin_alg::f32::Quaternion::new(q.w as f32, q.x as f32, q.y as f32, q.z as f32)
}

/// Utility function to linearly map an input value to an output
/// todo: Modified from how we usually do this, with out being a Vec.
pub fn map_linear(val: f64, range_in: (f64, f64), range_out: (Vec3, Vec3)) -> Vec3 {
    // todo: You may be able to optimize calls to this by having the ranges pre-store
    // todo the total range vals.
    let portion = (val - range_in.0) / (range_in.1 - range_in.0);

    (range_out.1 - range_out.0) * portion + range_out.0
}

/// Create a set of values in a given range, with a given number of values.
/// Similar to `numpy.linspace`.
/// The result terminates one step before the end of the range.
pub fn linspace(range: (f64, f64), num_vals: usize) -> Vec<f64> {
    let step = (range.1 - range.0) / num_vals as f64;

    let mut result = Vec::new();

    let mut val = range.0;
    for _ in 0..num_vals {
        result.push(val);
        val += step;
    }

    result
}
