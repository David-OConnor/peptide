//! Misc utility functionality

use rand;

use lin_alg2::f64::Quaternion;

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
    ).to_normalized()
}