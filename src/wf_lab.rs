//! Wave function lab.

use crate::util;

use lin_alg2::f64::{Vec3, Quaternion};

use rand;

#[derive(Clone, Default, Debug)]
pub struct Proton {
    pub position: Vec3,
    pub velocity: Vec3,
}

/// todo: Own module
/// An experimental rendering of the electron sampled at random locations according
/// to a wavefunction; then coulomb attraction applied.
#[derive(Debug)]
pub struct WaveFunctionState {
    pub protons: Vec<Proton>,
    /// Temp struct that models the wave function as a sphere with decreasing probabilities
    /// with distance. This is a crude wave function-like function we use for testing infrastructure
    /// and rendering.
    pub electron_centers: Vec<Vec3>,
    /// Used for rendering and charge calculations; the on-the-fly electron positions
    /// generated each iteration.
    pub electron_posits_dynamic: Vec<Vec3>,

}

impl WaveFunctionState {
    /// Update proton positions based on their velocity. Update electroncs based on their wavefunction.
    pub fn update_posits(&mut self, dt: f32) {
        for prot in &mut self.protons {
            prot.position += prot.velocity * dt as f64;
        }

        self.electron_posits_dynamic = Vec::new();
        for ctr_pt in &mut self.electron_centers {
            let rotation = util:: rand_orientation();

            let dist = rand::random::<f64>(); // todo: Concentrate in center?

            let posit = *ctr_pt + rotation.rotate_vec(Vec3::new(dist, 0., 0.));

            self.electron_posits_dynamic.push(posit);
        }
    }
}
