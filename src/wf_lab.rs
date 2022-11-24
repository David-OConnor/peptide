//! Wave function lab.

use core::f64::consts::{PI};

use crate::util;

use lin_alg2::f64::{Vec3, Quaternion};

use rand;

const A_0: f64 = 1.; // Bohr radius.
const Z_H: f64 = 1.; // Z is the atomic number.

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

/// https://chem.libretexts.org/Courses/University_of_California_Davis/UCD_Chem_107B%3A_Physical_Chemistry_for_Life_Scientists/Chapters/4%3A_Quantum_Theory/
/// 4.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom
fn h_wavefn(r: f64, n: u8, l: u8, m: u8) -> f64 {
    // todo: QC this.
    // todo: Hardcoded for n=1, l=m=0
    let ρ = Z_H * r / A_0;
    1. / PI.sqrt() * (Z_H / A_0).powf(3./2.) * (-ρ).exp() * r
}

/// todo: Temp sloppy way to map random uniform numbers to hydrogen distance
fn map_h_wf(v: f64) -> f64 {
    let x = if v < 0.1 {
        0.68
    } else if v < 0.2 {
        0.74
    } else if v < 0.3 {
        0.78
    } else if v < 0.4 {
        0.80
    } else if v < 0.5 {
        0.83
    } else if v < 0.6 {
        0.85
    } else if v < 0.7 {
        0.88
    } else if v < 0.8 {
        0.91
    } else if v < 0.9 {
        0.94
    } else {
        1.
    };

    x * 1.18
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

            // Radial component of hydrogen with n = 1, l = 0, m = 0
            let uniform_sample = rand::random::<f64>();
            let r = map_h_wf(uniform_sample);

            // let dist = h_wavefn( * dist, 1, 0, 0);

            // let r = HWaveFunction.sample(R);


            let posit = *ctr_pt + rotation.rotate_vec(Vec3::new(r, 0., 0.));

            self.electron_posits_dynamic.push(posit);
        }
    }
}

// struct HWaveFunction {
//
// }
//
// impl Distribution<T> for HWaveFunction {
//         fn sample<R>(&self, rng: &mut R) -> T
//     where
//         R: Rng + ?Sized {
//
//         }
//
//     fn sample_iter<R>(self, rng: R) -> DistIter<Self, R, T>
//     where
//         R: Rng,
//     {
//
//     }
//
//     fn map<F, S>(self, func: F) -> DistMap<Self, F, T, S>
//     where
//         F: Fn(T) -> S,
//     {
//
//     }
//
// }