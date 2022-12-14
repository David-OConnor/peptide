//! Wave function lab.

use core::f64::consts::PI;
use graphics::Entity;

use crate::{time_sim::SIM_BOX_DIST, util, wf_lab};

use lin_alg2::f64::{Quaternion, Vec3};

use crate::render_wgpu::Q_I;
use rand;

const A_0: f64 = 1.; // Bohr radius.
const Z_H: f64 = 1.; // Z is the atomic number.

const N_MOLECULES: usize = 200;
const VEL_SCALER: f64 = 0.0; // todo: increase A/R

// Render these electrons per the wave function, but don't apply them to
// force calculations; helps to better visualize the cloud.
pub const N_EXTRA_VISIBLE_ELECTRONS: usize = 30;

#[derive(Clone, Default, Debug)]
pub struct Proton {
    pub position: Vec3,
    pub velocity: Vec3,
}

#[derive(Clone, Copy, Debug)]
pub enum Spin {
    Up,
    Down,
}

impl Default for Spin {
    fn default() -> Self {
        Self::Up
    }
}

#[derive(Clone, Default, Debug)]
pub struct Electron {
    // todo: Is this OK for electron centers?
    pub position: Vec3,
    pub velocity: Vec3,
    pub spin: Spin,
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
    pub electron_centers: Vec<Electron>,
    /// Used for rendering and charge calculations; the on-the-fly electron positions
    /// generated each iteration.
    pub electron_posits_dynamic: Vec<Vec3>,
    pub extra_visible_elecs_dynamic: Vec<Vec3>,
    /// Map 0-1 uniform distribution values to radius, for Hydrogen wavefunction,
    /// n = 1, l=m=0
    pub h_pdf_map: Vec<(f64, f64)>,
}

impl WaveFunctionState {
    pub fn build() -> Self {
        let mut protons = Vec::new();
        let mut electron_centers = Vec::new();

        for _ in 0..N_MOLECULES {
            // todo: QC and clean up this logic.
            let position = Vec3::new(
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
            );

            let velocity = Vec3::new(
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
            );

            protons.push(Proton { position, velocity });
            electron_centers.push(Electron {
                position,
                velocity,
                spin: Spin::Up,
            });
        }

        let mut result = WaveFunctionState {
            protons,
            electron_centers,
            // protons: vec![
            //     Proton {
            //         position: Vec3::new(5., -10., 0.),
            //         velocity: Vec3::new_zero(),
            //     },
            //     Proton {
            //         position: Vec3::new(5., -14., 0.),
            //         velocity: Vec3::new_zero(),
            //     }
            // ],
            // electron_centers: vec![
            //     Vec3::new(5., -10., 0.),
            //     Vec3::new(5., -14., 0.),
            // ],
            electron_posits_dynamic: Vec::new(),
            extra_visible_elecs_dynamic: Vec::new(),
            h_pdf_map: generate_pdf_map(&h_wavefn, (0., 10.), 1_000),
        };

        result.update_posits(0.);

        result
    }

    /// Update proton positions based on their velocity. Update electroncs based on their wavefunction.
    pub fn update_posits(&mut self, dt: f64) {
        for prot in &mut self.protons {
            prot.position += prot.velocity * dt;
        }

        self.electron_posits_dynamic = Vec::new();
        self.extra_visible_elecs_dynamic = Vec::new();

        for electron in &mut self.electron_centers {
            self.electron_posits_dynamic
                .push(gen_electron_posit(electron.position, &self.h_pdf_map));

            // electron.position += electron.velocity * dt;

            for _ in 0..N_EXTRA_VISIBLE_ELECTRONS {
                self.extra_visible_elecs_dynamic
                    .push(gen_electron_posit(electron.position, &self.h_pdf_map));
            }
        }
    }
}

/// Create a set of values in a given range, with a given number of values.
/// Similar to `numpy.linspace`.
/// The result terminates one step before the end of the range.
fn linspace(range: (f64, f64), num_vals: usize) -> Vec<f64> {
    let step = (range.1 - range.0) / num_vals as f64;

    let mut result = Vec::new();

    let mut val = range.0;
    for _ in 0..num_vals {
        result.push(val);
        val += step;
    }

    result
}

/// Generate discrete mappings between a 0. - 1. uniform distribution
/// to the wave function's PDF.
/// Output is in the format `(uniform distro value, radius)`.
fn generate_pdf_map(
    wave_fn: &dyn Fn(f64, u8, u8, u8) -> f64,
    r_range: (f64, f64),
    num_r_vals: usize,
) -> Vec<(f64, f64)> {
    let r_vals = linspace(r_range, num_r_vals);

    let mut pdf_cum = 0.;

    let mut gates = Vec::new();
    for r in r_vals {
        // This PDF is psi^2.
        let pdf_at_r = wave_fn(r, 0, 1, 1).powi(2);

        pdf_cum += pdf_at_r;
        gates.push((pdf_cum, r));
    }

    // Now that we have our gates maping r to a cumulative PDF,
    // map this PDF to our 0-1 output range.

    let scale_factor = pdf_cum / (1.0 - 0.0); // Always works out to be pdf_cum.

    let mut result = Vec::new();
    for (pdf, r) in gates {
        // Implicit in this is that the output range starts at 0.
        // todo: This currently assumes pdf starts at 0 as well.
        // todo: You will need to change this for other functions.
        result.push((pdf / scale_factor, r));
    }

    result
}

/// https://chem.libretexts.org/Courses/University_of_California_Davis/UCD_Chem_107B%3A_Physical_Chemistry_for_Life_Scientists/Chapters/4%3A_Quantum_Theory/
/// 4.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom
fn h_wavefn(r: f64, n: u8, l: u8, m: u8) -> f64 {
    // todo: QC this.
    // todo: Hardcoded for n=1, l=m=0
    let ρ = Z_H * r / A_0;
    1. / PI.sqrt() * (Z_H / A_0).powf(3. / 2.) * (-ρ).exp() * r
}

/// Using a cumultive probability map, map a uniform RNG value
/// to a wavefunction value, eg radial position. Assumes increasing
/// PDF values in the map (it's 0 index.)
pub fn map_h_wf(v: f64, map: &Vec<(f64, f64)>) -> f64 {
    for (pdf, r) in map {
        if v < *pdf {
            return *r;
        }
    }
    // If it's the final value.
    map[map.len() - 1].1
}

/// Generate a random electron position, per a center point, and the wave function.
fn gen_electron_posit(ctr_pt: Vec3, map: &Vec<(f64, f64)>) -> Vec3 {
    let rotation = util::rand_orientation();

    // Radial component of hydrogen with n = 1, l = 0, m = 0
    let uniform_sample = rand::random::<f64>();
    let r = map_h_wf(uniform_sample, map);

    ctr_pt + rotation.rotate_vec(Vec3::new(r, 0., 0.))
}
