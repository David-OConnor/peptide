//! Wave function lab.

use core::f64::consts::PI;

use crate::{time_sim::SIM_BOX_DIST, util};

use lin_alg2::f64::Vec3;

use rand;

const A_0: f64 = 1.; // Bohr radius.
const Z_H: f64 = 1.; // Z is the atomic number.

const N_MOLECULES: usize = 2;
const VEL_SCALER: f64 = 0.0; // todo: increase A/R

// WF mapping precision, per dimension. Memory use and some parts of computation
// scale with the cube of this.
const WF_PRECISION: usize = 60;
// These dists are around each charge, for atomic orbitals. Smaller gives
// more precision around the value we care about; larger takes into account
// less-likely values(?) farther out.
const WF_DIST_MIN: f64 = -8.;
const WF_DIST_MAX: f64 = 8.;

// Render these electrons per the wave function, but don't apply them to
// force calculations; helps to better visualize the cloud.
pub const N_EXTRA_VISIBLE_ELECTRONS: usize = 30;

#[derive(Clone, Default, Debug)]
pub struct Nucleus {
    pub position: Vec3,
    pub velocity: Vec3,
    pub charge: f64,
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
    // pub position: Vec3,
    // pub velocity: Vec3,
    pub spin: Spin,
}

/// todo: Own module
/// An experimental rendering of the electron sampled at random locations according
/// to a wavefunction; then coulomb attraction applied.
#[derive(Debug)]
pub struct WaveFunctionState {
    pub nuclei: Vec<Nucleus>,
    // pub charges: Vec<Charge>,
    /// Temp struct that models the wave function as a sphere with decreasing probabilities
    /// with distance. This is a crude wave function-like function we use for
    /// testing infrastructure
    /// and rendering.
    pub electrons: Vec<Electron>,
    /// Used for rendering and charge calculations; the on-the-fly electron positions
    /// generated each iteration.
    pub electron_posits_dynamic: Vec<Vec3>,
    pub extra_visible_elecs_dynamic: Vec<Vec3>,
    /// Map 0-1 uniform distribution values to radius, for Hydrogen wavefunction,
    /// n = 1, l=m=0
    pub pdf_map: Vec<(f64, Vec3)>,
}

impl WaveFunctionState {
    pub fn build() -> Self {
        let mut nuclei = Vec::new();
        let mut electrons = Vec::new();

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

            nuclei.push(Nucleus {
                position,
                velocity,
                charge: 1.,
            });
            electrons.push(Electron {
                // position,
                // velocity,
                spin: Spin::Up,
            });
        }

        let wfs = vec![(&h_wf_100, 1.), (&h_wf_100, -1.)];

        let nuclei = vec![
            Nucleus {
                position: Vec3::new(5., -10., 0.),
                velocity: Vec3::new_zero(),
                charge: 1.,
            },
            Nucleus {
                position: Vec3::new(5., -12., 0.),
                velocity: Vec3::new_zero(),
                charge: 1.,
            },
        ];

        // We construct a box centered on our area of interest.
        // todo: Dynamic, instead of with fixed dimensions.

        // WF_DIST_MIN, WF_DIST_MAX

        // todo: Temp. Needs to be re-built or modified over time. Maybe by shifting the box
        // todo post map-gen? Not sure. Actually maybe we need to re-gen each time and here is fine?

        // todo: Hard-coded for 2 nuclei
        let ctr_pt = (nuclei[0].position + nuclei[1].position) / 2.;

        let (x_min, x_max) = (ctr_pt.x - WF_DIST_MAX, ctr_pt.x + WF_DIST_MAX);
        let (y_min, y_max) = (ctr_pt.y - WF_DIST_MAX, ctr_pt.y + WF_DIST_MAX);
        let (z_min, z_max) = (ctr_pt.z - WF_DIST_MAX, ctr_pt.z + WF_DIST_MAX);

        // todo: Put wfs back. Having trouble passing fns as args.
        // let pdf_map = generate_pdf_map(&nuclei, wfs, (WF_DIST_MIN, WF_DIST_MAX), WF_PRECISION);
        let pdf_map = generate_pdf_map(
            &nuclei,
            (x_min, x_max),
            (y_min, y_max),
            (z_min, z_max),
            WF_PRECISION,
        );

        let mut result = WaveFunctionState {
            // Charges,
            // electron_centers,
            nuclei,
            electrons: vec![
                Electron {
                    // position: Vec3::new(5., -10., 0.),
                    // velocity: Vec3::new_zero(),
                    spin: Spin::Up,
                },
                // Electron {
                //     // position: Vec3::new(5., -14., 0.),
                //     // velocity: Vec3::new_zero(),
                //     spin: Spin::Up,
                // },
            ],
            electron_posits_dynamic: Vec::new(),
            extra_visible_elecs_dynamic: Vec::new(),

            pdf_map,
        };

        result.update_posits(0.);

        result
    }

    /// Update Charge positions based on their velocity. Update electroncs based on
    /// their wavefunction.
    pub fn update_posits(&mut self, dt: f64) {
        // todo: You only want to update nuclei! Chagne this a/r once you
        // todo have electron charges.

        // todo: Regen map here??

        // todo: DRY this hard-coded 2-nuclei thing from init.
        let ctr_pt = (self.nuclei[0].position + self.nuclei[1].position) / 2.;

        let (x_min, x_max) = (ctr_pt.x - WF_DIST_MAX, ctr_pt.x + WF_DIST_MAX);
        let (y_min, y_max) = (ctr_pt.y - WF_DIST_MAX, ctr_pt.y + WF_DIST_MAX);
        let (z_min, z_max) = (ctr_pt.z - WF_DIST_MAX, ctr_pt.z + WF_DIST_MAX);

        // todo: Put wfs back. Having trouble passing fns as args.
        // let pdf_map = generate_pdf_map(&nuclei, wfs, (WF_DIST_MIN, WF_DIST_MAX), WF_PRECISION);
        self.pdf_map = generate_pdf_map(
            &self.nuclei,
            (x_min, x_max),
            (y_min, y_max),
            (z_min, z_max),
            WF_PRECISION,
        );

        for nuc in &mut self.nuclei {
            nuc.position += nuc.velocity * dt;
        }

        // todo: Temp

        self.electron_posits_dynamic = Vec::new();
        self.extra_visible_elecs_dynamic = Vec::new();

        for electron in &mut self.electrons {
            self.electron_posits_dynamic
                // .push(gen_electron_posit(Vec3::new(0., 0., 0.), &self.pdf_map));
                .push(gen_electron_posit(&self.pdf_map));

            // electron.position += electron.velocity * dt;

            for _ in 0..N_EXTRA_VISIBLE_ELECTRONS {
                self.extra_visible_elecs_dynamic
                    // .push(gen_electron_posit(Vec3::new(0., 0., 0.), &self.pdf_map));
                    .push(gen_electron_posit(&self.pdf_map));
            }
        }
    }
}

/// https://chem.libretexts.org/Courses/University_of_California_Davis/UCD_Chem_107B%3A_
/// Physical_Chemistry_for_Life_Scientists/Chapters/4%3A_Quantum_Theory/
/// 4.10%3A_The_Schr%C3%B6dinger_Wave_Equation_for_the_Hydrogen_Atom
/// Analytic solution for n=1, s orbital
fn h_wf_100(posit_nuc: Vec3, posit_sample: Vec3) -> f64 {
    let diff = posit_sample - posit_nuc;
    let r = (diff.x.powi(2) + diff.y.powi(2) + diff.z.powi(2)).sqrt();

    let ρ = Z_H * r / A_0;
    1. / PI.sqrt() * (Z_H / A_0).powf(3. / 2.) * (-ρ).exp()
    // 1. / sqrt(pi) * 1./ A_0.powf(3. / 2.) * (-ρ).exp()
}

/// Analytic solution for n=2, s orbital
fn h_wf_200(posit_nuc: Vec3, posit_sample: Vec3) -> f64 {
    let diff = posit_sample - posit_nuc;
    let r = (diff.x.powi(2) + diff.y.powi(2) + diff.z.powi(2)).sqrt();

    let ρ = Z_H * r / A_0;
    1. / (32. * PI).sqrt() * (Z_H / A_0).powf(3. / 2.) * (2. - ρ) * (-ρ / 2.).exp()
}

/// Generate discrete mappings between a 0. - 1. uniform distribution
/// to the wave function's PDF: Discretized through 3D space. ie, each
/// PDF value maps to a cube of space.
fn generate_pdf_map(
    // WF, weight
    // todo: Enforce same len of charges and wfs.
    nuclei: &Vec<Nucleus>,
    // todo: Put wfs back. Having trouble with fn passing
    // wave_fns: Vec<(dyn Fn(Vec3, Vec3) -> f64, f64)>,
    x_range: (f64, f64),
    y_range: (f64, f64),
    z_range: (f64, f64),
    // Of a cube, centered on... center-of-mass of system??
    vals_per_side: usize,
) -> Vec<(f64, Vec3)> {
    let x_vals = util::linspace(x_range, vals_per_side);
    let y_vals = util::linspace(y_range, vals_per_side);
    let z_vals = util::linspace(z_range, vals_per_side);

    let mut pdf_cum = 0.;

    let mut gates = Vec::new();

    for x in &x_vals {
        for y in &y_vals {
            for z in &z_vals {
                let posit_sample = Vec3::new(*x, *y, *z);
                // todo: Normalize.

                let mut pdf_this_cube = 0.;

                // todo: Put back! Having trouble with fn passing
                // for (i, (wf, weight)) in wave_fns.into_iter().enumerate() {
                //     // This PDF is psi^2.

                //     // todo: Take into account the charges... charge. (Here?)
                //     pdf_this_cube += (wf(nuclei[i].position, posit_sample) * weight).powi(2);
                // }
                // todo temp hardcoded fns/weights
                pdf_this_cube += h_wf_100(nuclei[0].position, posit_sample) * 1.;
                pdf_this_cube += h_wf_100(nuclei[1].position, posit_sample) * 1.;

                // We are interested in amplitude of wf**2.
                pdf_cum += pdf_this_cube.powi(2);
                gates.push((pdf_cum, posit_sample));
            }
        }
    }

    // Now that we have our gates maping r to a cumulative PDF,
    // map this PDF to our 0-1 output range.

    let scale_factor = pdf_cum / (1.0 - 0.0); // Always works out to be pdf_cum.

    let mut result = Vec::new();
    for (pdf, cube) in gates {
        // Implicit in this is that the output range starts at 0.
        // todo: This currently assumes pdf starts at 0 as well.
        // todo: You will need to change this for other functions.
        result.push((pdf / scale_factor, cube));
    }

    result
}

// /// Using a cumultive probability map, map a uniform RNG value
// /// to a wavefunction value. Assumes increasing
// /// PDF values in the map (it's 0 index.)
// pub fn map_wf(v: f64, map: &Vec<(f64, Vec3)>) -> Vec3 {
//     for (pdf, posit) in map {
//         if v < *pdf {
//             return *posit;
//         }
//     }
//     // If it's the final value.
//     map[map.len() - 1].1
// }

// todo: May need to combine above and below fns to turn this cartesian vice radial.

/// Using a cumultive probability map, map a uniform RNG value
/// to a wavefunction value. Assumes increasing
/// PDF values in the map (it's 0 index.)
///
/// Generate a random electron position, per a center reference point, and the wave
/// function.
// fn gen_electron_posit(ctr_pt: Vec3, map: &Vec<(f64, Vec3)>) -> Vec3 {
fn gen_electron_posit(map: &Vec<(f64, Vec3)>) -> Vec3 {
    let uniform_sample = rand::random::<f64>();

    // todo: we can't interpolate unless the grid mapping is continuous.
    // todo currently, it wraps.

    for (i, (pdf, posit)) in map.into_iter().enumerate() {
        if uniform_sample < *pdf {
            // Interpolate.
            let v = if i > 0 {
                // todo: QC this. If you're having trouble, TS by using just *posit,
                // todo as below.
                // util::map_linear(uniform_sample, (map[i - 1].0, *pdf), (map[i - 1].1, *posit))
                *posit
            } else {
                // todo: Map to 0 here?
                *posit
            };

            // return ctr_pt + v;
            return v;
        }
    }

    // If it's the final value, return this.
    // todo: Map lin on this too.?

    // ctr_pt + map[map.len() - 1].1
    map[map.len() - 1].1

    // center_pt
    //     + util::map_linear(
    //         uniform_sample,
    //         (map[map.len() - 2].0, map[map.len() - 1].0),
    //         (map[map.len() - 2].1, map[map.len() - 1].1),
    //     )
}
