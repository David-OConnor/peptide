//! Wave function lab.

use core::f64::consts::PI;

use lin_alg2::f64::Vec3;
use rand;
use wf_lab::{basis_wfs, complex_nums::Cplx, util::linspace, wf_ops::N, Arr3d, Arr3dReal};

use crate::{
    forces::{CHARGE_ELECTRON, CHARGE_PROTON},
    time_sim::SIM_BOX_DIST,
    util,
};

const A_0: f64 = 1.; // Bohr radius.
const Z_H: f64 = 1.; // Z is the atomic number.

// const N_MOLECULES: usize = 2;
const VEL_SCALER: f64 = 0.0; // todo: increase A/R

// These dists are around each charge, for atomic orbitals. Smaller gives
// more precision around the value we care about; larger takes into account
// less-likely values(?) farther out.

// Each electron will have this many point charges, generated per the PDF.
// Each charge value will be divided by this.
const POINT_CHARGES_PER_ELECTRON: usize = 1;
// Render this many electron points in a frame, for a given electron.
pub const EXTRA_VISIBILE_ELECTRONS: usize = 20;

// WF mapping precision, per dimension. Memory use and some parts of computation
// scale with the cube of this.
// const WF_N: usize = wf_lab::wf_ops::N; // number of values on a side. Computation scales with this ^3.

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
    // /// Temp struct that models the wave function as a sphere with decreasing probabilities
    // /// with distance. This is a crude wave function-like function we use for
    // /// testing infrastructure
    // /// and rendering.
    pub electrons: Vec<Electron>,
    /// Used for rendering and charge calculations; the on-the-fly electron positions
    /// generated each iteration.
    pub electron_posits_dynamic: Vec<Vec3>,
    pub extra_visible_elecs_dynamic: Vec<Vec3>,
    /// Map 0-1 uniform distribution values to radius, for Hydrogen wavefunction,
    /// n = 1, l=m=0
    // pub pdf_map: Vec<(f64, Vec3)>,
    /// Outer `Vec` iterates over electrons; inner one iterates over
    /// PDF gates. (ie cubes or grid points). This is set up as PDF map per electron so
    /// we can exclude an electron's own charge when calculating it's next position.
    pub pdf_maps: Vec<Vec<(f64, Vec3)>>,
    /// Grid min and max. For now, it's a cube centered on 0.
    /// todo: Consider a more flexible setup in the future.
    pub grid_range: (f64, f64),
}

impl WaveFunctionState {
    pub fn build() -> Self {
        // let mut nuclei = Vec::new();

        let n_electrons = 4;

        let mut electrons = Vec::new();
        for _ in 2..n_electrons {
            electrons.push(Electron {
                // position: Vec3::new(5., -10., 0.),
                // velocity: Vec3::new_zero(),
                spin: Spin::Up,
            });
        }
        //
        // for _ in 0..N_MOLECULES {
        //     // todo: QC and clean up this logic.
        //     let position = Vec3::new(
        //         rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
        //         rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
        //         rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
        //     );
        //
        //     let velocity = Vec3::new(
        //         rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
        //         rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
        //         rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
        //     );
        //
        //     nuclei.push(Nucleus {
        //         position,
        //         velocity,
        //         charge: 1.,
        //     });
        //     electrons.push(Electron {
        //         // position,
        //         // velocity,
        //         spin: Spin::Up,
        //     });
        // }

        let nuclei = vec![
            Nucleus {
                position: Vec3::new(5., -10., 0.),
                velocity: Vec3::new_zero(),
                charge: 2.,
            },
            // Nucleus {
            //     position: Vec3::new(5., -12., 0.),
            //     velocity: Vec3::new_zero(),
            //     charge: 1.,
            // },
        ];

        let mut result = WaveFunctionState {
            // Charges,
            // electron_centers,
            nuclei,
            electrons,
            electron_posits_dynamic: Vec::new(),
            extra_visible_elecs_dynamic: Vec::new(),

            pdf_maps: Vec::new(),
            grid_range: (-1., 1.),
        };

        result.update_posits(0.);

        result
    }

    /// Update Charge positions based on their velocity. Update electroncs based on
    /// their wavefunction.
    pub fn update_posits(&mut self, dt: f64) {
        // todo: Regen map here??

        // todo: DRY this hard-coded 2-nuclei thing from init.
        // let ctr_pt = (self.nuclei[0].position + self.nuclei[1].position) / 2.;
        let ctr_pt = self.nuclei[0].position;

        let (x_min, x_max) = (ctr_pt.x - self.grid_range.0, ctr_pt.x + self.grid_range.1);
        let (y_min, y_max) = (ctr_pt.y - self.grid_range.0, ctr_pt.y + self.grid_range.1);
        let (z_min, z_max) = (ctr_pt.z - self.grid_range.0, ctr_pt.z + self.grid_range.1);

        self.pdf_maps = Vec::new();

        /// We generate separate wave functions for each electron, so as to exclude
        /// an electron's charge from influencing its wavefunction.
        for (i, _electron) in self.electrons.iter().enumerate() {
            // Assemble charges from nuclei, and the previous iteration's generated electrons
            // shift by center point so center will be at 0 in the coords we pass to `wf_lab`.
            let mut charges = Vec::new();
            for nuc in &self.nuclei {
                charges.push((nuc.position - ctr_pt, CHARGE_PROTON));
            }

            // Charges from the previous iteration - exclude this electron's own charge.
            for (j, elec_posit) in self.electron_posits_dynamic.iter().enumerate() {
                if i == j {
                    // Don't let the electron interact with itself
                    continue;
                }
                // We divide by the number of point charges per electron, since each dynamic electron
                // generated represents a fraction of an electron charge.
                // More charges per electron is probably more accurate, but is more computationally intensive.
                charges.push((
                    *elec_posit - ctr_pt,
                    CHARGE_ELECTRON / POINT_CHARGES_PER_ELECTRON as f64,
                ));
            }

            // Perform wavefunction computations centered around 0, to avoid floating-point
            // precision issues.
            let psi = wf_lab::psi_from_pt_charges(&charges, &mut self.grid_range);

            let mut charge_density = wf_lab::wf_ops::new_data_real(N);
            wf_lab::wf_ops::charge_density_fm_psi(&psi, &mut charge_density, self.electrons.len());

            // todo: Shift center point from 0 to our local coords of interest.

            self.pdf_maps.push(generate_pdf_map(
                &charge_density,
                (x_min, x_max),
                (y_min, y_max),
                (z_min, z_max),
                N,
            ));
        }

        for nuc in &mut self.nuclei {
            nuc.position += nuc.velocity * dt;
        }

        self.electron_posits_dynamic = Vec::new();
        self.extra_visible_elecs_dynamic = Vec::new();

        for (i, _electron) in self.electrons.iter_mut().enumerate() {
            for _ in 0..POINT_CHARGES_PER_ELECTRON {
                self.electron_posits_dynamic
                    .push(gen_electron_posit(&self.pdf_maps[i]));
            }

            for _ in 0..EXTRA_VISIBILE_ELECTRONS {
                self.extra_visible_elecs_dynamic
                    .push(gen_electron_posit(&self.pdf_maps[i]));
            }
        }
    }
}

/// Generate discrete mappings between a 0. - 1. uniform distribution
/// to the wave function's PDF: Discretized through 3D space. ie, each
/// PDF value maps to a cube of space. <ψ|ψ> is normalized here.
/// This function generates a random number when called.
fn generate_pdf_map(
    charge_density: &Arr3dReal,
    x_range: (f64, f64),
    y_range: (f64, f64),
    z_range: (f64, f64),
    // Of a cube, centered on... center-of-mass of system??
    vals_per_side: usize,
) -> Vec<(f64, Vec3)> {
    let x_vals = linspace(x_range, vals_per_side);
    let y_vals = linspace(y_range, vals_per_side);
    let z_vals = linspace(z_range, vals_per_side);

    let mut pdf_cum = 0.;
    let mut gates = Vec::new();

    // Log the cumulative values of the PDF (charge density), at each point,
    // as we iterate through the array in a specific, but arbitrary order.
    // Note that pdf_cum will range from 0. to 1., since charge_density
    // is normalized. This maps well to a RNG of 0. to 1.
    for (i, x) in x_vals.iter().enumerate() {
        for (j, y) in y_vals.iter().enumerate() {
            for (k, z) in z_vals.iter().enumerate() {
                let posit_sample = Vec3::new(*x, *y, *z);

                // Note: If you end up with non-size-uniform chunks of space,
                // you'll need to incorporate a dVolume term.s

                pdf_cum += charge_density[i][j][k];
                gates.push((pdf_cum, posit_sample));
            }
        }
    }

    // Now that we have our gates maping r to a cumulative PDF,
    // map this PDF to our 0-1 output range.
    const RNG_RANGE: (f64, f64) = (0., 1.);

    let scale_factor = pdf_cum / (RNG_RANGE.1 - RNG_RANGE.0); // Always works out to be pdf_cum.

    let mut result = Vec::new();
    for (pdf, grid_pt) in gates {
        result.push((pdf / scale_factor, grid_pt));
    }

    result
}

// todo: May need to combine above and below fns to turn this cartesian vice radial.

/// Using a cumultive probability map, map a uniform RNG value
/// to a wavefunction value. Assumes increasing
/// PDF values in the map (it's 0 index.)
///
/// Generate a random electron position, per a center reference point, and the wave
/// function.
fn gen_electron_posit(map: &Vec<(f64, Vec3)>) -> Vec3 {
    let uniform_sample = rand::random::<f64>();

    // todo: we can't interpolate unless the grid mapping is continuous.
    // todo currently, it wraps.

    // todo: This approach will need a 3D map, and possibly 3 RNG values. (Or the RNG range
    // todo split in 3). You'll need to modify the PDF map to make this happen.
    // todo: Alternative approach using interpolation:
    // for (i, (pdf, posit)) in map.into_iter().enumerate() {
    //     wf_lab::util::interpolate_spline3pt(surface: Arr3d, val: f64, sfc_range: (f64, f64))
    // }

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

    map[map.len() - 1].1

    // center_pt
    //     + util::map_linear(
    //         uniform_sample,
    //         (map[map.len() - 2].0, map[map.len() - 1].0),
    //         (map[map.len() - 2].1, map[map.len() - 1].1),
    //     )
}
