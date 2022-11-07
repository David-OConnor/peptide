//! Code related to forces and interactions.
//! See [this Water Model Wikipedia article](https://en.wikipedia.org/wiki/Water_model) for an
//! overview of our initial approach.
//!
//! https://docs.lammps.org/Howto_tip4p.html
//! http://www.sklogwiki.org/SklogWiki/index.php/TIP4P_model_of_water

// todo: 104.52 vice 104.5 degrees

// todo: SPC/E's polarization correction?

// Note: You should seriously consider a flexible model. (eg flexible SPC)

// TIP4P/F?

use once_cell::sync::Lazy;

use crate::water::WaterMolecule;

use lin_alg2::f64::Vec3;

// todo: Temp. Get rid of these A/R
pub const MAG_O: f64 = 1.;
pub const MAG_N: f64 = 1.;

// Consts for TIP4P/2005:
pub const O_MASS: f64 = 15.9994; // Atomic mass units
pub const H_MASS: f64 = 1.008; // Atomic mass units

pub const H_CHARGE: f64 = 0.5564; // Coulombs?
                                  // The O char is placed at M in this model; no charge on O itself.
pub const O_CHARGE: f64 = 0.; // Coulombs?
pub const M_CHARGE: f64 = -2. * H_CHARGE; // Coulombs?
pub const R_OH_BOND: f64 = 0.9572;

// todo: Move these to the `water` module?

// Called in the creation of our bond vecs
pub const θ_HOH_ANGLE: f64 = 1.82421813; // 104.52 to radians.

// Distance between teh center of the oxygen atom, and the charge-center point M; along
// axis between the bonds to H.
pub const O_M_DIST: f64 = 0.1546;

// todo: Does this vary with temp? Ie the `sklogwiki` article above lists ε/k (K) as 78.0
const LJ_ε_OO: f64 = 0.1852;
const LJ_σ_OO: f64 = 3.1589;
const LJ_εσ_OH_HH: f64 = 0.;
const coulomb_cutoff: f64 = 8.5;

// todo: Is this the unit we want?
const K_C: f64 = 332.1; // Å·kcal/(mol·e²);
                        // const K_C: f64 = 14.3996; // eV·Å·e^−2

// `A` and `B` are Lennard-Jones parameters.
// todo: `lazy_static` as required.

// static ref COUNT: usize = HASHMAP.len();
// todo: If it works here, consider for `bond_vecs` too.

// https://www.researchgate.net/figure/The-parameters-of-the-Lennard-Jones-potential-in-water_tbl1_227066211
// const σ: f64 = 2.725; // Angstroms
// const ε: f64 = 4.9115; // *10^-21 J

static A: Lazy<f64> = Lazy::new(|| 4. * LJ_ε_OO * LJ_σ_OO.powi(12));
static B: Lazy<f64> = Lazy::new(|| 4. * LJ_ε_OO * LJ_σ_OO.powi(6));

/// Calculate the potential at a given point in space, using Coulomb's law for the electrostatic interaction, and the Lennard-Jones
/// potential for dispersion and repulsion forces. q is in Coulombs.
pub fn potential(point: Vec3, q: f64, water_molecules: &Vec<WaterMolecule>) -> f64 {
    let mut result = 0.;

    for water in water_molecules.iter() {
        let dist_o = (point - water.o_posit_world).magnitude();
        let dist_m = (point - water.m_posit_world).magnitude();
        let dist_h_a = (point - water.h_a_posit_world).magnitude();
        let dist_h_b = (point - water.h_b_posit_world).magnitude();

        // todo: Charge is from M, but is potential from O?

        let coulomb =
            K_C * ((q * M_CHARGE) / dist_m + (q * H_CHARGE) / dist_h_a + (q * H_CHARGE) / dist_h_b);

        let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        result += coulomb + lj;
    }

    result
}

/// Attempt at force; dup of above
pub fn force(point: Vec3, q: f64, water_molecules: &Vec<WaterMolecule>, i: usize) -> Vec3 {
    let mut result = Vec3::new_zero();

    for (j, water) in water_molecules.iter().enumerate() {
        if i == j {
            continue;
        }

        // Charge is from M, but LJ potential is from O?

        let dist_o = (point - water.o_posit_world).magnitude();
        let dist_m = (point - water.m_posit_world).magnitude();
        let dist_h_a = (point - water.h_a_posit_world).magnitude();
        let dist_h_b = (point - water.h_b_posit_world).magnitude();

        let coulomb =
            K_C * ((q * M_CHARGE) / dist_m + (q * H_CHARGE) / dist_h_a + (q * H_CHARGE) / dist_h_b);

        let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        result += (water.o_posit_world - point) * (coulomb + lj) / dist_o
    }

    result
}
