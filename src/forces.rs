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

use lin_alg2::f64::{Quaternion, Vec3};

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
        let dist_o = (point - water.o_posit).magnitude();
        let dist_h_a = (point - water.ha_posit).magnitude();
        let dist_h_b = (point - water.hb_posit).magnitude();
		let dist_m = (point - water.m_posit).magnitude();

        // todo: Charge is from M, but is potential from O?

        let coulomb =
            K_C * ((q * M_CHARGE) / dist_m + (q * H_CHARGE) / dist_h_a + (q * H_CHARGE) / dist_h_b);

        let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        result += coulomb + lj;
    }

    result
}

/// Calculate the force between 2 water molecules.
/// todo: This is probably too simplistic; would need more than just a point
/// todo to represent the atom being acted on; need a way to calculate it's rotation.
// todo: So, force, and torque as results?
pub fn _force(point: Vec3, q: f64, water_molecules: &Vec<WaterMolecule>, i: usize) -> Vec3 {
    let mut result = Vec3::new_zero();

    for (j, water) in water_molecules.iter().enumerate() {
        if i == j {
            continue;
        }

        // Charge is from M, but LJ potential is from O?

        let dist_o = (point - water.o_posit).magnitude();
        let dist_h_a = (point - water.ha_posit).magnitude();
        let dist_h_b = (point - water.hb_posit).magnitude();
		let dist_m = (point - water.m_posit).magnitude();

        let coulomb =
            K_C * ((q * M_CHARGE) / dist_m + (q * H_CHARGE) / dist_h_a + (q * H_CHARGE) / dist_h_b);

        let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        result += (water.o_posit - point) * (coulomb + lj) / dist_o
    }

    result
}

/// Attempt at force; dup of above
/// Output: Forces on O, H_A, H_B, M.
/// todo: Maybe output a single force for the O molecule, and a 3-axis torque (or a quaternion torque?)
pub fn force(acted_on: &WaterMolecule, water_molecules: &Vec<WaterMolecule>, i: usize) -> (Vec3, Vec3, Vec3, Vec3) {
    let mut f_o = Vec3::new_zero();
    let mut f_h_a = Vec3::new_zero();
    let mut f_h_b = Vec3::new_zero();
    let mut f_m = Vec3::new_zero();

    for (j, water) in water_molecules.iter().enumerate() {
        if i == j {
            continue;
        }

        // Charge is from M, but LJ potential is from O?
		
        let o_o = acted_on.o_posit - water.o_posit;
        let o_ha = acted_on.o_posit - water.ha_posit;
        let o_hb = acted_on.o_posit - water.hb_posit;
        let o_m = acted_on.o_posit - water.m_posit;
		
        let ha_o = acted_on.ha_posit - water.o_posit;
        let ha_ha = acted_on.ha - water.ha_posit;
        let ha_hb = acted_on.ha - water.hb_posit;
        let ha_m = acted_on.ha - water.m_posit;
		
		let hb_o = acted_on.hb_posit - water.o_posit;
        let hb_ha = acted_on.hb_posit - water.ha_posit;
        let hb_hb = acted_on.hb_posit - water.hb_posit;
        let hb_m = acted_on.hb_posit - water.o_posit;
		
		let m_o = acted_on.m_posit - water.o_posit;
        let m_ha = acted_on.m_posit - water.ha_posit;
        let m_hb = acted_on.m_posit - water.hb_posit;
        let m_m = acted_on.m_posit - water.m_posit;

        let r_o_ha = o_ha.magnitude();
        let r_o_hb = o_ha.magnitude();
        let r_o_ha = o_ha.magnitude();
        let r_o_m = o_ha.magnitude();

        let r_ha_ha = ha_ha.magnitude();
        let r_ha_hb = ha_hb.magnitude();
        let r_ha_m = ha_m.magnitude();

        let r_hb_hb = hb_hb.magnitude();
        let r_hb_m = hb_m.magnitude();

        let r_m_m = m_m.magnitude();

        // The first letter here is the molecule acted on.
        // Coulomb forces
        // todo: Consts for these.
        let q_h_h = H_CHARGE * H_CHARGE; // todo: Const?
        let q_h_m = H_CHARGE * M_CHARGE;

        let f_ha_ha = K_C * q_h_h * ha_ha / r_ha_ha.powi(2);
        let f_ha_hb = K_C * q_h_h * ha_hb / r_ha_hb.powi(2);
        let f_ha_m = K_C * q_h_m * ha_m / r_ha_m.powi(2);

        let f_hb_ha = K_C * q_h_h * hb_ha / r_ha_m.powi(2);
        let f_hb_hb = K_C * q_h_h * hb_hb / r_hb_hb.powi(2);
        let f_hb_m = K_C * H_CHARGE * M_CHARGE * hb_m / r_hb_m.powi(2);

        let f_m_ha = K_C * q_h_m * m_ha / r_ha_m.powi(2);
        let f_m_hb = K_C * q_h_m * m_hb / r_hb_m.powi(2);
        let f_m_m = K_C * M_CHARGE * M_CHARGE * m_m / r_m_m.powi(2);

        // todo: COme back to the Lennard Jones potential.
        // let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        // todo: You should factor out at least K_C above to save computation! Leaving as-is now
        // todo for clarity.
        let f_total = f_ha_ha + f_ha_hb + f_ha_m + f_hb_ha + f_hb_hb + f_hb_m + f_m_ha + f_m_hb + f_m_m;

        // No torque on O, since we're treating it as the system's rotation center.
        let torque_ha = Vec3::new();
        let torque_hb = Vec3::new();

        let mut torque = Quaternion::new_identity();

        // Calculate a torque on O, so we can rotate the (rigid for now) molecule.

    }

	(f_o, f_h_a, f_h_b, f_m)
}

