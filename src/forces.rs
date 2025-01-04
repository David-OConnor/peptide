//! Code related to Neutonian forces and interactions.
//!
//! https://docs.lammps.org/Howto_tip4p.html
//! http://www.sklogwiki.org/SklogWiki/index.php/TIP4P_model_of_water

// todo: 104.52 vice 104.5 degrees

// todo: SPC/E's polarization correction?

// Note: You should seriously consider a flexible model. (eg flexible SPC)

// TIP4P/F?

use lin_alg::f64::Vec3;
use once_cell::sync::Lazy;

use crate::{
    bond_vecs::{WATER_BOND_H_A, WATER_BOND_H_B, WATER_BOND_M},
    quantum::{Electron, Nucleus},
    water::{WaterMolecule, A, B, H_CHARGE, M_CHARGE, O_H_DIST, O_M_DIST},
};

// todo: Is this the unit we want?
pub const K_C: f64 = 332.1;

// todo: Temp. Get rid of these mag values. A/R
pub const MAG_O: f64 = 1.;
pub const MAG_N: f64 = 1.; // Atomic mass units
                           // Coulombs?

// todo?
pub const CHARGE_PROTON: f64 = 1.;
pub const CHARGE_ELECTRON: f64 = -1.;
pub const MASS_PROT: f64 = 1.; // todo: Wrong if using Hartree units.

// todo: Move these to the `water` module?
// 104.52 to radians.
// Å·kcal/(mol·e²);
// const K_C: f64 = 14.3996; // eV·Å·e^−2

// `A` and `B` are Lennard-Jones parameters.
// todo: `lazy_static` as required.

// static ref COUNT: usize = HASHMAP.len();
// todo: If it works here, consider for `bond_vecs` too.

// https://www.researchgate.net/figure/The-parameters-of-the-Lennard-Jones-potential-in-water_tbl1_227066211
// const σ: f64 = 2.725; // Angstroms
// const ε: f64 = 4.9115; // *10^-21 J

/// Calculate the potential at a given point in space, using Coulomb's law for the electrostatic interaction, and the Lennard-Jones
/// potential for dispersion and repulsion forces. q is in Coulombs.
pub fn _potential(point: Vec3, q: f64, water_molecules: &Vec<WaterMolecule>) -> f64 {
    let mut result = 0.;

    for water in water_molecules.iter() {
        let dist_o = (point - water.o_posit).magnitude();
        let dist_h_a = (point - water.ha_posit).magnitude();
        let dist_h_b = (point - water.hb_posit).magnitude();
        let dist_m = (point - water.m_posit).magnitude();

        // todo: Charge is from M, but is potential from O?

        let coulomb =
            K_C * ((q * M_CHARGE) / dist_m + (q * H_CHARGE) / dist_h_a + (q * H_CHARGE) / dist_h_b);

        // let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        // result += coulomb + lj;
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

        // let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        // result += (water.o_posit - point) * (coulomb + lj) / dist_o
    }

    result
}

/// Calculate the coulomb force on a particle, from other particles in a system. The arguments are
/// of type (position, charge). Note that we don't apply the coulomb const here; factored
/// out for computational reasons; apply to the result.
pub fn force_coulomb(particle: (Vec3, f64), others: &Vec<(Vec3, f64)>) -> Vec3 {
    let mut result = Vec3::new_zero();

    let (posit_this, charge_this) = particle;

    for (posit_other, charge_other) in others {
        let posit_diff = posit_this - *posit_other; // todo: Order
        let r = posit_diff.magnitude();

        // Without coulomb const
        let force = posit_diff / r * charge_this * *charge_other / r.powi(2);

        result += force;
    }

    result
}

// todo: Consider a compute shader for these iterate-over-water-mol operations.

/// Simulates, in a neutonian sense, the force on a nucleus from other nuclei, and
/// from electrons. Returns the total force on a specific nucleus.
pub fn atoms(
    nuc_acted_on: &Nucleus,
    protons: &Vec<Nucleus>,
    elec_posits_dynamic: &Vec<Vec3>,
    i: usize,
) -> Vec3 {
    let mut charges = Vec::new();

    for (j, prot_other) in protons.iter().enumerate() {
        if i == j {
            continue;
        }
        // todo: Come back to this
        charges.push((prot_other.position, CHARGE_PROTON));
    }

    // todo: DRY with forces::coulomb_force
    for elec in elec_posits_dynamic {
        charges.push((*elec, CHARGE_ELECTRON));
    }

    force_coulomb((nuc_acted_on.position, CHARGE_PROTON), &charges) * K_C
}
