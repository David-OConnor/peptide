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

use crate::{
    bond_vecs::{WATER_BOND_H_A, WATER_BOND_H_B, WATER_BOND_M},
    water::{WaterMolecule, A, B, H_CHARGE, M_CHARGE, O_H_DIST, O_M_DIST},
    wf_lab::{Electron, Nucleus},
};
use lin_alg2::f64::Vec3;

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
pub fn coulomb_force(particle: (Vec3, f64), others: &Vec<(Vec3, f64)>) -> Vec3 {
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

/// Attempt at force; dup of above
/// Output: Forces on O, H_A, H_B, M.
/// todo: Maybe output a single force for the O molecule, and a 3-axis torque (or a quaternion torque?)
/// This approach for coulomb force skips some of our optomizations, but is more general than the
/// hand-tuned code below, and is easier to inspect and debug.
/// This terser approach involves looping through the other molecules multiple times, which isn't great!
// pub fn force(acted_on: &WaterMolecule, water_molecules: &Vec<WaterMolecule>, i: usize) -> (Vec3, Vec3, Vec3, Vec3) {
pub fn water_tipt4(
    acted_on: &WaterMolecule,
    water_molecules: &Vec<WaterMolecule>,
    i: usize,
) -> (Vec3, Vec3) {
    // This approach for coulomb force skips some of our optomizations, but is more general than the
    // hand-tuned code below, and is easier to inspect and debug.
    // This terser approach involves looping through the other molecules multiple times, which isn't great!
    let mut coulomb_data = Vec::new();
    for (j, w) in water_molecules.iter().enumerate() {
        if i == j {
            continue;
        }
        coulomb_data.push((w.m_posit, M_CHARGE));
        coulomb_data.push((w.ha_posit, H_CHARGE));
        coulomb_data.push((w.hb_posit, H_CHARGE));
    }

    // todo: maybe handle this itereating in `coulomb_force()` by making the acted_on first arg a Vec.
    let mut force_coulomb = Vec3::new_zero();
    let mut torque_coulomb = Vec3::new_zero();

    let mut force_on_m = coulomb_force((acted_on.m_posit, M_CHARGE), &coulomb_data);
    let mut force_on_ha = coulomb_force((acted_on.ha_posit, H_CHARGE), &coulomb_data);
    let mut force_on_hb = coulomb_force((acted_on.hb_posit, H_CHARGE), &coulomb_data);

    // Apply K_C prior to calculating torques.
    force_on_m *= K_C; // Factor out for computational reasons.
    force_on_ha *= K_C; // Factor out for computational reasons.
    force_on_hb *= K_C; // Factor out for computational reasons.

    let torque_on_m =
        (acted_on.o_orientation.rotate_vec(unsafe { WATER_BOND_M }) * O_M_DIST).cross(force_on_m);

    let torque_on_ha =
        (acted_on.o_orientation.rotate_vec(WATER_BOND_H_A) * O_H_DIST).cross(force_on_ha);

    let torque_on_hb = (acted_on.o_orientation.rotate_vec(unsafe { WATER_BOND_H_B }) * O_H_DIST)
        .cross(force_on_hb);

    force_coulomb += force_on_m + force_on_ha + force_on_hb;
    torque_coulomb += torque_on_m + torque_on_ha + torque_on_hb;

    let mut force_lj = Vec3::new_zero();
    for (j, water_other) in water_molecules.iter().enumerate() {
        if i == j {
            continue;
        }

        let o_o = acted_on.o_posit - water_other.o_posit;
        let r_o_o = o_o.magnitude();
        force_lj += -(o_o / r_o_o) * 6. * (*B * r_o_o.powi(6) - 2. * *A) / r_o_o.powi(13);
    }

    (force_coulomb + force_lj, torque_coulomb)
}

/// Attempt at force; dup of above
/// Output: Forces on O, H_A, H_B, M.
///  todo: THis commented-out code is more efficient, but less flexible, and tougher to
/// todo reason about.
/// todo: Maybe output a single force for the O molecule, and a 3-axis torque (or a quaternion torque?)
// pub fn force(acted_on: &WaterMolecule, water_molecules: &Vec<WaterMolecule>, i: usize) -> (Vec3, Vec3, Vec3, Vec3) {
pub fn force_tipt4(
    acted_on: &WaterMolecule,
    water_molecules: &Vec<WaterMolecule>,
    i: usize,
) -> (Vec3, Vec3) {
    let mut f_total = Vec3::new_zero();
    let mut τ_total = Vec3::new_zero();

    // The first letter here is the molecule acted on.
    // Coulomb forces
    // todo: Consts for these.
    let q_h_h = H_CHARGE * H_CHARGE; // todo: Const?
    let q_h_m = H_CHARGE * M_CHARGE;

    for (j, water_other) in water_molecules.iter().enumerate() {
        if i == j {
            continue;
        }

        // Charge is from M, but LJ potential is from O.

        // todo: I think these subtraction directions are backwards by traditional definitions,
        // todo but we've compensated elsewhere.
        // These intermediate variables prevent preventcalculation repetitions.
        let o_o = acted_on.o_posit - water_other.o_posit;

        let ha_ha = acted_on.ha_posit - water_other.ha_posit;
        let ha_hb = acted_on.ha_posit - water_other.hb_posit;
        let ha_m = acted_on.ha_posit - water_other.m_posit;

        let hb_ha = acted_on.hb_posit - water_other.ha_posit;
        let hb_hb = acted_on.hb_posit - water_other.hb_posit;
        let hb_m = acted_on.hb_posit - water_other.m_posit;

        let m_ha = acted_on.m_posit - water_other.ha_posit;
        let m_hb = acted_on.m_posit - water_other.hb_posit;
        let m_m = acted_on.m_posit - water_other.m_posit;

        let r_o_o = o_o.magnitude();

        let r_ha_ha = ha_ha.magnitude();
        let r_ha_hb = ha_hb.magnitude();
        let r_ha_m = ha_m.magnitude();

        let r_hb_hb = hb_hb.magnitude();
        let r_hb_m = hb_m.magnitude();

        let r_m_m = m_m.magnitude();

        // Note: We don't include the coulomb const `K_C` in this section; we've factored
        // it out to save computations. Additionally, we factor out the -1 term. See below.
        let f_ha_ha = (ha_ha / r_ha_ha) * q_h_h / r_ha_ha.powi(2);
        let f_ha_hb = (ha_hb / r_ha_hb) * q_h_h / r_ha_hb.powi(2);
        let f_ha_m = (ha_m / r_ha_m) * q_h_m / r_ha_m.powi(2);

        let f_hb_ha = (hb_ha / r_ha_hb) * q_h_h / r_ha_hb.powi(2);
        let f_hb_hb = (hb_hb / r_hb_hb) * q_h_h / r_hb_hb.powi(2);
        let f_hb_m = (hb_m / r_hb_m) * q_h_m / r_hb_m.powi(2);

        let f_m_ha = (m_ha / r_ha_m) * q_h_m / r_ha_m.powi(2);
        let f_m_hb = (m_hb / r_hb_m) * q_h_m / r_hb_m.powi(2);
        let f_m_m = (m_m / r_m_m) * M_CHARGE * M_CHARGE / r_m_m.powi(2);

        // todo: COme back to the Lennard Jones potential.
        // [This article shows force](https://nznano.blogspot.com/2017/11/molecular-dynamics-in-python.html),
        // vice the more easily-discoverable potential.
        // let lj = *A / dist_o.powi(12) - *B / dist_o.powi(6);

        let f_ha = f_ha_ha + f_ha_hb + f_ha_m;
        let f_hb = f_hb_ha + f_hb_hb + f_hb_m;
        let f_m = f_m_ha + f_m_hb + f_m_m;

        let f_coulomb = (f_ha + f_hb + f_m) * K_C;

        // let σ_div_r = LJ_σ_OO / r_o_o; // computation simplifier
        // let f_lj = (o_o / r_o_o) * LJ_4ε * c.powi(12) - c.powi(6));

        // Differentiate the LJ potential to find its force, in conjunction with the vector
        // bewteen O molecules.
        let f_lj = -(o_o / r_o_o) * 6. * (*B * r_o_o.powi(6) - 2. * *A) / r_o_o.powi(13);
        // todo: Try the answer here:
        // https://math.stackexchange.com/questions/1742524/numerical-force-due-to-lennard-jones-potential
        // todo: May or may not be equiv to the A/B formulation above.
        // todo: Or here: https://towardsdatascience.com/the-lennard-jones-potential-35b2bae9446c

        // todo: This proably won't work, since the atoms have diff masses. (?)

        f_total += f_coulomb + f_lj;

        // I suppose Lennard-Jones forces don't affect torque, due to them interacting between O,
        // the center of our rotation. (?) // todo: This might not be right.
        let τ_m =
            (acted_on.o_orientation.rotate_vec(unsafe { WATER_BOND_M }) * O_M_DIST).cross(f_m);

        let τ_ha = (acted_on.o_orientation.rotate_vec(WATER_BOND_H_A) * O_H_DIST).cross(f_ha);

        let τ_hb =
            (acted_on.o_orientation.rotate_vec(unsafe { WATER_BOND_H_B }) * O_H_DIST).cross(f_hb);

        // let t_total =
        τ_total += τ_ha + τ_hb + τ_m;

        // todo: Is this right?

        // No torque on O, since we're treating it as the system's rotation center.
        // let mut torque_ha = Quaternion::from_axis_angle(t_ha, f_ha.magnitude());
        // let mut torque_hb = Quaternion::from_axis_angle(t_hb, f_hb.magnitude());
        // let mut torque_m = Quaternion::from_axis_angle(t_m, f_m.magnitude());

        // let torque = torque_ha * torque_hb * torque_m;

        // Calculate a torque on O, so we can rotate the (rigid for now) molecule.
    }

    // (f_o, f_h_a, f_h_b, f_m)
    (f_total, τ_total)
}

pub fn hydrogen_atoms(
    acted_on_prot: &Nucleus,
    acted_on_elec: &Electron,
    protons: &Vec<Nucleus>,
    elec_posits_dynamic: &Vec<Vec3>,
    i: usize,
) -> (Vec3, Vec3) {
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

    (
        coulomb_force((acted_on_prot.position, CHARGE_PROTON), &charges) * K_C,
        Vec3::new_zero(), // todo?

                          // todo: Posit dynamic?
                          // coulomb_force((acted_on_elec.position, CHARGE_ELECTRON), &charges) * K_C,
    )
}
