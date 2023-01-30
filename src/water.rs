//! This module contains state related to water molecule simulation.
//! See [this Water Model Wikipedia article](https://en.wikipedia.org/wiki/Water_model) for an
//! overview of our initial approach.

use rand;

use crate::{
    bond_vecs::{H_BOND_IN, H_BOND_OUT, WATER_BOND_H_A, WATER_BOND_H_B, WATER_BOND_M},
    forces, kinematics, util,
};

use crate::forces::K_C;
use crate::time_sim::SIM_BOX_DIST;
use lin_alg2::f64::{Quaternion, Vec3};
use once_cell::sync::Lazy;

const VEL_SCALER: f64 = 0.0; // todo: increase A/R
const ANG_VEL_SCALER: f64 = 0.0; // todo: increase A/R

pub const N_MOLECULES: usize = 0;

// Water model vars below.

// Consts for TIP4P/2005 (See LAMMPS docs above):
pub const O_MASS: f64 = 15.9994;
// Atomic mass units
pub const H_MASS: f64 = 1.008;

pub const H_CHARGE: f64 = 0.5564;
// Coulombs?
// The O char is placed at M in this model; no charge on O itself.
pub const _O_CHARGE: f64 = 0.;
// Coulombs?
pub const M_CHARGE: f64 = -2. * H_CHARGE;

// Called in the creation of our bond vecs
pub const θ_HOH_ANGLE: f64 = 1.82421813;

pub const O_H_DIST: f64 = 0.9572;

// Distance between the center of the oxygen atom, and the charge-center point M; along
// axis between the bonds to H.
pub const O_M_DIST: f64 = 0.1546;

// todo: Does this vary with temp? Ie the `sklogwiki` article above lists ε/k (K) as 78.0
const LJ_ε_OO: f64 = 0.1852;
const LJ_σ_OO: f64 = 3.1589;
const LJ_εσ_OH_HH: f64 = 0.;
const coulomb_cutoff: f64 = 8.5;

pub static A: Lazy<f64> = Lazy::new(|| 4. * LJ_ε_OO * LJ_σ_OO.powi(12));
pub static B: Lazy<f64> = Lazy::new(|| 4. * LJ_ε_OO * LJ_σ_OO.powi(6));

// Saves computation?
const LJ_4ε: f64 = 4. * LJ_ε_OO;

// #[derive(Clone, Copy, Debug)]
// pub enum WaterAtomType {
//     O,
//     H,
// }

// #[derive(Debug)]
// // todo: Consider if you want this to be a struct, a const of some other struct etc.
// pub struct WaterAtom {
//     pub atom_type: WaterAtomType,
//     pub position: Vec3,
// }

/// Describes a water molecule. These aren't directly part of a protein, but may play a role in its
/// folding, among other potential roles. All positions are in world space.
#[derive(Clone, Debug)]
// todo: Consider if you want this to be a struct, a const of some other struct etc.
pub struct WaterMolecule {
    /// Worldspace coordinates of the O atom.
    pub o_posit: Vec3,
    /// Using the same orientation ref as protein atoms.
    pub o_orientation: Quaternion,
    pub velocity: Vec3,
    // todo: Is this right? Seems reasonable if we assume it means apply this rotation
    // todo per second?
    // pub angular_velocity: Quaternion,
    pub angular_velocity: Vec3, // magnitude is rate.
    /// Generated from above vars.
    pub ha_posit: Vec3,
    /// Generated from above vars.
    pub hb_posit: Vec3,
    /// TIP4P M position. Generated from above vars.
    pub m_posit: Vec3,
}

impl WaterMolecule {
    /// Update h and m posits.
    pub fn update_posits(&mut self) {
        let (h_a_posit, _) = kinematics::find_atom_placement(
            self.o_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            0.,
            self.o_posit,
            Vec3::new_zero(), // todo?
            WATER_BOND_H_A,
            O_H_DIST,
        );

        let (h_b_posit, _) = kinematics::find_atom_placement(
            self.o_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            0.,
            self.o_posit,
            Vec3::new_zero(), // todo?
            unsafe { WATER_BOND_H_B },
            O_H_DIST,
        );

        // todo: At least at the start, visualize this by drawing a small sphere on it.
        let (m_posit, _) = kinematics::find_atom_placement(
            self.o_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            0.,
            self.o_posit,
            Vec3::new_zero(), // todo?
            unsafe { WATER_BOND_M },
            O_M_DIST,
        );
        self.ha_posit = h_a_posit;
        self.hb_posit = h_b_posit;
        self.m_posit = m_posit;
    }
}

/// State for the water molecule environment surrounding a protein or other molecule collection
/// of interest.
pub struct WaterEnvironment {
    /// Persistent state
    pub water_molecules: Vec<WaterMolecule>,
    // /// Generated from `water_molecules`.
    // pub atom_positions: Vec<WaterAtom>,
}

impl WaterEnvironment {
    /// Temp is in kelvin
    pub fn build(num_molecules: usize, temp: f64) -> Self {
        let mut molecules = Vec::new();

        // todo: Consts for these probably.

        for _ in 0..num_molecules {
            // todo: QC and clean up this logic.
            let position_o_world = Vec3::new(
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST / 2.,
            );

            let orientation = util::rand_orientation();

            let velocity = Vec3::new(
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER / 2.,
            );

            let angular_velocity = Vec3::new(
                rand::random::<f64>() * ANG_VEL_SCALER - ANG_VEL_SCALER / 2.,
                rand::random::<f64>() * ANG_VEL_SCALER - ANG_VEL_SCALER / 2.,
                rand::random::<f64>() * ANG_VEL_SCALER - ANG_VEL_SCALER / 2.,
            );

            molecules.push(WaterMolecule {
                o_posit: position_o_world,
                o_orientation: orientation,
                velocity,
                angular_velocity,
                ha_posit: Vec3::new_zero(),
                hb_posit: Vec3::new_zero(),
                m_posit: Vec3::new_zero(),
            });
        }

        // todo: Temp
        // let molecules = vec![
        //     WaterMolecule {
        //         o_posit: Vec3::new(10., -10., 0.),
        //         o_orientation: Quaternion::new_identity(),
        //         velocity: Vec3::new_zero(),
        //         angular_velocity: Vec3::new_zero(),
        //         ha_posit: Vec3::new_zero(),
        //         hb_posit: Vec3::new_zero(),
        //         m_posit: Vec3::new_zero(),
        //     },
        //     WaterMolecule {
        //         o_posit: Vec3::new(16., -10., 0.),
        //         o_orientation: Quaternion::from_axis_angle(
        //             Vec3::new(0., 1., 1.).to_normalized(),
        //             5.,
        //         ),
        //         velocity: Vec3::new_zero(),
        //         angular_velocity: Vec3::new_zero(),
        //         ha_posit: Vec3::new_zero(),
        //         hb_posit: Vec3::new_zero(),
        //         m_posit: Vec3::new_zero(),
        //     },
        // ];

        let mut result = Self {
            water_molecules: molecules,
            // atom_positions: Vec::new(),
        };

        result.update_atom_posits();

        result
    }

    pub fn update_atom_posits(&mut self) {
        for water_molecule in &mut self.water_molecules {
            water_molecule.update_posits()
        }
    }
}

/// Experimental Born-Oppenheimer model of water, where positive charges are point charges at the
/// center of each atom, and negative charges are spread out over the electron cloud, using
/// experimental data. They sum to equal the full negative charge. More points is more computationally
/// expensive, but more accurate.
///
/// The coordinate system is with regards to the 2 bond vectors we use for H, with the oxygen atom
/// at its origin.
struct WaterElecCloud {
    /// Our negative point charges, representing electron density. Sums to (10? 8 O + 2 H)
    electron_charges: Vec<(Vec3, f64)>,
}

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

    let mut force_on_m = forces::force_coulomb((acted_on.m_posit, M_CHARGE), &coulomb_data);
    let mut force_on_ha = forces::force_coulomb((acted_on.ha_posit, H_CHARGE), &coulomb_data);
    let mut force_on_hb = forces::force_coulomb((acted_on.hb_posit, H_CHARGE), &coulomb_data);

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
