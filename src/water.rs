//! This module contains state related to water molecule simulation.

use rand;

use crate::{
    bond_vecs::{H_BOND_IN, H_BOND_OUT, WATER_BOND_H_A, WATER_BOND_H_B, WATER_BOND_M},
    forces::{self},
    kinematics,
};

use lin_alg2::f64::{Quaternion, Vec3};
use once_cell::sync::Lazy;

// Distance from the origin.
const SIM_BOX_DIST: f64 = 20.;
const VEL_SCALER: f64 = 0.0; // todo: increase A/R
const ANG_VEL_SCALER: f64 = 0.0; // todo: increase A/R

pub const N_MOLECULES: usize = 200;

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

// todo: Is this the unit we want?
pub const K_C: f64 = 332.1;

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

            let orientation = Quaternion::new(
                rand::random::<f64>(),
                rand::random::<f64>(),
                rand::random::<f64>(),
                rand::random::<f64>(),
            )
            .to_normalized();

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
        //         o_orientation: Quaternion::new_identity(),
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
