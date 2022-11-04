//! This module contains state related to water molecule simulation.

use rand;

use crate::{
    bond_vecs::{H_BOND_IN, H_BOND_OUT, LEN_O_H_WATER, PLANAR3_C, WATER_BOND_A, WATER_BOND_B},
    kinematics,
};

use crate::atom_coords::AtomCoords;
use crate::chem_definitions::AtomType;
use egui::Key::Q;
use lin_alg2::f64::{Quaternion, Vec3};

// Distance from the origin.
const SIM_BOX_DIST: f64 = 100.;
const VEL_SCALER: f64 = 1.;

pub const N_MOLECULES: usize = 1_000;

#[derive(Clone, Copy, Debug)]
pub enum WaterAtomType {
    O,
    H,
}

#[derive(Debug)]
// todo: Consider if you want this to be a struct, a const of some other struct etc.
pub struct WaterAtom {
    pub atom_type: WaterAtomType,
    pub position: Vec3,
}

/// Describes a water molecule. These aren't directly part of a protein, but may play a role in its
/// folding, among other potential roles.
#[derive(Debug)]
// todo: Consider if you want this to be a struct, a const of some other struct etc.
pub struct WaterMolecule {
    /// Worldspace coordinates of the O atom.
    pub position_o_world: Vec3,
    /// Using the same orientation ref as protein atoms.
    pub orientation: Quaternion,
    pub velocity: Vec3,
}

/// State for the water molecule environment surrounding a protein or other molecule collection
/// of interest.
pub struct WaterEnvironment {
    /// Persistent state
    pub water_molecules: Vec<WaterMolecule>,
    /// Generated from `water_molecules`.
    pub atom_positions: Vec<WaterAtom>,
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

            molecules.push(WaterMolecule {
                position_o_world,
                orientation,
                velocity,
            });
        }

        let mut result = Self {
            water_molecules: molecules,
            atom_positions: Vec::new(),
        };

        result.update_atom_posits();

        result
    }

    pub fn update_atom_posits(&mut self) {
        self.atom_positions = Vec::new();

        for water_molecule in &self.water_molecules {
            // Oxygen
            self.atom_positions.push(WaterAtom {
                atom_type: WaterAtomType::O,
                position: water_molecule.position_o_world,
            });

            let (h_a_posit, _) = kinematics::find_atom_placement(
                water_molecule.orientation,
                H_BOND_IN,
                unsafe { H_BOND_OUT },
                0.,
                water_molecule.position_o_world,
                Vec3::new_zero(), // todo?
                WATER_BOND_A,
                LEN_O_H_WATER,
            );

            let (h_b_posit, _) = kinematics::find_atom_placement(
                water_molecule.orientation,
                H_BOND_IN,
                unsafe { H_BOND_OUT },
                0.,
                water_molecule.position_o_world,
                Vec3::new_zero(), // todo?
                unsafe { WATER_BOND_B },
                LEN_O_H_WATER,
            );

            self.atom_positions.push(WaterAtom {
                atom_type: WaterAtomType::H,
                position: h_a_posit,
            });

            self.atom_positions.push(WaterAtom {
                atom_type: WaterAtomType::H,
                position: h_b_posit,
            });
        }
    }
}
