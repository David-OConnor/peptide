//! This module contains state related to water molecule simulation.

use rand;

use egui::Key::Q;
use lin_alg2::f64::{Quaternion, Vec3};

// Distance from the origin.
const SIM_BOX_DIST: f64 = 100.;
const VEL_SCALER: f64 = 1.;

pub const N_MOLECULES: usize = 1_000;

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
    pub water_molecules: Vec<WaterMolecule>,
}

impl WaterEnvironment {
    /// Temp is in kelvin
    pub fn build(num_molecules: usize, temp: f64) -> Self {
        let mut molecules = Vec::new();

        // todo: Consts for these probably.

        for _ in 0..num_molecules {
            // todo: QC and clean up this logic.
            let position_o_world = Vec3::new(
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST/2.,
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST/2.,
                rand::random::<f64>() * SIM_BOX_DIST - SIM_BOX_DIST/2.,
            );

            let orientation = Quaternion::new(
                rand::random::<f64>(),
                rand::random::<f64>(),
                rand::random::<f64>(),
                rand::random::<f64>(),
            )
            .to_normalized();

            let velocity = Vec3::new(
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER/2.,
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER/2.,
                rand::random::<f64>() * VEL_SCALER - VEL_SCALER/2.,
            );

            molecules.push(WaterMolecule {
                position_o_world,
                orientation,
                velocity,
            });
        }

        Self {
            water_molecules: molecules,
        }
    }
}
