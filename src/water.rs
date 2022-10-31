//! This module contains state related to water molecule simulation.

use rand::Rand;

use egui::Key::Q;
use lin_alg2::f64::{Quaternion, Vec3};

// Centered around the origin.
const SIM_BOX_WIDTH: f64 = 10_000.;
const SIM_BOX_HEIGHT: f64 = 10_000.;

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
        for i in 0..num_molecules {
            molecules.push(WaterMolecule {
                position_o_world: Vec3::new(0., 0., 0.).to_normalized(), // todo: Random
                orientation: Quaternion::new(0., 0., 0., 0.).to_normalized(), // todo random
                velocity: Vec3::new(0., 0., 0.).to_normalized(),         // todo: Random
            });
        }

        Self {
            water_molecules: molecules,
        }
    }
}
