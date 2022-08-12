//! Model, render, and predict protein structure. Uses the peptide bond
//! as an immutable basis for structure. Attempts to find foldin patterns
//! that may lead to an ultimate structure.

// todo: Switch from Bevy to Ash or similar. Or perhaps Godot's rust bindings.
// todo: Temperature sensitivity. Surroundin water molecules.
// Initially, focus on modeling the bond angles and backbone. Both
// data to describe, and a 3d render

// Doe sfolding begin starting at the end extruded?

use std::f64::consts::TAU;

mod lin_alg;

mod chem_definitions;
mod coord_gen;
mod render;
mod render_bevy;
// mod graphics_wgpu;

use chem_definitions::AminoAcidType;
use coord_gen::{ProteinCoords, ProteinDescription, Residue};
use lin_alg::{Quaternion, Vec3};
use render::Camera;

// todo: model the oxygen double-bounded to Cp next.

const ROTATION_SPEED: f64 = 1.; // radians/s

// #[derive(Clone, Copy)]
// enum CarbonBond {
//     // todo: Come back to what these mean for the 3 backbone atoms
//     A,
//     B,
//     C,
//     D,
// }

/// Store our atom descriptions here, for global state the renderer can access.
struct State {
    /// Descriptions of each amino acid, including its name, and bond angles.
    pub protein_descrip: ProteinDescription,
    /// Stored coordinates, calculated in `coord_gen`.
    pub protein_coords: ProteinCoords,
    /// Residue id that's selected for rotation.
    pub active_residue: usize,
    /// Camera position and orientation
    pub cam: Camera,
}

impl Default for State {
    // Required for Bevy init
    fn default() -> Self {
        Self {
            protein_descrip: ProteinDescription {
                residues: Vec::new(),
            },
            protein_coords: ProteinCoords {
                atoms_backbone: Vec::new(),
            },
            active_residue: 0,
            cam: Camera {
                position: Vec3::new(0., 0., -7.),
                orientation: Quaternion::new_identity(),
            },
        }
    }
}

/// Set up our protein; passed to our initial render state.
fn init_protein() -> ProteinDescription {
    let a = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let b = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let c = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 0.2 * TAU,
        ψ: 0. * TAU,
    };

    ProteinDescription {
        // residues: vec![a, b, c],
        residues: vec![a],
    }
}

fn main() {
    // init_protein(); // todo temp
    render_bevy::run();
    // graphics::run();
}
