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
                position: Vec3::new(0., 0., 7.),
                orientation: Quaternion::new_identity(),
            },
        }
    }
}

/// Set up our protein; passed to our initial render state.
fn init_protein() -> ProteinDescription {
    let φ_helix = -0.715584993317675;
    let ψ_helix = -0.715584993317675;

    let mut residues = Vec::new();
    for i in 0..10 {
        residues.push(
            Residue {
                aa: AminoAcidType::A,
                ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
                φ: φ_helix,
                ψ: ψ_helix,
            }
        );
    }

    let r0 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r1 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r2 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r3 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r4 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r5 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r6 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r7 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r8 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r9 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    let r10 = Residue {
        aa: AminoAcidType::A,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
    };

    ProteinDescription {
        residues,
    }
}

fn main() {
    coord_gen::init_local_bond_vecs();
    render_bevy::run();
    // graphics::run();
}
