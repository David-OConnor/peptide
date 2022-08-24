//! Model, render, and predict protein structure. Uses the peptide bond
//! as an immutable basis for structure. Attempts to find foldin patterns
//! that may lead to an ultimate structure.
//!
//! [A paper on modelling proteins](https://cnx.org/contents/9cMfjngH@6.3:WjXbYFJI@15/Representing-Proteins-in-Silico-and-Protein-Forward-Kinematics)
//! [Article on hydrogen bond modelling](https://www.nature.com/articles/ncomms6803)

// todo: Instead of fixed DT,use a dynamic DT based on frame time subtraction

// todo: Switch from Bevy to Ash or similar. Or perhaps Godot's rust bindings.
// todo: Temperature sensitivity. Surroundin water molecules.
// Initially, focus on modeling the bond angles and backbone. Both
// data to describe, and a 3d render

// Doe sfolding begin starting at the end extruded?

use std::f64::consts::TAU;

mod lin_alg;

mod atom_coords;
mod chem_definitions;
mod kinematics;
mod render;
mod render_bevy;
mod render_wgpu;
mod sidechain;

use atom_coords::ProteinCoords;
use chem_definitions::AminoAcidType;
use kinematics::{ProteinDescription, Residue};
use lin_alg::{Quaternion, Vec3};
use render::Camera;
use sidechain::Sidechain;

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

    let φ_sheet = -140. * TAU / 360.;
    let ψ_sheet = 135. * TAU / 360.;

    let mut residues = Vec::new();
    for i in 0..4 {
        residues.push(Residue {
            // aa_type: AminoAcidType::Ala,
            ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
            // φ: φ_helix,
            φ: φ_sheet,
            // ψ: ψ_helix,
            ψ: ψ_sheet,
            sidechain: Sidechain::Asp(sidechain::Asp {
                χ_1: 1. / 2. * TAU,
                χ_2: 1. / 2. * TAU,
                // χ_3: 1. / 2. * TAU,
                // χ_4: 1. / 2. * TAU,
                // χ_5: 1. / 2. * TAU,
            }),
        });
    }

    let r0 = Residue {
        // aa_type: AminoAcidType::Ala,
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r1 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r2 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r3 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r4 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r5 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r6 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r7 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r8 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r9 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    let r10 = Residue {
        ω: 1. / 2. * TAU, // ω Assumed to be TAU/2 for most cases
        φ: 1. / 2. * TAU,
        ψ: 1. / 2. * TAU,
        sidechain: Sidechain::Arg(Default::default()),
    };

    ProteinDescription { residues }
}

fn main() {
    kinematics::init_local_bond_vecs();
    render_bevy::run();
    // render_wgpu::run();
}
