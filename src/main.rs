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
// use graphics;

use chem_definitions::AminoAcidType;
use coord_gen::{ProteinDescription, Residue};

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

/// Set up our protein; passed to our initial render state.
fn init_protein() -> ProteinDescription {
    let a = Residue {
        aa: AminoAcidType::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.6 * TAU,
        ψ: 0.6 * TAU,
    };

    let b = Residue {
        aa: AminoAcidType::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.75 * TAU,
        ψ: 0.25 * TAU,
    };

    let c = Residue {
        aa: AminoAcidType::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.2 * TAU,
        ψ: 0. * TAU,
    };

    let a = ProteinDescription { aas: vec![a, b, c] };

    // println!("COORSD: {:?}",
    //     ProteinCoords::from_descrip(&a)
    // );

    a
}

fn main() {
    // init_protein(); // todo temp
    render::run();
    // graphics::run();
}
