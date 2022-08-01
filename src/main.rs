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

mod render;
mod chem_definitions;
// use graphics;

use lin_alg::{Quaternion, Vec3};
use chem_definitions::{AtomType, AminoAcidType, BackboneRole};

// todo temp while we figure out how to pass val to bevy

// Double bond len of C' to N.
// todo: 1.089000
const LEN_CP_N: f64 = 1.32; // angstrom
const LEN_N_CALPHA: f64 = 1.; // angstrom // todo
const LEN_CALPHA_CP: f64 = 1.; // angstrom // todo

// todo: These may change under diff circumstances, and be diff for the diff atoms
// and diff things on the atom.
const BACKBONE_BOND_ANGLES: f64 = 1.911; // radians

const ROTATION_SPEED: f64 = 1.; // radians/s

#[derive(Debug)]
struct _Protein {
    pub aas: Vec<AaDihedralAngles>,
}

// #[derive(Clone, Copy)]
// enum CarbonBond {
//     // todo: Come back to what these mean for the 3 backbone atoms
//     A,
//     B,
//     C,
//     D,
// }



/// Set up our protein; passed to our initial render state.
fn init_protein() -> Vec<AaDihedralAngles> {
    let a = AaDihedralAngles {
        aa: AminoAcidType::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.6 * TAU,
        ψ: 0.6 * TAU,
    };

    let b = AaDihedralAngles {
        aa: AminoAcidType::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.75 * TAU,
        ψ: 0.25 * TAU,
    };

    let c = AaDihedralAngles {
        aa: AminoAcidType::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.2 * TAU,
        ψ: 0. * TAU,
    };

    vec![a, b, c]
}

fn main() {
    render::run();
    // graphics::run();
}
