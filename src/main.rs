//! Model, render, and predict protein structure. Uses the peptide bond
//! as an immutable basis for structure. Attempts to find foldin patterns
//! that may lead to an ultimate structure.
//!
//! [A paper on modelling proteins](https://cnx.org/contents/9cMfjngH@6.3:WjXbYFJI@15/Representing-Proteins-in-Silico-and-Protein-Forward-Kinematics)
//! [Article on hydrogen bond modelling](https://www.nature.com/articles/ncomms6803)

#![allow(non_upper_case_globals)]

// todo: Consider starting with a known folded protein, applying heat, and seeing how
// todo it unfolds; this may provide some insight. Perhaps it's a reversal??

// todo: Examine this approach as a basic force model?:
// "We use the hydrophobic/polar (HP)
// model.19 In the HP lattice model, a protein is represented
// as a self-avoiding chain of hydrophobic (P) and hydro-
// philic (H) residues living on a two-dimensional square
// lattice.19 For a given conformation, each pair of non-
// bonded hydrophobic residues in contact contributes one
// unit of favorable energy, e. The HP model is exactly enu-
// merable, and therefore allows us to know for sure
// whether a search method reaches the true global mini-
// mum, or just a local minimum, and yet the model
// presents a search problem that has the same challenges
// a protein has—the folding process seeks a single lowest
// energy native state in a conformational space that grows
// exponentially with chain length, and its complexity
// arises from steric constraints due to chain connectivity
// and excluded volume, and energetic roughness."

// Pymol cheater
// `fetch 1l2y` - Download and display a protein from its PDB identifier
// `get_angle 4/n,4/c,4/ca` - Find the angle between 3 bonds.

// todo: Instead of fixed DT,use a dynamic DT based on frame time subtraction

// todo: Switch from Bevy to Ash or similar. Or perhaps Godot's rust bindings.
// todo: Temperature sensitivity. Surroundin water molecules.
// Initially, focus on modeling the bond angles and backbone. Both
// data to describe, and a 3d render

// Is there value in modeling bonds as something other than straight lines? Curves?
// Something else?

// Doe sfolding begin starting at the end extruded?

extern crate graphics;

use types::ProteinDescription;

mod atom_coords;
mod bond_vecs;
mod chem_definitions;
mod forces;
mod gui;
mod kinematics;
mod proteins;
mod render;
mod render_wgpu;
mod save_load;
mod sc_atom_placement;
mod sidechain;
mod simulate;
mod types;
mod water;

use crate::types::State;

use crate::chem_definitions::AminoAcidType;
use crate::sidechain::{PRO_PHI_MAX, PRO_PHI_MIN};
use std::f64::consts::TAU;

// todo: model the oxygen double-bounded to Cp next.

const BOND_ROTATION_SPEED: f64 = 1.; // radians/s

/// Set up our protein; passed to our initial render state.
fn init_protein() -> ProteinDescription {
    // let mut residues = Vec::new();
    // for i in 0..4 {
    //     residues.push(
    //         // ω, φ, ψ
    //         Residue::new(
    //             1. / 2. * TAU,
    //             φ_sheet,
    //             ψ_sheet,
    //             Sidechain::Gly(sidechain::Gly {}),
    //             // Sidechain::Arg(sidechain::Arg {
    //             //     χ_1: 1. / 2. * TAU,
    //             //     χ_2: 1. / 2. * TAU,
    //             //     χ_3: 1. / 2. * TAU,
    //             //     χ_4: 1. / 2. * TAU,
    //             //     χ_5: 1. / 2. * TAU,
    //             // }),
    //         ),
    //     );
    // }

    let mut prot_test = proteins::make_trp_cage();

    // Clamp to within 0, TAU. Otherwise, the GUI sliders will glitch out.
    // todo: Would this make sense elsewhere, given it's only for the GUI?

    for res in &mut prot_test.residues {
        clamp_angle(&mut res.φ, res.sidechain.aa_type() == AminoAcidType::Pro);
        clamp_angle(&mut res.ψ, false);
        clamp_angle(&mut res.ω, false);

        if let Some(χ) = res.sidechain.get_mut_χ1() {
            clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ2() {
            clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ3() {
            clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ4() {
            clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ5() {
            clamp_angle(χ, false);
        }
    }

    // prot_test.save("trp_cage.prot");
    // let prot_test_load = ProteinDescription::load("trp_cage.prot");

    prot_test
}

// todo: util mod?

pub fn clamp_angle(angle: &mut f64, pro_φ: bool) {
    // Clamp to > 0
    if *angle < 0. {
        // Note that this assumes it's not lower than -TAU.
        *angle += TAU;
    // Clamp to < TAU
    } else if *angle > TAU {
        *angle = *angle % TAU;
    }

    if pro_φ {
        if *angle < PRO_PHI_MIN {
            *angle = PRO_PHI_MIN;
            // *angle = PRO_PHI_MIN + (PRO_PHI_MIN - *angle);
        }
        if *angle > PRO_PHI_MAX {
            *angle = PRO_PHI_MAX;
        }
    }
}

fn main() {
    bond_vecs::init_local_bond_vecs();

    let state = State::new(crate::init_protein());

    render_wgpu::run(state);
}
