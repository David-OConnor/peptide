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

// Is there value in modeling bonds as something other than straight lines? Curves?
// Something else?

// Doe sfolding begin starting at the end extruded?

use std::{f64::consts::TAU, thread};

mod atom_coords;
mod bond_vecs;
mod chem_definitions;
mod gui;
mod kinematics;
mod proteins;
mod render;
mod render_wgpu;
mod save_load;
mod sidechain;

use atom_coords::ProteinCoords;
use kinematics::{ProteinDescription, Residue};
use render::Camera;
use sidechain::Sidechain;

use lin_alg2::f64::{Quaternion, Vec3};

use eframe::egui;

// todo: model the oxygen double-bounded to Cp next.

const BOND_ROTATION_SPEED: f64 = 1.; // radians/s

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
    /// Residue id that's selected for rotation. Starts at 1.
    pub active_residue: usize,
    /// Camera position and orientation
    /// todo: DO we want this? Probably not.
    pub cam: Camera,
}

impl Default for State {
    // Required for Bevy init
    fn default() -> Self {
        Self {
            protein_descrip: ProteinDescription {
                name: "".to_owned(),
                pdb_ident: "".to_owned(),
                residues: Vec::new(),
            },
            protein_coords: ProteinCoords {
                atoms_backbone: Vec::new(),
            },
            active_residue: 1,
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
        residues.push(
            // ω, φ, ψ
            Residue::new(
                1. / 2. * TAU,
                φ_sheet,
                ψ_sheet,
                Sidechain::Gly(sidechain::Gly {}),
                // Sidechain::Arg(sidechain::Arg {
                //     χ_1: 1. / 2. * TAU,
                //     χ_2: 1. / 2. * TAU,
                //     χ_3: 1. / 2. * TAU,
                //     χ_4: 1. / 2. * TAU,
                //     χ_5: 1. / 2. * TAU,
                // }),
            ),
        );
    }

    // ProteinDescription { residues }
    proteins::make_trp_cage()
}

struct MyApp {
    name: String,
    age: u32,
}

impl Default for MyApp {
    fn default() -> Self {
        Self {
            name: "Arthur".to_owned(),
            age: 42,
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("My egui Application");
            ui.horizontal(|ui| {
                ui.label("Your name: ");
                ui.text_edit_singleline(&mut self.name);
            });
            ui.add(egui::Slider::new(&mut self.age, 0..=120).text("age"));
            if ui.button("Click each year").clicked() {
                self.age += 1;
            }
            ui.label(format!("Hello '{}', age {}", self.name, self.age));
        });
    }
}

fn main() {
    bond_vecs::init_local_bond_vecs();

    // thread::spawn(|| {
    //     gui::setup();
    // });

    // todo: unsafe here is temp due to not getting Fn closure support working.
    unsafe {
        render_wgpu::run();
    }

    //
    // let options = eframe::NativeOptions::default();
    // eframe::run_native(
    //     "My egui App",
    //     options,
    //     Box::new(|_cc| Box::new(MyApp::default())),
    // );
}
