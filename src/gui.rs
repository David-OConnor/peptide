use core::f64::consts::TAU;

use egui;

use crate::{render_wgpu::STATE, types::State};

const WINDOW_TITLE: &str = "Peptide info";
const WINDOW_SIZE_X: f32 = 900.0;
const WINDOW_SIZE_Y: f32 = 600.0;

// Note: This is draggable.
const SIDE_PANEL_SIZE: f32 = 400.;

// todo: Quick and dirty here, just like in `render_gpu`. Probably avoidable
// todo by using Fn traits instead of `fn` pointers.
pub static mut PROT_NAME: &'static str = "";
pub static mut PDB_IDENT: &'static str = "";

pub static mut ACTIVE_RES_ID: usize = 1;
pub static mut ACTIVE_RES_AA_NAME: &'static str = "";

pub static mut ACTIVE_RES_PHI: f64 = 0.;
pub static mut ACTIVE_RES_PSI: f64 = 0.;
pub static mut ACTIVE_RES_OMEGA: f64 = 0.;

pub static mut ACTIVE_RES_XI_1: Option<f64> = None;
pub static mut ACTIVE_RES_XI_2: Option<f64> = None;
pub static mut ACTIVE_RES_XI_3: Option<f64> = None;
pub static mut ACTIVE_RES_XI_4: Option<f64> = None;
pub static mut ACTIVE_RES_XI_5: Option<f64> = None;

// todo: This is giving us lifetime issues.
/// Initialize our static muts, which we currently (temporarily...) use to manage UI state.
pub unsafe fn init_statics(state: &State) {
    // todo: Temp code here for UI due to temp use of static muts
    //
    //     // let state = state_.clone(); // dodges lifetime issues.
    //
    //     PROT_NAME = &state.protein_descrip.name;
    //     PDB_IDENT = &state.protein_descrip.pdb_ident;
    //
    //     ACTIVE_RES_ID = state.active_residue;
    //
    //     ACTIVE_RES_AA_NAME = state.protein_descrip.residues[state.active_residue - 1]
    //         .sidechain
    //         .aa_name();
    //
    //     ACTIVE_RES_PSI = state.protein_descrip.residues[state.active_residue - 1].φ;
    //     ACTIVE_RES_PHI = state.protein_descrip.residues[state.active_residue - 1].ψ;
    //     ACTIVE_RES_OMEGA = state.protein_descrip.residues[state.active_residue - 1].ω;
}

/// This function draws the (immediate-mode) GUI. We're currently editing it
/// from the main program by modifying the `static mut` variables above.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn draw_ui(ctx: &egui::Context) {
    unsafe {
        let panel = egui::SidePanel::left(0) // ID must be unique among panels.
            .default_width(SIDE_PANEL_SIZE);

        panel.show(ctx, |ui| {
            // println!("{:?}", ui.spacing());
            ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!("Protein: {PROT_NAME}. PDB: {PDB_IDENT}"));

            ui.label(format!("Active Residue: {ACTIVE_RES_ID}"));

            ui.label(ACTIVE_RES_AA_NAME);

            ui.horizontal(|ui| {
                // ui.text_edit_singleline(&mut aa_name);
            });

            ui.label("Backbone dihedral angles:");

            // todo: Put the display in terms of Tau.

            // todo: Ask or otherwise investigate about `from_get_set` vs `new`,
            // todo with not creating a get/set loop.
            // ui.add(
            //     egui::Slider::from_get_set(0.0..=TAU, |v| {
            //         // if v.is_none() {
            //         // //     return 0., // todo???
            //         // } else
            //         //
            //         // {
            //         let ar_i = STATE.active_residue - 1;
            //         let mut active_res = &mut STATE.protein_descrip.residues[ar_i];
            //
            //         // todo: This line is causing problems... might be STATIC_MUT-induced UB.
            //         // active_res.ψ = v.unwrap_or(0.2); // todo: What causes this to be `None`?
            //
            //         // ACTIVE_RES_PSI = active_res.ψ;
            //
            //         // ACTIVE_RES_PSI
            //         0.
            //     })
            //     .text("ψ"),
            // );
            ui.add(egui::Slider::new(&mut ACTIVE_RES_PSI, 0.0..=TAU).text("ψ"));

            ui.add(egui::Slider::new(&mut ACTIVE_RES_PHI, 0.0..=TAU).text("φ"));
            ui.add(egui::Slider::new(&mut ACTIVE_RES_OMEGA, 0.0..=TAU).text("ω"));

            ui.label("Sidechain dihedral angles:");

            if let Some(mut χ) = ACTIVE_RES_XI_1 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ1"));
            }
            if let Some(mut χ) = ACTIVE_RES_XI_2 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ2"));
            }
            if let Some(mut χ) = ACTIVE_RES_XI_3 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ3"));
            }
            if let Some(mut χ) = ACTIVE_RES_XI_4 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ4"));
            }
            if let Some(mut χ) = ACTIVE_RES_XI_5 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ5"));
            }

            // if ui.button("Click each year").clicked() {
                // Perform action here.
            // }
        });
    }
}
