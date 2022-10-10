use core::f64::consts::TAU;

use egui;

use crate::{render_wgpu::STATE, types::State};

const WINDOW_TITLE: &str = "Peptide info";
const WINDOW_SIZE_X: f32 = 900.0;
const WINDOW_SIZE_Y: f32 = 600.0;

// Note: This is draggable.
const SIDE_PANEL_SIZE: f32 = 400.;

#[derive(Clone, Default)] // todo: Remove `Clone` a/r
pub struct UiState {
    pub prot_name: String,
    pub pdb_ident: String,
    pub active_res_id: usize,
    pub active_res_aa_name: String,
    pub active_res_φ: f64,
    pub active_res_ψ: f64,
    pub active_res_ω: f64,
    pub active_res_χ1: Option<f64>,
    pub active_res_χ2: Option<f64>,
    pub active_res_χ3: Option<f64>,
    pub active_res_χ4: Option<f64>,
    pub active_res_χ5: Option<f64>,
}

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn run(mut state: UiState) -> Box<dyn Fn(&egui::Context)> {
    Box::new(move |ctx: &egui::Context| {
        let panel = egui::SidePanel::left(0) // ID must be unique among panels.
            .default_width(SIDE_PANEL_SIZE);

        panel.show(ctx, |ui| {
            // println!("{:?}", ui.spacing());
            ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!(
                "Protein: {state.prot_name}. PDB: {state.pdb_ident}"
            ));

            ui.label(format!("Active Residue: {state.active_res_id}"));

            ui.label(state.active_res_aa_name);

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
            ui.add(egui::Slider::new(&mut state.active_res_ψ, 0.0..=TAU).text("ψ"));

            ui.add(egui::Slider::new(&mut state.active_res_φ, 0.0..=TAU).text("φ"));
            ui.add(egui::Slider::new(&mut state.active_res_ω, 0.0..=TAU).text("ω"));

            ui.label("Sidechain dihedral angles:");

            if let Some(mut χ) = state.active_res_χ1 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ1"));
            }
            if let Some(mut χ) = state.active_res_χ2 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ2"));
            }
            if let Some(mut χ) = state.active_res_χ3 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ3"));
            }
            if let Some(mut χ) = state.active_res_χ4 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ4"));
            }
            if let Some(mut χ) = state.active_res_χ5 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ5"));
            }

            // if ui.button("Click each year").clicked() {
            // Perform action here.
            // }
        });
    })
}
