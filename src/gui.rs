use core::f64::consts::TAU;

use egui;

use crate::types::State;

// Note: This is draggable.
const SIDE_PANEL_SIZE: f32 = 400.;

#[derive(Default)]
pub struct StateUi {
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

//fn make_event_handler(
//     state: &mut State,
// ) -> impl FnMut(&mut State, DeviceEvent, &mut Scene, f32) -> bool {
//     // todo: Higher level api from winit or otherwise instead of scancode?
//     let mut coords_changed = false;
//     let mut active_res_backbone_changed = false;
//     let mut active_res_sidechain_changed = false;
//     let mut active_res_changed = false;
//
//     // We count starting at 1, per chem conventions.
//     let ar_i = state.active_residue - 1;
//     // code shortener
//
//     // Box::new(move |event: DeviceEvent, scene: &mut Scene, dt: f32| {
//     move |state: &mut State, event: DeviceEvent, scene: &mut Scene, dt: f32| {

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn run() -> impl FnMut(&mut State, &egui::Context) {
    move |state: &mut State, ctx: &egui::Context| {
        // pub fn run(mut state: StateUi) -> dyn FnMut(&egui::Context) {
        // Box::new(move |ctx: &egui::Context| {
        let panel = egui::SidePanel::left(0) // ID must be unique among panels.
            .default_width(SIDE_PANEL_SIZE);

        panel.show(ctx, |ui| {
            // println!("{:?}", ui.spacing());
            ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!(
                "Protein: {}. PDB: {}",
                state.ui.prot_name, state.ui.pdb_ident
            ));

            ui.label(format!("Active Residue: {}", state.ui.active_res_id));

            ui.label(&state.ui.active_res_aa_name);

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
            //         let ar_i = state.ui.active_residue - 1;
            //         let mut active_res = &mut state.ui.protein_descrip.residues[ar_i];
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
            ui.add(egui::Slider::new(&mut state.ui.active_res_ψ, 0.0..=TAU).text("ψ"));

            ui.add(egui::Slider::new(&mut state.ui.active_res_φ, 0.0..=TAU).text("φ"));
            ui.add(egui::Slider::new(&mut state.ui.active_res_ω, 0.0..=TAU).text("ω"));

            ui.label("Sidechain dihedral angles:");

            if let Some(mut χ) = state.ui.active_res_χ1 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ1"));
            }
            if let Some(mut χ) = state.ui.active_res_χ2 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ2"));
            }
            if let Some(mut χ) = state.ui.active_res_χ3 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ3"));
            }
            if let Some(mut χ) = state.ui.active_res_χ4 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ4"));
            }
            if let Some(mut χ) = state.ui.active_res_χ5 {
                ui.add(egui::Slider::new(&mut χ, 0.0..=TAU).text("χ5"));
            }

            // if ui.button("Click each year").clicked() {
            // Perform action here.
        });
    }
}
