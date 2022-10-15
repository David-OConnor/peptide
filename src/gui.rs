use core::f64::consts::TAU;

use egui;

use graphics::Scene;


use crate::{types::State, chem_definitions::BackboneRole, render_wgpu, atom_coords::ProteinCoords};

// Note: This is draggable.
const SIDE_PANEL_SIZE: f32 = 400.;

use lin_alg2::f32::Vec3;

// todo: Unused
#[derive(Default)]
pub struct _StateUi {
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

/// Helper function to add a slider.
fn angle_slider(val: &mut f64, label: &str, scene_changed: &mut bool, ui: &mut egui::Ui) {
    ui.add(
        // todo: Make the sliders wider.
        egui::Slider::from_get_set(0.0..=TAU, |v| {
            if let Some(v_) = v {
                *val = v_;
                *scene_changed = true;
            }

            *val
        })
            .text(label),
    );
}

/// Updating the position of the light source on the active residue.
pub fn change_lit_res(state: &State, scene: &mut Scene) {
    let active_n_posit = state
        .protein_coords
        .atoms_backbone
        .iter()
        .find(|a| a.residue_id == state.active_residue && a.role == BackboneRole::N)
        .unwrap()
        .position;

    scene.lighting.point_lights[1].position = Vec3::new(
        active_n_posit.x as f32,
        active_n_posit.y as f32,
        active_n_posit.z as f32,
    );
}

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn run() -> impl FnMut(&mut State, &egui::Context, &mut Scene) -> (bool, bool) {
    move |state: &mut State, ctx: &egui::Context, scene: &mut Scene| {
        // pub fn run(mut state: StateUi) -> dyn FnMut(&egui::Context) {
        // Box::new(move |ctx: &egui::Context| {

        let mut scene_changed = false;
        let mut lighting_changed = false;

        let panel = egui::SidePanel::left(0) // ID must be unique among panels.
            .default_width(SIDE_PANEL_SIZE);

        panel.show(ctx, |ui| {
            // println!("{:?}", ui.spacing());
            ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

            // todo: Wider values on larger windows?
            // todo: Delegate the numbers etc here to consts etc.
            ui.spacing_mut().slider_width = 160.0;

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!(
                "Protein: {}. PDB: {}",
                state.protein_descrip.name, state.protein_descrip.pdb_ident
            ));

            let ar_i = state.active_residue - 1;

            // todo: Is this the right way for text input?
            // ui.add(egui::TextEdit::singleline(&mut my_string))
            let mut res_entry = state.active_residue.to_string();
            // ui.text_edit_singleline(&mut res_entry);

            ui.label(format!("Active Residue: {}", state.active_residue));

            ui.horizontal(|ui| {
                ui.label("Change active: res:");

                let response = ui.add(egui::TextEdit::singleline(&mut res_entry)
                    .desired_width(30.));
                if response.changed() {
                    let num = res_entry.parse::<usize>().unwrap_or(1);
                    if num < state.protein_descrip.residues.len() + 1 {
                        state.active_residue = num;

                        change_lit_res(state, scene);
                        scene_changed = true;
                        lighting_changed = true;
                    }
                }

                // todo: Set button width.
                if ui.button("+").clicked() {
                    if state.active_residue != state.protein_descrip.residues.len() {
                        state.active_residue += 1;

                        change_lit_res(state, scene);
                        scene_changed = true;
                        lighting_changed = true;
                    }
                }
                if ui.button("-").clicked() {
                    if state.active_residue != 1 {
                        state.active_residue -= 1;

                        change_lit_res(state, scene);
                        scene_changed = true;
                        lighting_changed = true;
                    }
                }
            });

            // if response.lost_focus() && ui.input().key_pressed(egui::Key::Enter) {
            //     // …
            // }

            ui.label(state.protein_descrip.residues[ar_i]
                .sidechain
                .aa_name());

            ui.horizontal(|ui| {
                // ui.text_edit_singleline(&mut aa_name);
            });

            ui.label("Backbone dihedral angles:");

            // todo: Put the display in terms of Tau.

            let mut active_res = &mut state.protein_descrip.residues[ar_i];

            // We use this syntax instead of the more concise `new` syntax, so we know
            // if we need to change the scene.
            angle_slider(&mut active_res.ψ, "ψ", &mut scene_changed, ui);
            angle_slider(&mut active_res.φ, "φ", &mut scene_changed, ui);
            angle_slider(&mut active_res.ω, "ω", &mut scene_changed, ui);

            ui.label("Sidechain dihedral angles:");

            let mut sc = &mut active_res.sidechain;

            if let Some(mut χ) = sc.get_mut_χ1() {
                angle_slider(χ, "χ1", &mut scene_changed, ui);
            }
            if let Some(mut χ) = sc.get_mut_χ2() {
                angle_slider(χ, "χ2", &mut scene_changed, ui);
            }
            if let Some(mut χ) = sc.get_mut_χ3() {
                angle_slider(χ, "χ3", &mut scene_changed, ui);
            }
            if let Some(mut χ) = sc.get_mut_χ4() {
                angle_slider(χ, "χ4", &mut scene_changed, ui);
            }
            if let Some(mut χ) = sc.get_mut_χ5() {
                angle_slider(χ, "χ5", &mut scene_changed, ui);
            }
        });

        if scene_changed {
            // Recalculate coordinates now that we've updated our bond angles
            state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
            scene.entities = render_wgpu::generate_entities(&state);
        }

        // todo: ONly rue if you've changed a slider.
        (scene_changed, lighting_changed)
    }
}
