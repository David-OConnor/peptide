use core::f64::consts::TAU;

use egui;

use graphics::Scene;

use crate::{
    atom_coords::ProteinCoords,
    chem_definitions::{AminoAcidType, BackboneRole},
    render_wgpu,
    sidechain::Sidechain,
    types::{Residue, State},
};

// Note: This is draggable.
const SIDE_PANEL_SIZE: f32 = 400.;
const SLIDER_WIDTH: f32 = 160.;

const SPACE_BETWEEN_SECTIONS: f32 = 20.;

const WINDOW_MARGIN: egui::style::Margin = egui::style::Margin { // todo not working
    left: 10.,
    right: 10.,
    top: 10.,
    bottom: 30.,
};

use lin_alg2::f32::Vec3;

// todo: Unused
#[derive(Default)]
pub struct _StateUi {}

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
fn add_angle_slider(val: &mut f64, label: &str, scene_changed: &mut bool, ui: &mut egui::Ui) {
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

/// Adds a ui area for editing the primary (AA) sequence.
fn add_active_aa_editor(
    ui: &mut egui::Ui,
    state: &mut State,
    scene: &mut Scene,
    scene_changed: &mut bool,
    lighting_changed: &mut bool,
) {
    // todo: Is this the right way for text input?
    let ar_i = state.active_residue - 1;
    let mut res_entry = state.active_residue.to_string();

    // ui.label(format!("Active Residue: {}", state.active_residue));

    ui.horizontal(|ui| {
        ui.label("Active residue:");

        let response = ui.add(egui::TextEdit::singleline(&mut res_entry).desired_width(30.));
        if response.changed() {
            let num = res_entry.parse::<usize>().unwrap_or(1);
            if num < state.protein_descrip.residues.len() + 1 {
                state.active_residue = num;

                change_lit_res(state, scene);
                *scene_changed = true;
                *lighting_changed = true;
            }
        }

        // todo: Set button width.
        if ui.button("+").clicked() {
            if state.active_residue != state.protein_descrip.residues.len() {
                state.active_residue += 1;

                change_lit_res(state, scene);
                *scene_changed = true;
                *lighting_changed = true;
            }
        }
        if ui.button("-").clicked() {
            if state.active_residue != 1 {
                state.active_residue -= 1;

                change_lit_res(state, scene);
                *scene_changed = true;
                *lighting_changed = true;
            }
        }
    });

    ui.label(state.protein_descrip.residues[ar_i].sidechain.aa_name());

    let active_aa_type = state.protein_descrip.residues[ar_i].sidechain.aa_type();
    add_aa_selector(ui, state, scene_changed, ar_i, active_aa_type, 0);

    // if response.lost_focus() && ui.input().key_pressed(egui::Key::Enter) {
    //     // …
    // }

    ui.add_space(SPACE_BETWEEN_SECTIONS);

    ui.horizontal(|ui| {
        // ui.text_edit_singleline(&mut aa_name);
    });

    ui.label("Backbone dihedral angles:");

    // todo: Put the display in terms of Tau.

    let mut active_res = &mut state.protein_descrip.residues[ar_i];

    // We use this syntax instead of the more concise `new` syntax, so we know
    // if we need to change the scene.
    add_angle_slider(&mut active_res.ψ, "ψ", scene_changed, ui);
    add_angle_slider(&mut active_res.φ, "φ", scene_changed, ui);
    add_angle_slider(&mut active_res.ω, "ω", scene_changed, ui);

    ui.label("Sidechain dihedral angles:");

    let mut sc = &mut active_res.sidechain;

    if let Some(mut χ) = sc.get_mut_χ1() {
        add_angle_slider(χ, "χ1", scene_changed, ui);
    }
    if let Some(mut χ) = sc.get_mut_χ2() {
        add_angle_slider(χ, "χ2", scene_changed, ui);
    }
    if let Some(mut χ) = sc.get_mut_χ3() {
        add_angle_slider(χ, "χ3", scene_changed, ui);
    }
    if let Some(mut χ) = sc.get_mut_χ4() {
        add_angle_slider(χ, "χ4", scene_changed, ui);
    }
    if let Some(mut χ) = sc.get_mut_χ5() {
        add_angle_slider(χ, "χ5", scene_changed, ui);
    }
}

/// Adds a ui area for editing the primary (AA) sequence.
fn add_aa_selector(
    ui: &mut egui::Ui,
    state: &mut State,
    scene_changed: &mut bool,
    ar_i: usize,
    initial_type: AminoAcidType,
    sel_id: usize,
) {
    // todo: Maybe a pure enum?
    let mut selected = initial_type;

    // todo: A bit of your API mismatch between AminoAcidType and SideChain
    // todo come out here. Consider just using Sidechain, and making Partial Eq
    // todo match top level only?
    // todo: Change width
    egui::ComboBox::from_id_source(sel_id)
        .width(100.)
        .selected_text(format!("{:?}", selected))
        .show_ui(ui, |ui| {
            // todo: Code shortener
            ui.selectable_value(&mut selected, AminoAcidType::Arg, "Arg");
            ui.selectable_value(&mut selected, AminoAcidType::His, "His");
            ui.selectable_value(&mut selected, AminoAcidType::Lys, "Lys");
            ui.selectable_value(&mut selected, AminoAcidType::Asp, "Asp");
            ui.selectable_value(&mut selected, AminoAcidType::Glu, "Glu");
            ui.selectable_value(&mut selected, AminoAcidType::Ser, "Ser");
            ui.selectable_value(&mut selected, AminoAcidType::Thr, "Thr");
            ui.selectable_value(&mut selected, AminoAcidType::Asn, "Asn");
            ui.selectable_value(&mut selected, AminoAcidType::Gln, "Gln");
            ui.selectable_value(&mut selected, AminoAcidType::Cys, "Cys");
            ui.selectable_value(&mut selected, AminoAcidType::Sec, "Sec");
            ui.selectable_value(&mut selected, AminoAcidType::Gly, "Gly");
            ui.selectable_value(&mut selected, AminoAcidType::Pro, "Pro");
            ui.selectable_value(&mut selected, AminoAcidType::Ala, "Ala");
            ui.selectable_value(&mut selected, AminoAcidType::Val, "Val");
            ui.selectable_value(&mut selected, AminoAcidType::Ile, "Ile");
            ui.selectable_value(&mut selected, AminoAcidType::Leu, "Leu");
            ui.selectable_value(&mut selected, AminoAcidType::Met, "Met");
            ui.selectable_value(&mut selected, AminoAcidType::Phe, "Phe");
            ui.selectable_value(&mut selected, AminoAcidType::Tyr, "Tyr");
            ui.selectable_value(&mut selected, AminoAcidType::Trp, "Trp");
        });

    if selected != state.protein_descrip.residues[ar_i].sidechain.aa_type() {
        let sc = &mut state.protein_descrip.residues[ar_i].sidechain;

        *sc = match selected {
            AminoAcidType::Arg => Sidechain::Arg(Default::default()),
            AminoAcidType::His => Sidechain::His(Default::default()),
            AminoAcidType::Lys => Sidechain::Lys(Default::default()),
            AminoAcidType::Asp => Sidechain::Asp(Default::default()),
            AminoAcidType::Glu => Sidechain::Glu(Default::default()),
            AminoAcidType::Ser => Sidechain::Ser(Default::default()),
            AminoAcidType::Thr => Sidechain::Thr(Default::default()),
            AminoAcidType::Asn => Sidechain::Asn(Default::default()),
            AminoAcidType::Gln => Sidechain::Gln(Default::default()),
            AminoAcidType::Cys => Sidechain::Cys(Default::default()),
            AminoAcidType::Sec => Sidechain::Sec(Default::default()),
            AminoAcidType::Gly => Sidechain::Gly(Default::default()),
            AminoAcidType::Pro => Sidechain::Pro(Default::default()),
            AminoAcidType::Ala => Sidechain::Ala(Default::default()),
            AminoAcidType::Val => Sidechain::Val(Default::default()),
            AminoAcidType::Ile => Sidechain::Ile(Default::default()),
            AminoAcidType::Leu => Sidechain::Leu(Default::default()),
            AminoAcidType::Met => Sidechain::Met(Default::default()),
            AminoAcidType::Phe => Sidechain::Phe(Default::default()),
            AminoAcidType::Tyr => Sidechain::Tyr(Default::default()),
            AminoAcidType::Trp => Sidechain::Trp(Default::default()),
        };

        *scene_changed = true;
    }
}

/// Adds a ui area for editing the primary (AA) sequence.
fn add_sequence_editor(ui: &mut egui::Ui, state: &mut State, scene_changed: &mut bool) {
    ui.horizontal(|ui| {
        ui.label("Add a new residue");
        // todo: Set button width.
        if ui.button("+").clicked() {
            state.protein_descrip.residues.push(Default::default());
            *scene_changed = true;
        }
    });

    egui::containers::ScrollArea::vertical().show(ui, |ui| {
        // Add a lot of widgets here.
        // for (i, res) in state.protein_descrip.residues.iter().enumerate() {
        for ar_i in 0..state.protein_descrip.residues.len() {
            // let response =
            //     ui.add(egui::TextEdit::singleline(&mut res_entry).desired_width(30.));
            // if response.changed() {
            // }
            // Each selector (combo box) needs a unique id.
            let sel_id = 100 + ar_i; // todo?
            let aa_type = state.protein_descrip.residues[ar_i].sidechain.aa_type();

            ui.horizontal(|ui| {
                ui.label((ar_i + 1).to_string());
                add_aa_selector(ui, state, scene_changed, ar_i, aa_type, sel_id);
            });
        }

        ui.spacing_mut().window_margin = WINDOW_MARGIN;
    });
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

        // ctx.set_mar

        panel.show(ctx, |ui| {
            // println!("{:?}", ui.spacing());
            ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

            // todo: Wider values on larger windows?
            // todo: Delegate the numbers etc here to consts etc.
            ui.spacing_mut().slider_width = SLIDER_WIDTH;


            // ui.spacing_mut().window_margin = WINDOW_MARGIN; // todo not working

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!(
                "Protein: {}. PDB: {}",
                state.protein_descrip.name, state.protein_descrip.pdb_ident
            ));

            ui.horizontal(|ui| {
                if ui.button("Save").clicked() {
                    scene_changed = true;
                    lighting_changed = true;
                }
                if ui.button("Load").clicked() {
                    scene_changed = true;
                    lighting_changed = true;
                }
            });

            add_active_aa_editor(ui, state, scene, &mut scene_changed, &mut lighting_changed);

            ui.add_space(SPACE_BETWEEN_SECTIONS);

            add_sequence_editor(ui, state, &mut scene_changed);

            // ui.add_space(SPACE_BETWEEN_SECTIONS); // todo: Not working to add margin.
        });

        if scene_changed {
            // Recalculate coordinates, eg if we've changed an angle, or AA.
            state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
            scene.entities = render_wgpu::generate_entities(&state);
        }

        // todo: Sequence editor, where you can easily see and edit the whole sequence.

        // todo: Only true if you've changed a slider.
        (scene_changed, lighting_changed)
    }
}
