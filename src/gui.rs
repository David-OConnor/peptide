use core::f64::consts::TAU;

use egui::{self, RichText};

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
const AA_SEL_WIDTH: f32 = 60.;
// const KEYBOARD_HELP_WIDTH: f32 = 200.;

const SPACE_BETWEEN_SECTIONS: f32 = 20.;
pub const UI_WIDTH: f32 = 260.;

// When focusing on the active residue, move to this distance.
pub const FOCUS_TARGET_DIST: f32 = 15.;

const ACTIVE_MODE_COLOR: egui::Color32 = egui::Color32::LIGHT_BLUE;
const INACTIVE_MODE_COLOR: egui::Color32 = egui::Color32::LIGHT_GRAY;

const WINDOW_MARGIN: egui::style::Margin = egui::style::Margin {
    // todo not working
    left: 10.,
    right: 10.,
    top: 10.,
    bottom: 30.,
};

const SIM_TIME_SCALE_MIN: f64 = 0.;
const SIM_TIME_SCALE_MAX: f64 = 0.1;

// Water freezing vs boiling
const TEMPERATURE_MIN: f64 = 273.15;
const TEMPERATURE_MAX: f64 = 373.15;

use lin_alg2::f32::Vec3;

// todo: Unused
#[derive(Default)]
pub struct _StateUi {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum UiMode {
    ActiveAaEditor,
    SeqEditor,
    MotionSim,
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

/// Helper function to add a slider.
fn add_angle_slider(val: &mut f64, label: &str, entities_changed: &mut bool, ui: &mut egui::Ui) {
    ui.add(
        egui::Slider::from_get_set(0.0..=TAU, |mut v| {
            if let Some(v_) = v {
                *val = v_;
                *entities_changed = true;
            }

            *val
        })
        .text(label),
    );

    // ui.add(egui::Slider::new(val, 0.0..=TAU).text(label));
    // *entities_changed = true;
}

/// Look at the anchor N atom of the desired res.
fn add_focus_btn(ui: &mut egui::Ui, state: &mut State, scene: &mut Scene, res_id: usize) {
    if ui.button("Focus").clicked() {
        // todo: DRY here between this and the lighting update.
        let active_n_posit = state
            .protein_coords
            .atoms_backbone
            .iter()
            .find(|a| a.residue_id == res_id && a.role == BackboneRole::N)
            .unwrap()
            .position;

        let active_n_posit = Vec3::new(
            active_n_posit.x as f32,
            active_n_posit.y as f32,
            active_n_posit.z as f32,
        );

        render_wgpu::look_at(&mut scene.camera, active_n_posit, FOCUS_TARGET_DIST);
    }
}

/// Adds a ui area for editing the primary (AA) sequence.
fn add_active_aa_editor(
    ui: &mut egui::Ui,
    state: &mut State,
    scene: &mut Scene,
    entities_changed: &mut bool,
    lighting_changed: &mut bool,
) {
    // todo: Is this the right way for text input?
    let ar_i = state.active_residue - 1;
    ui.label(state.protein_descrip.residues[ar_i].sidechain.aa_name());

    let active_aa_type = state.protein_descrip.residues[ar_i].sidechain.aa_type();
    add_aa_selector(
        ui,
        state,
        scene,
        entities_changed,
        lighting_changed,
        ar_i,
        active_aa_type,
        0,
    );

    // if response.lost_focus() && ui.input().key_pressed(egui::Key::Enter) {
    //     // …
    // }

    ui.add_space(SPACE_BETWEEN_SECTIONS);

    ui.label("Backbone dihedral angles:");

    // todo: Put the display in terms of Tau.

    // Update `ar_i` since it may have changed.
    let ar_i = state.active_residue - 1;
    let mut active_res = &mut state.protein_descrip.residues[ar_i];

    // We use this syntax instead of the more concise `new` syntax, so we know
    // if we need to change the scene.
    // todo: The backbone sliders aren't initializing correctly
    add_angle_slider(&mut active_res.ψ, "ψ", entities_changed, ui);
    let mut active_res = &mut state.protein_descrip.residues[ar_i];
    add_angle_slider(&mut active_res.φ, "φ", entities_changed, ui);
    let mut active_res = &mut state.protein_descrip.residues[ar_i];
    add_angle_slider(&mut active_res.ω, "ω", entities_changed, ui);

    ui.label("Sidechain dihedral angles:");

    let mut sc = &mut active_res.sidechain;

    if let Some(χ) = sc.get_mut_χ1() {
        add_angle_slider(χ, "χ1", entities_changed, ui);
    }
    if let Some(χ) = sc.get_mut_χ2() {
        add_angle_slider(χ, "χ2", entities_changed, ui);
    }
    if let Some(χ) = sc.get_mut_χ3() {
        add_angle_slider(χ, "χ3", entities_changed, ui);
    }
    if let Some(χ) = sc.get_mut_χ4() {
        add_angle_slider(χ, "χ4", entities_changed, ui);
    }
    if let Some(χ) = sc.get_mut_χ5() {
        add_angle_slider(χ, "χ5", entities_changed, ui);
    }
}

/// Adds a ui area for editing the primary (AA) sequence.
fn add_aa_selector(
    ui: &mut egui::Ui,
    state: &mut State,
    scene: &mut Scene,
    entities_changed: &mut bool,
    lighting_changed: &mut bool,
    ar_i: usize,
    initial_type: AminoAcidType,
    sel_id: usize,
) {
    // todo: Maybe a pure enum?
    let mut selected = initial_type;

    // todo: A bit of your API mismatch between AminoAcidType and SideChain
    // todo come out here. Consider just using Sidechain, and making Partial Eq
    // todo match top level only?

    // todo: Click one to make it active. Perhaps with a button next to each.

    ui.horizontal(|ui| {
        egui::ComboBox::from_id_source(sel_id)
            .width(AA_SEL_WIDTH)
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

            *entities_changed = true;
        }

        // Click this button to change the active residue to this.
        // todo: button .fill() method to change BG color
        if ui
            .button(
                RichText::new("Sel").color(if state.active_residue == ar_i + 1 {
                    ACTIVE_MODE_COLOR
                } else {
                    INACTIVE_MODE_COLOR
                }),
            )
            .clicked()
        {
            state.active_residue = ar_i + 1;
            change_lit_res(state, scene);

            *entities_changed = true; // to change entity color.
            *lighting_changed = true;
        }

        add_focus_btn(ui, state, scene, ar_i + 1);
    });
}

/// Adds a ui area for editing the primary (AA) sequence.
fn add_sequence_editor(
    ui: &mut egui::Ui,
    state: &mut State,
    scene: &mut Scene,
    entities_changed: &mut bool,
    lighting_changed: &mut bool,
) {
    ui.horizontal(|ui| {
        ui.label("Add a new residue");
        // todo: Set button width.
        if ui.button("+").clicked() {
            state.protein_descrip.residues.push(Default::default());
            *entities_changed = true;
        }
    });

    egui::containers::ScrollArea::vertical().show(ui, |ui| {
        for ar_i in 0..state.protein_descrip.residues.len() {
            // Each selector (combo box) needs a unique id.
            let sel_id = 100 + ar_i; // todo?
            let aa_type = state.protein_descrip.residues[ar_i].sidechain.aa_type();

            ui.horizontal(|ui| {
                ui.label((ar_i + 1).to_string());
                add_aa_selector(
                    ui,
                    state,
                    scene,
                    entities_changed,
                    lighting_changed,
                    ar_i,
                    aa_type,
                    sel_id,
                );
            });
        }

        ui.spacing_mut().window_margin = WINDOW_MARGIN;
    });
}

/// Adds a ui area for editing the primary (AA) sequence.
fn add_motion_sim(
    ui: &mut egui::Ui,
    state: &mut State,
    // scene: &mut Scene,
    // entities_changed: &mut bool,
) {
    ui.horizontal(|ui| {
        if ui.button("Start thermal sim").clicked() {
            state.sim_running = true;
        }
        if ui.button("End thermal sim").clicked() {
            state.sim_running = false;
        }
    });

    ui.add_space(SPACE_BETWEEN_SECTIONS);

    ui.label("Simulation speed");

    ui.add(
        egui::Slider::new(
            &mut state.sim_time_scale,
            SIM_TIME_SCALE_MIN..=SIM_TIME_SCALE_MAX,
        )
        .logarithmic(true)
        .text(""),
    );

    ui.label("Temperature (K)");

    ui.add(
        egui::Slider::new(&mut state.temperature, TEMPERATURE_MIN..=TEMPERATURE_MAX)
            .fixed_decimals(0) // Or it bounces between 0 and 1 rapidly
            .text(""),
    );
}

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn run() -> impl FnMut(&mut State, &egui::Context, &mut Scene) -> (bool, bool) {
    move |state: &mut State, ctx: &egui::Context, scene: &mut Scene| {
        let mut entities_changed = false;
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

            ui.set_max_width(UI_WIDTH);


            // ui.spacing_mut().window_margin = WINDOW_MARGIN; // todo not working

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!(
                "Protein: {}. PDB: {}",
                state.protein_descrip.name, state.protein_descrip.pdb_ident
            ));

            ui.horizontal(|ui| {
                if ui.button("Save").clicked() {
                }
                if ui.button("Load").clicked() {
                }
            });

            ui.horizontal(|ui| {
                ui.label("Active residue:");

                let mut res_entry = state.active_residue.to_string();

                let response = ui.add(egui::TextEdit::singleline(&mut res_entry).desired_width(30.));
                if response.changed() {
                    let num = res_entry.parse::<usize>().unwrap_or(1);
                    if num < state.protein_descrip.residues.len() + 1 {
                        state.active_residue = num;

                        change_lit_res(state, scene);
                        entities_changed = true;
                        lighting_changed = true;
                    }
                }

                // todo: Set button width.
                if ui.button("+").clicked() {
                    if state.active_residue != state.protein_descrip.residues.len() {
                        state.active_residue += 1;

                        change_lit_res(state, scene);
                        entities_changed = true;
                        lighting_changed = true;
                    }
                }
                if ui.button("-").clicked() {
                    if state.active_residue != 1 {
                        state.active_residue -= 1;

                        change_lit_res(state, scene);
                        entities_changed = true;
                        lighting_changed = true;
                    }
                }

                add_focus_btn(ui, state, scene, state.active_residue);
            });

            ui.add_space(SPACE_BETWEEN_SECTIONS);

            ui.horizontal(|ui| {
                // todo: BG color instead of text color?

                // todo: Something to let you know current mode; highlight etc.
                // if ui.selectable_label(state.ui_mode == UiMode::ActiveAaEditor, "Aa editor").clicked() {
                if ui.button(RichText::new("AA editor").color(
                    if state.ui_mode == UiMode::ActiveAaEditor { ACTIVE_MODE_COLOR } else { INACTIVE_MODE_COLOR }
                )).clicked() {
                    state.ui_mode = UiMode::ActiveAaEditor
                }
                // if ui.selectable_label(state.ui_mode == UiMode::SeqEditor, "Seq editor").clicked() {
                if ui.button(RichText::new("Seq editor").color(
                    if state.ui_mode == UiMode::SeqEditor { ACTIVE_MODE_COLOR } else { INACTIVE_MODE_COLOR }
                )).clicked() {
                    state.ui_mode = UiMode::SeqEditor;
                }
                // if ui.selectable_label(state.ui_mode == UiMode::MotionSim, "Motion sim").clicked() {
                if ui.button(RichText::new("Motion sim").color(
                    if state.ui_mode == UiMode::MotionSim { ACTIVE_MODE_COLOR } else { INACTIVE_MODE_COLOR }
                )).clicked() {
                    state.ui_mode = UiMode::MotionSim;
                }
            });

            ui.add_space(SPACE_BETWEEN_SECTIONS);

            match state.ui_mode {
                UiMode::ActiveAaEditor => {
                    add_active_aa_editor(ui, state, scene, &mut entities_changed, &mut lighting_changed);
                }
                UiMode::SeqEditor => {
                    add_sequence_editor(ui, state, scene, &mut entities_changed, &mut lighting_changed);
                }
                UiMode::MotionSim => {
                    // add_motion_sim(ui, state, &mut entities_changed);
                    add_motion_sim(ui, state);
                }
            }

            ui.add_space(SPACE_BETWEEN_SECTIONS * 2.);

            ui.add(egui::Label::new("Keyboard commands:").wrap(true));
            ui.add(egui::Label::new("W/S/A/D: Move forward, back, left, right. C: down. Space: up. Q/E: Roll.").wrap(true));
            ui.add(egui::Label::new("Mouse left click + drag: Pitch and yaw. Up and down arrows: change active residue.").wrap(true));
            ui.add(egui::Label::new("Hold shift to move faster").wrap(true));
            // ui.label("Keyboard commands:").wrap(true);
            // ui.label("WSAD: Move forward, back, left, right. C: down. Space: up. Q/E: Roll.");
            // ui.label("Mouse left click + drag: Pitch and yaw. Up and down arrows: change active residue.");
        });

        if entities_changed {
            // Recalculate coordinates, eg if we've changed an angle, or AA.
            state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
            scene.entities = render_wgpu::generate_entities(&state);
        }

        // todo: Sequence editor, where you can easily see and edit the whole sequence.

        // todo: Only true if you've changed a slider.
        (entities_changed, lighting_changed)
    }
}
