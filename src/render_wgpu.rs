//! This module contains code for use with our custom renderer.

use std::f64::consts::TAU;

use egui::Color32;
use graphics::{
    self, event::WindowEvent, winit::platform::scancode::PhysicalKeyExtScancode, Camera,
    ControlScheme, DeviceEvent, ElementState, EngineUpdates, Entity, InputSettings, LightType,
    Lighting, Mesh, PointLight, Scene, UiSettings, FWD_VEC,
};
use lin_alg::{
    self,
    f32::{Quaternion as QuatF32, Vec3 as Vec3F32},
    f64::{Quaternion, Vec3},
};
use na_seq::AminoAcid;

use crate::{
    atom_coords::{AtomCoords, ProteinCoords},
    bond_vecs::{LEN_CALPHA_CP, LEN_CALPHA_H, LEN_CP_N, LEN_CP_O, LEN_N_CALPHA, LEN_N_H},
    chem_definitions::AtomRole,
    gui,
    // quantum::{self, EXTRA_VISIBILE_ELECTRONS},
    render::{
        self, ACTIVE_CALPHA_COLOR, ATOM_SHINYNESS, BACKGROUND_COLOR, BOND_COLOR_BACKBONE,
        BOND_COLOR_SIDECHAIN, BOND_RADIUS_BACKBONE, BOND_RADIUS_SIDECHAIN, BOND_SHINYNESS,
        CALPHA_COLOR, H_COLOR, H_SCALE, M_SCALE, O_COLOR, RENDER_DIST, WINDOW_SIZE_X,
        WINDOW_SIZE_Y, WINDOW_TITLE,
    },
    sidechain::LEN_SC,
    time_sim,
    types::State,
    util,
    util::{quat_to_f32, vec3_to_f32},
};

// The length-wise axis of our graphics engine's cylinder mesh.
const BOND_MODEL_AXIS: Vec3 = Vec3 {
    x: 0.,
    y: 1.,
    z: 0.,
};

pub const Q_I: QuatF32 = QuatF32 {
    w: 1.,
    x: 0.,
    y: 0.,
    z: 0.,
};

/// Adjust the camera to focus on a new point. Does so by rotating in the shortest direction
/// to point at the new point, then moving forward or back to get to the requested distance.
pub fn look_at(cam: &mut Camera, focus_pt: Vec3F32, dist: f32) {
    // todo: A smooth animation later. For now, just calculate the resulting camera
    // todo position and orientation.

    // todo: Add this to your quaternion article.

    let vec_currently_looking_at = cam.orientation.rotate_vec(FWD_VEC);
    let dir_to_pt = (focus_pt - cam.position).to_normalized();

    let rotation = QuatF32::from_unit_vecs(vec_currently_looking_at, dir_to_pt);

    cam.orientation = rotation * cam.orientation;

    // todo: Why so slow?

    // Now that the orientation is set, move forward or backwards in relation
    // to the new focus point until at the target distance.
    cam.position = focus_pt - dir_to_pt * dist;
}

fn event_dev_handler(
    state: &mut State,
    event: DeviceEvent,
    scene: &mut Scene,
    dt: f32,
) -> EngineUpdates {
    // todo: Higher level api from winit or otherwise instead of scancode?
    let mut entities_changed = false;
    let mut active_res_backbone_changed = false;
    let mut active_res_sidechain_changed = false;
    let mut active_res_changed = false;
    let mut lighting_changed = false;

    // We count starting at 1, per chem conventions.
    let ar_i = state.ui.active_residue - 1;
    // code shortener

    let mut active_res = &mut state.protein.descrip.residues[ar_i];
    let rotation_amt = crate::BOND_ROTATION_SPEED * dt as f64;

    match event {
        DeviceEvent::Key(key) => {
            if key.state == ElementState::Pressed {
                // todo: These should probably be every 10 residues.
                match key.physical_key.to_scancode().unwrap() {
                    // todo: Why are these scan codes for up/down so high??
                    57_416 => {
                        // Up arrow
                        if state.ui.active_residue != state.protein.descrip.residues.len() {
                            state.ui.active_residue += 1;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                    }
                    57_424 => {
                        // Down arrow
                        if state.ui.active_residue != 1 {
                            state.ui.active_residue -= 1;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                    }
                    20 => {
                        // T
                        active_res.φ += rotation_amt;
                        active_res_backbone_changed = true;
                        entities_changed = true;
                    }
                    34 => {
                        // G
                        active_res.φ -= rotation_amt;
                        active_res_backbone_changed = true;
                        entities_changed = true;
                    }
                    21 => {
                        // Y
                        active_res.ψ += rotation_amt;
                        active_res_backbone_changed = true;
                        entities_changed = true;
                    }
                    35 => {
                        // H
                        active_res.ψ -= rotation_amt;
                        active_res_backbone_changed = true;
                        entities_changed = true;
                    }
                    22 => {
                        // U
                        active_res.ω += rotation_amt;
                        active_res_backbone_changed = true;
                        entities_changed = true;
                    }
                    36 => {
                        // J
                        active_res.ω -= rotation_amt;
                        active_res_backbone_changed = true;
                        entities_changed = true;
                    }
                    23 => {
                        // I
                        active_res.sidechain.add_to_χ1(rotation_amt);
                        active_res_sidechain_changed = true;
                        entities_changed = true;
                    }
                    37 => {
                        // K
                        active_res.sidechain.add_to_χ1(-rotation_amt);
                        active_res_sidechain_changed = true;
                        entities_changed = true;
                    }
                    24 => {
                        // O
                        active_res.sidechain.add_to_χ2(rotation_amt);
                        active_res_sidechain_changed = true;
                        entities_changed = true;
                    }
                    38 => {
                        // L
                        active_res.sidechain.add_to_χ2(-rotation_amt);
                        active_res_sidechain_changed = true;
                        entities_changed = true;
                    }
                    36 => {
                        // P
                        active_res.sidechain.add_to_χ3(rotation_amt);
                        active_res_sidechain_changed = true;
                        entities_changed = true;
                    }
                    39 => {
                        // ;
                        active_res.sidechain.add_to_χ3(-rotation_amt);
                        active_res_sidechain_changed = true;
                        entities_changed = true;
                    }
                    // 29 => {
                    //     // Left ctrl
                    // }
                    // todo: Sidechain dihedral angles
                    _ => {}
                }
            }
        }
        // todo: Gamepad not yet supported by winit, but WIP.
        // GamepadEvent::Axis { axis_id, axis, value, stick } => {
        //     println!("Axis id: {:?}, axis: {:?}, value: {:?}, stick: {:?}", axis_id, axis, value, stick)
        // }
        // // GamepadEvent::Stick => {
        // GamepadEvent::Button{ button_id, button, state } => {
        //     // todo: Println etc
        // }
        _ => {}
    }

    if entities_changed {
        // Recalculate coordinates now that we've updated our bond angles
        state.protein.coords = ProteinCoords::from_descrip(&state.protein.descrip);
        scene.entities = generate_entities(&state);
    }

    if active_res_changed {
        // state.ui.active_res_id = state.ui.active_residue;
        //
        // // let aa_name = format!("{}", state.protein.descrip.residues[state.ui.active_residue].sidechain);
        // // let aa_name =
        // state.ui.active_res_aa_name = state.protein.descrip.residues[ar_i]
        //     .sidechain
        //     .aa_name()
        //     .to_owned();

        gui::change_lit_res(state, scene);
        lighting_changed = true;
    }

    if active_res_changed || active_res_backbone_changed {
        // todo: Break this out by psi, phi etc instead of always updating all?
        let res = &mut state.protein.descrip.residues[state.ui.active_residue - 1];

        // todo: Only do this for backbone changd; not acive res

        crate::clamp_angle(&mut res.φ, res.sidechain.aa_type() == AminoAcid::Pro);
        crate::clamp_angle(&mut res.ψ, false);
        crate::clamp_angle(&mut res.ω, false);
    }

    // if active_res_changed || active_res_sidechain_changed {
    //     // todo: Break this out by psi, phi etc instead of always updating all?
    //     // let sc = &state.protein.descrip.residues[state.ui.active_residue - 1].sidechain;
    //     //
    //     // state.ui.active_res_χ1 = sc.get_χ1();
    //     // state.ui.active_res_χ2 = sc.get_χ2();
    //     // state.ui.active_res_χ3 = sc.get_χ3();
    //     // state.ui.active_res_χ4 = sc.get_χ4();
    //     // state.ui.active_res_χ5 = sc.get_χ5();
    // }

    EngineUpdates {
        entities: entities_changed,
        lighting: lighting_changed,
        ..Default::default()
    }
}

fn event_win_handler(
    _state: &mut State,
    _event: WindowEvent,
    _scene: &mut Scene,
    _dt: f32,
) -> EngineUpdates {
    EngineUpdates::default()
}

/// This runs each frame. Update our time-based simulation here.
fn render_handler(state: &mut State, scene: &mut Scene, dt: f32) -> EngineUpdates {
    let mut entities_changed = false;

    // delegate to a sim fn/module

    if state.ui.sim_running {
        time_sim::run_frame(state, dt);

        scene.entities = generate_entities(&state);

        entities_changed = true;
    }

    EngineUpdates {
        entities: entities_changed,
        ..Default::default()
    }
    // todo: This may be where you need to update the render after changing a slider
}

/// todo: Use a vec with vec ops instead?
fn avg_colors(color1: (f32, f32, f32), color2: (f32, f32, f32)) -> (f32, f32, f32) {
    let a = color1.0 + color2.0 / 2.;
    let b = color1.1 + color2.1 / 2.;
    let c = color1.2 + color2.2 / 2.;

    let mag = (a.powi(2) + b.powi(2) + c.powi(2)).sqrt();
    // todo: Normalize
    (a / mag, b / mag, c / mag)
}

/// Add bonds from a given atom. Usually 1.
fn add_bond(
    atom: &AtomCoords,
    atoms_backbone: &Vec<AtomCoords>,
    atom_id: usize,
    // n_id: &mut usize,
    // cα_id: &mut usize,
    // cp_id: &mut usize,
    atom_prev_id: usize,
    entities: &mut Vec<Entity>,
) {
    // Calculate the position of the bond mesh: This is the cylinder's z point,
    // half way between the 2 atoms it connects.

    let atom_prev = &atoms_backbone[atom_prev_id];

    let bond_center_position = (atom.position + atom_prev.position) * 0.5;
    let bond_dir = (atom.position - atom_prev.position).to_normalized();

    let bond_orientation = Quaternion::from_unit_vecs(BOND_MODEL_AXIS, bond_dir);

    let color_a = match atom.role {
        AtomRole::CSidechain | AtomRole::OSidechain | AtomRole::NSidechain => BOND_COLOR_SIDECHAIN,
        _ => BOND_COLOR_BACKBONE,
    };

    // let color = Color::rgb(color_a.0, color_a.1, color_a.2).into();

    let (bond_len, bond_mesh, bond_color) = match atom.role {
        AtomRole::N => (LEN_CP_N as f32, 2, BOND_COLOR_BACKBONE),
        AtomRole::Cα => (LEN_N_CALPHA as f32, 2, BOND_COLOR_BACKBONE),
        AtomRole::Cp => (LEN_CALPHA_CP as f32, 2, BOND_COLOR_BACKBONE),
        AtomRole::O => (LEN_CP_O as f32, 2, BOND_COLOR_BACKBONE),
        AtomRole::HCα => (LEN_CALPHA_H as f32, 2, BOND_COLOR_BACKBONE),
        AtomRole::HN => (LEN_N_H as f32, 2, BOND_COLOR_BACKBONE),
        AtomRole::CSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
        AtomRole::OSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
        AtomRole::NSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
        AtomRole::HSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
        AtomRole::SSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
        AtomRole::SeSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
    };

    // todo: Sidechain mesh and color.
    // The bond
    entities.push(Entity::new(
        bond_mesh,
        vec3_to_f32(bond_center_position),
        quat_to_f32(bond_orientation),
        1.,
        bond_color,
        BOND_SHINYNESS,
    ));
}

/// Generates entities from protein coordinates.
/// todo: don't take both state and atoms_backbone, since ab is part of state.
// fn generate_entities(state: &State, atoms_backbone: &Vec<AtomCoords>) -> Vec<Entity> {
pub fn generate_entities(state: &State) -> Vec<Entity> {
    let mut result = Vec::new();

    // Store cα and c' so we can properly assign bonds after sidechains.
    let mut n_id = 0; // Residue 0, index 0 for first N
    let mut cα_id = 1; // Residue 0, index 1 for first Cα
    let mut cp_id = 0;

    // Atom id is used for station-keeping here.
    for (atom_id, atom) in state.protein.coords.atoms_backbone.iter().enumerate() {
        // let atom_color = if state.ui.active_residue == atom.residue_id {
        //     avg_colors(ACTIVE_COLOR_ATOM, atom.role.render_color())
        // } else {
        //     atom.role.render_color()
        // };

        let mut atom_color = atom.role.render_color();

        // todo, until we can find a better way to highlight the active atoms.
        // Highlight the active Nitrogen in an eye-catching color.
        if state.ui.active_residue == atom.residue_id {
            match atom.role {
                AtomRole::Cα => {
                    atom_color = ACTIVE_CALPHA_COLOR;
                }
                _ => (),
            }
        }

        let atom_mesh = match atom.role {
            AtomRole::N | AtomRole::Cα | AtomRole::Cp => 0, // Cube
            _ => 1,                                         // Sphere
        };

        let scale = match atom.role {
            AtomRole::HCα | AtomRole::HN | AtomRole::HSidechain => render::H_SCALE,
            _ => 1.,
        };

        if !state.ui.show_hydrogens {
            match atom.role {
                AtomRole::HCα | AtomRole::HN | AtomRole::HSidechain => continue,
                _ => (),
            }
        }

        //    CSidechain,
        //     OSidechain,
        //     NSidechain,
        //     SSidechain,
        //     SeSidechain,
        //     HSidechain,
        if !state.ui.show_sidechains {
            match atom.role {
                AtomRole::CSidechain
                | AtomRole::OSidechain
                | AtomRole::NSidechain
                | AtomRole::SSidechain
                | AtomRole::SeSidechain
                | AtomRole::HSidechain => continue,
                _ => (),
            }
        }

        // The atom
        result.push(Entity::new(
            atom_mesh,
            vec3_to_f32(atom.position),
            quat_to_f32(atom.orientation),
            scale,
            atom_color,
            ATOM_SHINYNESS,
        ));

        // Anchor N at position=0; we don't have a bond connected to it.
        if atom_id != 0 {
            // Find the previous atom in the chain: The one that connects to this.
            let atom_prev_id = match atom.role {
                AtomRole::N => {
                    n_id = atom_id;
                    cp_id
                }
                AtomRole::Cα => {
                    cα_id = atom_id;
                    n_id
                }
                AtomRole::CSidechain
                | AtomRole::OSidechain
                | AtomRole::NSidechain
                | AtomRole::SSidechain
                | AtomRole::SeSidechain
                | AtomRole::HSidechain => {
                    // This assumes the prev atom added before the sidechain was Cα.
                    atom_id - atom.sidechain_bond_step
                }
                AtomRole::Cp => {
                    cp_id = atom_id;
                    cα_id
                }
                AtomRole::O => cp_id,
                AtomRole::HCα => cα_id,
                AtomRole::HN => n_id,
            };

            add_bond(
                &atom,
                &state.protein.coords.atoms_backbone,
                atom_id,
                atom_prev_id,
                &mut result,
            );

            // todo: DRY between these calls.
            if let Some(second_bond_step) = atom.second_bond_step {
                let atom_prev_id = match atom.role {
                    AtomRole::N => {
                        n_id = atom_id;
                        cp_id
                    }
                    AtomRole::Cα => {
                        cα_id = atom_id;
                        n_id
                    }
                    AtomRole::CSidechain
                    | AtomRole::OSidechain
                    | AtomRole::NSidechain
                    | AtomRole::SSidechain
                    | AtomRole::SeSidechain
                    | AtomRole::HSidechain => {
                        // This assumes the prev atom added before the sidechain was Cα.
                        atom_id - second_bond_step
                    }
                    AtomRole::Cp => {
                        cp_id = atom_id;
                        cα_id
                    }
                    AtomRole::O => cp_id,
                    AtomRole::HCα => cα_id,
                    AtomRole::HN => n_id,
                };

                add_bond(
                    &atom,
                    &state.protein.coords.atoms_backbone,
                    atom_id,
                    atom_prev_id,
                    &mut result,
                );
            }
        }
    }

    // Genererate entities for water molecules.
    // todo: We shouldn't create this here; this should be stored somewhere else, perhaps.

    if state.ui.show_water_molecules {
        for water in &state.water_env.water_molecules {
            // Oxygen
            result.push(Entity::new(
                1, // sphere
                vec3_to_f32(water.o_posit),
                Q_I, // todo: Is this OK? Probably.
                1.,
                O_COLOR,
                ATOM_SHINYNESS,
            ));

            // H a
            result.push(Entity::new(
                1, // sphere
                vec3_to_f32(water.ha_posit),
                Q_I, // todo: Is this OK? Probably.
                H_SCALE,
                H_COLOR,
                ATOM_SHINYNESS,
            ));

            // H b
            result.push(Entity::new(
                1, // sphere
                vec3_to_f32(water.hb_posit),
                Q_I, // todo: Is this OK? Probably.
                H_SCALE,
                H_COLOR,
                ATOM_SHINYNESS,
            ));

            // M (artificial point)
            result.push(Entity::new(
                1, // sphere
                vec3_to_f32(water.m_posit),
                Q_I, // todo: Is this OK? Probably.
                M_SCALE,
                CALPHA_COLOR,
                ATOM_SHINYNESS,
            ));
            // todo: bonds.
        }
    }

    for nuc in &state.wavefunction_lab.nuclei {
        result.push(Entity::new(
            1,
            vec3_to_f32(nuc.position),
            Q_I,
            0.6,
            O_COLOR,
            ATOM_SHINYNESS,
        ));
    }

    for electron in &state.wavefunction_lab.electron_posits_dynamic {
        result.push(Entity::new(
            1,
            vec3_to_f32(*electron),
            Q_I,
            0.2,
            render::ELECTRON_COLOR,
            1.,
        ));
    }

    // Non-interacting electrons to better visualize the cloud.
    // Only show these when the sim is running, so we only see a single electron when stopped.
    if state.ui.sim_running {
        for electron in &state.wavefunction_lab.extra_visible_elecs_dynamic {
            result.push(Entity::new(
                1,
                util::vec3_to_f32(*electron),
                Q_I,
                0.2,
                render::ELECTRON_COLOR,
                1.,
            ));
        }
    }

    result
}

/// The entry point for our renderer.
pub fn run(state: State) {
    // Render our atoms.
    let entities = generate_entities(&state);

    let scene = Scene {
        meshes: vec![
            // Mesh::new_tetrahedron(render::SIDE_LEN),
            Mesh::new_box(
                render::SIDE_LEN * 0.7,
                render::SIDE_LEN * 0.7,
                render::SIDE_LEN * 0.7,
            ),
            Mesh::new_sphere(render::SIDE_LEN / 2., 20, 20),
            // todo: Temp bond len. You prob need a mesh per possible len.
            Mesh::new_cylinder(1.2, BOND_RADIUS_BACKBONE, render::BOND_N_SIDES),
            Mesh::new_cylinder(1.2, BOND_RADIUS_SIDECHAIN, render::BOND_N_SIDES),
        ],
        entities,
        camera: Camera {
            fov_y: TAU as f32 / 7.,
            position: Vec3F32::new(0., 0., -30.),
            far: RENDER_DIST,
            // orientation: QuatF32::from
            ..Default::default()
        },
        lighting: Lighting {
            ambient_color: [1., 1., 1., 0.5],
            ambient_intensity: 0.05,
            point_lights: vec![
                // Light from above. The sun?
                PointLight {
                    type_: LightType::Omnidirectional,
                    position: Vec3F32::new(0., 100., 0.),
                    diffuse_color: [0.6, 0.4, 0.3, 1.],
                    specular_color: [0.6, 0.4, 0.3, 1.],
                    diffuse_intensity: 15_000.,
                    specular_intensity: 15_000.,
                },
                // The light on the active residue. Moves
                PointLight {
                    type_: LightType::Omnidirectional,
                    position: Vec3F32::new_zero(),
                    diffuse_color: [0.3, 0.3, 0.3, 0.5],
                    specular_color: [1., 1., 0.7, 0.5],
                    diffuse_intensity: 100.,
                    specular_intensity: 100.,
                },
            ],
        },
        background_color: BACKGROUND_COLOR,
        window_size: (WINDOW_SIZE_X, WINDOW_SIZE_Y),
        window_title: WINDOW_TITLE.to_owned(),
        ..Default::default()
    };

    let input_settings = InputSettings {
        initial_controls: ControlScheme::FreeCamera,
        ..Default::default()
    };
    let ui_settings = UiSettings {
        // todo: How to handle this? For blocking keyboard and moues inputs when over the UI.
        // width: gui::UI_WIDTH as f64, // todo: Not working correctly.
        icon_path: Some("./resources/icon.png".to_owned()),
        ..Default::default()
    };

    // Of note, our event_handler and gui_handler functions could be used directly,
    // vice as closures.
    // Leaving them as closure-creators now for flexibility.

    graphics::run(
        state,
        scene,
        input_settings,
        ui_settings,
        Default::default(),
        render_handler,
        event_dev_handler,
        event_win_handler,
        gui::run(),
    );
}

// todo: Use state cam? Should it be in `Scene`?
