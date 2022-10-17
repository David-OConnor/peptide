//! This module contains code for use with our custom renderer.

use std::f64::consts::TAU;

use rand;

use graphics::{
    self, Camera, ControlScheme, DeviceEvent, ElementState, EngineUpdates, Entity, InputSettings,
    LightType, Lighting, Mesh, PointLight, Scene, UiSettings,
};

use lin_alg2::{
    self,
    f32::Quaternion as QuatF32,
    f32::Vec3 as Vec3F32,
    f64::{Quaternion, Vec3},
};

use crate::{
    atom_coords::ProteinCoords,
    bond_vecs::{LEN_CALPHA_CP, LEN_CALPHA_H, LEN_CP_N, LEN_CP_O, LEN_N_CALPHA, LEN_N_H},
    chem_definitions::AtomRole,
    gui,
    render::{
        self, ACTIVE_COLOR_ATOM, ATOM_SHINYNESS, BOND_COLOR_BACKBONE, BOND_COLOR_SIDECHAIN,
        BOND_RADIUS_BACKBONE, BOND_RADIUS_SIDECHAIN, BOND_SHINYNESS,
    },
    sidechain::LEN_SC,
    types::State,
};

const WINDOW_TITLE: &str = "Peptide";
const WINDOW_SIZE_X: f32 = 900.0;
const WINDOW_SIZE_Y: f32 = 600.0;

// Changes the far end of the frustrum; try to have this shortly past the farthest
// view distance you expect. Possible to make this dynamic?
const RENDER_DIST: f32 = 200.;

// The length-wise axis of our graphics engine's cylinder mesh.
const BOND_MODEL_AXIS: Vec3 = Vec3 {
    x: 0.,
    y: 1.,
    z: 0.,
};

fn vec3_to_f32(v: Vec3) -> lin_alg2::f32::Vec3 {
    lin_alg2::f32::Vec3::new(v.x as f32, v.y as f32, v.z as f32)
}

fn quat_to_f32(q: Quaternion) -> lin_alg2::f32::Quaternion {
    lin_alg2::f32::Quaternion::new(q.w as f32, q.x as f32, q.y as f32, q.z as f32)
}

/// Adjust the camera to focus on a new point. Does so by rotating in the shortest direction
/// to point at the new point, then moving forward or back to get to the requested distance.
pub fn look_at(cam: &mut Camera, focus_pt: Vec3F32, dist: f32) {
    // todo: A smooth animation later. For now, just calculate the resulting camera
    // todo position and orientation.

    // todo: Add this to your quaternion article.

    let vec_currently_looking_at = cam.orientation.rotate_vec(Vec3F32::new(0., 0., 1.)); // fwd vec
    let dir_to_pt = (focus_pt - cam.position).to_normalized();

    let rotation = QuatF32::from_unit_vecs(vec_currently_looking_at, dir_to_pt);

    cam.orientation = rotation * cam.orientation;

    // todo: Why so slow?

    // Now that the orientation is set, move forward or backwards in relation
    // to the new focus point until at the target distance.
    cam.position = focus_pt - dir_to_pt * dist;
}

fn make_event_handler() -> impl FnMut(&mut State, DeviceEvent, &mut Scene, f32) -> EngineUpdates {
    // Box::new(move |event: DeviceEvent, scene: &mut Scene, dt: f32| {
    move |state: &mut State, event: DeviceEvent, scene: &mut Scene, dt: f32| {
        // todo: Higher level api from winit or otherwise instead of scancode?
        let mut entities_changed = false;
        let mut active_res_backbone_changed = false;
        let mut active_res_sidechain_changed = false;
        let mut active_res_changed = false;
        let mut lighting_changed = false;

        // We count starting at 1, per chem conventions.
        let ar_i = state.active_residue - 1;
        // code shortener

        let mut active_res = &mut state.protein_descrip.residues[ar_i];
        let rotation_amt = crate::BOND_ROTATION_SPEED * dt as f64;

        match event {
            DeviceEvent::Key(key) => {
                if key.state == ElementState::Pressed {
                    // todo: These should probably be every 10 residues.
                    match key.scancode {
                        2 => {
                            state.active_residue = 1;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        3 => {
                            state.active_residue = 2;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        4 => {
                            state.active_residue = 3;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        5 => {
                            state.active_residue = 4;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        6 => {
                            state.active_residue = 5;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        7 => {
                            state.active_residue = 6;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        8 => {
                            state.active_residue = 7;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        9 => {
                            state.active_residue = 8;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        10 => {
                            state.active_residue = 9;
                            entities_changed = true;
                            active_res_changed = true;
                        }
                        // todo: Why are these scan codes for up/down so high??
                        57_416 => {
                            // Up arrow
                            if state.active_residue != state.protein_descrip.residues.len() {
                                state.active_residue += 1;
                                entities_changed = true;
                                active_res_changed = true;
                            }
                        }
                        57_424 => {
                            // Down arrow
                            if state.active_residue != 1 {
                                state.active_residue -= 1;
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
            _ => {}
        }

        if entities_changed {
            // Recalculate coordinates now that we've updated our bond angles
            state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
            scene.entities = generate_entities(&state);
        }

        if active_res_changed {
            // state.ui.active_res_id = state.active_residue;
            //
            // // let aa_name = format!("{}", state.protein_descrip.residues[state.active_residue].sidechain);
            // // let aa_name =
            // state.ui.active_res_aa_name = state.protein_descrip.residues[ar_i]
            //     .sidechain
            //     .aa_name()
            //     .to_owned();

            gui::change_lit_res(state, scene);
            lighting_changed = true;
        }

        if active_res_changed || active_res_backbone_changed {
            // todo: Break this out by psi, phi etc instead of always updating all?
            let mut res = &mut state.protein_descrip.residues[state.active_residue - 1];

            // todo: Only do this for backbone changd; not acive res

            crate::clamp_angle(&mut res.φ);
            crate::clamp_angle(&mut res.ψ);
            crate::clamp_angle(&mut res.ω);
        }

        // if active_res_changed || active_res_sidechain_changed {
        //     // todo: Break this out by psi, phi etc instead of always updating all?
        //     // let sc = &state.protein_descrip.residues[state.active_residue - 1].sidechain;
        //     //
        //     // state.ui.active_res_χ1 = sc.get_χ1();
        //     // state.ui.active_res_χ2 = sc.get_χ2();
        //     // state.ui.active_res_χ3 = sc.get_χ3();
        //     // state.ui.active_res_χ4 = sc.get_χ4();
        //     // state.ui.active_res_χ5 = sc.get_χ5();
        // }

        EngineUpdates {
            entities: entities_changed,
            camera: false,
            lighting: lighting_changed,
        }
    }
}

// todo: Util?
/// Get a random value from -0.5 to 0.5
fn rng() -> f64 {
    rand::random::<f64>() - 0.5
}

fn render_handler(state: &mut State, scene: &mut Scene, dt: f32) -> EngineUpdates {
    let mut entities_changed = false;

    if state.sim_running {
        // Constant term in our noise.
        let c = state.temperature * state.sim_time_scale * dt as f64;
        // Crude approach, where we add random noise to the angles.

        for res in &mut state.protein_descrip.residues {
            res.ψ += rng() * c;
            res.φ += rng() * c;

            if let Some(χ) = res.sidechain.get_mut_χ1() {
                *χ += rng() * c;
            }
            if let Some(χ) = res.sidechain.get_mut_χ2() {
                *χ += rng() * c;
            }
            if let Some(χ) = res.sidechain.get_mut_χ3() {
                *χ += rng() * c;
            }
            if let Some(χ) = res.sidechain.get_mut_χ4() {
                *χ += rng() * c;
            }
            if let Some(χ) = res.sidechain.get_mut_χ5() {
                *χ += rng() * c;
            }
        }

        state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
        scene.entities = generate_entities(&state);

        entities_changed = true;
    }

    EngineUpdates {
        entities: entities_changed,
        camera: false,
        lighting: false,
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

/// Generates entities from protein coordinates.
/// todo: don't take both state and atoms_backbone, since ab is part of state.
// fn generate_entities(state: &State, atoms_backbone: &Vec<AtomCoords>) -> Vec<Entity> {
pub fn generate_entities(state: &State) -> Vec<Entity> {
    let mut result = Vec::new();

    // Store cα and c' so we can properly assign bonds after sidechains.
    let mut n_id = 0; // Residue 0, index 0 for first N
    let mut cα_id = 1; // Residue 0, index 1 for first Cα
    let mut cp_id = 0;

    for (id, atom) in state.protein_coords.atoms_backbone.iter().enumerate() {
        let atom_color = if state.active_residue == atom.residue_id {
            avg_colors(ACTIVE_COLOR_ATOM, atom.role.render_color())
        } else {
            atom.role.render_color()
        };

        let atom_mesh = match atom.role {
            AtomRole::N | AtomRole::Cα | AtomRole::Cp | AtomRole::O => 0,
            _ => 1,
        };

        // The atom
        result.push(Entity::new(
            atom_mesh,
            vec3_to_f32(atom.position),
            quat_to_f32(atom.orientation),
            1.,
            atom_color,
            ATOM_SHINYNESS,
        ));

        // Anchor N at position=0; we don't have a bond connected to it.
        if id != 0 {
            // Find the previous atom in the chain: The one that connects to this.
            let atom_prev_id = match atom.role {
                AtomRole::N => {
                    n_id = id;
                    cp_id
                }
                AtomRole::Cα => {
                    cα_id = id;
                    n_id
                }
                AtomRole::CSidechain
                | AtomRole::OSidechain
                | AtomRole::NSidechain
                | AtomRole::SSidechain
                | AtomRole::SeSidechain => {
                    // This assumes the prev atom added before the sidechain was Cα.
                    id - atom.sidechain_bond_step
                }
                AtomRole::Cp => {
                    cp_id = id;
                    cα_id
                }
                AtomRole::O => cp_id,
                AtomRole::HCα => cα_id,
                AtomRole::HN => n_id,
            };

            // Calculate the position of the bond mesh: This is the cylinder's z point,
            // half way between the 2 atoms it connects.

            let atom_prev = &state.protein_coords.atoms_backbone[atom_prev_id];

            let bond_center_position = (atom.position + atom_prev.position) * 0.5;
            let bond_dir = (atom.position - atom_prev.position).to_normalized();

            let bond_orientation = Quaternion::from_unit_vecs(BOND_MODEL_AXIS, bond_dir);

            let color_a = match atom.role {
                AtomRole::CSidechain | AtomRole::OSidechain | AtomRole::NSidechain => {
                    BOND_COLOR_SIDECHAIN
                }
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
                AtomRole::SSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
                AtomRole::SeSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
            };

            // todo: Sidechain mesh and color.
            // The bond
            result.push(Entity::new(
                bond_mesh,
                vec3_to_f32(bond_center_position),
                quat_to_f32(bond_orientation),
                1.,
                bond_color,
                BOND_SHINYNESS,
            ));
        }
    }

    result
}

/// The entry point for our renderer.
pub fn run(mut state: State) {
    let res = &state.protein_descrip.residues[state.active_residue - 1];

    // Initialize the GUI state here.
    // state.ui = gui::StateUi {
    //     // prot_name: state.protein_descrip.name.clone(),
    //     // pdb_ident: state.protein_descrip.pdb_ident.clone(),
    //     // active_res_id: state.active_residue,
    //     // active_res_aa_name: res.sidechain.aa_name().to_owned(),
    //     // active_res_ψ: res.ψ,
    //     // active_res_φ: res.φ,
    //     // active_res_ω: res.ω,
    //     active_res_χ1: res.sidechain.get_χ1(),
    //     active_res_χ2: res.sidechain.get_χ2(),
    //     active_res_χ3: res.sidechain.get_χ3(),
    //     active_res_χ4: res.sidechain.get_χ4(),
    //     active_res_χ5: res.sidechain.get_χ5(),
    // };

    // Render our atoms.
    // let entities = generate_entities(&state, &state.protein_coords.atoms_backbone);
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
            Mesh::new_cylinder(1.2, render::BOND_RADIUS_BACKBONE, render::BOND_N_SIDES),
            Mesh::new_cylinder(1.2, render::BOND_RADIUS_SIDECHAIN, render::BOND_N_SIDES),
        ],
        entities,
        camera: Camera {
            position: Vec3F32::new(0., 0., -30.),
            far: RENDER_DIST,
            // orientation: QuatF32::from
            ..Default::default()
        },
        lighting: Lighting {
            ambient_color: [1., 1., 1., 0.5],
            ambient_intensity: 0.05,
            // Light from above. The sun?
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
        background_color: render::BACKGROUND_COLOR,
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
        width: 500.,
    };

    // Of note, these functions could be used directly, vice as closures.
    // Leaving them as closure-creators now for flexibility.
    let event_handler = make_event_handler();
    let gui_handler = gui::run();

    graphics::run(
        state,
        scene,
        input_settings,
        ui_settings,
        render_handler,
        event_handler,
        gui_handler,
    );
}

// todo: Use state cam? Should it be in `Scene`?
