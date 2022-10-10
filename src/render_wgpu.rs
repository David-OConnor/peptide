//! This module contains code for use with our custom renderer.

use std::{
    f64::consts::TAU,
    // boxed::Box,
    mem,
    rc::Rc,
    sync::Mutex,
};

use graphics::{
    self, Camera, DeviceEvent, ElementState, Entity, InputSettings, Lighting, Mesh, Scene,
    UiSettings,
};
use lin_alg2::{
    self,
    f64::{Quaternion, Vec3},
};

use crate::types::{ProteinDescription, State};
use crate::{
    atom_coords::{AtomCoords, ProteinCoords},
    bond_vecs::{LEN_CALPHA_CP, LEN_CALPHA_H, LEN_CP_N, LEN_CP_O, LEN_N_CALPHA, LEN_N_H},
    chem_definitions::BackboneRole,
    gui,
    render::{
        self, ACTIVE_COLOR_ATOM, ATOM_SHINYNESS, BOND_COLOR_BACKBONE, BOND_COLOR_SIDECHAIN,
        BOND_RADIUS_BACKBONE, BOND_RADIUS_SIDECHAIN, BOND_SHINYNESS,
    },
    sidechain::LEN_SC,
};

// todo: Use a mutex instead of static mut here.
// todo: pub for access in `gui`.
pub static mut STATE: State = unsafe {
    State {
        protein_descrip: ProteinDescription {
            name: String::new(),
            pdb_ident: String::new(),
            residues: Vec::new(),
        },
        protein_coords: ProteinCoords {
            atoms_backbone: Vec::new(),
        },
        active_residue: 1,
        cam: render::Camera {
            position: Vec3 {
                x: 0.,
                y: 0.,
                z: 0.,
            },
            orientation: Quaternion {
                w: 0.,
                x: 0.,
                y: 0.,
                z: 0.,
            },
        },
    }
};

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

// todo: Don't make dependency import winit? rexport etc?
fn device_event_handler(
    event: DeviceEvent,
    scene: &mut Scene,
    dt: f32, // in seconds
) -> bool {
    // todo: Higher level api from winit or otherwise instead of scancode?
    let mut coords_changed = false;
    let mut active_res_backbone_changed = false;
    let mut active_res_sidechain_changed = false;
    let mut active_res_changed = false;

    // We count starting at 1, per chem conventions.
    let ar_i = unsafe { STATE.active_residue } - 1;
    // code shortener
    let mut active_res = unsafe { &mut STATE.protein_descrip.residues[ar_i] };

    let rotation_amt = crate::BOND_ROTATION_SPEED * dt as f64;

    unsafe {
        match event {
            DeviceEvent::Key(key) => {
                if key.state == ElementState::Pressed {
                    // todo: These should probably be every 10 residues.
                    match key.scancode {
                        2 => {
                            STATE.active_residue = 1;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        3 => {
                            STATE.active_residue = 2;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        4 => {
                            STATE.active_residue = 3;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        5 => {
                            STATE.active_residue = 4;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        6 => {
                            STATE.active_residue = 5;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        7 => {
                            STATE.active_residue = 6;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        8 => {
                            STATE.active_residue = 7;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        9 => {
                            STATE.active_residue = 8;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        10 => {
                            STATE.active_residue = 9;
                            coords_changed = true;
                            active_res_changed = true;
                        }
                        // todo: Why are these scan codes for up/down so high??
                        57_416 => {
                            // Up arrow
                            if STATE.active_residue != STATE.protein_descrip.residues.len() {
                                STATE.active_residue += 1;
                                coords_changed = true;
                                active_res_changed = true;
                            }
                        }
                        57_424 => {
                            // Down arrow
                            if STATE.active_residue != 1 {
                                STATE.active_residue -= 1;
                                coords_changed = true;
                                active_res_changed = true;
                            }
                        }
                        20 => {
                            // T
                            active_res.φ += rotation_amt;
                            active_res_backbone_changed = true;
                            coords_changed = true;
                        }
                        34 => {
                            // G
                            active_res.φ -= rotation_amt;
                            active_res_backbone_changed = true;
                            coords_changed = true;
                        }
                        21 => {
                            // Y
                            active_res.ψ += rotation_amt;
                            active_res_backbone_changed = true;
                            coords_changed = true;
                        }
                        35 => {
                            // H
                            active_res.ψ -= rotation_amt;
                            active_res_backbone_changed = true;
                            coords_changed = true;
                        }
                        22 => {
                            // U
                            active_res.ω += rotation_amt;
                            active_res_backbone_changed = true;
                            coords_changed = true;
                        }
                        36 => {
                            // J
                            active_res.ω -= rotation_amt;
                            active_res_backbone_changed = true;
                            coords_changed = true;
                        }
                        23 => {
                            // I
                            active_res.sidechain.add_to_χ1(rotation_amt);
                            active_res_sidechain_changed = true;
                            coords_changed = true;
                        }
                        37 => {
                            // K
                            active_res.sidechain.add_to_χ1(-rotation_amt);
                            active_res_sidechain_changed = true;
                            coords_changed = true;
                        }
                        24 => {
                            // O
                            active_res.sidechain.add_to_χ2(rotation_amt);
                            active_res_sidechain_changed = true;
                            coords_changed = true;
                        }
                        38 => {
                            // L
                            active_res.sidechain.add_to_χ2(-rotation_amt);
                            active_res_sidechain_changed = true;
                            coords_changed = true;
                        }
                        36 => {
                            // P
                            active_res.sidechain.add_to_χ3(rotation_amt);
                            active_res_sidechain_changed = true;
                            coords_changed = true;
                        }
                        39 => {
                            // ;
                            active_res.sidechain.add_to_χ3(-rotation_amt);
                            active_res_sidechain_changed = true;
                            coords_changed = true;
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
    }

    unsafe {
        if coords_changed {
            // Recalculate coordinates now that we've updated our bond angles
            STATE.protein_coords = ProteinCoords::from_descrip(&STATE.protein_descrip);
            scene.entities = generate_entities(&STATE.protein_coords.atoms_backbone);
        }

        if active_res_changed {
            gui::ACTIVE_RES_ID = STATE.active_residue;

            // let aa_name = format!("{}", state.protein_descrip.residues[state.active_residue].sidechain);
            // let aa_name =
            gui::ACTIVE_RES_AA_NAME = &STATE.protein_descrip.residues[STATE.active_residue - 1]
                .sidechain
                .aa_name();
        }

        if active_res_changed || active_res_backbone_changed {
            // todo: Break this out by psi, phi etc instead of always updating all?
            let res = &STATE.protein_descrip.residues[STATE.active_residue - 1];

            gui::ACTIVE_RES_PSI = res.φ;
            gui::ACTIVE_RES_PHI = res.ψ;
            gui::ACTIVE_RES_OMEGA = res.ω;

            // todo: Only do this for backbone changd; not acive res
            // Note: We use modulus here to make integrating with the GUI
            // easier, ie clamping the range between 0 and TAU.

            // Clamp to > 0
            if active_res.φ < 0. {
                active_res.φ += TAU;
            }
            if active_res.ψ < 0. {
                active_res.ψ += TAU;
            }
            if active_res.ω < 0. {
                active_res.ω += TAU;
            }

            // Clamp to < TAU
            active_res.φ = (active_res.φ) % TAU;
            active_res.ψ = (active_res.ψ) % TAU;
            active_res.ω = (active_res.ω) % TAU;
        }

        if active_res_changed || active_res_sidechain_changed {
            // todo: Break this out by psi, phi etc instead of always updating all?
            let sc = &STATE.protein_descrip.residues[STATE.active_residue - 1].sidechain;

            gui::ACTIVE_RES_XI_1 = sc.get_χ1();
            gui::ACTIVE_RES_XI_2 = sc.get_χ2();
            gui::ACTIVE_RES_XI_3 = sc.get_χ3();
            gui::ACTIVE_RES_XI_4 = sc.get_χ4();
            gui::ACTIVE_RES_XI_5 = sc.get_χ5();
        }
    }

    coords_changed
}

fn render_handler(scene: &mut Scene) -> bool {
    // todo: This may be where you need to update the render after changing a slider
    false
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
fn generate_entities(atoms_backbone: &Vec<AtomCoords>) -> Vec<Entity> {
    let mut result = Vec::new();

    // Store cα and c' so we can properly assign bonds after sidechains.
    let mut n_id = 0; // Residue 0, index 0 for first N
    let mut cα_id = 1; // Residue 0, index 1 for first Cα
    let mut cp_id = 0;

    for (id, atom) in atoms_backbone.iter().enumerate() {
        let atom_color = if unsafe { STATE.active_residue } == atom.residue_id {
            avg_colors(ACTIVE_COLOR_ATOM, atom.role.render_color())
        } else {
            atom.role.render_color()
        };

        let atom_mesh = match atom.role {
            BackboneRole::N | BackboneRole::Cα | BackboneRole::Cp | BackboneRole::O => 0,
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
                BackboneRole::N => {
                    n_id = id;
                    cp_id
                }
                BackboneRole::Cα => {
                    cα_id = id;
                    n_id
                }
                BackboneRole::CSidechain | BackboneRole::OSidechain | BackboneRole::NSidechain => {
                    // This assumes the prev atom added before the sidechain was Cα.
                    id - atom.sidechain_bond_step
                }
                BackboneRole::Cp => {
                    cp_id = id;
                    cα_id
                }
                BackboneRole::O => cp_id,
                BackboneRole::HCα => cα_id,
                BackboneRole::HN => n_id,
            };

            // Calculate the position of the bond mesh: This is the cylinder's z point,
            // half way between the 2 atoms it connects.

            let atom_prev = unsafe { &STATE.protein_coords.atoms_backbone[atom_prev_id] };

            let bond_center_position = (atom.position + atom_prev.position) * 0.5;
            let bond_dir = (atom.position - atom_prev.position).to_normalized();

            let bond_orientation = Quaternion::from_unit_vecs(BOND_MODEL_AXIS, bond_dir);

            let color_a = match atom.role {
                BackboneRole::CSidechain | BackboneRole::OSidechain | BackboneRole::NSidechain => {
                    BOND_COLOR_SIDECHAIN
                }
                _ => BOND_COLOR_BACKBONE,
            };

            // let color = Color::rgb(color_a.0, color_a.1, color_a.2).into();

            let (bond_len, bond_mesh, bond_color) = match atom.role {
                BackboneRole::N => (LEN_CP_N as f32, 2, BOND_COLOR_BACKBONE),
                BackboneRole::Cα => (LEN_N_CALPHA as f32, 2, BOND_COLOR_BACKBONE),
                BackboneRole::Cp => (LEN_CALPHA_CP as f32, 2, BOND_COLOR_BACKBONE),
                BackboneRole::O => (LEN_CP_O as f32, 2, BOND_COLOR_BACKBONE),
                BackboneRole::HCα => (LEN_CALPHA_H as f32, 2, BOND_COLOR_BACKBONE),
                BackboneRole::HN => (LEN_N_H as f32, 2, BOND_COLOR_BACKBONE),
                BackboneRole::CSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
                BackboneRole::OSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
                BackboneRole::NSidechain => (LEN_SC as f32, 3, BOND_COLOR_SIDECHAIN),
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
pub unsafe fn run() {
    // let mut state = State::default();
    // let state = Rc::new(Mutex::new(State::default());
    // todo: Initialize our static mut; bit of a hack
    STATE = State::default();

    // Initialize the render state here.

    STATE.protein_descrip = crate::init_protein();
    STATE.protein_coords = ProteinCoords::from_descrip(&STATE.protein_descrip);

    // let state = gui::init_statics(&state); // todo: lifetime issues with this

    gui::PROT_NAME = &STATE.protein_descrip.name;
    gui::PDB_IDENT = &STATE.protein_descrip.pdb_ident;

    gui::ACTIVE_RES_ID = STATE.active_residue;

    let res = &STATE.protein_descrip.residues[STATE.active_residue - 1];

    gui::ACTIVE_RES_AA_NAME = res.sidechain.aa_name();

    gui::ACTIVE_RES_PSI = res.φ;
    gui::ACTIVE_RES_PHI = res.ψ;
    gui::ACTIVE_RES_OMEGA = res.ω;

    gui::ACTIVE_RES_XI_1 = res.sidechain.get_χ1();
    gui::ACTIVE_RES_XI_2 = res.sidechain.get_χ2();
    gui::ACTIVE_RES_XI_3 = res.sidechain.get_χ3();
    gui::ACTIVE_RES_XI_4 = res.sidechain.get_χ4();
    gui::ACTIVE_RES_XI_5 = res.sidechain.get_χ5();

    // Render our atoms.
    let entities = generate_entities(&STATE.protein_coords.atoms_backbone);

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
        camera: Camera::default(),
        lighting: Lighting::default(),
        background_color: render::BACKGROUND_COLOR,
        window_size: (900., 600.),
        window_title: "Peptide".to_owned(),
    };

    let input_settings = InputSettings::default();
    let ui_settings = UiSettings::default();

    // let deh = |device_event, scene: &mut _| {
    //     device_event_handler(device_event, &mut state, &mut scene);
    //     println!("TEST");
    //     None
    // };

    graphics::run(
        scene,
        input_settings,
        ui_settings,
        Box::new(render_handler),
        Box::new(device_event_handler),
        Box::new(gui::draw_ui),
    );
}

// todo: Use state cam? Should it be in `Scene`?
