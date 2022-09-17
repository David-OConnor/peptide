//! This module contains code for use with our custom renderer.

use std::{
    // boxed::Box,
    mem,
    rc::Rc,
    sync::Mutex,
};

use graphics::{
    self, lighting, Camera, DeviceEvent, ElementState, Entity, InputSettings, Lighting, Mesh, Scene,
};

use lin_alg2::{
    self,
    f64::{Quaternion, Vec3},
};

use crate::{
    atom_coords::{AtomCoords, ProteinCoords},
    bond_vecs::{LEN_CALPHA_CP, LEN_CP_N, LEN_CP_O, LEN_N_CALPHA},
    chem_definitions::BackboneRole,
    kinematics::ProteinDescription,
    render::{
        self, ACTIVE_COLOR_ATOM, ATOM_SHINYNESS, BOND_COLOR_BACKBONE, BOND_COLOR_SIDECHAIN,
        BOND_RADIUS_BACKBONE, BOND_RADIUS_SIDECHAIN, BOND_SHINYNESS,
    },
    sidechain::LEN_SC,
    State,
};

// todo: Temp until we sort out how to use Fn traits vice function pointers.
// static mut state: State = unsafe { mem::zeroed() };

// todo: Use a mutex instead of static mut here.
static mut state: State = unsafe {
    State {
        protein_descrip: ProteinDescription {
            name: String::new(),
            pdb_ident: String::new(),
            residues: Vec::new(),
        },
        protein_coords: ProteinCoords {
            atoms_backbone: Vec::new(),
        },
        active_residue: 0,
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
    // state: &mut State,
    scene: &mut Scene,
    dt: f32, // in seconds
) -> bool {
    // todo: Higher level api from winit or otherwise instead of scancode?
    let mut changed = false;

    // We count starting at 1, per chem conventions.
    let ar_i = unsafe { state.active_residue } - 1;
    // code shortener
    let mut active_res = unsafe { &mut state.protein_descrip.residues[ar_i] };

    let mut rotation_amt = crate::BOND_ROTATION_SPEED * dt as f64;

    unsafe {
        match event {
            DeviceEvent::Key(key) => {
                if key.state == ElementState::Pressed {
                    // todo: These should probably be every 10 residues.
                    match key.scancode {
                        2 => {
                            state.active_residue = 1;
                            changed = true;
                        }
                        3 => {
                            state.active_residue = 2;
                            changed = true;
                        }
                        4 => {
                            state.active_residue = 3;
                            changed = true;
                        }
                        5 => {
                            state.active_residue = 4;
                            changed = true;
                        }
                        6 => {
                            state.active_residue = 5;
                            changed = true;
                        }
                        7 => {
                            state.active_residue = 6;
                            changed = true;
                        }
                        8 => {
                            state.active_residue = 7;
                            changed = true;
                        }
                        9 => {
                            state.active_residue = 8;
                            changed = true;
                        }
                        10 => {
                            state.active_residue = 9;
                            changed = true;
                        }
                        // todo: Why are these scan codes for up/down so high??
                        57_416 => {
                            // Up arrow
                            if state.active_residue != state.protein_descrip.residues.len() {
                                state.active_residue += 1;
                            }
                            changed = true;
                        }
                        57_424 => {
                            // Down arrow
                            if state.active_residue != 1 {
                                state.active_residue -= 1;
                            }
                            changed = true;
                        }
                        20 => {
                            // T
                            active_res.φ += rotation_amt;
                            changed = true;
                        }
                        34 => {
                            // G
                            active_res.φ -= rotation_amt;
                            changed = true;
                        }
                        21 => {
                            // Y
                            active_res.ψ += rotation_amt;
                            changed = true;
                        }
                        35 => {
                            // H
                            active_res.ψ -= rotation_amt;
                            changed = true;
                        }
                        22 => {
                            // U
                            active_res.ω += rotation_amt;
                            changed = true;
                        }
                        36 => {
                            // J
                            active_res.ω -= rotation_amt;
                            changed = true;
                        }
                        23 => {
                            // I
                            active_res.sidechain.add_to_χ1(rotation_amt);
                            changed = true;
                        }
                        37 => {
                            // K
                            active_res.sidechain.add_to_χ1(-rotation_amt);
                            changed = true;
                        }
                        24 => {
                            // O
                            active_res.sidechain.add_to_χ2(rotation_amt);
                            changed = true;
                        }
                        38 => {
                            // L
                            active_res.sidechain.add_to_χ2(-rotation_amt);
                            changed = true;
                        }
                        36 => {
                            // P
                            active_res.sidechain.add_to_χ3(rotation_amt);
                            changed = true;
                        }
                        39 => {
                            // ;
                            active_res.sidechain.add_to_χ3(-rotation_amt);
                            changed = true;
                        }
                        29 => {
                            // Left ctrl
                            println!("Active res id: {}\n{}", state.active_residue, active_res);
                        }
                        // todo: Sidechain dihedral angles
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }
    if changed {
        // Recalculate coordinates now that we've updated our bond angles
        unsafe {
            state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
            scene.entities = generate_entities(&state.protein_coords.atoms_backbone);
        }
    }

    changed
}

fn render_handler() -> Option<Vec<Entity>> {
    // fn render_handler(scene: &mut Scene) -> Option<Vec<Entity>> {
    None
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
        let atom_color = if unsafe { state.active_residue } == atom.residue_id {
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
            };

            // Calculate the position of the bond mesh: This is the cylinder's z point,
            // half way between the 2 atoms it connects.

            let atom_prev = unsafe { &state.protein_coords.atoms_backbone[atom_prev_id] };

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

            let (bond_len, bond_mesh) = match atom.role {
                BackboneRole::N => (LEN_CP_N as f32, 2),
                BackboneRole::Cα => (LEN_N_CALPHA as f32, 2),
                BackboneRole::Cp => (LEN_CALPHA_CP as f32, 2),
                BackboneRole::O => (LEN_CP_O as f32, 2),
                BackboneRole::CSidechain => (LEN_SC as f32, 3),
                BackboneRole::OSidechain => (LEN_SC as f32, 3),
                BackboneRole::NSidechain => (LEN_SC as f32, 3),
            };

            // todo: Sidechain mesh and color.
            // The bond
            result.push(Entity::new(
                bond_mesh,
                vec3_to_f32(bond_center_position),
                quat_to_f32(bond_orientation),
                1.,
                render::BOND_COLOR_BACKBONE,
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
    state = State::default();

    // Initialize the render state here.

    state.protein_descrip = crate::init_protein();
    state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);

    // Render our atoms.
    let entities = generate_entities(&state.protein_coords.atoms_backbone);

    let mut scene = Scene {
        meshes: vec![
            // Mesh::new_tetrahedron(render::SIDE_LEN),
            Mesh::new_box(render::SIDE_LEN, render::SIDE_LEN, render::SIDE_LEN),
            Mesh::new_sphere(render::SIDE_LEN / 2., 20, 20),
            // todo: Temp bond len. You prob need a mesh per possible len.
            Mesh::new_cylinder(1.2, render::BOND_RADIUS_BACKBONE, render::BOND_N_SIDES),
            Mesh::new_cylinder(1.2, render::BOND_RADIUS_SIDECHAIN, render::BOND_N_SIDES),
        ],
        entities,
        camera: Camera::default(),
        lighting: Lighting::default(),
    };

    let input_settings = InputSettings::default();

    // let deh = |device_event, scene: &mut _| {
    //     device_event_handler(device_event, &mut state, &mut scene);
    //     println!("TEST");
    //     None
    // };

    graphics::run(scene, input_settings, render_handler, device_event_handler);
    // graphics::run(scene, input_settings, render_handler, Box::new(device_event_handler));
}

// todo: Use state cam? Should it be in `Scene`?
