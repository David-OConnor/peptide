//! This module contains code for use with our custom renderer.

use crate::{
    atom_coords::{AtomCoords, ProteinCoords},
    chem_definitions::BackboneRole,
    kinematics::{ProteinDescription, LEN_N_CALPHA, LEN_CP_O, LEN_CP_N, LEN_CALPHA_CP},
    render::{self, BOND_RADIUS_SIDECHAIN, BOND_RADIUS_BACKBONE,
             BOND_COLOR_BACKBONE, BOND_COLOR_SIDECHAIN},
    sidechain::{LEN_SC},
};

use graphics::{self, lighting, Entity, InputSettings, Lighting, Mesh, Scene};

use lin_alg2::{self, f64::{Vec3, Quaternion}};

// The length-wise axis of our graphics engine's cylinder mesh.
const BOND_MODEL_AXIS: Vec3 = Vec3 { x: 0., y: 1., z: 0. };

fn vec3_to_f32(v: Vec3) -> lin_alg2::f32::Vec3 {
    lin_alg2::f32::Vec3::new(v.x as f32, v.y as f32, v.z as f32)
}

fn quat_to_f32(q: Quaternion) -> lin_alg2::f32::Quaternion {
    lin_alg2::f32::Quaternion::new(q.w as f32, q.x as f32, q.y as f32, q.z as f32)
}

/// The entry point for our renderer.
pub fn run() {
    let mut state = crate::State::default();

    // Initialize the render state here.
    state.protein_descrip = crate::init_protein();
    state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);

    let mut entities = Vec::new();

    // Store cα and c' so we can properly assign bonds after sidechains.
    let mut n_id = 0; // Residue 0, index 0 for first N
    let mut cα_id = 1; // Residue 0, index 1 for first Cα
    let mut cp_id = 0;

    // Render our atoms.
    for (id, atom) in state.protein_coords.atoms_backbone.iter().enumerate() {

        entities.push(Entity::new(
            0,
            vec3_to_f32(atom.position),
            quat_to_f32(atom.orientation),
            1.,
            atom.role.render_color(),
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

            let atom_prev = &state.protein_coords.atoms_backbone[atom_prev_id];

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
                BackboneRole::N => (LEN_CP_N as f32, 1),
                BackboneRole::Cα => (LEN_N_CALPHA as f32, 1),
                BackboneRole::Cp => (LEN_CALPHA_CP as f32, 1),
                BackboneRole::O => (LEN_CP_O as f32, 1),
                BackboneRole::CSidechain => (LEN_SC as f32, 2),
                BackboneRole::OSidechain => (LEN_SC as f32, 2),
                BackboneRole::NSidechain => (LEN_SC as f32, 2),
            };

            // todo: Sidechain mesh and color.
            entities.push(Entity::new(
                bond_mesh,
                vec3_to_f32(bond_center_position),
                quat_to_f32(bond_orientation),
                1.,
                render::BOND_COLOR_BACKBONE,
            ));
        }
    }

    let scene = Scene {
        meshes: vec![
            Mesh::new_tetrahedron(render::SIDE_LEN),
            // todo: Temp bond len. You prob need a mesh per possible len.
            Mesh::new_cylinder(1.0, render::BOND_RADIUS_BACKBONE, render::BOND_N_SIDES),
            Mesh::new_cylinder(1.0, render::BOND_RADIUS_SIDECHAIN, render::BOND_N_SIDES),
        ],
        entities,
        lighting: Lighting::default(),
    };

    let input_settings = InputSettings::default();

    graphics::run(scene, input_settings);
}

// todo: Use state cam? Should it be in `Scene`?
