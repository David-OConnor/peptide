use crate::{
    atom_coords::{AtomCoords, ProteinCoords},
    kinematics::{ProteinDescription},
    render,
};

use graphics::{self, lighting, Entity, InputSettings, Lighting, Mesh, Scene};

use lin_alg2;

/// The entry point for our renderer.
pub fn run() {
    let mut state = crate::State::default();

    // Initialize the render state here.
    state.protein_descrip = crate::init_protein();
    state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);

    let mut entities = Vec::new();

    // Render our atoms.
    for (id, atom) in state.protein_coords.atoms_backbone.iter().enumerate() {
        entities.push(Entity::new(
            0,
            // convert from f64 to f32. Find a cleaner way, like a fn.
            lin_alg2::f32::Vec3::new(
                atom.position.x as f32,
                atom.position.y as f32,
                atom.position.z as f32,
            ),
            lin_alg2::f32::Quaternion::new(
                atom.orientation.w as f32,
                atom.orientation.x as f32,
                atom.orientation.y as f32,
                atom.orientation.z as f32,
            ),
            1.,
            atom.role.render_color(),
        ));
    }

    let scene = Scene {
        meshes: vec![Mesh::new_tetrahedron(render::SIDE_LEN)],
        entities,
        lighting: Lighting::default(),
    };

    let input_settings = InputSettings::default();

    graphics::run(scene, input_settings);
}

// todo: Use state cam? Should it be in `Scene`?
