//! Render, using Bevy. Bevy's API and docs are a mess; ideally transition to a different
//! engine, or a custom one using a Vulkan API like WGPU or Ash. It works for
//! now as a quick+dirty solution.

use bevy::{
    prelude::shape::{Capsule, CapsuleUvProfile, UVSphere},
    prelude::*,
    time::FixedTimestep,
};

// use bevy::{core::FixedTimestep, render::primitives::Sphere};

use crate::{
    chem_definitions::BackboneRole,
    coord_gen::{ProteinCoords, LEN_CP_N},
    lin_alg, State, ROTATION_SPEED,
};
// Reference: https://bevyengine.org/examples/games/alien-cake-addict/

const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);

const LIGHT_INTENSITY: f32 = 1500.0;

const BOND_RADIUS: f32 = 0.03;
const BOND_SIDES: usize = 10;

const BOND_COLOR: Color = Color::rgb(0.2, 0.2, 0.2);

const DT: f64 = 1. / 60.;

#[derive(Clone, Component, Debug)]
struct AtomRender {
    pub id: usize,
}

#[derive(Component, Debug)]
struct BondRender {
    pub id: usize,
}
impl BackboneRole {
    pub fn render_color(&self) -> Color {
        match self {
            Self::Cα => Color::rgb(1., 0., 1.),
            Self::Cp => Color::rgb(0., 1., 1.),
            Self::N => Color::rgb(0., 0., 1.),
            Self::O => Color::rgb(1., 0., 0.),
        }
    }
}

/// set up a simple 3D scene
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut state: ResMut<State>,
) {
    // Initialize the render state here.
    state.protein_descrip = crate::init_protein();
    let coords = ProteinCoords::from_descrip(&state.protein_descrip);

    // Render our atoms.
    for (id, atom) in coords.atoms_backbone.iter().enumerate() {
        let atom_render = AtomRender { id };

        let bond_renders = [
            BondRender { id: id * 4 },
            BondRender { id: id * 4 + 1 },
            BondRender { id: id * 4 + 2 },
            BondRender { id: id * 4 + 3 },
        ];

        // todo: DRY from coord_gen
        let bond_vecs = [
            lin_alg::Vec3::new(-1., 1., 1.).to_normalized(),
            lin_alg::Vec3::new(1., 1., -1.).to_normalized(),
            lin_alg::Vec3::new(1., -1., 1.).to_normalized(),
            lin_alg::Vec3::new(-1., -1., -1.).to_normalized(),
        ];

        commands
            .spawn()
            .insert(atom_render.clone())
            .insert_bundle(PbrBundle {
                //UVSphere {
                //     radius: f32,
                //     sectors: usize,
                //     stacks: usize,
                // }
                mesh: meshes.add(Mesh::from(shape::Cube { size: 0.25 })),
                material: materials.add(atom.role.render_color().into()),
                transform: Transform::from_xyz(
                    atom.position.x as f32,
                    atom.position.y as f32,
                    atom.position.z as f32,
                )
                .with_rotation(Quat::from_xyzw(
                    atom.orientation.x as f32,
                    atom.orientation.y as f32,
                    atom.orientation.z as f32,
                    atom.orientation.w as f32,
                )),
                ..Default::default()
            });

        for (i, bond_render) in bond_renders.into_iter().enumerate() {
            let bond_orientation = atom.orientation.rotate_vec(bond_vecs[i]);

            commands
                .spawn()
                .insert(bond_render)
                .insert_bundle(PbrBundle {
                    mesh: meshes.add(Mesh::from(Capsule {
                        radius: BOND_RADIUS,
                        rings: 1,
                        depth: LEN_CP_N as f32,
                        latitudes: 6, // todo what is this?
                        longitudes: BOND_SIDES,
                        uv_profile: CapsuleUvProfile::Uniform, // todo
                    })),
                    material: materials.add(BOND_COLOR.into()),
                    transform: Transform::from_xyz(
                        atom.position.x as f32,
                        atom.position.y as f32,
                        atom.position.z as f32,
                    )
                    .with_rotation(Quat::from_xyzw(
                        bond_orientation.x as f32,
                        atom.orientation.y as f32,
                        atom.orientation.z as f32,
                        atom.orientation.w as f32,
                    )),
                    ..Default::default()
                });
        }
    }

    // todo: Also render bonds, perhaps as a cylinder or line for now.

    // light
    commands.spawn_bundle(PointLightBundle {
        point_light: PointLight {
            intensity: LIGHT_INTENSITY,
            shadows_enabled: true,
            ..default()
        },

        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });

    commands.spawn_bundle(Camera3dBundle {
        transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
}

fn change_dihedral_angle(
    keyboard_input: Res<Input<KeyCode>>,
    // mut query: Query<&mut Transform, With<AtomRender>>,
    mut query: Query<(&mut AtomRender, &mut Transform), With<AtomRender>>,
    mut state: ResMut<State>,
) {
    // Update dihedral angles.
    // todo: Move keyboard logic elsewhere.
    let rotation_amt = ROTATION_SPEED * DT;

    // Change active residue.
    // todo: You need to debounce this for it to work
    if keyboard_input.pressed(KeyCode::Up) {
        if state.active_residue != state.protein_descrip.residues.len() - 1 {
            state.active_residue += 1;
        }
    } else if keyboard_input.pressed(KeyCode::Down) && state.active_residue != 0 {
        state.active_residue -= 1;
    }

    // todo: Lots of DRY here

    let max_i = state.protein_descrip.residues.len(); // Code shortener

    if keyboard_input.pressed(KeyCode::Key0) && max_i >= 1 {
        state.active_residue = 0;
    } else if keyboard_input.pressed(KeyCode::Key1) && max_i >= 2 {
        state.active_residue = 1;
    } else if keyboard_input.pressed(KeyCode::Key2) && max_i >= 3 {
        state.active_residue = 2;
    } else if keyboard_input.pressed(KeyCode::Key3) && max_i >= 4 {
        state.active_residue = 3;
    } else if keyboard_input.pressed(KeyCode::Key4) && max_i >= 5 {
        state.active_residue = 4;
    } else if keyboard_input.pressed(KeyCode::Key5) && max_i >= 6 {
        state.active_residue = 5;
    } else if keyboard_input.pressed(KeyCode::Key6) && max_i >= 7 {
        state.active_residue = 6;
    } else if keyboard_input.pressed(KeyCode::Key7) && max_i >= 8 {
        state.active_residue = 7;
    } else if keyboard_input.pressed(KeyCode::Key8) && max_i >= 9 {
        state.active_residue = 8;
    } else if keyboard_input.pressed(KeyCode::Key9) && max_i >= 10 {
        state.active_residue = 9;
    }

    let ar = state.active_residue; // code shortener.

    if keyboard_input.pressed(KeyCode::Q) {
        state.protein_descrip.residues[ar].ω += rotation_amt;
    } else if keyboard_input.pressed(KeyCode::A) {
        state.protein_descrip.residues[ar].ω -= rotation_amt;
    }
    if keyboard_input.pressed(KeyCode::W) {
        state.protein_descrip.residues[ar].φ += rotation_amt;
    } else if keyboard_input.pressed(KeyCode::S) {
        state.protein_descrip.residues[ar].φ -= rotation_amt;
    }
    if keyboard_input.pressed(KeyCode::E) {
        state.protein_descrip.residues[ar].ψ += rotation_amt;
    } else if keyboard_input.pressed(KeyCode::D) {
        state.protein_descrip.residues[ar].ψ -= rotation_amt;
        // println!("Changing ψ: {:?}", state.protein_descrip.residues[0].ψ);
    }

    // todo: return if no key is pressed, so as not to update coordinates
    // and the transform.
    // else {
    //     return // Exit if no key changed
    // }

    for (atom_render, mut transform) in query.iter_mut() {
        // Note: We are using an id field on `AtomRender`, but from initial tests, this
        // appears to be in sync with enumerating the query. (But isn't guaranteed
        // to be by Bevy?)

        // Recalculate coordinates now that we've updated our bond angles
        let coords = ProteinCoords::from_descrip(&state.protein_descrip).atoms_backbone;

        let atom = &coords[atom_render.id];

        // println!("Atom: {:?}\n\n", atom);

        // Convert from our vector and quaternion types to Bevy's.
        let position = Vec3::new(
            atom.position.x as f32,
            atom.position.y as f32,
            atom.position.z as f32,
        );

        let orientation = Quat::from_xyzw(
            atom.orientation.x as f32,
            atom.orientation.y as f32,
            atom.orientation.z as f32,
            atom.orientation.w as f32,
        );

        transform.translation = position;
        transform.rotation = orientation;
        // }
    }
}

pub fn run() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .init_resource::<State>()
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup)
        .add_system(change_dihedral_angle)
        .insert_resource(ClearColor(BACKGROUND_COLOR))
        .add_system_set(SystemSet::new().with_run_criteria(FixedTimestep::step(DT as f64)))
        .run();
}
