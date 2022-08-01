use bevy::prelude::*;
use bevy::render::primitives::Sphere;

use crate::{
    // Don't import `Vec3` here since Bevy has its own
    lin_alg::{self, Quaternion},
    ROTATION_SPEED,
};

// Reference: https://bevyengine.org/examples/games/alien-cake-addict/

// use gdnative::prelude::*;

const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);

const TIME_STEP: f32 = 1.0 / 60.0;

#[derive(Default, Component, Debug)]
struct AaRender {
    pub atom_cα: AtomRender,
    pub atom_cp: AtomRender,
    pub atom_n_next: AtomRender,
}

#[derive(Default, Component, Debug)]
struct AtomRender {
    // pub type_: crate::AtomType,
    pub entity: Option<Entity>, // todo: Do we want this?
    pub position: lin_alg::Vec3,
    pub orientation: Quaternion,
}

/// Store our atom coordinates and orientations here, for use
/// with the renderer.
#[derive(Default, Component)] // Bevy resources must implement `Default` or `FromWorld`.
struct RenderState {
    /// Descriptions of each amino acid, including its name, and bond angles.
    pub amino_acids: Vec<crate::AaDihedralAngles>,
    /// Geometry of each atom, including position, and orientation. Updated whenever
    /// any bond angle in the protein changes.
    pub aa_renders: Vec<AaRender>,
}

/// set up a simple 3D scene
fn setup_render(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut render_state: ResMut<RenderState>,
) {
    // Initialize the render state here.
    render_state.amino_acids = crate::init_protein();

    let mut aas = Vec::new();

    for aa in &render_state.amino_acids {
        let coords =
            aa.backbone_cart_coords(lin_alg::Vec3::new(0., 0., 0.), Quaternion::new_identity());

        aas.push(AaRender {
            atom_cα: AtomRender {
                entity: Some(commands.spawn().id()),
                position: coords.cα,
                orientation: coords.cα_orientation,
            },
            atom_cp: AtomRender {
                entity: Some(commands.spawn().id()),
                position: coords.cp,
                orientation: coords.cp_orientation,
            },
            atom_n_next: AtomRender {
                entity: Some(commands.spawn().id()),
                position: coords.n_next,
                orientation: coords.n_next_orientation,
            },
        });
    }

    render_state.aa_renders = aas;

    // Render our atoms.
    for aa in &render_state.aa_renders {
        let atom = &aa.atom_cα; // todo temp

        println!("atom {:?}", atom);

        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(1., 0., 0.).into()),
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

        let atom = &aa.atom_cp; // todo temp

        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(0., 1., 0.).into()),
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

        let atom = &aa.atom_n_next; // todo temp

        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(0., 0., 1.).into()),
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
    }

    // todo: Also render bonds, perhaps as a cylinder or line for now.

    // light
    commands.spawn_bundle(PointLightBundle {
        point_light: PointLight {
            intensity: 1500.0,
            shadows_enabled: true,
            ..default()
        },

        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });

    // camera
    commands.spawn_bundle(PerspectiveCameraBundle {
        transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
}

fn change_dihedral_angle(
    keyboard_input: Res<Input<KeyCode>>,
    mut query: Query<&mut Transform, With<RenderState>>,
    mut render_state: ResMut<RenderState>,
) {
    // Update ω
    if keyboard_input.pressed(KeyCode::Left) {
        println!("Changing ω");
        render_state.amino_acids[0].ω += ROTATION_SPEED;
    }

    // todo: Dry here with init
    // Recalculate coordinates now that we've updated our bond angles
    let coords = render_state.amino_acids[0]
        .backbone_cart_coords(lin_alg::Vec3::new(0., 0., 0.), Quaternion::new_identity());

    // We've now recalculated the positions and orientations. Will this trigger a re-render
    // by updating the state as a resource here?
    // render_state.aa_renders.push(atom_cα.position = coords.cα);
    // render_state.aa_renders[0].position = coords.cp;
    // render_state.aa_renders[0].position = coords.n_next;

    // render_state.aa_renders.push(atom_cα.orientation = coords.cα_orientation);
    // render_state.aa_renders[0].orientation = coords.cp_orientation;
    // render_state.aa_renders[0].orientation = coords.n_next_orientation;

    for aa in &render_state.aa_renders {
        // for atom in &[aa.atom_cα, aa.atom_cp, aa.atom_n_next] {
        let mut atom = &aa.atom_cα; // todo temp.
        if let Ok(mut transform) = query.get(atom.entity.unwrap()) {
            // transform.translation = Vec3::new(
            //     atom.position.x as f32,
            //     atom.position.y as f32,
            //     atom.position.z as f32,
            // );
            // transform.rotation = Quat::from_xyzw(
            //     atom.orientation.x as f32,
            //     atom.orientation.y as f32,
            //     atom.orientation.z as f32,
            //     atom.orientation.w as f32,
            // );
        }
        // }
    }
}

pub fn run() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .init_resource::<RenderState>()
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup_render)
        .add_system(change_dihedral_angle)
        .insert_resource(ClearColor(BACKGROUND_COLOR))
        // .add_system_set(
        //     SystemSet::new()
        //         .with_run_criteria(FixedTimestep::step(render::TIME_STEP as f64))
        // )
        .run();
}
