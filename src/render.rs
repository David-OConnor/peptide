use bevy::prelude::*;
use bevy::render::primitives::Sphere;

use crate::{
    chem_definitions::{AtomType, BackboneRole},
    coord_gen::{ProteinCoords, ProteinDescription},
    // Don't import `Vec3` here since Bevy has its own
    lin_alg::{self, Quaternion},
    ROTATION_SPEED,
};
// Reference: https://bevyengine.org/examples/games/alien-cake-addict/

// use gdnative::prelude::*;

const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);

const TIME_STEP: f32 = 1.0 / 60.0;
const LIGHT_INTENSITY: f32 = 1500.0;

const DT: f64 = 1. / 60.;

#[derive(Clone, Component, Debug)]
struct AtomRender {
    // pub entity: Option<Entity>, // todo: Do we want this?
    pub id: usize,
    pub role: BackboneRole,
    pub position: lin_alg::Vec3,
    pub orientation: Quaternion,
}

impl Default for AtomRender {
    // Required for Bevy init
    fn default() -> Self {
        Self {
            // entity: None,
            id: 0,
            role: BackboneRole::N,
            position: Default::default(),
            orientation: Quaternion::new_identity(),
        }
    }
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

/// Store our atom coordinates and orientations here, for use
/// with the renderer.
#[derive(Component)] // Bevy resources must implement `Default` or `FromWorld`.
struct RenderState {
    /// Descriptions of each amino acid, including its name, and bond angles.
    pub protein_descrip: ProteinDescription,
}

impl Default for RenderState {
    // Required for Bevy init
    fn default() -> Self {
        Self {
            protein_descrip: ProteinDescription { aas: Vec::new() },
        }
    }
}

/// set up a simple 3D scene
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut render_state: ResMut<RenderState>,
) {
    // Initialize the render state here.
    render_state.protein_descrip = crate::init_protein();
    let coords = ProteinCoords::from_descrip(&render_state.protein_descrip);

    // Render our atoms.
    for (id, coords) in coords.atoms_backbone.iter().enumerate() {
        let atom = AtomRender {
            id,
            role: coords.role,
            position: coords.position,
            orientation: coords.orientation,
        };
        commands
            .spawn()
            .insert(atom.clone())
            .insert_bundle(PbrBundle {
                mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
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

    // camera
    commands.spawn_bundle(PerspectiveCameraBundle {
        transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
}

fn change_dihedral_angle(
    keyboard_input: Res<Input<KeyCode>>,
    // mut query: Query<&mut Transform, With<AtomRender>>,
    mut query: Query<(&mut AtomRender, &mut Transform), With<AtomRender>>,
    mut state: ResMut<RenderState>,
) {
    // Update ω
    if keyboard_input.pressed(KeyCode::Left) {
        state.protein_descrip.aas[0].ω += ROTATION_SPEED * DT;
        println!("Changing ω: {:?}", state.protein_descrip.aas[0].ω);

        // Recalculate coordinates now that we've updated our bond angles
        let coords = ProteinCoords::from_descrip(&state.protein_descrip).atoms_backbone;

        for (_id, (mut atom_render, mut transform)) in query.iter_mut().enumerate() {
            // Note: We are using an id field on `AtomRender`, but from initial tests, this
            // appears to be in sync with enumerating the query. (But isn't guaranteed
            // to be by Bevy?)

            let atom = &coords[atom_render.id];

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
        }
    }
}

pub fn run() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .init_resource::<RenderState>()
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup)
        .add_system(change_dihedral_angle)
        .insert_resource(ClearColor(BACKGROUND_COLOR))
        .add_system_set(
            SystemSet::new()
                .with_run_criteria(FixedTimestep::step(render::TIME_STEP as f64))
        )
        .run();
}
