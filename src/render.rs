use bevy::prelude::*;
use bevy::render::primitives::Sphere;

use crate::lin_alg::{Vec3, Quaternion};

// Reference: https://bevyengine.org/examples/games/alien-cake-addict/

// use gdnative::prelude::*;

pub const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);

pub const TIME_STEP: f32 = 1.0 / 60.0;

/// set up a simple 3D scene
pub fn setup_render(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut render_state: Res<crate::RenderState>
) {
    // todo: Unsafe - temp for using static mut pts
    unsafe {
        // sphere

        for atom in render_state.atoms {
            commands.spawn_bundle(PbrBundle {
                mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
                material: materials.add(Color::rgb(1., 0., 0.).into()),
                transform: Transform::from_xyz(
                    atom.0.x as f32,
                    atom.0.y as f32,
                    atom.0.z as f32,
                ),
                ..default()
            });

            // todo: Also render bonds, perhaps as a cylinder or line for now.

        }
    }

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

// #[derive(NativeClass)]
// #[inherit(Node)]
// pub struct HelloWorld;
//
// #[methods]
// impl HelloWorld {
//     fn new(_owner: &Node) -> Self {
//         HelloWorld
//     }
//
//     #[export]
//     fn _ready(&self, _owner: &Node) {
//         godot_print!("Hello, world.");
//     }
// }
//
// fn init(handle: InitHandle) {
//     handle.add_class::<HelloWorld>();
// }
//
// godot_init!(init);
