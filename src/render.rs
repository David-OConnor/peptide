use bevy::prelude::*;
use bevy::render::primitives::Sphere;

/// set up a simple 3D scene
pub fn setup_render(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // todo: Unsafe - temp for using static mut pts
    unsafe {
        // sphere
        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(1., 0., 0.).into()),
            transform: Transform::from_xyz(
                crate::PT1.x as f32,
                crate::PT1.y as f32,
                crate::PT1.z as f32,
            ),
            ..default()
        });

        // sphere
        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(0., 1., 0.).into()),
            transform: Transform::from_xyz(
                crate::PT2.x as f32,
                crate::PT2.y as f32,
                crate::PT2.z as f32,
            ),
            ..default()
        });

        // sphere
        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(0., 0., 1.).into()),
            transform: Transform::from_xyz(
                crate::PT3.x as f32,
                crate::PT3.y as f32,
                crate::PT3.z as f32,
            ),
            ..default()
        });

        // sphere
        commands.spawn_bundle(PbrBundle {
            mesh: meshes.add(Mesh::from(shape::Cube { size: 0.2 })),
            material: materials.add(Color::rgb(1., 1., 0.).into()),
            transform: Transform::from_xyz(
                crate::PT4.x as f32,
                crate::PT4.y as f32,
                crate::PT4.z as f32,
            ),
            ..default()
        });
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
