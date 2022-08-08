//! Render, using Bevy. Bevy's API and docs are a mess; ideally transition to a different
//! engine, or a custom one using a Vulkan API like WGPU or Ash. It works for
//! now as a quick+dirty solution.

use bevy::{
    input::mouse::MouseMotion,
    prelude::shape::{Capsule, CapsuleUvProfile, UVSphere},
    prelude::*,
    time::FixedTimestep,
};

// use bevy::{core::FixedTimestep, render::primitives::Sphere};

use crate::{
    chem_definitions::BackboneRole,
    coord_gen::{ProteinCoords, LEN_CP_N},
    lin_alg::{self, Quaternion},
    render::{
        BACKGROUND_COLOR, BOND_COLOR, BOND_N_SIDES, BOND_RADIUS, CAM_MOVE_SENS,
        CAM_ROTATE_KEY_SENS, CAM_ROTATE_SENS, DT, FWD_VEC, LIGHT_INTENSITY, RIGHT_VEC, UP_VEC,
    },
    State, ROTATION_SPEED,
};

impl lin_alg::Vec3 {
    pub fn to_bevy(self) -> Vec3 {
        Vec3 {
            x: self.x as f32,
            y: self.y as f32,
            z: self.z as f32,
        }
    }
}

impl Quaternion {
    pub fn to_bevy(self) -> Quat {
        Quat::from_xyzw(self.x as f32, self.y as f32, self.z as f32, self.w as f32)
    }
}

#[derive(Component)]
struct Cam;

#[derive(Clone, Component)]
struct AtomRender {
    pub id: usize,
}

#[derive(Component)]
struct BondRender {
    pub bond_id: usize, // Note: Currently not a unique ID; only unique per atom.
    pub atom_id: usize,
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
    state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);

    // Render our atoms.
    for (id, atom) in state.protein_coords.atoms_backbone.iter().enumerate() {
        let atom_render = AtomRender { id };

        commands
            .spawn()
            .insert(atom_render.clone())
            .insert_bundle(PbrBundle {
                //UVSphere {
                //     radius: f32,
                //     sectors: usize,
                //     stacks: usize,
                // }
                mesh: meshes.add(Mesh::from(shape::Cube { size: 0.3 })),
                material: materials.add(atom.role.render_color().into()),
                ..Default::default()
            });

        let bond_renders = [
            BondRender {
                bond_id: 0,
                atom_id: id,
            },
            BondRender {
                bond_id: 1,
                atom_id: id,
            },
            BondRender {
                bond_id: 2,
                atom_id: id,
            },
            BondRender {
                bond_id: 3,
                atom_id: id,
            },
        ];

        for (i, bond_render) in bond_renders.into_iter().enumerate() {
            commands
                .spawn()
                .insert(bond_render)
                .insert_bundle(PbrBundle {
                    mesh: meshes.add(Mesh::from(Capsule {
                        radius: BOND_RADIUS,
                        rings: 1,
                        depth: LEN_CP_N as f32,
                        latitudes: 4,
                        longitudes: BOND_N_SIDES,
                        uv_profile: CapsuleUvProfile::Uniform, // todo
                    })),
                    material: materials
                        .add(Color::rgb(BOND_COLOR.0, BOND_COLOR.1, BOND_COLOR.2).into()),
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

    let mut cam_transform = Transform::identity();
    cam_transform.translation = state.cam.position.to_bevy();
    cam_transform.rotation = state.cam.orientation.to_bevy();

    commands
        .spawn_bundle(Camera3dBundle {
            transform: cam_transform.looking_at(Vec3::Z, Vec3::Y),
            // transform: cam_transform,
            // transform: Transform::from_xyz(
            //     state.cam.position.x as f32,
            //     state.cam.position.y as f32,
            //     state.cam.position.z as f32,
            // )
            //     // todo: Default orientation from state
            //     .looking_at(Vec3::Z, Vec3::Y),
            ..default()
        })
        .insert(Cam);
}

fn change_dihedral_angle(
    keyboard_input: Res<Input<KeyCode>>,
    // mut query: Query<(&mut AtomRender, &mut Transform), With<AtomRender>>,
    // mut query_bond: Query<(&mut BondRender, &mut Transform), With<BondRender>>,
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

    let mut changed = false;

    if keyboard_input.pressed(KeyCode::Key0) && max_i >= 1 {
        state.active_residue = 0;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key1) && max_i >= 2 {
        state.active_residue = 1;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key2) && max_i >= 3 {
        state.active_residue = 2;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key3) && max_i >= 4 {
        state.active_residue = 3;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key4) && max_i >= 5 {
        state.active_residue = 4;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key5) && max_i >= 6 {
        state.active_residue = 5;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key6) && max_i >= 7 {
        state.active_residue = 6;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key7) && max_i >= 8 {
        state.active_residue = 7;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key8) && max_i >= 9 {
        state.active_residue = 8;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::Key9) && max_i >= 10 {
        state.active_residue = 9;
        changed = true;
    }

    let ar = state.active_residue; // code shortener.

    if keyboard_input.pressed(KeyCode::T) {
        state.protein_descrip.residues[ar].φ += rotation_amt;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::G) {
        state.protein_descrip.residues[ar].φ -= rotation_amt;
        changed = true;
    }
    if keyboard_input.pressed(KeyCode::Y) {
        state.protein_descrip.residues[ar].ψ += rotation_amt;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::H) {
        state.protein_descrip.residues[ar].ψ -= rotation_amt;
        changed = true;
    }
    if keyboard_input.pressed(KeyCode::U) {
        state.protein_descrip.residues[ar].ω += rotation_amt;
        changed = true;
    } else if keyboard_input.pressed(KeyCode::J) {
        state.protein_descrip.residues[ar].ω -= rotation_amt;
        changed = true;
    }

    if changed {
        // Recalculate coordinates now that we've updated our bond angles
        state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
    }
}

/// Re-render all atoms, based on the latest calculated positions and orientations.
fn render_atoms(
    mut query: Query<(&mut AtomRender, &mut Transform), With<AtomRender>>,
    state: Res<State>,
) {
    for (atom_render, mut transform) in query.iter_mut() {
        // Note: We are using an id field on `AtomRender`, but from initial tests, this
        // appears to be in sync with enumerating the query. (But isn't guaranteed
        // to be by Bevy?)

        let atom = &state.protein_coords.atoms_backbone[atom_render.id];

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

/// Re-render all bonds, based on the latest calculated positions and orientations.
fn render_bonds(
    mut query: Query<(&mut BondRender, &mut Transform), With<BondRender>>,
    state: Res<State>,
) {
    // todo: DRY from coord_gen
    let bond_vecs = [
        lin_alg::Vec3::new(-1., 1., 1.).to_normalized(),
        lin_alg::Vec3::new(1., 1., -1.).to_normalized(),
        lin_alg::Vec3::new(1., -1., 1.).to_normalized(),
        lin_alg::Vec3::new(-1., -1., -1.).to_normalized(),
    ];

    for (bond_render, mut transform) in query.iter_mut() {
        let atom = &state.protein_coords.atoms_backbone[bond_render.atom_id];

        let bond_worldspace = atom.orientation.rotate_vec(bond_vecs[bond_render.bond_id]);

        // Angle doesn't matter, since it's radially symetric.
        // let bond_orientation = Quaternion::from_axis_angle(bond_worldspace, 0.);

        // todo: Figure this out.
        // let a = Quaternion::from_unit_vecs(lin_alg::Vec3::new(0., 0., 1.), bond_worldspace);
        // let b = Quaternion::new_identity().rotate_vec(bond_worldspace);

        // in: Vec, and Q_I. Out: Q

        // let bond_orientation = a * atom.orientation;
        // let bond_orientation = atom.orientation;
        // let bond_orientation = Quaternion::from_axis_angle(bond_worldspace, 0.);

        //
        // let objectRay = atom_
        // normalize3(objectRay);
        // angleDif = acos(dotProduct(targetRay,objectRay));
        // if (angleDif!=0) {
        //     orthoRay = crossProduct(objectRay,targetRay);
        //     normalize3(orthoRay);
        //     deltaQ = quaternionFromAxisAngle(orthoRay,angleDif);
        //     rotationQuaternion = deltaQ*rotationQuaternion;
        //     normalize4(rotationQuaternion);
        // }

        // let bond_orientation =
        //     Quaternion::from_euler(bond_worldspace.x, bond_worldspace.y, bond_worldspace.z);

        let bond_orientation =
            Quaternion::from_vec_direction(bond_worldspace, lin_alg::Vec3::new(0., 1., 0.));

        transform.translation = atom.position.to_bevy();
        transform.rotation = bond_orientation.to_bevy();
    }
}

/// Re-render all bonds, based on the latest calculated positions and orientations.
fn adjust_camera(
    mut query: Query<(&Cam, &mut Transform)>,
    keyboard_input: Res<Input<KeyCode>>,
    mut mouse_motion_events: EventReader<MouseMotion>,
    mut state: ResMut<State>,
) {
    const MOVE_AMT: f64 = CAM_MOVE_SENS * DT;
    const ROTATE_AMT: f64 = CAM_ROTATE_SENS * DT;
    const ROTATE_KEY_AMT: f64 = CAM_ROTATE_KEY_SENS * DT;

    let mut cam_moved = false;
    let mut cam_rotated = false;

    if keyboard_input.pressed(KeyCode::W) {
        state.cam.position.z += MOVE_AMT;
        cam_moved = true;
    } else if keyboard_input.pressed(KeyCode::S) {
        state.cam.position.z -= MOVE_AMT;
        cam_moved = true;
    }

    if keyboard_input.pressed(KeyCode::D) {
        state.cam.position.x += MOVE_AMT;
        cam_moved = true;
    } else if keyboard_input.pressed(KeyCode::A) {
        state.cam.position.x -= MOVE_AMT;
        cam_moved = true;
    }

    let mut rotation = lin_alg::Quaternion::new_identity();

    if keyboard_input.pressed(KeyCode::E) {
        rotation = Quaternion::from_axis_angle(FWD_VEC, ROTATE_KEY_AMT);
        cam_rotated = true;
    } else if keyboard_input.pressed(KeyCode::Q) {
        rotation = Quaternion::from_axis_angle(FWD_VEC, -ROTATE_KEY_AMT);
        cam_rotated = true;
    }

    // len 1, but Bevy requires iter.
    for mouse_motion in mouse_motion_events.iter() {
        let rotation =
            Quaternion::from_axis_angle(UP_VEC, mouse_motion.delta.x as f64 * ROTATE_AMT)
                * Quaternion::from_axis_angle(RIGHT_VEC, mouse_motion.delta.y as f64 * ROTATE_AMT)
                * rotation;

        state.cam.orientation = rotation * state.cam.orientation;
        cam_rotated = true;
    }

    if cam_moved {
        for (_camera, mut transform) in query.iter_mut() {
            transform.translation = state.cam.position.to_bevy();
        }
    }

    if cam_rotated {
        for (_camera, mut transform) in query.iter_mut() {
            transform.rotation = state.cam.orientation.to_bevy();
        }
    }
}

/// The entry point for our renderer.
pub fn run() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .init_resource::<State>()
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup)
        .add_system(adjust_camera)
        .add_system(change_dihedral_angle)
        .add_system(render_atoms)
        .add_system(render_bonds)
        .insert_resource(ClearColor(
            Color::rgb(BACKGROUND_COLOR.0, BACKGROUND_COLOR.1, BACKGROUND_COLOR.2).into(),
        ))
        .add_system_set(SystemSet::new().with_run_criteria(FixedTimestep::step(DT as f64)))
        .run();
}
