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
    lin_alg::{self, Quaternion},
    State, ROTATION_SPEED,
};
// Reference: https://bevyengine.org/examples/games/alien-cake-addict/

const BACKGROUND_COLOR: Color = Color::rgb(0.9, 0.9, 0.9);

const LIGHT_INTENSITY: f32 = 1500.0;

const BOND_RADIUS: f32 = 0.03;
const BOND_SIDES: usize = 10;

const BOND_COLOR: Color = Color::rgb(0.2, 0.2, 0.2);

const CAM_MOVE_AMT: f64 = 1.;
const DT: f64 = 1. / 60.;

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
    let coords = ProteinCoords::from_descrip(&state.protein_descrip);

    // Render our atoms.
    for (id, atom) in coords.atoms_backbone.iter().enumerate() {
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
                        longitudes: BOND_SIDES,
                        uv_profile: CapsuleUvProfile::Uniform, // todo
                    })),
                    material: materials.add(BOND_COLOR.into()),
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
        transform: Transform::from_xyz(0.0, 0.0, -7.0).looking_at(Vec3::Z, Vec3::Y),
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
        state.protein_descrip.residues[ar].φ += rotation_amt;
    } else if keyboard_input.pressed(KeyCode::A) {
        state.protein_descrip.residues[ar].φ -= rotation_amt;
    }
    if keyboard_input.pressed(KeyCode::W) {
        state.protein_descrip.residues[ar].ψ += rotation_amt;
    } else if keyboard_input.pressed(KeyCode::S) {
        state.protein_descrip.residues[ar].ψ -= rotation_amt;
    }
    if keyboard_input.pressed(KeyCode::E) {
        state.protein_descrip.residues[ar].ω += rotation_amt;
    } else if keyboard_input.pressed(KeyCode::D) {
        state.protein_descrip.residues[ar].ω -= rotation_amt;
    }

    // todo: return if no key is pressed, so as not to update coordinates
    // and the transform.
    // else {
    //     return // Exit if no key changed
    // }

    // Recalculate coordinates now that we've updated our bond angles
    state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);
}

/// Re-render all atoms, based on the latest calculated positions and orientations.
fn render_atoms(
    mut query: Query<(&mut AtomRender, &mut Transform), With<AtomRender>>,
    state: Res<State>,
) {
    if state.protein_coords.atoms_backbone.len() == 0 {
        return; // This occurs once on init.
    }

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
    if state.protein_coords.atoms_backbone.len() == 0 {
        return; // This occurs once on init.
    }

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

        let bond_orientation = Quaternion::from_vec_direction(bond_worldspace, lin_alg::Vec3::new(0., 1., 0.));

        transform.translation = atom.position.to_bevy();
        transform.rotation = bond_orientation.to_bevy();
    }
}

/// Re-render all bonds, based on the latest calculated positions and orientations.
fn adjust_camera(
    mut query: Query<(&Cam, &mut Transform)>,
    keyboard_input: Res<Input<KeyCode>>,
) {

    let mut translation_y = 0.;
    if keyboard_input.pressed(KeyCode::Up) {
        translation_y += CAM_MOVE_AMT * DT;
    } else if keyboard_input.pressed(Dn) {
        translation_y -= CAM_MOVE_AMT * DT;
    }

    for (_camera, mut transform) in query.iter_mut() {
        transform.translation = // ...;
        // transform.rotation = // ...;
    }
}

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
        .insert_resource(ClearColor(BACKGROUND_COLOR))
        .add_system_set(SystemSet::new().with_run_criteria(FixedTimestep::step(DT as f64)))
        .run();
}
