//! Model, render, and predict protein structure. Uses the peptide bond
//! as an immutable basis for structure. Attempts to find foldin patterns
//! that may lead to an ultimate structure.

// todo: Switch from Bevy to Ash or similar.
// todo: Temperature sensitivity. Surroundin water molecules.
// Initially, focus on modeling the bond angles and backbone. Both
// data to describe, and a 3d render

// Doe sfolding begin starting at the end extruded?

use std::f64::consts::TAU;

// todo: Don' glob import

use bevy::prelude::*;

/// Location of an atom in 3d space, using the AA it's part of's
/// coordinate system. Pt 0 is defined as the __ atom.
#[derive(Debug)]
struct AtomLocation {
    atom: Atom,
    point: PtCart,
}

impl AtomLocation {
    pub fn new(atom: Atom, point: PtCart) -> Self {
        Self { atom, point }
    }
}

#[derive(Clone, Copy, Debug)]
enum Atom {
    C,
    N,
    H,
    O,
    P,
    S,
}

impl Atom {
    // todo: What is the significance of this? It's a bit nebulous
    pub fn charge(&self) -> f64 {
        match self {
            Self::C => 4.,
            Self::N => -3.,
            Self::H => 1.,
            Self::O => -2.,
            Self::P => 0., // (5., 3., -3.)
            Self::S => 0., // (-2., 2., 4., 6.)
        }
    }
}

/// A point in cartesian coordinates
#[derive(Debug)]
struct PtCart {
    x: f64,
    y: f64,
    z: f64,
}

impl PtCart {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

#[derive(Clone, Copy, Debug)]
enum AminoAcid {
    A,
    R,
    N,
    D,
    C,
    Q,
    E,
    G,
    H,
    I,
    L,
    K,
    M,
    F,
    P,
    S,
    T,
    W,
    Y,
    V,
}

impl AminoAcid {
    pub fn symbol(&self) -> &str {
        match self {
            Self::A => "Ala",
            Self::R => "Arg",
            Self::N => "Asn",
            Self::C => "Cys",
            Self::Q => "Gln",
            Self::E => "Glu",
            Self::G => "Gly",
            Self::H => "His",
            Self::I => "Ile",
            Self::L => "Leu",
            Self::K => "Lys",
            Self::M => "Met",
            Self::F => "Phe",
            Self::P => "Pro",
            Self::S => "Ser",
            Self::T => "Thr",
            Self::W => "Trp",
            Self::Y => "Tyr",
            Self::V => "Val",
        }
    }

    pub fn structure(&self) -> Vec<AtomLocation> {
        match self {
            Self::A => {
                vec![AtomLocation::new(Atom::C, PtCart::new(0., 0., 0.))]
            }

            _ => Vec::new(), // todo
        }
    }
}

/// An amino acid in a protein structure, including position information.

#[derive(Debug)]
struct AaInProtein {
    aa: AminoAcid,
    sequence_num: usize, // todo: Instead of this, use an array?
    // omega: f64,  // ω Assumed to be TAU/2
    /// φ bond, between C^α and N
    phi: f64,
    /// ψ bond, between C^α and C'
    psi: f64,
}

/// set up a simple 3D scene

fn setup_render(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // plane
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Plane { size: 5.0 })),
        material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
        ..default()
    });

    // cube
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(Mesh::from(shape::Cube { size: 1.0 })),
        material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
        transform: Transform::from_xyz(0.0, 0.5, 0.0),
        ..default()
    });

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

#[derive(Debug)]
struct Protein {
    aas: Vec<AaInProtein>,
}

fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup_render)
        .run();
}
