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

mod lin_alg;
mod render;

use lin_alg::{Quaternion, Vec3, Pt3};

// Double bond len of C' to N.
const LEN_CP_N: f64 = 1.32; // angstrom
const LEN_N_CALPHA: f64 = 1.; // angstrom // todo
const LEN_CALPHA_CP: f64 = 1.; // angstrom // todo

// todo: These may change under diff circumstances, and be diff for the diff atoms
// and diff things on the atom.
const BACKBONE_BOND_ANGLES: f64 = 1.911; // radians

#[derive(Clone, Copy)]
enum CarbonBond {
    // todo: Come back to what these mean for the 3 backbone atoms
    A,
    B,
    C,
    D,
}

/// Location of an atom in 3d space, using the AA it's part of's
/// coordinate system. Pt 0 is defined as the α atom.
#[derive(Debug)]
struct AtomLocation {
    atom: Atom,
    point: Pt3,
}

impl AtomLocation {
    pub fn new(atom: Atom, point: Pt3) -> Self {
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
    pub fn _charge(&self) -> f64 {
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
            Self::D => "Asp",
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

    /// See note on `BackBoneCoords` for direction and related conventions.
    pub fn structure(&self) -> Vec<AtomLocation> {
        match self {
            Self::A => {
                vec![AtomLocation::new(Atom::C, Pt3::new(0., 0., 0.))]
            }

            _ => Vec::new(), // todo
        }
    }
}

/// A struct holding backbone atom coordinates, relative to the alpha carbon.
/// We will use the convention of the chain starting at NH2 (amine) end (N terminus), and ending in
/// COOH (carboxyl) (C terminus).
/// Flow is N -> C_alpha -> C'.
/// Coordinates are in relation to the Nitrogen molecule at the start of the segement, using the
/// directional convention described above. (todo: How do we orient axes?)
/// todo: To start, convention is the previous C' to this AA's starting A is on the positive X axis.
#[derive(Debug)]
struct BackBoneCoords {
    // pub c_alpha: f64, // This field may be unecessary if we define it to be 0.
    // /// Nitrogen atom bonded to C_alpha // defined as 0.
    // pub n: f64,

    // note: We currently don't store rotations (eg quaternions) here for each atom,
    // since the relation of the point to each other contains this information.
    /// Cα
    pub c_alpha: Pt3,
    /// Carbon' atom bound to C_alpha
    pub c_p: Pt3,
    /// Nitrogen atom of the next module.
    pub n_next: Pt3,
}

/// An amino acid in a protein structure, including position information.
#[derive(Debug)]
struct AaInProtein {
    aa: AminoAcid,
    // sequence_num: usize, // todo: Instead of this, use an array?
    ω: f64, // ω Assumed to be TAU/2 for most cases
    /// φ bond, between Cα and N
    ϕ: f64,
    /// ψ bond, between Cα and C'
    ψ: f64,
    // todo: Include bond lengths here if they're not constant.
}

// todo?
// fn bond_angle(orientation: Quaternion, bond: CarbonBond) -> Vec3 {
//
// }

/// Calculate the orientation, as a quaternion, of a backbone atom, given the orientation of the
/// previous atom, and the bond angle. `torsion_angle` is the vector representing the bond to
/// this atom from the previous atom's orientation.
fn find_backbone_atom_orientation(q_prev: Quaternion, torsion_angle_prev: Vec3, bond_angle: f64) {
    // We split up determining the orientation of each backbone atom into 3 steps:

    // #1: Find the bond axis. This is the world-space vector of the bond by rotating the
    // torsion angle by the previous atom's orientation quaternion.
    let bond_axis = q_prev.rotate_vec(torsion_angle_prev);

    // #2: Calculate the orientation quaternion of the current atom so its bond to the previous atom
    // is aligned with the bond axis found above.
    let current_atom_orientation = Quaternion::from_axis_angle(bond_axis);

    // #3: Rotate the current atom's orientation by the bond's rotation. (eg ϕ, or ψ).
    // let todo bond_angle
    current_atom_orientation * bond_angle
}

impl AaInProtein {
    /// Generate cartesian coordinates of points from diahedral angles and bond lengths. Starts with
    /// N, and ends with C'.
    /// Accepts position, and orientation of the N atom that starts this segment.
    /// Also returns orientation of the N atom, for use when calculating coordinates for the next
    /// AA in the chain.
    // pub fn backbone_cart_coords(&self, cp_n_bond_angle_prev: Vec3, ω_prev: f64) -> BackBoneCoords {
    pub fn backbone_cart_coords(
        &self,
        pos_n: Pt3,
        q_n: Quaternion,
    ) -> (BackBoneCoords, Quaternion) {
        let torsion_angle_next_atom = Vec3::new(0.577, 0.577, 0.577); // todo temp
        let torsion_angle_a = Vec3::new(-0.577, 0.577, 0.577); // todo temp
        let torsion_angle_b = Vec3::new(0.577, -0.577, 0.577); // todo temp
        let torsion_angle_c = Vec3::new(0., 0., -1.); // todo temp

        let q_calpha = find_backbone_atom_orientation(q_prev: q_n, torsion_angle_next_atom, self.ψ);
        let q_cp = find_backbone_atom_orientation(q_prev: q_calpha, torsion_angle_next_atom, self.ϕ);
        let q_n_next = find_backbone_atom_orientation(q_prev: q_cp, torsion_angle_next_atom, self.ω);

        // Orientations for the various atoms, Z up. todo ?

        // Calculate the orientation of each atom from the orientation, and torsion angle of the
        // bond from the previous atom, and the dihedral angle between the two.

        // Set up atom-oriented vectors encoding the angle and len of each bond.
        let bond_n_calpha = torsion_angle_next_atom * LEN_N_CALPHA;
        let bond_calpha_cp = torsion_angle_next_atom * LEN_CALPHA_CP;
        let bond_cp_n = torsion_angle_next_atom * LEN_CP_N;

        // Calculate the position of each atom from the position, orientation and torsion angle of the
        // bond from the previous atom.
        // todo: Start by describing each term in terms of its dependencies.
        let c_alpha = pos_n + q_n.rotate_vec(bond_n_c_alpha);
        let c_p = c_alpha + q_calpha.rotate_vec(bond_c_alpha_cp);
        let n_next = c_p + q_cp.rotate_vec(bond_cp_n);

        (
            BackBoneCoords {
                c_alpha,
                c_p,
                n_next,
            },
            q_n_next,
        )
    }
}

#[derive(Debug)]
struct Protein {
    pub aas: Vec<AaInProtein>,
}

fn main() {
    let a = AaInProtein {
        aa: AminoAcid::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        ϕ: 0. * TAU,
        ψ: 0. * TAU,
    };

    let b = AaInProtein {
        aa: AminoAcid::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        ϕ: 0. * TAU,
        ψ: 0. * TAU,
    };

    let c = AaInProtein {
        aa: AminoAcid::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        ϕ: 0. * TAU,
        ψ: 0. * TAU,
    };

    let prot = Protein { aas: vec![a, b, c] };

    let coords = a.backbone_cart_coords(
        Pt3 {
            x: 0.,
            y: 0.,
            z: 0.,
        },
        Quaternion::new_identity(),
    );

    println!("Coords: {:?}", coords);

    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_startup_system(render::setup_render)
        .run();
}
