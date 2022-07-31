//! Model, render, and predict protein structure. Uses the peptide bond
//! as an immutable basis for structure. Attempts to find foldin patterns
//! that may lead to an ultimate structure.

// todo: Switch from Bevy to Ash or similar. Or perhaps Godot's rust bindings.
// todo: Temperature sensitivity. Surroundin water molecules.
// Initially, focus on modeling the bond angles and backbone. Both
// data to describe, and a 3d render

// Doe sfolding begin starting at the end extruded?

use std::f64::consts::TAU;

// use graphics;
// use gdnative::prelude::*;

mod lin_alg;
mod render;

use lin_alg::{Quaternion, Vec3};

// todo temp while we figure out how to pass val to bevy

// Double bond len of C' to N.
// todo: 1.089000
const LEN_CP_N: f64 = 1.32; // angstrom
const LEN_N_CALPHA: f64 = 1.; // angstrom // todo
const LEN_CALPHA_CP: f64 = 1.; // angstrom // todo

// todo: These may change under diff circumstances, and be diff for the diff atoms
// and diff things on the atom.
const BACKBONE_BOND_ANGLES: f64 = 1.911; // radians

const ROTATION_SPEED: f64 = 1.; // radians/s

#[derive(Debug)]
struct _Protein {
    pub aas: Vec<AaInProtein>,
}

#[derive(Clone, Copy)]
enum CarbonBond {
    // todo: Come back to what these mean for the 3 backbone atoms
    A,
    B,
    C,
    D,
}

// /// Location of an atom in 3d space, using the AA it's part of's
// /// coordinate system. Pt 0 is defined as the α atom.
// #[derive(Debug)]
// struct AtomLocation {
//     atom: Atom,
//     point: Vec3,
// }
//
// impl AtomLocation {
//     pub fn new(atom: Atom, point: Vec3) -> Self {
//         Self { atom, point }
//     }
// }

#[derive(Clone, Copy, Debug)]
enum AtomType {
    C,
    N,
    H,
    O,
    P,
    S,
}

impl AtomType {
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

    // /// See note on `BackBoneCoords` for direction and related conventions.
    // pub fn structure(&self) -> Vec<AtomLocation> {
    //     match self {
    //         Self::A => {
    //             vec![AtomLocation::new(Atom::C, Vec3::new(0., 0., 0.))]
    //         }
    //
    //         _ => Vec::new(), // todo
    //     }
    // }
}

/// A struct holding backbone atom coordinates, relative to the alpha carbon.
/// Also hols orientations.
/// We will use the convention of the chain starting at NH2 (amine) end (N terminus), and ending in
/// COOH (carboxyl) (C terminus).
/// Flow is N -> C_alpha -> C'.
/// Coordinates are in relation to the Nitrogen molecule at the start of the segement, using the
/// directional convention described above. (todo: How do we orient axes?)
/// todo: To start, convention is the previous C' to this AA's starting A is on the positive X axis.
#[derive(Debug)]
struct BackBoneCoords {
    // /// Nitrogen atom bonded to C_alpha // defined as 0.
    // pub n: f64,
    /// Cα
    pub cα: Vec3,
    /// Carbon' atom bound to C_alpha
    pub cp: Vec3,
    /// Nitrogen atom of the next module.
    pub n_next: Vec3,
    pub cα_orientation: Quaternion,
    pub cp_orientation: Quaternion,
    pub n_next_orientation: Quaternion,
}

// todo?
// fn bond_angle(orientation: Quaternion, bond: CarbonBond) -> Vec3 {
//
// }

/// Calculate the orientation, as a quaternion, of a backbone atom, given the orientation of the
/// previous atom, and the bond angle. `bond_angle` is the vector representing the bond to
/// this atom from the previous atom's orientation.
fn find_backbone_atom_orientation(
    q_prev: Quaternion,
    bond_prev: Vec3,
    dihedral_angle: f64,
) -> Quaternion {
    // We split up determining the orientation of each backbone atom into 3 steps:

    // #1: Find the bond axis. This is the world-space vector of the backbone bond by rotating the
    // bond angle with the previous atom's orientation quaternion.
    let bond_axis = q_prev.rotate_vec(bond_prev);

    // todo: Rotate from prev bond to next bond.

    // #2: Calculate the orientation quaternion of the current atom so its bond to the previous atom
    // is aligned with the bond axis found above, and rotated per the dihedral angle.
    let current_atom_orientation = Quaternion::from_axis_angle(bond_axis, dihedral_angle);

    // #3: Rotate the current atom's orientation by the dihedral angle between the two atoms. (eg φ, or ψ).

    // todo: Instead of actual dihedral angles, we are using a somewhat arbitrary rotation. todo: Change
    // todo this later to be actual dihedral angle. To do so, you probably need to pass the plane formed.
    // todo: Consider changing this variable name to `bond_rotation` in the meanwhile.

    // by the prev 2 angles, and subtract it from the plane formed by this and the next.

    // current_atom_orientation * dihedral_angle

    current_atom_orientation
}

/// Used to represent one atom in a system built of atoms.
struct _ZmatrixItem {
    atomic_number: u8,
    bond_len: f64,       // Angstrom
    bond_angle: f64,     // radians
    dihedral_angle: f64, // radians
}

/// An amino acid in a protein structure, including position information.
#[derive(Debug)]
struct AaInProtein {
    pub aa: AminoAcid,
    /// Dihedral angle between C' and N
    /// Tor (Cα, C, N, Cα) is the ω torsion angle
    /// Assumed to be TAU/2 for most cases
    pub ω: f64,
    /// Dihedral angle between Cα and N.
    /// Tor (C, N, Cα, C) is the φ torsion angle
    pub φ: f64,
    /// Dihedral angle, between Cα and C'
    ///  Tor (N, Cα, C, N) is the ψ torsion angle
    pub ψ: f64,
    // todo: Include bond lengths here if they're not constant.
}

impl AaInProtein {
    /// Generate cartesian coordinates of points from diahedral angles and bond lengths. Starts with
    /// N, and ends with C'.
    /// Accepts position, and orientation of the N atom that starts this segment.
    /// Also returns orientation of the N atom, for use when calculating coordinates for the next
    /// AA in the chain.
    pub fn backbone_cart_coords(&self, pos_n: Vec3, n_orientation: Quaternion) -> BackBoneCoords {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.
        let bond_prev = Vec3::new(-1., 1., 1.).to_normalized();
        let bond_next = Vec3::new(1., 1., -1.).to_normalized();

        // The other bonds, ie 2xH for C' and N, and H + side chain for Cα.
        let _bond_a = Vec3::new(1., -1., 1.).to_normalized();
        let _bond_b = Vec3::new(-1., -1., -1.).to_normalized();

        // todo: You must anchor the dihedral angle in terms of the plane formed by the atom pairs
        // todo on each side of the angle; that's how it's defined.
        let cα_orientation = find_backbone_atom_orientation(n_orientation, bond_prev, self.ψ);
        let cp_orientation = find_backbone_atom_orientation(cα_orientation, bond_prev, self.φ);
        let n_next_orientation = find_backbone_atom_orientation(cp_orientation, bond_prev, self.ω);

        // Orientations for the various atoms, Z up. todo ?

        // Calculate the orientation of each atom from the orientation, the angle of the
        // bond from the previous atom, and the dihedral angle between the two.

        // Set up atom-oriented vectors encoding the angle and len of each bond.
        let bond_n_cα = bond_next * LEN_N_CALPHA;
        let bond_cα_cp = bond_next * LEN_CALPHA_CP;
        let bond_cp_n = bond_next * LEN_CP_N;

        // Calculate the position of each atom from the position, orientation and bond angle
        // from the previous atom.
        // todo: Start by describing each term in terms of its dependencies.
        let cα = pos_n + n_orientation.rotate_vec(bond_n_cα);
        let cp = cα + cα_orientation.rotate_vec(bond_cα_cp);
        let n_next = cp + cp_orientation.rotate_vec(bond_cp_n);

        BackBoneCoords {
            cα,
            cp,
            n_next,
            cα_orientation,
            cp_orientation,
            n_next_orientation,
        }
    }
}

/// Set up our protein; passed to our initial render state.
fn init_protein() -> Vec<AaInProtein> {
    let a = AaInProtein {
        aa: AminoAcid::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.6 * TAU,
        ψ: 0.6 * TAU,
    };

    let b = AaInProtein {
        aa: AminoAcid::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.75 * TAU,
        ψ: 0.25 * TAU,
    };

    let c = AaInProtein {
        aa: AminoAcid::A,
        ω: TAU / 2., // ω Assumed to be TAU/2 for most cases
        φ: 0.2 * TAU,
        ψ: 0. * TAU,
    };

    vec![a]
}

fn main() {
    render::run();
}
