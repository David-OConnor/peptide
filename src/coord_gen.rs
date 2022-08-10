//! This module contains code for calculating atom coordinates.

// todo: QC quaternion error creep, and re-normalize A/R

use core::f64::consts::TAU;

use crate::{
    chem_definitions::{AminoAcidType, BackboneRole},
    lin_alg::{Quaternion, Vec3},
};

// Double bond len of C' to N.
pub const LEN_CP_N: f64 = 1.33; // angstrom
pub const LEN_N_CALPHA: f64 = 1.46; // angstrom
pub const LEN_CALPHA_CP: f64 = 1.53; // angstrom

pub const LEN_CP_O: f64 = 1.0; // angstrom // todo!

// Ideal bond angles. There are an approximation; from averages. Consider replacing with something
// more robust later. All angles are in radians. We use degrees with math to match common sources.
// R indicates the side chain.
const BOND_ANGLE_N: f64 = 121.7 * 360. / TAU; // This is to Calpha and C'.

// Bond from the Calpha atom
const BOND_ANGLE_CALPHA_N_R: f64 = 110.6 * 360. / TAU;
const BOND_ANGLE_CALPHA_R_CP: f64 = 110.6 * 360. / TAU;
const BOND_ANGLE_CALPHA_N_CP: f64 = 111.0 * 360. / TAU;

// Bonds from the C' atom
const BOND_ANGLE_CP_CALPHA_O: f64 = 120.1 * 360. / TAU;
const BOND_ANGLE_CP_CALPHA_N: f64 = 117.2 * 360. / TAU;
const BOND_ANGLE_CP_O_N: f64 = 122.7 * 360. / TAU;

// todo temp
const V: f64 = 0.57735027;

const CALPHA_N_BOND: Vec3 = Vec3 { x: V, y: -V, z: -V };

const CALPHA_CP_BOND: Vec3 = Vec3 { x: V, y: V, z: V };

const CALPHA_R_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

const CP_CALPHA_BOND: Vec3 = Vec3 { x: V, y: -V, z: -V };

const CP_N_BOND: Vec3 = Vec3 { x: V, y: V, z: V };

const CP_O_BOND: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 0.,
};

const N_CALPHA_BOND: Vec3 = Vec3 { x: V, y: V, z: V };
const N_CP_BOND: Vec3 = Vec3 { x: V, y: -V, z: -V };

// 360 degrees in tau rad
// rad = 360 /tau degrees

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub residues: Vec<Residue>,
}

#[derive(Debug)]
/// Describes the sequence of atoms that make up a protein backbone, with worldspace coordinates.
/// this is what is needed for the render, and spacial manipulations.
pub struct ProteinCoords {
    pub atoms_backbone: Vec<AtomCoords>,
}

impl ProteinCoords {
    pub fn from_descrip(descrip: &ProteinDescription) -> Self {
        let mut backbone = Vec::new();

        let mut id = 0;

        // N-terminus nitrogen, at the *start* of our chain.
        let starting_n = AtomCoords {
            residue_id: id,
            role: BackboneRole::N,
            position: Vec3::new(0., 0., 0.),
            orientation: Quaternion::new_identity(),
        };

        let mut prev_position = starting_n.position;
        // todo: This may need adjustment. to match physical reality.
        // this position affects the first dihedral angle.
        let prev_2_position = Vec3::new(1., 0., 0.);
        let mut prev_orientation = starting_n.orientation;

        backbone.push(starting_n);
        id += 1;

        for aa in &descrip.residues {
            let aa_coords =
                aa.backbone_cart_coords(prev_position, prev_orientation, prev_2_position);

            backbone.push(AtomCoords {
                residue_id: id,
                role: BackboneRole::Cα,
                position: aa_coords.cα,
                orientation: aa_coords.cα_orientation,
            });
            id += 1;

            backbone.push(AtomCoords {
                residue_id: id,
                role: BackboneRole::Cp,
                position: aa_coords.cp,
                orientation: aa_coords.cp_orientation,
            });
            id += 1;

            let n = AtomCoords {
                residue_id: id,
                role: BackboneRole::N,
                position: aa_coords.n_next,
                orientation: aa_coords.n_next_orientation,
            };

            prev_position = n.position;
            prev_orientation = n.orientation;

            backbone.push(n);
            id += 1;

            let o = AtomCoords {
                residue_id: id,
                role: BackboneRole::O,
                position: aa_coords.o,
                orientation: Quaternion::new_identity(), // todo temp?
            };

            backbone.push(o);
            id += 1;
        }

        Self {
            atoms_backbone: backbone,
        }
    }
}

/// Location of an atom, in the worldspace coordinate system.
#[derive(Debug)]
pub struct AtomCoords {
    /// id of the Amino Acid this atom is part of
    pub residue_id: usize, // todo: Do we want this id, or use an index?
    pub role: BackboneRole,
    pub position: Vec3,
    pub orientation: Quaternion,
}

/// Holds backbone atom coordinates and orientations, relative to the alpha carbon, for
/// a single amino acid. Generated from an AA's bond angles.
///
/// We will use the convention of the chain starting at NH2 (amine) end (N terminus), and ending in
/// COOH (carboxyl) (C terminus).
/// Flow is N -> C_alpha -> C'.
/// Coordinates are in relation to the Nitrogen molecule at the start of the segement, using the
/// directional convention described above. (todo: How do we orient axes?)
/// todo: To start, convention is the previous C' to this AA's starting A is on the positive X axis.
#[derive(Debug)]
pub struct BackboneCoordsAa {
    /// Cα
    pub cα: Vec3,
    /// Carbon' atom bound to C_alpha
    pub cp: Vec3,
    /// Nitrogen atom of the next module.
    pub n_next: Vec3,
    /// Oxygen atom bonded to c'
    pub o: Vec3,
    pub cα_orientation: Quaternion,
    pub cp_orientation: Quaternion,
    pub n_next_orientation: Quaternion,
}

/// Calculate the dihedral angle between 4 atoms.
fn calc_dihedral_angle(bond_middle: Vec3, bond_adjacent1: Vec3, bond_adjacent2: Vec3) -> f64 {
    // Project the next and previous bonds onto the plane that has this bond as its normal.
    let bond1_on_plane = bond_adjacent1.project_to_plane(bond_middle).to_normalized();
    let bond2_on_plane = bond_adjacent2.project_to_plane(bond_middle).to_normalized();

    println!(
        "\nPlane: {:?}\n1: {:?}, \n2: {:?}\n",
        bond_middle, bond_adjacent1, bond_adjacent2
    );

    (bond1_on_plane.dot(bond2_on_plane)).acos();
}

/// Calculate the orientation, as a quaternion, of a backbone atom, given the orientation of the
/// previous atom, and the bond angle. `bond_angle` is the vector representing the bond to
/// this atom from the previous atom's orientation. `bond_prev` and `bond_next` are in the atom's
/// coordinates; not worldspace.
pub fn find_backbone_atom_orientation(
    q_prev: Quaternion,
    bond_to_prev: Vec3,         // Local space
    bond_to_next: Vec3,         // Local space
    prev_bond_worldspace: Vec3, // World space
    dihedral_angle: f64,
) -> Quaternion {
    // #1: Align the prev atom's bond vector to world space based on the prev atom's orientation.
    let bond_to_this_worldspace = q_prev.rotate_vec(bond_to_next);

    // // #2: Align this bond's vec to the prev to worldspace, as above
    // let bond_to_prev_worldspace = q_prev.rotate_vec(bond_to_prev); // todo figure out if you wnat this step!

    // #2: Find the rotation quaternion that aligns the (inverse of) the local(world?)-space bond to
    // the prev atom with the world-space "to" bond of the previous atom. This is also the
    // orientation of our atom, without applying the dihedral angle.
    let bond_alignment_rotation =
        Quaternion::from_unit_vecs(bond_to_prev * -1., bond_to_this_worldspace);
    // Quaternion::from_unit_vecs(bond_to_prev_worldspace * -1., bond_to_this_worldspace); // todo local or world for first arg?

    // #3: Rotate the orientation around the dihedral angle. We must do this so that our
    // dihedral angle is in relation to the previous and next bonds.
    // Adjust the dihedral angle to be in reference to the previous 2 atoms, per the convention.
    let next_bond_worldspace = bond_alignment_rotation.rotate_vec(bond_to_next);

    let dihedral_angle_current = calc_dihedral_angle(
        bond_to_this_worldspace,
        prev_bond_worldspace,
        next_bond_worldspace,
    );

    let dihedral_rotation = Quaternion::from_axis_angle(
        bond_to_this_worldspace,
        // dihedral_angle - dihedral_angle_current,
        // dihedral_angle, // todo: Temp while TS.
        dihedral_angle_current - dihedral_angle,
    );

    let result = dihedral_rotation * bond_alignment_rotation;

    // todo start checking diihedral angle

    let prev_bond_worldspace = result.rotate_vec(bond_to_prev);
    let next_bond_worldspace = result.rotate_vec(bond_to_next);

    let prev_bond_on_plane = prev_bond_worldspace
        .project_to_plane(bond_to_this_worldspace)
        .to_normalized();
    let next_bond_on_plane = next_bond_worldspace
        .project_to_plane(bond_to_this_worldspace)
        .to_normalized();
    let dihedral_angle_current = (next_bond_on_plane.dot(prev_bond_on_plane)).acos(); // todo temp for testing if our algo worked

    println!("Dihedral angle expected: {dihedral_angle}, Actual: {dihedral_angle_current}");
    // todo end checking dihedral angle

    result
}

/// In this approach we model atoms with a vector indicating their position, and we model
/// each bond as a quaternion.
pub fn find_backbone_atom_orientation2(
    bond_to_this: Quaternion,
    // bond_to_this_worldspace: Vec3,
    bond_angle: f64, // The angle between the incoming and outgoing bond.
    dihedral_angle: f64,
) -> Quaternion {
    let to_this_vec = bond_to_this.to_vec().to_normalized();

    let v = Vec3::new(1., 0., 0.); // todo: Arbitrary, fornow.
    let normal_to_bond = to_this_vec.cross(v);
    let rotation_to_next_bond = Quaternion::from_axis_angle(normal_to_bond, bond_angle);

    println!("BTT: {:?}, This to vec: {:?}", bond_to_this, to_this_vec);

    let angle = dihedral_angle; // todo
    let bond_to_this_worldspace = to_this_vec; // todo is this right?

    let rotation_around_axis = Quaternion::from_axis_angle(bond_to_this_worldspace, angle);

    println!(
        "TEST {:?}, {:?}, {:?}",
        rotation_around_axis, rotation_to_next_bond, bond_to_this
    );

    rotation_around_axis * rotation_to_next_bond * bond_to_this
}

/// An amino acid in a protein structure, including position information.
#[derive(Debug)]
pub struct Residue {
    pub aa: AminoAcidType,
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

impl Residue {
    /// Generate cartesian coordinates of points from diahedral angles and bond lengths. Starts with
    /// N, and ends with C'.
    /// Accepts position, and orientation of the N atom that starts this segment.
    /// Also returns orientation of the N atom, for use when calculating coordinates for the next
    /// AA in the chain.
    /// `prev_cp_pos` is usd to anchor the dihedral angle properly, since it's defined by planes
    ///  of 3 atoms.
    pub fn backbone_cart_coords(
        &self,
        n_pos: Vec3,
        n_orientation: Quaternion,
        prev_cp_pos: Vec3,
    ) -> BackboneCoordsAa {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.

        // todo: You must anchor the dihedral angle in terms of the plane formed by the atom pairs
        // todo on each side of the angle; that's how it's defined.
        let cα_orientation = find_backbone_atom_orientation(
            n_orientation,
            CALPHA_N_BOND,
            CALPHA_CP_BOND,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            n_pos - prev_cp_pos,
            self.φ,
        );

        let cα = n_pos + n_orientation.rotate_vec(N_CALPHA_BOND);

        let cp_orientation = find_backbone_atom_orientation(
            cα_orientation,
            CP_CALPHA_BOND,
            CP_N_BOND,
            cα - n_pos,
            self.ψ,
        );

        let cp = cα + cα_orientation.rotate_vec(CP_CALPHA_BOND);

        let n_next_orientation = find_backbone_atom_orientation(
            cp_orientation,
            N_CP_BOND,
            N_CALPHA_BOND,
            cp - cα,
            self.ω,
        );

        let n_next = cp + cp_orientation.rotate_vec(N_CP_BOND);

        let o = cp + cp_orientation.rotate_vec(CP_O_BOND);

        // Calculate the position of each atom from the orientation, the angle of the
        // bond from the previous atom, and the dihedral angle between the two.

        // Calculate the position of each atom from the position, orientation and bond angle
        // from the previous atom.
        // todo: Start by describing each term in terms of its dependencies.

        BackboneCoordsAa {
            cα,
            cp,
            n_next,
            o,
            cα_orientation,
            cp_orientation,
            n_next_orientation,
        }
    }

    // pub fn backbone_cart_coords_(
    //     &self,
    //     n_pos: Vec3,
    //     n_orientation: Quaternion,
    //     prev_cp_pos: Vec3,
    // ) -> BackboneCoordsAa {
    //     // Set up atom-oriented vectors encoding the angle and len of each bond.
    //     let bond_to_cα =
    //         find_backbone_atom_orientation2(n_orientation, BOND_ANGLE_CALPHA_N_CP, self.φ);
    //
    //     let bond_to_cp =
    //         find_backbone_atom_orientation2(bond_to_cα, BOND_ANGLE_CP_CALPHA_N, self.ψ);
    //
    //     let bond_to_n_next =
    //         find_backbone_atom_orientation2(bond_to_cp, BOND_ANGLE_CALPHA_N_CP, self.ω);
    //
    //     let cα = n_pos + bond_to_cα.to_vec().to_normalized() * LEN_N_CALPHA;
    //     let cp = cα + bond_to_cp.to_vec().to_normalized() * LEN_CALPHA_CP;
    //     let n_next = cp + bond_to_n_next.to_vec().to_normalized() * LEN_CP_N;
    //     let o = cp + bond_to_n_next.to_vec().to_normalized() * LEN_CP_O;
    //
    //     let cα_orientation = Quaternion::new_identity();
    //     let cp_orientation = Quaternion::new_identity();
    //     let n_next_orientation = Quaternion::new_identity();
    //
    //     BackboneCoordsAa {
    //         cα,
    //         cp,
    //         n_next,
    //         0,
    //         cα_orientation,
    //         cp_orientation,
    //         n_next_orientation,
    //     }
    // }
}

/// Used to represent one atom in a system built of atoms.
pub struct _ZmatrixItem {
    atomic_number: u8,
    bond_len: f64,       // Angstrom
    bond_angle: f64,     // radians
    dihedral_angle: f64, // radians
}
