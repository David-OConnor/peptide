//! This module contains structs that describe our proteins and render states.

use std::fmt;

use crate::{atom_coords::ProteinCoords, kinematics::Residue, render::Camera};

use lin_alg2::f64::{Quaternion, Vec3};

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub name: String,
    pub pdb_ident: String,
    pub residues: Vec<Residue>,
}

impl fmt::Display for ProteinDescription {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "\nName: {} PDB ident: {}, Residues:\n\n",
            self.name, self.pdb_ident
        )?;

        for (i, residue) in self.residues.iter().enumerate() {
            write!(f, "Id {} - {}\n", i + 1, residue)?;
        }

        Ok(())
    }
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
pub struct BackboneCoords {
    /// Cα
    pub cα: Vec3,
    /// Carbon' atom bound to Cα
    pub cp: Vec3,
    /// Nitrogen atom of the next module.
    pub n_next: Vec3,
    /// Oxygen atom bonded to C'
    pub o: Vec3,
    /// HYdrogen atom bonded to N.
    pub h_cα: Vec3,
    pub h_n: Vec3,
    pub cα_orientation: Quaternion,
    pub cp_orientation: Quaternion,
    pub n_next_orientation: Quaternion,
    // todo: Do we want h and o orientations/
    pub o_orientation: Quaternion,
    pub h_cα_orientation: Quaternion,
    pub h_n_orientation: Quaternion,
}

/// Describes a water molecule. These aren't directly part of a protein, but may play a role in its
/// folding, among other potential roles.
#[derive(Debug)]
// todo: Consider if you want this to be a struct, a const of some other struct etc.
pub struct _WaterMolecule {
    /// Worldspace coordinates of the O atom.
    position_o_world: Vec3,
    /// Using the same orientation ref as protein atoms.
    orientation: Quaternion,
}

/// Store our atom descriptions here, for global state the renderer can access.
pub struct State {
    /// Descriptions of each amino acid, including its name, and bond angles.
    pub protein_descrip: ProteinDescription,
    /// Stored coordinates, calculated in `coord_gen`.
    pub protein_coords: ProteinCoords,
    /// Residue id that's selected for rotation. Starts at 1.
    pub active_residue: usize,
    /// Camera position and orientation
    /// todo: DO we want this? Probably not.
    pub cam: Camera,
}

impl Default for State {
    fn default() -> Self {
        Self {
            protein_descrip: ProteinDescription {
                name: "".to_owned(),
                pdb_ident: "".to_owned(),
                residues: Vec::new(),
            },
            protein_coords: ProteinCoords {
                atoms_backbone: Vec::new(),
            },
            active_residue: 1,
            cam: Camera {
                position: Vec3::new(0., 0., 7.),
                orientation: Quaternion::new_identity(),
            },
        }
    }
}

impl State {
    /// Change the amino acid sequence.
    pub fn change_sequence(
        &mut self,
        name: Option<String>,
        pdb_ident: Option<String>,
        residues: Option<Vec<Residue>>,
    ) {
        if let Some(n) = name {
            self.protein_descrip.name = n;
        }

        if let Some(p) = pdb_ident {
            self.protein_descrip.pdb_ident = p
        }

        if let Some(r) = residues {
            self.protein_descrip.residues = r;
        }
    }
}