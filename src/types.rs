//! This module contains structs that describe our proteins and render states.

use std::{f64::consts::TAU, fmt};

use crate::{
    atom_coords::ProteinCoords, gui::UiMode, proteins, sidechain::Sidechain,
    water::WaterEnvironment,
};

use lin_alg2::f64::{Quaternion, Vec3};
use crate::wf_lab::{Proton, WaveFunctionState};

#[derive(Debug)]
/// A protein defined by AminoAcids: Name and bond angle.
pub struct ProteinDescription {
    pub name: String,
    pub pdb_ident: String,
    pub residues: Vec<Residue>,
    /// We may need to adjust anchor position and orientation in response to forces on
    /// atoms in the chain. Eg on the anchor atom itself, or propogated torque from forces
    /// anywhere on the chain.
    pub anchor_n_posit: Vec3,
    pub anchor_n_orientation: Quaternion,
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

/// Information about a single protein.
pub struct ProteinState {
    /// Descriptions of each amino acid, including its name, and bond angles.
    pub descrip: ProteinDescription,
    /// Cartesian coordinates; generated from the description. Calculated in `coord_gen`.
    pub coords: ProteinCoords,
}

impl ProteinState {
    /// Change the amino acid sequence.
    pub fn change_sequence(
        &mut self,
        name: Option<String>,
        pdb_ident: Option<String>,
        residues: Option<Vec<Residue>>,
    ) {
        if let Some(n) = name {
            self.descrip.name = n;
        }

        if let Some(p) = pdb_ident {
            self.descrip.pdb_ident = p
        }

        if let Some(r) = residues {
            self.descrip.residues = r;
        }
    }
}

pub struct UiState {
    /// Residue id that's selected for rotation. Starts at 1. (Proteins)
    pub active_residue: usize,
    pub ui_mode: UiMode,
    // pub ui: StateUi,
    /// Simulation temperature. Nominally K, but currently unhinged from real time etc.
    pub temperature: f64,
    /// Ratio of sim time to real time. A higher value (closer to 1) is faster.
    pub sim_time_scale: f64,
    pub sim_running: bool,
    pub show_sidechains: bool,
    pub show_hydrogens: bool,
    pub show_water_molecules: bool,
}

/// Store our atom descriptions here, for global state the renderer can access.
pub struct State {
    pub protein: ProteinState, // todo: As required, allow multiple proteins.
    pub ui: UiState,
    pub water_env: WaterEnvironment,
    pub wavefunction_lab: WaveFunctionState,
}

impl State {
    /// Create state from a protein description
    pub fn new(protein_descrip: ProteinDescription) -> Self {
        let protein_coords = ProteinCoords::from_descrip(&protein_descrip);

        let protein = ProteinState {
            descrip: protein_descrip,
            coords: protein_coords,
        };

        let ui = UiState {
            active_residue: 1,
            ui_mode: UiMode::ActiveAaEditor,
            // ui: Default::default(),
            temperature: 293.,
            sim_time_scale: 0.001,
            sim_running: false,
            show_sidechains: true,
            show_hydrogens: true,
            show_water_molecules: true,
        };

        let mut wavefunction_lab = WaveFunctionState {
            protons: vec![
                Proton {
                    position: Vec3::new(5., -10., 0.),
                    velocity: Vec3::new_zero(),
                },
                Proton {
                    position: Vec3::new(5., -14., 0.),
                    velocity: Vec3::new_zero(),
                }
            ],
            electron_centers: vec![
                Vec3::new(5., -10., 0.),
                Vec3::new(5., -14., 0.),
            ],
            electron_posits_dynamic: Vec::new(),
        };
        wavefunction_lab.update_posits(0.);

        Self {
            protein,
            ui,
            water_env: WaterEnvironment::build(crate::water::N_MOLECULES, 308.15),
            wavefunction_lab,
        }
    }
}

/// An amino acid in a protein structure, including all dihedral angles required to determine
/// the conformation. Includes backbone and side chain dihedral angles. Doesn't store coordinates,
/// but coordinates can be generated using forward kinematics from the angles.
#[derive(Debug)]
pub struct Residue {
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
    /// Contains the χ angles that define t
    pub sidechain: Sidechain,
    // pub dipole: Vec3,
}

impl Default for Residue {
    /// This is what new AAs are added as.
    fn default() -> Self {
        Self {
            ω: TAU / 2.,
            // φ: proteins::PHI_SHEET,
            // ψ: proteins::PSI_SHEET,
            φ: proteins::PHI_HELIX,
            ψ: proteins::PSI_HELIX,
            sidechain: Sidechain::Arg(Default::default()),
            // dipole: Vec3::new_zero(),
        }
    }
}

impl fmt::Display for Residue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\nω: {:.2}τ, φ: {:.2}τ ψ: {:.2}τ\n",
            self.sidechain,
            self.ω / TAU,
            self.φ / TAU,
            self.ψ / TAU
        )
    }
}

impl Residue {
    pub fn new(ω: f64, φ: f64, ψ: f64, sidechain: Sidechain) -> Self {
        Self {
            ω,
            φ,
            ψ,
            sidechain,
            // dipole: Vec3::new_zero(),
        }
    }
}
