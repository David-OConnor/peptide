//! This module contains respresentations of atom coords.

use crate::{
    chem_definitions::BackboneRole,
    kinematics::ProteinDescription,
    lin_alg::{Quaternion, Vec3},
    sidechain::Sidechain,
};

/// Location of an atom, in the worldspace coordinate system.
#[derive(Debug)]
pub struct AtomCoords {
    /// id of the Amino Acid this atom is part of
    pub residue_id: usize, // todo: Do we want this id, or use an index?
    pub role: BackboneRole,
    pub position: Vec3,
    pub orientation: Quaternion,
}

#[derive(Debug)]
/// Describes the sequence of atoms that make up a protein backbone, with worldspace coordinates.
/// this is what is needed for the render, and spacial manipulations.
pub struct ProteinCoords {
    pub atoms_backbone: Vec<AtomCoords>,
}

/// Helper function to reduce code repetition
fn add_atom(
    role: BackboneRole,
    position: Vec3,
    orientation: Quaternion,
    backbone: &mut Vec<AtomCoords>,
    id: &mut usize,
) {
    backbone.push(AtomCoords {
        residue_id: *id,
        role,
        position,
        orientation,
    });

    *id += 1;
}

impl ProteinCoords {
    /// Creates coordinates from bond angles etc.
    pub fn from_descrip(descrip: &ProteinDescription) -> Self {
        // todo: Only update downstream atoms for a given rotation.
        let mut backbone = Vec::new();

        let mut id = 0;

        // N-terminus nitrogen, at the *start* of our chain. This is our anchor atom, with 0 position,
        // and an identity-quaternion orientation.
        let n_anchor = AtomCoords {
            residue_id: id,
            role: BackboneRole::N,
            position: Vec3::new(0., 0., 0.),
            orientation: Quaternion::new_identity(),
        };

        // Store these values, to anchor each successive residue to the previous. We update them
        // after each residue.
        let mut prev_n_posit = n_anchor.position;
        let mut prev_n_or = n_anchor.orientation;
        // todo: This may need adjustment. to match physical reality.
        // this position affects the first dihedral angle.
        let mut prev_cp_posit = Vec3::new(1., 1., 0.).to_normalized();

        // We store c_alpha posit and orientation for anchoring the sidechain
        let mut c_alpha_posit = Vec3::new_zero();
        let mut c_alpha_or = Quaternion::new_identity();

        backbone.push(n_anchor);
        id += 1;

        for res in &descrip.residues {
            let bb_coords =
                res.backbone_cart_coords(prev_n_posit, prev_n_or, prev_cp_posit);

            add_atom(
                BackboneRole::Cα,
                bb_coords.cα,
                bb_coords.cα_orientation,
                &mut backbone,
                &mut id,
            );

            c_alpha_posit = backbone[id - 1].position;
            c_alpha_or = backbone[id - 1].orientation;

            add_atom(
                BackboneRole::Cp,
                bb_coords.cp,
                bb_coords.cp_orientation,
                &mut backbone,
                &mut id,
            );

            prev_cp_posit = backbone[id - 1].position;

            add_atom(
                BackboneRole::N,
                bb_coords.n_next,
                bb_coords.n_next_orientation,
                &mut backbone,
                &mut id,
            );

            prev_n_posit = backbone[id - 1].position;
            prev_n_or = backbone[id - 1].orientation;

            add_atom(
                BackboneRole::O,
                bb_coords.o,
                bb_coords.o_orientation,
                &mut backbone,
                &mut id,
            );

            // todo: Consider moving this elsewhere, since it's verbose.
            // Add sidechains
            match &res.sidechain {
                Sidechain::Arg(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        BackboneRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        &mut id,
                    );
                    add_atom(
                        BackboneRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        &mut id,
                    );
                    add_atom(
                        BackboneRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        &mut id,
                    );
                    add_atom(
                        BackboneRole::NSidechain,
                        sc_coords.n_eps,
                        sc_coords.n_eps_orientation,
                        &mut backbone,
                        &mut id,
                    );
                    add_atom(
                        BackboneRole::CSidechain,
                        sc_coords.c_zeta,
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        &mut id,
                    );
                    add_atom(
                        BackboneRole::NSidechain,
                        sc_coords.n_eta1,
                        Quaternion::new_identity(),
                        &mut backbone,
                        &mut id,
                    );
                    add_atom(
                        BackboneRole::NSidechain,
                        sc_coords.n_eta2,
                        Quaternion::new_identity(),
                        &mut backbone,
                        &mut id,
                    );
                }
                _ => (), // todo; rest are not implemented yet in `sidechains` module
            }
        }

        Self {
            atoms_backbone: backbone,
        }
    }
}
