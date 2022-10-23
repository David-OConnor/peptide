//! This module contains respresentations of atom coords.

use lin_alg2::f64::{Quaternion, Vec3};

use crate::types::ProteinDescription;
use crate::{chem_definitions::AtomRole, sidechain::Sidechain};

const Q_I: Quaternion = Quaternion {
    w: 1.,
    x: 0.,
    y: 0.,
    z: 0.,
};

/// Location of an atom, in the worldspace coordinate system.
#[derive(Debug)]
pub struct AtomCoords {
    /// id of the Amino Acid this atom is part of
    pub residue_id: usize,
    pub role: AtomRole,
    pub position: Vec3,
    pub orientation: Quaternion,
    // Normally 1. Set this to 2 to go back 2 in the sidechain to find its partner,
    // eg ASPs second O.
    pub sidechain_bond_step: usize,
}

#[derive(Debug)]
/// Describes the sequence of atoms that make up a protein backbone, with worldspace coordinates.
/// this is what is needed for the render, and spacial manipulations.
pub struct ProteinCoords {
    pub atoms_backbone: Vec<AtomCoords>,
}

/// Helper function to reduce code repetition
fn add_atom(
    role: AtomRole,
    position: Vec3,
    orientation: Quaternion,
    backbone: &mut Vec<AtomCoords>,
    step: usize,
    residue_id: usize,
    atom_id: &mut usize,
) {
    backbone.push(AtomCoords {
        residue_id,
        role,
        position,
        orientation,
        sidechain_bond_step: step,
    });

    *atom_id += 1;
}

impl ProteinCoords {
    /// Creates coordinates from bond angles etc.
    /// The order we add them matters in terms of rendering indices, including bond arrangement.
    pub fn from_descrip(descrip: &ProteinDescription) -> Self {
        // todo: Only update downstream atoms for a given rotation.
        let mut backbone = Vec::new();

        let mut residue_id = 0;
        let mut atom_id = 0;

        // N-terminus nitrogen, at the *start* of our chain. This is our anchor atom, with 0 position,
        // and an identity-quaternion orientation.
        add_atom(
            AtomRole::N,
            Vec3::new(0., 0., 0.),
            Q_I,
            &mut backbone,
            1,
            residue_id,
            &mut atom_id,
        );

        // Start at residue id = 1 after the anchor N.
        residue_id += 1;

        // Store these values, to anchor each successive residue to the previous. We update them
        // after each residue.
        let mut prev_n_posit = backbone[atom_id - 1].position;
        let mut prev_n_or = backbone[atom_id - 1].orientation;

        // todo: This may need adjustment. to match physical reality.
        // this position affects the first dihedral angle.
        let mut prev_cp_posit = Vec3::new(1., 1., 0.).to_normalized();

        // We store c_alpha posit and orientation for anchoring the sidechain
        // let mut c_alpha_posit = Vec3::new_zero();
        // let mut c_alpha_or = Q_I;

        for res in &descrip.residues {
            let bb_coords = res.backbone_cart_coords(prev_n_posit, prev_n_or, prev_cp_posit);

            add_atom(
                AtomRole::HN,
                bb_coords.h_n,
                bb_coords.h_n_orientation,
                &mut backbone,
                1,
                residue_id,
                &mut atom_id,
            );

            add_atom(
                AtomRole::Cα,
                bb_coords.cα,
                bb_coords.cα_orientation,
                &mut backbone,
                2, // Back 2 to skip the H bonded to N.
                residue_id,
                &mut atom_id,
            );

            // c_alpha_posit = backbone[atom_id - 1].position;
            // c_alpha_or = backbone[atom_id - 1].orientation;

            add_atom(
                AtomRole::HCα,
                bb_coords.h_cα,
                bb_coords.h_cα_orientation,
                &mut backbone,
                1,
                residue_id,
                &mut atom_id,
            );

            // Add sidechains
            match &res.sidechain {
                // The first atom of each sidechain is back 2, to skip over the H bonded
                // to calpha.
                Sidechain::Arg(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eps,
                        sc_coords.n_eps_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta,
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eta1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eta2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Lys(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps,
                        sc_coords.c_eps_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_zeta,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Asp(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_delta1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_delta2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Ser(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_gamma,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Thr(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma2,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_gamma1,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Asn(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_delta1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_delta2,
                        sc_coords.n_delta2_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_amine1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_amine2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Gln(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eps1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eps2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Gly(_angles) => {} // No sidechain on Gly
                Sidechain::Pro(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Ile(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma2,
                        sc_coords.c_gamma2_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Leu(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Phe(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta2,
                        sc_coords.c_delta2_orientation,
                        &mut backbone,
                        2, // to gamma
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps1,
                        sc_coords.c_eps1_orientation,
                        &mut backbone,
                        2, // to delta1
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        sc_coords.c_eps2_orientation,
                        &mut backbone,
                        2, // to delta2
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta,
                        Q_I,
                        &mut backbone,
                        2, // to eps1
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Tyr(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta2,
                        sc_coords.c_delta2_orientation,
                        &mut backbone,
                        2, // to gamma
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps1,
                        sc_coords.c_eps1_orientation,
                        &mut backbone,
                        2, // to delta1
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        sc_coords.c_eps2_orientation,
                        &mut backbone,
                        2, // to delta2
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta,
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        2, // to eps1
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eta,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Trp(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_delta2,
                        sc_coords.n_delta2_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps1,
                        sc_coords.c_eps1_orientation,
                        &mut backbone,
                        3, // Back to gamma.
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        sc_coords.c_eps2_orientation,
                        &mut backbone,
                        1, // Connected to eps1.
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta1,
                        sc_coords.c_zeta1_orientation,
                        &mut backbone,
                        2, // Back to eps1.
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta2,
                        sc_coords.c_zeta2_orientation,
                        &mut backbone,
                        2, // Back to eps2.
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eta1,
                        Q_I,
                        &mut backbone,
                        2, // Back to zeta1
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eta2,
                        Q_I,
                        &mut backbone,
                        2, // Back to zeta2
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Cys(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::SSidechain,
                        sc_coords.s_gamma,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Sec(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::SeSidechain,
                        sc_coords.se_gamma,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Met(angles) => {
                    // todo: Why are there no kinks in the Met chain?
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::SSidechain,
                        sc_coords.s_delta,
                        sc_coords.s_delta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::His(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.n_delta2,
                        sc_coords.n_delta2_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.n_eps1,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Glu(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eps1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eps2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Ala(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
                Sidechain::Val(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_beta,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma1,
                        Q_I,
                        &mut backbone,
                        1,
                        residue_id,
                        &mut atom_id,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma2,
                        Q_I,
                        &mut backbone,
                        2,
                        residue_id,
                        &mut atom_id,
                    );
                }
            }

            add_atom(
                AtomRole::Cp,
                bb_coords.cp,
                bb_coords.cp_orientation,
                &mut backbone,
                1,
                residue_id,
                &mut atom_id,
            );

            prev_cp_posit = backbone[atom_id - 1].position;

            add_atom(
                AtomRole::O,
                bb_coords.o,
                bb_coords.o_orientation,
                &mut backbone,
                1,
                residue_id,
                &mut atom_id,
            );

            add_atom(
                AtomRole::N,
                bb_coords.n_next,
                bb_coords.n_next_orientation,
                &mut backbone,
                1,
                residue_id,
                &mut atom_id,
            );

            prev_n_posit = backbone[atom_id - 1].position;
            prev_n_or = backbone[atom_id - 1].orientation;

            residue_id += 1;
        }

        Self {
            atoms_backbone: backbone,
        }
    }
}
