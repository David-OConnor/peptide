//! This module contains respresentations of atom coords.

use lin_alg::f64::{Quaternion, Vec3};

use crate::{
    chem_definitions::{AminoAcidType, AtomRole},
    forces::{MAG_N, MAG_O},
    sidechain::Sidechain,
    types::ProteinDescription,
};

// todo: For chain termination AA, you may need to replace the final backbone N with an O.

// todo: Where should this go? interactions mo

pub const Q_I: Quaternion = Quaternion {
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
    /// Use use this step to find the previous atom to create the bond to. Normally 1. Set this to 2 to go back 2
    /// in the sidechain to find its partner, eg ASPs second O.
    pub sidechain_bond_step: usize,
    /// Used when this atom needs to initiate 2 bonds; eg for terminating rings, or maybe
    /// for double-bonds?
    pub second_bond_step: Option<usize>,
    /// Hydrogen-bond dipole. Length is proportional to strenght. Perhaps an Option allows us to easily
    /// skip computations.
    pub dipole: Option<Vec3>,
    /// How an atom interacts with a dipole.
    pub mag_force: Option<f64>, // todo?
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
    bond_step: usize,
    second_bond_step: Option<usize>,
    residue_id: usize,
    atom_id: &mut usize,
    dipole: Option<Vec3>,
    mag_force: Option<f64>,
) {
    // todo: Can we infer magnetic force and dipole here? Probably not without nearby H?
    backbone.push(AtomCoords {
        residue_id,
        role,
        position,
        orientation,
        sidechain_bond_step: bond_step,
        second_bond_step,
        dipole,
        mag_force,
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
            descrip.anchor_n_posit,
            descrip.anchor_n_orientation,
            &mut backbone,
            1,
            None,
            residue_id,
            &mut atom_id,
            None,
            Some(MAG_N),
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

            // Proline sees the nitrogen bonded to one of the sidechain carbons (the one that
            // completes the ring) vice a hydrogen.
            let mut cα_bond_step = 2;
            if res.sidechain.aa_type() != AminoAcidType::Pro {
                add_atom(
                    AtomRole::HN,
                    bb_coords.h_n,
                    bb_coords.h_n_orientation,
                    &mut backbone,
                    1,
                    None,
                    residue_id,
                    &mut atom_id,
                    None,
                    None,
                );
                cα_bond_step = 1;
            }

            add_atom(
                AtomRole::Cα,
                bb_coords.cα,
                bb_coords.cα_orientation,
                &mut backbone,
                cα_bond_step, // Back 2 to skip the H bonded to N.
                None,
                residue_id,
                &mut atom_id,
                None,
                None,
            );

            // c_alpha_posit = backbone[atom_id - 1].position;
            // c_alpha_or = backbone[atom_id - 1].orientation;

            add_atom(
                AtomRole::HCα,
                bb_coords.h_cα,
                bb_coords.h_cα_orientation,
                &mut backbone,
                1,
                None,
                residue_id,
                &mut atom_id,
                None,
                None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eps,
                        sc_coords.n_eps_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta,
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eta1,
                        sc_coords.n_eta1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eta2,
                        sc_coords.n_eta2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eta1_a,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eps,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eta1_b,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eta2_a,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eta2_b,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        12,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        13,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_a,
                        Q_I,
                        &mut backbone,
                        13,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_b,
                        Q_I,
                        &mut backbone,
                        14,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_a,
                        Q_I,
                        &mut backbone,
                        14,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_b,
                        Q_I,
                        &mut backbone,
                        15,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps,
                        sc_coords.c_eps_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_zeta,
                        sc_coords.n_zeta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_zeta_a,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_zeta_b,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    // todo: Amonium group with 3 Hs for Lys?? H3N+?
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_a,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_b,
                        Q_I,
                        &mut backbone,
                        9,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_a,
                        Q_I,
                        &mut backbone,
                        9,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_b,
                        Q_I,
                        &mut backbone,
                        10,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps_a,
                        Q_I,
                        &mut backbone,
                        10,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps_b,
                        Q_I,
                        &mut backbone,
                        11,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_delta1,
                        sc_coords.o_delta1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_delta2,
                        sc_coords.o_delta2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    // todo: H on Os?
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_gamma,
                        sc_coords.o_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_o_gamma,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma2,
                        sc_coords.c_gamma2_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_gamma1,
                        sc_coords.o_gamma1_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_o,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma2,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma3,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_delta1,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_delta2,
                        sc_coords.n_delta2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_delta_a,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_delta_b,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eps1,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );

                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eps2,
                        sc_coords.n_eps2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eps_a,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eps_b,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_a,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_b,
                        Q_I,
                        &mut backbone,
                        9,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                }
                Sidechain::Gly(angles) => {
                    let sc_coords = angles.sidechain_cart_coords(
                        bb_coords.cα,
                        bb_coords.cα_orientation,
                        prev_n_posit,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                }
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        Some(5), // to backbone N.
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_a,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_b,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_a,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_b,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma1,
                        sc_coords.c_gamma1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma2,
                        sc_coords.c_gamma2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1_a,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1_b,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1_c,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma2_a,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma2_b,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_a,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_b,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta_c,
                        Q_I,
                        &mut backbone,
                        9,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta2,
                        sc_coords.c_delta2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta1_a,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta1_b,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta1_c,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta2_a,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta2_b,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta2_c,
                        Q_I,
                        &mut backbone,
                        9,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta2,
                        sc_coords.c_delta2_orientation,
                        &mut backbone,
                        2, // to gamma
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps1,
                        sc_coords.c_eps1_orientation,
                        &mut backbone,
                        2, // to delta1
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        sc_coords.c_eps2_orientation,
                        &mut backbone,
                        2, // to delta2
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta,
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        2,       // to eps1
                        Some(1), // ep eps2)
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta1,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta2,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps1,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps2,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_zeta,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta2,
                        sc_coords.c_delta2_orientation,
                        &mut backbone,
                        2, // to gamma
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps1,
                        sc_coords.c_eps1_orientation,
                        &mut backbone,
                        2, // to delta1
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        sc_coords.c_eps2_orientation,
                        &mut backbone,
                        2, // to delta2
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta,
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        2,       // to eps1
                        Some(1), // to eps2.
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eta,
                        sc_coords.o_eta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        9,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta1,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta2,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps1,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps2,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_o_eta,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eps,
                        sc_coords.n_eps_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_zeta, // Border between rings
                        sc_coords.c_zeta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eta, // Border between rings
                        sc_coords.c_eta_orientation,
                        &mut backbone,
                        1,
                        Some(4), // Connect to C gamma.
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_theta,
                        sc_coords.c_theta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_iota,
                        sc_coords.c_iota_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_kappa,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_lambda,
                        Q_I,
                        &mut backbone,
                        1,
                        // Connect to the hinge carbon that's bonded to N.
                        Some(5),
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        10,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        11,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta,
                        Q_I,
                        &mut backbone,
                        10,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_eps,
                        Q_I,
                        &mut backbone,
                        10,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_theta,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_iota,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_kappa,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_lambda,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::SSidechain,
                        sc_coords.s_gamma,
                        sc_coords.s_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_s_gamma,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::SeSidechain,
                        sc_coords.se_gamma,
                        sc_coords.se_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::SSidechain,
                        sc_coords.s_delta,
                        sc_coords.s_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps,
                        sc_coords.c_eps_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_a,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_b,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps_a,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps_b,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps_c,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta1,
                        sc_coords.c_delta1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_delta2,
                        sc_coords.n_delta2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::NSidechain,
                        sc_coords.n_eps1,
                        sc_coords.n_eps1_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_N),
                    );
                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_eps2,
                        sc_coords.c_eps2_orientation,
                        &mut backbone,
                        2,
                        Some(1),
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_n_delta,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                    // add_atom(
                    //     AtomRole::HSidechain,
                    //     sc_coords.h_n_eps,
                    //     Q_I,
                    //     &mut backbone,
                    //     3,
                    //     None, residue_id,
                    //     &mut atom_id,
                    // );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        8,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_delta1,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_eps2,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma,
                        sc_coords.c_gamma_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_delta,
                        sc_coords.c_delta_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eps1,
                        sc_coords.o_eps1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );

                    add_atom(
                        AtomRole::OSidechain,
                        sc_coords.o_eps2,
                        sc_coords.o_eps2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        Some(MAG_O),
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_a,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma_b,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_a,
                        Q_I,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_b,
                        Q_I,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta_c,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
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
                        sc_coords.c_beta_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma1,
                        sc_coords.c_gamma1_orientation,
                        &mut backbone,
                        1,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::CSidechain,
                        sc_coords.c_gamma2,
                        sc_coords.c_gamma2_orientation,
                        &mut backbone,
                        2,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_beta,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1_a,
                        Q_I,
                        &mut backbone,
                        3,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1_b,
                        Q_I,
                        &mut backbone,
                        4,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma1_c,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma2_a,
                        Q_I,
                        &mut backbone,
                        5,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma2_b,
                        Q_I,
                        &mut backbone,
                        6,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );

                    add_atom(
                        AtomRole::HSidechain,
                        sc_coords.h_c_gamma2_c,
                        Q_I,
                        &mut backbone,
                        7,
                        None,
                        residue_id,
                        &mut atom_id,
                        None,
                        None,
                    );
                }
            }

            add_atom(
                AtomRole::Cp,
                bb_coords.cp,
                bb_coords.cp_orientation,
                &mut backbone,
                1,
                None,
                residue_id,
                &mut atom_id,
                None,
                None,
            );

            prev_cp_posit = backbone[atom_id - 1].position;

            add_atom(
                AtomRole::O,
                bb_coords.o,
                bb_coords.o_orientation,
                &mut backbone,
                1,
                None,
                residue_id,
                &mut atom_id,
                None,
                Some(MAG_O),
            );

            add_atom(
                AtomRole::N,
                bb_coords.n_next,
                bb_coords.n_next_orientation,
                &mut backbone,
                1,
                None,
                residue_id,
                &mut atom_id,
                None,
                Some(MAG_O),
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
