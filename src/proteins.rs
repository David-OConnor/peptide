//! This module contains descriptions of known proteins.

use core::f64::consts::TAU;

use crate::{
    kinematics::{ProteinDescription, Residue},
    sidechain::Sidechain,
};

fn to_rad(deg: f64) -> f64 {
    deg * TAU / 360.
}

// Pymol, to get dihedral angles:

// omega: get_dihedral 1/ca,1/c,2/n,2/ca

// phi: get_dihedral 1/c,2/n,2/ca,2/c
// psi: get_dihedral 1/n,1/ca,1/c,2/n

// Get all backbone dihedral angles: phi_psi 1l2y

/// https://www.rcsb.org/sequence/1L2Y
pub fn make_trp_cage() -> ProteinDescription {
    let dihed = 1. / 2. * TAU; // todo temp

    let r1 = Residue::new(
        dihed,
        to_rad(0.),
        to_rad(0.),
        Sidechain::Asn(Default::default()),
    );

    let r2 = Residue::new(
        dihed,
        to_rad(-44.),
        to_rad(-51.3),
        Sidechain::Leu(Default::default()),
    );

    let r3 = Residue::new(
        dihed,
        to_rad(-66.5),
        to_rad(-30.9),
        Sidechain::Tyr(Default::default()),
    );

    let r4 = Residue::new(
        dihed,
        to_rad(-65.2),
        to_rad(-45.9),
        Sidechain::Ile(Default::default()),
    );

    let r5 = Residue::new(
        dihed,
        to_rad(-64.7),
        to_rad(-30.3),
        Sidechain::Gln(Default::default()),
    );

    let r6 = Residue::new(
        dihed,
        to_rad(-73.1),
        to_rad(-43.4),
        Sidechain::Trp(Default::default()), // todo: Impl this.
    );

    let r7 = Residue::new(
        dihed,
        to_rad(-64.9),
        to_rad(-43.3),
        Sidechain::Leu(Default::default()),
    );

    let r8 = Residue::new(
        dihed,
        to_rad(-59.5),
        to_rad(-25.7),
        Sidechain::Lys(Default::default()),
    );

    let r9 = Residue::new(
        dihed,
        to_rad(-78.),
        to_rad(-8.8),
        Sidechain::Asp(Default::default()),
    );

    let r10 = Residue::new(
        dihed,
        to_rad(110.8),
        to_rad(8.1),
        Sidechain::Gly(Default::default()),
    );

    let r11 = Residue::new(
        dihed,
        to_rad(55.2),
        to_rad(-124.4),
        Sidechain::Gly(Default::default()),
    );

    let r12 = Residue::new(
        dihed,
        to_rad(-58.),
        to_rad(-28.8),
        Sidechain::Pro(Default::default()),
    );

    let r13 = Residue::new(
        dihed,
        to_rad(-81.8),
        to_rad(19.1),
        Sidechain::Ser(Default::default()),
    );

    let r14 = Residue::new(
        dihed,
        to_rad(-124.1),
        to_rad(13.4),
        Sidechain::Ser(Default::default()),
    );

    let r15 = Residue::new(
        dihed,
        to_rad(67.9),
        to_rad(25.2),
        Sidechain::Gly(Default::default()),
    );

    let r16 = Residue::new(
        dihed,
        to_rad(-144.),
        to_rad(131.3),
        Sidechain::Arg(Default::default()),
    );

    let r17 = Residue::new(
        dihed,
        to_rad(-70.1),
        to_rad(160.1),
        Sidechain::Pro(Default::default()),
    );

    let r18 = Residue::new(
        dihed,
        to_rad(-69.5),
        to_rad(145.7),
        Sidechain::Pro(Default::default()),
    );

    let r19 = Residue::new(
        dihed,
        to_rad(-77.3),
        to_rad(124.2),
        Sidechain::Pro(Default::default()),
    );

    let r20 = Residue::new(
        dihed,
        to_rad(0.),
        to_rad(0.),
        Sidechain::Ser(Default::default()),
    );

    ProteinDescription {
        name: "Trp-Cage".to_owned(),
        pdb_ident: "1L2Y".to_owned(),
        residues: vec![
            r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19,
            r20,
        ],
    }
}
