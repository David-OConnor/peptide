//! This module contains descriptions of known proteins.

use core::f64::consts::TAU;

use crate::{
    kinematics::{ProteinDescription, Residue},
    sidechain::Sidechain,
};

///https://www.rcsb.org/sequence/1L2Y
pub fn make_trp_cage() -> ProteinDescription {
    let dihed = 1. / 2. * TAU; // todo temp

    let r1 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Asn(Default::default())
    );

    let r2 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Leu(Default::default())
    );

    let r3 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Tyr(Default::default())
    );

    let r4 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Ile(Default::default())
    );

    let r5 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Gln(Default::default())
    );

    let r6 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Trp(Default::default()) // todo: Impl this.
    );

    let r7 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Leu(Default::default())
    );

    let r8 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Lys(Default::default())
    );

    let r9 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Asp(Default::default())
    );

    let r10 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Gly(Default::default())
    );

    let r11 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Gly(Default::default())
    );

    let r12 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Pro(Default::default())
    );

    let r13 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Ser(Default::default())
    );

    let r14 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Ser(Default::default())
    );

    let r15 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Gly(Default::default())
    );

    let r16 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Arg(Default::default())
    );

    let r17 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Pro(Default::default())
    );

    let r18 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Pro(Default::default())
    );

    let r19 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Pro(Default::default())
    );

    let r20 = Residue::new(
        dihed,
        dihed,
        dihed,
        Sidechain::Ser(Default::default())
    );


    ProteinDescription {
        residues: vec![
            r1, r2, r3, r4, r5, r6, r7, r8, r9, r10,
            r11, r12, r13, r14, r15, r16, r17, r18, r19, r20
        ]
    }
}