//! This module contains info related to side chains, including their geometry

// Don't show warnings for un`
use std::{f64::consts::TAU, fmt};

use crate::{
    bond_vecs::*,
    chem_definitions::AminoAcidType,
    // todo: Find actual bonds / angles; these are temp!
    kinematics::find_atom_placement,
};

use lin_alg2::f64::{Quaternion, Vec3};

const TAU_DIV2: f64 = TAU / 2.;

// todo: These are temp
pub const LEN_SC: f64 = 1.53;

// Dummy bonds, for semantic clarity. They're the same, but we need both in our forward-kinematics API.
const H_BOND_IN: Vec3 = Vec3 { x: 1., y: 0., z: 0. };
const H_BOND_OUT: Vec3 = Vec3 { x: 1., y: 0., z: 0. };

// As a convention, we use generic tetrahedral and planar geometry in this module.
// This is a stopgap. Note that we are treating the generic-geometry `BOND_A` as to
// the previous atom in a chain, and `BOND_B` to the next one. For branching, we incorporate
// `BOND_C` and `BOND_D` as required. We use planar geometry as required vice tetrahedral
// when appropriate, eg for Nitrogen atoms. For bonds from hydrogens, we use the anchor bond directly.

// todo: Clean up tetra vs planar

// todo: Gauche+ and trans etc for beta C. EG opposite C' or opposite N?

#[derive(Debug, PartialEq)]
pub enum Sidechain {
    Arg(Arg),
    His(His),
    Lys(Lys),
    Asp(Asp),
    Glu(Glu),
    Ser(Ser),
    Thr(Thr),
    Asn(Asn),
    Gln(Gln),
    Cys(Cys),
    Gly(Gly),
    Pro(Pro),
    Ala(Ala),
    Val(Val),
    Ile(Ile),
    Leu(Leu),
    Met(Met),
    Phe(Phe),
    Tyr(Tyr),
    Trp(Trp),
    /// Sec is not one of the most-common 20.
    Sec(Sec),
}

impl fmt::Display for Sidechain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Arg(aa) => {
                write!(
                    f,
                    "Arg (R)\nχ1: {:.2}τ χ2: {:.2}τ χ3: {:.2}τ χ4: {:.2}τ χ5: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU,
                    aa.χ_3 / TAU,
                    aa.χ_4 / TAU,
                    aa.χ_5
                )
            }
            Self::His(aa) => {
                write!(
                    f,
                    "His (H)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Lys(aa) => {
                write!(
                    f,
                    "Lys (K)\nχ1: {:.2}τ χ2: {:.2}τ χ3: {:.2}τ χ4: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU,
                    aa.χ_3 / TAU,
                    aa.χ_4 / TAU
                )
            }
            Self::Asp(aa) => {
                write!(
                    f,
                    "Asp (D)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Glu(aa) => {
                write!(
                    f,
                    "Glu\nχ1: {:.2}τ χ2: {:.2}τ χ3: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU,
                    aa.χ_3 / TAU
                )
            }
            Self::Ser(aa) => {
                write!(f, "Ser(S)\nχ1: {:.2}τ", aa.χ_1 / TAU)
            }
            Self::Thr(aa) => {
                write!(f, "Thr (T)\nχ1: {:.2}τ", aa.χ_1 / TAU)
            }
            Self::Asn(aa) => {
                write!(
                    f,
                    "Asn (N)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Gln(aa) => {
                write!(
                    f,
                    "Gln (Q)\nχ1: {:.2}τ χ2: {:.2}τ χ3: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU,
                    aa.χ_3 / TAU
                )
            }
            Self::Cys(aa) => {
                write!(f, "Cys (C)\nχ1: {:.2}τ", aa.χ_1 / TAU)
            }
            Self::Sec(aa) => {
                write!(f, "Sec (U)\nχ1: {:.2}τ", aa.χ_1 / TAU)
            }
            Self::Gly(_aa) => {
                write!(f, "Gly (G)")
            }
            Self::Pro(aa) => {
                write!(f, "Pro (P)")
            }
            Self::Ala(_aa) => {
                write!(f, "Ala (A)")
            }
            Self::Val(aa) => {
                write!(f, "Val (V)\nχ1: {:.2}τ", aa.χ_1 / TAU)
            }
            Self::Ile(aa) => {
                write!(
                    f,
                    "Ile (I)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Leu(aa) => {
                write!(
                    f,
                    "Leu (L)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Met(aa) => {
                write!(
                    f,
                    "Met (M)\nχ1: {:.2}τ χ2: {:.2}τ χ3: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU,
                    aa.χ_3 / TAU
                )
            }
            Self::Phe(aa) => {
                write!(
                    f,
                    "Phe (F)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Tyr(aa) => {
                write!(
                    f,
                    "Tyr (Y)\nχ1: {:.2}τ χ2: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU
                )
            }
            Self::Trp(aa) => {
                write!(
                    f,
                    "Trp (W)\nχ1: {:.2}τ χ2: {:.2}τ χ3: {:.2}τ",
                    aa.χ_1 / TAU,
                    aa.χ_2 / TAU,
                    aa.χ_3 / TAU
                )
            }
        }
    }
}

impl Sidechain {
    /// Construct an AA (with default dihedral angles) from an amino acid type.
    pub fn from_aa_type(aa_type: AminoAcidType) -> Self {
        match aa_type {
            AminoAcidType::Arg => Self::Arg(Default::default()),
            AminoAcidType::His => Self::His(Default::default()),
            AminoAcidType::Lys => Self::Lys(Default::default()),
            AminoAcidType::Asp => Self::Asp(Default::default()),
            AminoAcidType::Glu => Self::Glu(Default::default()),
            AminoAcidType::Ser => Self::Ser(Default::default()),
            AminoAcidType::Thr => Self::Thr(Default::default()),
            AminoAcidType::Asn => Self::Asn(Default::default()),
            AminoAcidType::Gln => Self::Gln(Default::default()),
            AminoAcidType::Cys => Self::Cys(Default::default()),
            AminoAcidType::Sec => Self::Sec(Default::default()),
            AminoAcidType::Gly => Self::Gly(Default::default()),
            AminoAcidType::Pro => Self::Pro(Default::default()),
            AminoAcidType::Ala => Self::Ala(Default::default()),
            AminoAcidType::Val => Self::Val(Default::default()),
            AminoAcidType::Ile => Self::Ile(Default::default()),
            AminoAcidType::Leu => Self::Leu(Default::default()),
            AminoAcidType::Met => Self::Met(Default::default()),
            AminoAcidType::Phe => Self::Phe(Default::default()),
            AminoAcidType::Tyr => Self::Tyr(Default::default()),
            AminoAcidType::Trp => Self::Trp(Default::default()),
        }
    }

    /// Construct an AA (with default dihedral angles) from a single-letter identifier.
    /// Returns `None` if an invalid letter is passed.
    pub fn from_ident_single_letter(ident: &str) -> Option<Self> {
        match ident {
            "R" => Some(Self::Arg(Default::default())),
            "H" => Some(Self::His(Default::default())),
            "K" => Some(Self::Lys(Default::default())),
            "D" => Some(Self::Asp(Default::default())),
            "E" => Some(Self::Glu(Default::default())),
            "S" => Some(Self::Ser(Default::default())),
            "T" => Some(Self::Thr(Default::default())),
            "N" => Some(Self::Asn(Default::default())),
            "Q" => Some(Self::Gln(Default::default())),
            "C" => Some(Self::Cys(Default::default())),
            "U" => Some(Self::Sec(Default::default())),
            "G" => Some(Self::Gly(Default::default())),
            "P" => Some(Self::Pro(Default::default())),
            "A" => Some(Self::Ala(Default::default())),
            "V" => Some(Self::Val(Default::default())),
            "I" => Some(Self::Ile(Default::default())),
            "L" => Some(Self::Leu(Default::default())),
            "M" => Some(Self::Met(Default::default())),
            "F" => Some(Self::Phe(Default::default())),
            "Y" => Some(Self::Tyr(Default::default())),
            "W" => Some(Self::Trp(Default::default())),
            _ => None,
        }
    }

    /// todo: Instead of this, many impl partial eq in a way that makes sense?
    pub fn aa_type(&self) -> AminoAcidType {
        match self {
            Self::Arg(_) => AminoAcidType::Arg,
            Self::His(_) => AminoAcidType::His,
            Self::Lys(_) => AminoAcidType::Lys,
            Self::Asp(_) => AminoAcidType::Asp,
            Self::Glu(_) => AminoAcidType::Glu,
            Self::Ser(_) => AminoAcidType::Ser,
            Self::Thr(_) => AminoAcidType::Thr,
            Self::Asn(_) => AminoAcidType::Asn,
            Self::Gln(_) => AminoAcidType::Gln,
            Self::Cys(_) => AminoAcidType::Cys,
            Self::Sec(_) => AminoAcidType::Sec,
            Self::Gly(_) => AminoAcidType::Gly,
            Self::Pro(_) => AminoAcidType::Pro,
            Self::Ala(_) => AminoAcidType::Ala,
            Self::Val(_) => AminoAcidType::Val,
            Self::Ile(_) => AminoAcidType::Ile,
            Self::Leu(_) => AminoAcidType::Leu,
            Self::Met(_) => AminoAcidType::Met,
            Self::Phe(_) => AminoAcidType::Phe,
            Self::Tyr(_) => AminoAcidType::Tyr,
            Self::Trp(_) => AminoAcidType::Trp,
        }
    }

    pub fn aa_name(&self) -> &str {
        match self {
            Self::Arg(_) => "Arg (R)",
            Self::His(_) => "His (H)",
            Self::Lys(_) => "Lys (K)",
            Self::Asp(_) => "Asp (D)",
            Self::Glu(_) => "Glu (E)",
            Self::Ser(_) => "Ser (S)",
            Self::Thr(_) => "Thr (T)",
            Self::Asn(_) => "Asn (N)",
            Self::Gln(_) => "Gln (Q)",
            Self::Cys(_) => "Cys (C)",
            Self::Sec(_) => "Sec (U)",
            Self::Gly(_) => "Gly (G)",
            Self::Pro(_) => "Pro (P)",
            Self::Ala(_) => "Ala (A)",
            Self::Val(_) => "Val (V)",
            Self::Ile(_) => "Ile (I)",
            Self::Leu(_) => "Leu (L)",
            Self::Met(_) => "Met (M)",
            Self::Phe(_) => "Phe (F)",
            Self::Tyr(_) => "Tyr (Y)",
            Self::Trp(_) => "Trp (W)",
        }
    }

    pub fn aa_ident_single_letter(&self) -> &str {
        match self {
            Self::Arg(_) => "R",
            Self::His(_) => "H",
            Self::Lys(_) => "K",
            Self::Asp(_) => "D",
            Self::Glu(_) => "E",
            Self::Ser(_) => "S",
            Self::Thr(_) => "T",
            Self::Asn(_) => "N",
            Self::Gln(_) => "Q",
            Self::Cys(_) => "C",
            Self::Sec(_) => "U",
            Self::Gly(_) => "G",
            Self::Pro(_) => "P",
            Self::Ala(_) => "A",
            Self::Val(_) => "V",
            Self::Ile(_) => "I",
            Self::Leu(_) => "L",
            Self::Met(_) => "M",
            Self::Phe(_) => "F",
            Self::Tyr(_) => "Y",
            Self::Trp(_) => "W",
        }
    }

    pub fn get_χ1(&self) -> Option<f64> {
        match self {
            Self::Arg(aa) => Some(aa.χ_1),
            Self::His(aa) => Some(aa.χ_1),
            Self::Lys(aa) => Some(aa.χ_1),
            Self::Asp(aa) => Some(aa.χ_1),
            Self::Glu(aa) => Some(aa.χ_1),
            Self::Ser(aa) => Some(aa.χ_1),
            Self::Thr(aa) => Some(aa.χ_1),
            Self::Asn(aa) => Some(aa.χ_1),
            Self::Gln(aa) => Some(aa.χ_1),
            Self::Cys(aa) => Some(aa.χ_1),
            Self::Sec(aa) => Some(aa.χ_1),
            Self::Val(aa) => Some(aa.χ_1),
            Self::Ile(aa) => Some(aa.χ_1),
            Self::Leu(aa) => Some(aa.χ_1),
            Self::Met(aa) => Some(aa.χ_1),
            Self::Phe(aa) => Some(aa.χ_1),
            Self::Tyr(aa) => Some(aa.χ_1),
            Self::Trp(aa) => Some(aa.χ_1),
            _ => None,
        }
    }

    pub fn get_χ2(&self) -> Option<f64> {
        match self {
            Self::Arg(aa) => Some(aa.χ_2),
            Self::His(aa) => Some(aa.χ_2),
            Self::Lys(aa) => Some(aa.χ_2),
            Self::Asp(aa) => Some(aa.χ_2),
            Self::Glu(aa) => Some(aa.χ_2),
            Self::Asn(aa) => Some(aa.χ_2),
            Self::Gln(aa) => Some(aa.χ_2),
            Self::Ile(aa) => Some(aa.χ_2),
            Self::Leu(aa) => Some(aa.χ_2),
            Self::Met(aa) => Some(aa.χ_2),
            Self::Phe(aa) => Some(aa.χ_2),
            Self::Tyr(aa) => Some(aa.χ_2),
            Self::Trp(aa) => Some(aa.χ_2),
            _ => None,
        }
    }
    pub fn get_χ3(&self) -> Option<f64> {
        match self {
            Self::Arg(aa) => Some(aa.χ_3),
            Self::Lys(aa) => Some(aa.χ_3),
            Self::Glu(aa) => Some(aa.χ_3),
            Self::Gln(aa) => Some(aa.χ_3),
            Self::Met(aa) => Some(aa.χ_3),
            _ => None,
        }
    }
    pub fn get_χ4(&self) -> Option<f64> {
        match self {
            Self::Arg(aa) => Some(aa.χ_4),
            Self::Lys(aa) => Some(aa.χ_4),
            _ => None,
        }
    }
    pub fn get_χ5(&self) -> Option<f64> {
        match self {
            Self::Arg(aa) => Some(aa.χ_5),
            _ => None,
        }
    }

    pub fn get_mut_χ1(&mut self) -> Option<&mut f64> {
        match self {
            Self::Arg(aa) => Some(&mut aa.χ_1),
            Self::His(aa) => Some(&mut aa.χ_1),
            Self::Lys(aa) => Some(&mut aa.χ_1),
            Self::Asp(aa) => Some(&mut aa.χ_1),
            Self::Glu(aa) => Some(&mut aa.χ_1),
            Self::Ser(aa) => Some(&mut aa.χ_1),
            Self::Thr(aa) => Some(&mut aa.χ_1),
            Self::Asn(aa) => Some(&mut aa.χ_1),
            Self::Gln(aa) => Some(&mut aa.χ_1),
            Self::Cys(aa) => Some(&mut aa.χ_1),
            Self::Sec(aa) => Some(&mut aa.χ_1),
            Self::Val(aa) => Some(&mut aa.χ_1),
            Self::Ile(aa) => Some(&mut aa.χ_1),
            Self::Leu(aa) => Some(&mut aa.χ_1),
            Self::Met(aa) => Some(&mut aa.χ_1),
            Self::Phe(aa) => Some(&mut aa.χ_1),
            Self::Tyr(aa) => Some(&mut aa.χ_1),
            Self::Trp(aa) => Some(&mut aa.χ_1),
            _ => None,
        }
    }

    pub fn get_mut_χ2(&mut self) -> Option<&mut f64> {
        match self {
            Self::Arg(aa) => Some(&mut aa.χ_2),
            Self::His(aa) => Some(&mut aa.χ_2),
            Self::Lys(aa) => Some(&mut aa.χ_2),
            Self::Asp(aa) => Some(&mut aa.χ_2),
            Self::Glu(aa) => Some(&mut aa.χ_2),
            Self::Asn(aa) => Some(&mut aa.χ_2),
            Self::Gln(aa) => Some(&mut aa.χ_2),
            Self::Ile(aa) => Some(&mut aa.χ_2),
            Self::Leu(aa) => Some(&mut aa.χ_2),
            Self::Met(aa) => Some(&mut aa.χ_2),
            Self::Phe(aa) => Some(&mut aa.χ_2),
            Self::Tyr(aa) => Some(&mut aa.χ_2),
            Self::Trp(aa) => Some(&mut aa.χ_2),
            _ => None,
        }
    }
    pub fn get_mut_χ3(&mut self) -> Option<&mut f64> {
        match self {
            Self::Arg(aa) => Some(&mut aa.χ_3),
            Self::Lys(aa) => Some(&mut aa.χ_3),
            Self::Glu(aa) => Some(&mut aa.χ_3),
            Self::Gln(aa) => Some(&mut aa.χ_3),
            Self::Met(aa) => Some(&mut aa.χ_3),
            _ => None,
        }
    }
    pub fn get_mut_χ4(&mut self) -> Option<&mut f64> {
        match self {
            Self::Arg(aa) => Some(&mut aa.χ_4),
            Self::Lys(aa) => Some(&mut aa.χ_4),
            _ => None,
        }
    }
    pub fn get_mut_χ5(&mut self) -> Option<&mut f64> {
        match self {
            Self::Arg(aa) => Some(&mut aa.χ_5),
            _ => None,
        }
    }

    pub fn add_to_χ1(&mut self, val: f64) {
        match self {
            Self::Arg(aa) => aa.χ_1 += val,
            Self::His(aa) => aa.χ_1 += val,
            Self::Lys(aa) => aa.χ_1 += val,
            Self::Asp(aa) => aa.χ_1 += val,
            Self::Glu(aa) => aa.χ_1 += val,
            Self::Ser(aa) => aa.χ_1 += val,
            Self::Thr(aa) => aa.χ_1 += val,
            Self::Asn(aa) => aa.χ_1 += val,
            Self::Gln(aa) => aa.χ_1 += val,
            Self::Cys(aa) => aa.χ_1 += val,
            Self::Sec(aa) => aa.χ_1 += val,
            Self::Val(aa) => aa.χ_1 += val,
            Self::Ile(aa) => aa.χ_1 += val,
            Self::Leu(aa) => aa.χ_1 += val,
            Self::Met(aa) => aa.χ_1 += val,
            Self::Phe(aa) => aa.χ_1 += val,
            Self::Tyr(aa) => aa.χ_1 += val,
            Self::Trp(aa) => aa.χ_1 += val,
            _ => (),
        }
    }

    pub fn add_to_χ2(&mut self, val: f64) {
        match self {
            Self::Arg(aa) => aa.χ_2 += val,
            Self::His(aa) => aa.χ_2 += val,
            Self::Lys(aa) => aa.χ_2 += val,
            Self::Asp(aa) => aa.χ_2 += val,
            Self::Glu(aa) => aa.χ_2 += val,
            Self::Asn(aa) => aa.χ_2 += val,
            Self::Gln(aa) => aa.χ_2 += val,
            Self::Ile(aa) => aa.χ_2 += val,
            Self::Leu(aa) => aa.χ_2 += val,
            Self::Met(aa) => aa.χ_2 += val,
            Self::Phe(aa) => aa.χ_2 += val,
            Self::Tyr(aa) => aa.χ_2 += val,
            Self::Trp(aa) => aa.χ_2 += val,
            _ => (),
        }
    }
    pub fn add_to_χ3(&mut self, val: f64) {
        match self {
            Self::Arg(aa) => aa.χ_3 += val,
            Self::Lys(aa) => aa.χ_3 += val,
            Self::Glu(aa) => aa.χ_3 += val,
            Self::Gln(aa) => aa.χ_3 += val,
            Self::Met(aa) => aa.χ_3 += val,
            _ => (),
        }
    }
    pub fn add_to_χ4(&mut self, val: f64) {
        match self {
            Self::Arg(aa) => aa.χ_4 += val,
            Self::Lys(aa) => aa.χ_4 += val,
            _ => (),
        }
    }
    pub fn add_to_χ5(&mut self, val: f64) {
        match self {
            Self::Arg(aa) => aa.χ_5 += val,
            _ => (),
        }
    }

    // pub fn set_χ1(&mut self, val: f64) {
    //     match self {
    //         Self::Arg(aa) => aa.χ_1 = val,
    //         Self::His(aa) => aa.χ_1 = val,
    //         Self::Lys(aa) => aa.χ_1 = val,
    //         Self::Asp(aa) => aa.χ_1 = val,
    //         Self::Glu(aa) => aa.χ_1 = val,
    //         Self::Ser(aa) => aa.χ_1 = val,
    //         Self::Thr(aa) => aa.χ_1 = val,
    //         Self::Asn(aa) => aa.χ_1 = val,
    //         Self::Gln(aa) => aa.χ_1 = val,
    //         Self::Cys(aa) => aa.χ_1 = val,
    //         Self::Sec(aa) => aa.χ_1 = val,
    //         Self::Val(aa) => aa.χ_1 = val,
    //         Self::Ile(aa) => aa.χ_1 = val,
    //         Self::Leu(aa) => aa.χ_1 = val,
    //         Self::Met(aa) => aa.χ_1 = val,
    //         Self::Phe(aa) => aa.χ_1 = val,
    //         Self::Tyr(aa) => aa.χ_1 = val,
    //         Self::Trp(aa) => aa.χ_1 = val,
    //         _ => (),
    //     }
    // }
    //
    // pub fn set_χ2(&mut self, val: f64) {
    //     match self {
    //         Self::Arg(aa) => aa.χ_2 = val,
    //         Self::His(aa) => aa.χ_2 = val,
    //         Self::Lys(aa) => aa.χ_2 = val,
    //         Self::Asp(aa) => aa.χ_2 = val,
    //         Self::Glu(aa) => aa.χ_2 = val,
    //         Self::Asn(aa) => aa.χ_2 = val,
    //         Self::Gln(aa) => aa.χ_2 = val,
    //         Self::Ile(aa) => aa.χ_2 = val,
    //         Self::Leu(aa) => aa.χ_2 = val,
    //         Self::Met(aa) => aa.χ_2 = val,
    //         Self::Phe(aa) => aa.χ_2 = val,
    //         Self::Tyr(aa) => aa.χ_2 = val,
    //         Self::Trp(aa) => aa.χ_2 = val,
    //         _ => (),
    //     }
    // }
    // pub fn set_χ3(&mut self, val: f64) {
    //     match self {
    //         Self::Arg(aa) => aa.χ_3 = val,
    //         Self::Lys(aa) => aa.χ_3 = val,
    //         Self::Glu(aa) => aa.χ_3 = val,
    //         Self::Gln(aa) => aa.χ_3 = val,
    //         Self::Met(aa) => aa.χ_3 = val,
    //         _ => (),
    //     }
    // }
    // pub fn set_χ4(&mut self, val: f64) {
    //     match self {
    //         Self::Arg(aa) => aa.χ_4 = val,
    //         Self::Lys(aa) => aa.χ_4 = val,
    //         _ => (),
    //     }
    // }
    // pub fn set_χ5(&mut self, val: f64) {
    //     match self {
    //         Self::Arg(aa) => aa.χ_5 = val,
    //         _ => (),
    //     }
    // }
}

impl Arg {
    // todo: Equiv from `backbone_cart_coords`.
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        // todo: Do we want our prev bond anchor to be n-calpha?
        n_pos: Vec3,
    ) -> CoordsArg {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.

        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            self.χ_3,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (n_eps, n_eps_orientation) = find_atom_placement(
            c_delta_orientation,
            unsafe { PLANAR3_A },
            unsafe { PLANAR3_B },
            self.χ_4,
            c_delta,
            c_gamma,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_zeta, c_zeta_orientation) = find_atom_placement(
            n_eps_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            self.χ_5,
            n_eps,
            c_delta,
            unsafe { PLANAR3_B },
            LEN_SC,
        );

        let (n_eta1, n_eta1_orientation) = find_atom_placement(
            c_zeta_orientation,
            unsafe { PLANAR3_A },
            unsafe { PLANAR3_B },
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (n_eta2, n_eta2_orientation) = find_atom_placement(
            c_zeta_orientation,
            unsafe { PLANAR3_A },
            unsafe { PLANAR3_B },
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (h_amine_eta1a, _) = find_atom_placement(
            c_zeta_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { PLANAR3_B },
            LEN_SC,
        );

        let (h_amine_eta1b, _) = find_atom_placement(
            c_zeta_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { PLANAR3_C },
            LEN_SC,
        );

        let (h_amine_eta2a, _) = find_atom_placement(
            c_zeta_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { PLANAR3_B },
            LEN_SC,
        );

        let (h_amine_eta2b, _) = find_atom_placement(
            c_zeta_orientation,
            H_BOND_IN,
            H_BOND_OUT,
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { PLANAR3_C },
            LEN_SC,
        );

        CoordsArg {
            c_beta,
            c_gamma,
            c_delta,
            n_eps,
            c_zeta,
            n_eta1,
            n_eta2,
            h_amine_eta1a,
            h_amine_eta1b,
            h_amine_eta2a,
            h_amine_eta2b,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
            n_eps_orientation,
            c_zeta_orientation,
            n_eta1_orientation,
            n_eta2_orientation,
        }
    }
}

impl His {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsHis {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (n_delta2, n_delta2_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { PLANAR3_A },
            unsafe { PLANAR3_B },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (n_eps1, n_eps1_orientation) = find_atom_placement(
            c_delta1_orientation,
            unsafe { PLANAR3_A },
            unsafe { PLANAR3_B },
            0.,
            c_delta1,
            c_gamma,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_eps2, _) = find_atom_placement(
            n_delta2_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            0.,
            n_delta2,
            c_gamma,
            unsafe { PLANAR3_B},
            LEN_SC,
        );

        let (h_amine_delta, _) = find_atom_placement(
            n_delta2_orientation,
            unsafe { H_BOND_IN },
            unsafe { H_BOND_OUT },
            TAU_DIV2,
            n_delta2,
            c_gamma,
            unsafe { PLANAR3_C },
            LEN_N_H,
        );
        let (h_amine_eps, _) = find_atom_placement(
            n_eps1_orientation,
            unsafe { H_BOND_IN },
            unsafe { H_BOND_OUT },
            TAU_DIV2,
            n_eps1,
            c_delta1,
            unsafe { PLANAR3_C },
            LEN_N_H,
        );

        // todo: These bond vecs are wrong! Needs to be tighter angles
        // todo due to there only being 5 atoms in the ring.

        CoordsHis {
            c_beta,
            c_gamma,
            c_delta1,
            n_delta2,
            n_eps1,
            c_eps2,
            h_amine_delta,
            h_amine_eps,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta1_orientation,
            n_delta2_orientation,
            n_eps1_orientation,
            // c_eps2_orientation,
        }
    }
}

// Calpha R bond is tetra C
// SC prev is planar C

impl Lys {
    // todo: Equiv from `backbone_cart_coords`.
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        // todo: Do we want our prev bond anchor to be n-calpha?
        n_pos: Vec3,
    ) -> CoordsLys {
        // These are the angles between each of 2 4 equally-spaced atoms on a tetrahedron,
        // with center of (0., 0., 0.). They are the angle formed between 3 atoms.
        // We have chosen the two angles to describe the backbone. We have chosen these arbitrarily.

        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            self.χ_3,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_eps, c_eps_orientation) = find_atom_placement(
            c_delta_orientation,
            unsafe { TETRA_A },
            unsafe { TETRA_B },
            self.χ_4,
            c_delta,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (n_zeta, n_zeta_orientation) = find_atom_placement(
            c_eps_orientation,
            unsafe {PLANAR3_A },
            unsafe { PLANAR3_B },
            TAU_DIV2,
            c_eps,
            c_delta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (h_amine1, _) = find_atom_placement(
            n_zeta_orientation,
            ANCHOR_BOND_VEC,
            H_BOND_OUT,
            TAU_DIV2,
            n_zeta,
            c_eps,
            unsafe { PLANAR3_B },
            LEN_N_H,
        );

        let (h_amine2, _) = find_atom_placement(
            n_zeta_orientation,
            ANCHOR_BOND_VEC,
            H_BOND_OUT,
            TAU_DIV2,
            n_zeta,
            c_eps,
            unsafe { PLANAR3_C },
            LEN_N_H,
        );

        CoordsLys {
            c_beta,
            c_gamma,
            c_delta,
            c_eps,
            n_zeta,
            h_amine1,
            h_amine2,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
            c_eps_orientation,
            n_zeta_orientation,
        }
    }
}

impl Asp {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsAsp {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        CoordsAsp {
            c_beta,
            c_gamma,
            o_delta1,
            o_delta2,

            c_beta_orientation,
            c_gamma_orientation,
        }
    }
}

impl Glu {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsGlu {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_3,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_eps1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_delta,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_eps2, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_delta,
            c_gamma,
            unsafe { TETRA_B },
            LEN_SC,
        );

        CoordsGlu {
            c_beta,
            c_gamma,
            c_delta,
            o_eps1,
            o_eps2,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
        }
    }
}

impl Ser {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsSer {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (o_gamma, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsSer {
            c_beta,
            o_gamma,

            c_beta_orientation,
        }
    }
}

impl Thr {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsThr {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma2, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (o_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsThr {
            c_beta,
            c_gamma2,
            o_gamma1,

            c_beta_orientation,
        }
    }
}

impl Asn {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsAsn {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (n_delta2, n_delta2_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (h_amine1, _) = find_atom_placement(
            n_delta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            n_delta2,
            c_gamma,
            unsafe { TETRA_A },
            LEN_N_H,
        );

        let (h_amine2, _) = find_atom_placement(
            n_delta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            n_delta2,
            c_gamma,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            LEN_N_H,
        );

        CoordsAsn {
            c_beta,
            c_gamma,
            o_delta1,
            n_delta2,
            h_amine1,
            h_amine2,

            c_beta_orientation,
            c_gamma_orientation,
            n_delta2_orientation,
        }
    }
}

impl Gln {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsGln {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_3,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_eps1, _) = find_atom_placement(
            c_delta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_delta,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (n_eps2, n_eps2_orientation) = find_atom_placement(
            c_delta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_delta,
            c_gamma,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (h_amine1, _) = find_atom_placement(
            n_eps2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            n_eps2,
            c_delta,
            unsafe { TETRA_A },
            LEN_N_H,
        );

        let (h_amine2, _) = find_atom_placement(
            n_eps2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            n_eps2,
            c_delta,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            LEN_N_H,
        );

        CoordsGln {
            c_beta,
            c_gamma,
            c_delta,
            o_eps1,
            n_eps2,
            h_amine1,
            h_amine2,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
            n_eps2_orientation,
        }
    }
}

impl Cys {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsCys {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (s_gamma, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsCys {
            c_beta,
            s_gamma,

            c_beta_orientation,
        }
    }
}

impl Sec {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsSec {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (se_gamma, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsSec {
            c_beta,
            se_gamma,

            c_beta_orientation,
        }
    }
}

impl Gly {
    pub fn sidechain_cart_coords(
        &self,
        _c_alpha: Vec3,
        _c_alpha_orientation: Quaternion,
        _n_pos: Vec3,
    ) -> CoordsGly {
        CoordsGly {}
    }
}

impl Pro {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsPro {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );
        CoordsPro {
            c_beta,
            c_gamma,
            c_delta,

            c_beta_orientation,
            c_gamma_orientation,
        }
    }
}

impl Ala {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsAla {
        let (c_beta, _c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        CoordsAla { c_beta }
    }
}

impl Val {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsVal {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_gamma2, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_B },
            LEN_SC,
        );

        CoordsVal {
            c_beta,
            c_gamma1,
            c_gamma2,

            c_beta_orientation,
        }
    }
}

impl Ile {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsIle {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_gamma2, c_gamma2_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_B }, // todo?
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_delta, _) = find_atom_placement(
            c_gamma2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma2,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsIle {
            c_beta,
            c_gamma1,
            c_gamma2,
            c_delta, // off gamma2

            c_beta_orientation,
            c_gamma2_orientation,
        }
    }
}

impl Leu {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsLeu {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { TETRA_D },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        CoordsLeu {
            c_beta,
            c_gamma,
            c_delta1,
            c_delta2,

            c_beta_orientation,
            c_gamma_orientation,
        }
    }
}

impl Met {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsMet {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (s_delta, s_delta_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_3,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_eps, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            s_delta,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsMet {
            c_beta,
            c_gamma,
            s_delta,
            c_eps,

            c_beta_orientation,
            c_gamma_orientation,
            s_delta_orientation,
        }
    }
}

impl Phe {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsPhe {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta2, c_delta2_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_eps1, c_eps1_orientation) = find_atom_placement(
            c_delta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_delta1,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_eps2, c_eps2_orientation) = find_atom_placement(
            c_delta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_delta2,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        // We anchor c_zeta off eps1.
        let (c_zeta, _c_zeta_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_eps1,
            c_delta1,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsPhe {
            c_beta,
            c_gamma,
            c_delta1,
            c_delta2,
            c_eps1,
            c_eps2,
            c_zeta,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta1_orientation,
            c_delta2_orientation,
            c_eps1_orientation,
            c_eps2_orientation,
        }
    }
}

impl Tyr {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsTyr {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_delta2, c_delta2_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_eps1, c_eps1_orientation) = find_atom_placement(
            c_delta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_delta1,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_eps2, c_eps2_orientation) = find_atom_placement(
            c_delta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_delta2,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        // We anchor c_zeta off eps1.
        let (c_zeta, c_zeta_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_eps1,
            c_delta1,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (o_eta, _) = find_atom_placement(
            c_zeta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_zeta,
            c_eps2,
            unsafe { TETRA_A },
            LEN_SC,
        );

        CoordsTyr {
            c_beta,
            c_gamma,
            c_delta1,
            c_delta2,
            c_eps1,
            c_eps2,
            c_zeta,
            o_eta,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta1_orientation,
            c_delta2_orientation,
            c_eps1_orientation,
            c_eps2_orientation,
            c_zeta_orientation,
        }
    }
}

impl Trp {
    pub fn sidechain_cart_coords(
        &self,
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
    ) -> CoordsTrp {
        let (c_beta, c_beta_orientation) = find_atom_placement(
            c_alpha_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { TETRA_A },
            LEN_SC,
        );

        // Trp's delta 1 and 2 are the ones not connecting to the second
        // ring.
        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (n_delta2, n_delta2_orientation) = find_atom_placement(
            c_delta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_delta1,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        // eps1 and eps2 are shared between the 2 rings.
        let (c_eps1, c_eps1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { TETRA_B },
            LEN_SC,
        );

        // Note: We've chosen to attach eps2 to eps1; another option
        // would be attaching it to delta2.
        let (c_eps2, c_eps2_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_eps1,
            c_gamma,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_zeta1, c_zeta1_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            // TAU_DIV2,
            self.χ_3, // todo: How do you sort out this hinge?
            c_eps1,
            c_gamma,
            unsafe { TETRA_B },
            LEN_SC,
        );

        let (c_zeta2, c_zeta2_orientation) = find_atom_placement(
            c_eps2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            self.χ_3, // todo: How do you sort out this hinge?
            // 0., // todo: How do you sort out this hinge?
            c_eps2,
            c_eps1,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_eta1, _) = find_atom_placement(
            c_zeta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            0.,
            c_zeta1,
            c_eps1,
            unsafe { TETRA_A },
            LEN_SC,
        );

        let (c_eta2, _) = find_atom_placement(
            c_zeta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            c_zeta2,
            c_eps2,
            unsafe { TETRA_A },
            LEN_SC,
        );
        let (h_amine, _) = find_atom_placement(
            n_delta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { TETRA_A },
            TAU_DIV2,
            n_delta2,
            c_gamma,
            unsafe { TETRA_A },
            LEN_N_H,
        );

        CoordsTrp {
            c_beta,
            c_gamma,
            c_delta1,
            n_delta2,
            c_eps1,
            c_eps2,
            c_zeta1,
            c_zeta2,
            c_eta1,
            c_eta2,
            h_amine,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta1_orientation,
            n_delta2_orientation,
            c_eps1_orientation,
            c_eps2_orientation,
            c_zeta1_orientation,
            c_zeta2_orientation,
        }
    }
}

/// These are global coordinates. Analog to `BackboneCoordsAa`
#[derive(Debug, Default)]
pub struct CoordsArg {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,
    pub n_eps: Vec3,
    pub c_zeta: Vec3,
    pub n_eta1: Vec3,
    pub n_eta2: Vec3,
    pub h_amine_eta1a: Vec3,
    pub h_amine_eta1b: Vec3,
    pub h_amine_eta2a: Vec3,
    pub h_amine_eta2b: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
    pub n_eps_orientation: Quaternion,
    pub c_zeta_orientation: Quaternion,
    pub n_eta1_orientation: Quaternion,
    pub n_eta2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsHis {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta1: Vec3,
    pub n_delta2: Vec3,
    pub n_eps1: Vec3,
    pub c_eps2: Vec3,
    pub h_amine_delta: Vec3,
    pub h_amine_eps: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta1_orientation: Quaternion,
    pub n_delta2_orientation: Quaternion,
    pub n_eps1_orientation: Quaternion,
    // pub n_eps1_orientation: Quaternion,
    // pub c_eps2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsLys {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,
    pub c_eps: Vec3,
    pub n_zeta: Vec3,
    pub h_amine1: Vec3,
    pub h_amine2: Vec3,
    // todo: Third N here?
    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
    pub c_eps_orientation: Quaternion,
    pub n_zeta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsAsp {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub o_delta1: Vec3,
    pub o_delta2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsGlu {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,
    pub o_eps1: Vec3,
    pub o_eps2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsSer {
    pub c_beta: Vec3,
    pub o_gamma: Vec3,

    pub c_beta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsThr {
    pub c_beta: Vec3,
    pub o_gamma1: Vec3,
    pub c_gamma2: Vec3,

    pub c_beta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsAsn {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub o_delta1: Vec3,
    pub n_delta2: Vec3,
    pub h_amine1: Vec3,
    pub h_amine2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub n_delta2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsGln {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,
    pub o_eps1: Vec3,
    pub n_eps2: Vec3,
    pub h_amine1: Vec3,
    pub h_amine2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
    pub n_eps2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsAla {
    pub c_beta: Vec3,
}

#[derive(Debug, Default)]
pub struct CoordsVal {
    pub c_beta: Vec3,
    pub c_gamma1: Vec3,
    pub c_gamma2: Vec3,

    c_beta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsIle {
    pub c_beta: Vec3,
    pub c_gamma1: Vec3,
    pub c_gamma2: Vec3,
    pub c_delta: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsLeu {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta1: Vec3,
    pub c_delta2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsCys {
    pub c_beta: Vec3,
    pub s_gamma: Vec3,

    pub c_beta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsSec {
    pub c_beta: Vec3,
    pub se_gamma: Vec3,

    pub c_beta_orientation: Quaternion,
}

pub struct CoordsMet {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub s_delta: Vec3,
    pub c_eps: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub s_delta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsGly {}

#[derive(Debug, Default)]
pub struct CoordsPro {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsPhe {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta1: Vec3,
    pub c_delta2: Vec3,
    pub c_eps1: Vec3,
    pub c_eps2: Vec3,
    pub c_zeta: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta1_orientation: Quaternion,
    pub c_delta2_orientation: Quaternion,
    pub c_eps1_orientation: Quaternion,
    pub c_eps2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsTyr {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta1: Vec3,
    pub c_delta2: Vec3,
    pub c_eps1: Vec3,
    pub c_eps2: Vec3,
    pub c_zeta: Vec3,
    pub o_eta: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta1_orientation: Quaternion,
    pub c_delta2_orientation: Quaternion,
    pub c_eps1_orientation: Quaternion,
    pub c_eps2_orientation: Quaternion,
    pub c_zeta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsTrp {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta1: Vec3,
    pub n_delta2: Vec3,
    /// eps 1 and 2 are shared between the two rings.
    pub c_eps1: Vec3,
    pub c_eps2: Vec3,
    /// zeta 1 is connected to eps1. zeta 2 is connected to eps2.
    pub c_zeta1: Vec3,
    pub c_zeta2: Vec3,
    pub c_eta1: Vec3,
    pub c_eta2: Vec3,
    pub h_amine: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta1_orientation: Quaternion,
    pub n_delta2_orientation: Quaternion,
    pub c_eps1_orientation: Quaternion,
    pub c_eps2_orientation: Quaternion,
    pub c_zeta1_orientation: Quaternion,
    pub c_zeta2_orientation: Quaternion,
}

// todo: Coord structs for the remaining AAs.

// the AA-specific structs below specify dihedral angles for each AA instance
// `χ_1` for each is for the bond between the c_alpha, and the first atom in the
// sidechain (eg c_bravo)

#[derive(Debug, PartialEq)]
pub struct Arg {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
    pub χ_4: f64,
    pub χ_5: f64,
}

impl Default for Arg {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
            χ_3: TAU_DIV2,
            χ_4: TAU_DIV2,
            χ_5: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct His {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for His {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Lys {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
    pub χ_4: f64,
}

impl Default for Lys {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
            χ_3: TAU_DIV2,
            χ_4: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Asp {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for Asp {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Glu {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

impl Default for Glu {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
            χ_3: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Ser {
    pub χ_1: f64,
}

impl Default for Ser {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug, PartialEq)]
pub struct Thr {
    pub χ_1: f64,
}

impl Default for Thr {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug, PartialEq)]
pub struct Asn {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for Asn {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Gln {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

impl Default for Gln {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
            χ_3: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Cys {
    pub χ_1: f64,
}

impl Default for Cys {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug, PartialEq)]
pub struct Sec {
    pub χ_1: f64,
}

impl Default for Sec {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug, Default, PartialEq)]
pub struct Gly {}

#[derive(Debug, Default, PartialEq)]
pub struct Pro {}

#[derive(Debug, Default, PartialEq)]
pub struct Ala {}

#[derive(Debug, PartialEq)]
pub struct Val {
    pub χ_1: f64,
}

impl Default for Val {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug, PartialEq)]
pub struct Ile {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for Ile {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Leu {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for Leu {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Met {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

impl Default for Met {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
            χ_3: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Phe {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for Phe {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Tyr {
    pub χ_1: f64,
    pub χ_2: f64,
}

impl Default for Tyr {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Trp {
    pub χ_1: f64,
    pub χ_2: f64,
    /// χ_3 is the rotation angle between the 2 rings.
    pub χ_3: f64,
}

impl Default for Trp {
    fn default() -> Self {
        Self {
            χ_1: TAU_DIV2,
            χ_2: TAU_DIV2,
            // todo: Make sure this is anchored correctly; it may
            // todo behave differently from single-atom angles.
            χ_3: TAU_DIV2,
        }
    }
}
