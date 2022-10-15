//! This module contains info related to side chains, including their geometry

// Don't show warnings for unused AAs etc.
#![allow(dead_code)]

use std::{f64::consts::TAU, fmt};

use crate::{
    bond_vecs::*,
    // todo: Find actual bonds / angles; these are temp!
    kinematics::find_atom_placement,
};

use lin_alg2::f64::{Quaternion, Vec3};

const TAU_DIV2: f64 = TAU / 2.;

// todo: These are temp
pub const LEN_SC: f64 = 1.53;

// pub const BOND_IN: Vec3 = Vec3 {
//     x: 1.,
//     y: 0.,
//     z: 0.,
// };
// pub const BOND_OUT1: Vec3 = Vec3 {
//     x: 0.,
//     y: 1.,
//     z: 0.,
// };
// pub const BOND_OUT2: Vec3 = Vec3 {
//     x: 0.,
//     y: 0.,
//     z: 1.,
// };

// todo: Gauche+ and trans etc for beta C. EG opposite C' or opposite N?

#[derive(Debug)]
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
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_eps, n_eps_orientation) = find_atom_placement(
            c_delta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_4,
            c_delta,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_zeta, c_zeta_orientation) = find_atom_placement(
            n_eps_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_5,
            n_eps,
            c_delta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_eta1, _) = find_atom_placement(
            c_zeta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_zeta,
            n_eps,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_eta2, _) = find_atom_placement(
            c_zeta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_zeta,
            n_eps,
            unsafe { SIDECHAIN_BOND_OUT2 },
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

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
            n_eps_orientation,
            c_zeta_orientation,
        }
    }
}

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
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_eps, c_eps_orientation) = find_atom_placement(
            c_delta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_4,
            c_delta,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_zeta, _) = find_atom_placement(
            c_eps_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_eps,
            c_delta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        CoordsLys {
            c_beta,
            c_gamma,
            c_delta,
            c_eps,
            n_zeta,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
            c_eps_orientation,
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
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (o_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (o_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { SIDECHAIN_BOND_OUT2 },
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (o_gamma, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma2, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_beta,
            c_alpha,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        let (o_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (o_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        CoordsAsn {
            c_beta,
            c_gamma,
            o_delta1,
            n_delta2,

            c_beta_orientation,
            c_gamma_orientation,
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (o_eps1, _) = find_atom_placement(
            c_delta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_delta,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_eps2, _) = find_atom_placement(
            c_delta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_delta,
            c_gamma,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        CoordsGln {
            c_beta,
            c_gamma,
            c_delta,
            o_eps1,
            n_eps2,

            c_beta_orientation,
            c_gamma_orientation,
            c_delta_orientation,
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
            SIDECHAIN_BOND_OUT1,
            0.,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            0.,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_gamma2, c_gamma2_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            unsafe { SIDECHAIN_BOND_OUT2 }, // todo?
            self.χ_2,
            c_beta,
            c_alpha,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        let (c_delta, _) = find_atom_placement(
            c_gamma2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma2,
            c_beta,
            SIDECHAIN_BOND_OUT1,
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
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { SIDECHAIN_BOND_OUT2 },
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_delta2, c_delta2_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        let (c_eps1, c_eps1_orientation) = find_atom_placement(
            c_delta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            0.,
            c_delta1,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_eps2, c_eps2_orientation) = find_atom_placement(
            c_delta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            0.,
            c_delta2,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        // We anchor c_zeta off eps1.
        let (c_zeta, c_zeta_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_eps1,
            c_delta1,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (o_eta, _) = find_atom_placement(
            c_zeta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_zeta,
            c_eps2,
            SIDECHAIN_BOND_OUT1,
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
            SIDECHAIN_BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        // Trp's delta 1 and 2 are the ones not connecting to the second
        // ring.
        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (n_delta2, n_delta2_orientation) = find_atom_placement(
            c_delta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_delta1,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        // eps1 and eps2 are shared between the 2 rings.
        let (c_eps1, c_eps1_orientation) = find_atom_placement(
            c_gamma_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_gamma,
            c_beta,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        // Note: We've chosen to attach eps2 to eps1; another option
        // would be attaching it to delta2.
        let (c_eps2, c_eps2_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_eps1,
            c_gamma,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_zeta1, c_zeta1_orientation) = find_atom_placement(
            c_eps1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            // TAU_DIV2,
            self.χ_3, // todo: How do you sort out this hinge?
            c_eps1,
            c_gamma,
            unsafe { SIDECHAIN_BOND_OUT2 },
            LEN_SC,
        );

        let (c_zeta2, c_zeta2_orientation) = find_atom_placement(
            c_eps2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            self.χ_3, // todo: How do you sort out this hinge?
            // 0., // todo: How do you sort out this hinge?
            c_eps2,
            c_eps1,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_eta1, _) = find_atom_placement(
            c_zeta1_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            0.,
            c_zeta1,
            c_eps1,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
        );

        let (c_eta2, _) = find_atom_placement(
            c_zeta2_orientation,
            unsafe { SIDECHAIN_BOND_TO_PREV },
            SIDECHAIN_BOND_OUT1,
            TAU_DIV2,
            c_zeta2,
            c_eps2,
            SIDECHAIN_BOND_OUT1,
            LEN_SC,
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

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
    pub n_eps_orientation: Quaternion,
    pub c_zeta_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsHis {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta1: Vec3,
    pub n_delta2: Vec3,
    pub n_eps1: Vec3,
    pub c_eps2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta1_orientation: Quaternion,
    pub n_delta2_orientation: Quaternion,
    pub n_eps1_orientation: Quaternion,
    pub c_eps2_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsLys {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,
    pub c_eps: Vec3,
    pub n_zeta: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
    pub c_eps_orientation: Quaternion,
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

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
}

#[derive(Debug, Default)]
pub struct CoordsGln {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,
    pub o_eps1: Vec3,
    pub n_eps2: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
    pub c_delta_orientation: Quaternion,
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
pub struct Ser {
    pub χ_1: f64,
}

impl Default for Ser {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug)]
pub struct Thr {
    pub χ_1: f64,
}

impl Default for Thr {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
pub struct Cys {
    pub χ_1: f64,
}

impl Default for Cys {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug)]
pub struct Sec {
    pub χ_1: f64,
}

impl Default for Sec {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug, Default)]
pub struct Gly {}

#[derive(Debug, Default)]
pub struct Pro {}

#[derive(Debug, Default)]
pub struct Ala {}

#[derive(Debug)]
pub struct Val {
    pub χ_1: f64,
}

impl Default for Val {
    fn default() -> Self {
        Self { χ_1: TAU_DIV2 }
    }
}

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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

#[derive(Debug)]
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
