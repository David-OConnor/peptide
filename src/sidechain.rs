//! This module contains info related to side chains, including their geometry

use crate::{
    chem_definitions::{AminoAcidType, AtomType},
    // todo: Find actual bonds / angles; these are temp!
    kinematics::{find_atom_placement, CALPHA_R_BOND},
    lin_alg::{Quaternion, Vec3},
};

// /// An atom, in a sidechain
// #[derive(Debug)]
// pub struct AtomSidechain {
//     /// type of Atom, eg Carbon, Oxygen etc
//     pub role: AtomType,
//     /// Local position, anchored to Calpha = 0, 0, 0
//     /// todo: Orientation of side chain rel to Calpha?
//     pub position: Vec3,
// }

// Note: Carbon alpha here isn't the same as
// the backbone carbon alpha this is connected to!

// todo: These are temp
pub const LEN_SC: f64 = 1.53;
pub const BOND_IN: Vec3 = Vec3 {
    x: 1.,
    y: 0.,
    z: 0.,
};
pub const BOND_OUT1: Vec3 = Vec3 {
    x: 0.,
    y: 1.,
    z: 0.,
};
pub const BOND_OUT2: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 1.,
};

// todo: Consider moving AminoAcidType here, and making some of these methods of it?

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
    Sec(Sec),
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
            BOND_IN,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (n_eps, n_eps_orientation) = find_atom_placement(
            c_delta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_4,
            c_delta,
            c_gamma,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_zeta, c_zeta_orientation) = find_atom_placement(
            n_eps_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_5,
            n_eps,
            c_delta,
            BOND_OUT1,
            LEN_SC,
        );

        let (n_eta1, _) = find_atom_placement(
            c_zeta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_zeta,
            n_eps,
            BOND_OUT1,
            LEN_SC,
        );

        let (n_eta2, _) = find_atom_placement(
            c_zeta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_zeta,
            n_eps,
            BOND_OUT2,
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
            BOND_IN,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_eps, c_eps_orientation) = find_atom_placement(
            c_delta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_4,
            c_delta,
            c_gamma,
            BOND_OUT1,
            LEN_SC,
        );

        let (n_zeta, _) = find_atom_placement(
            c_eps_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_eps,
            c_delta,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            // Use our info about the previous 2 atoms so we can define the dihedral angle properly.
            // (world space)
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (o_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (o_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            BOND_OUT2,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (o_gamma, _) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_beta,
            c_alpha,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma2, _) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_beta,
            c_alpha,
            BOND_OUT2,
            LEN_SC,
        );

        let (o_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_beta,
            c_alpha,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (o_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (n_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            BOND_OUT2,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, c_delta_orientation) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (o_eps1, _) = find_atom_placement(
            c_delta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_delta,
            c_gamma,
            BOND_OUT1,
            LEN_SC,
        );

        let (n_eps2, _) = find_atom_placement(
            c_delta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_delta,
            c_gamma,
            BOND_OUT2,
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
        c_alpha: Vec3,
        c_alpha_orientation: Quaternion,
        n_pos: Vec3,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta, _) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_3,
            c_gamma,
            c_beta,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma1, _) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_gamma2, c_gamma2_orientation) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT2,
            LEN_SC,
        );

        let (c_delta, _) = find_atom_placement(
            c_gamma2_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma2,
            c_beta,
            BOND_OUT1,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta1, _) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta2, _) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_gamma,
            c_beta,
            BOND_OUT2,
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
            BOND_IN,
            BOND_OUT1,
            self.χ_1,
            c_alpha,
            n_pos,
            unsafe { CALPHA_R_BOND },
            LEN_SC,
        );

        let (c_gamma, c_gamma_orientation) = find_atom_placement(
            c_beta_orientation,
            BOND_IN,
            BOND_OUT1,
            self.χ_2,
            c_beta,
            c_alpha,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta1, c_delta1_orientation) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0., // todo
            c_gamma,
            c_beta,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_delta2, c_delta2_orientation) = find_atom_placement(
            c_gamma_orientation,
            BOND_IN,
            BOND_OUT1,
            0., // todo
            c_gamma,
            c_beta,
            BOND_OUT2,
            LEN_SC,
        );

        let (c_eps1, c_eps1_orientation) = find_atom_placement(
            c_delta1_orientation,
            BOND_IN,
            BOND_OUT1,
            0., // todo
            c_delta1,
            c_gamma,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_eps2, c_eps2_orientation) = find_atom_placement(
            c_delta2_orientation,
            BOND_IN,
            BOND_OUT1,
            0., // todo
            c_delta2,
            c_gamma,
            BOND_OUT1,
            LEN_SC,
        );

        let (c_zeta, c_zeta_orientation) = find_atom_placement(
            c_eps2_orientation,
            BOND_IN,
            BOND_OUT1,
            0., // todo
            c_eps2,
            c_delta2,
            BOND_OUT1,
            LEN_SC,
        );

        let (o_eta, _) = find_atom_placement(
            c_zeta_orientation,
            BOND_IN,
            BOND_OUT1,
            0.,
            c_zeta,
            c_eps2,
            BOND_OUT1,
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
pub struct CoordsGly {}

#[derive(Debug, Default)]
pub struct CoordsPro {
    pub c_beta: Vec3,
    pub c_gamma: Vec3,
    pub c_delta: Vec3,

    pub c_beta_orientation: Quaternion,
    pub c_gamma_orientation: Quaternion,
}

// todo: Coord structs for the remaining AAs.

// the AA-specific structs below specify dihedral angles for each AA instance
// `χ_1` for each is for the bond between the c_alpha, and the first atom in the
// sidechain (eg c_bravo)

#[derive(Debug, Default)]
pub struct Arg {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
    pub χ_4: f64,
    pub χ_5: f64,
}

#[derive(Debug, Default)]
pub struct His {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Lys {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
    pub χ_4: f64,
}

#[derive(Debug, Default)]
pub struct Asp {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Glu {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

#[derive(Debug, Default)]
pub struct Ser {
    pub χ_1: f64,
}

#[derive(Debug, Default)]
pub struct Thr {
    pub χ_1: f64,
}

#[derive(Debug, Default)]
pub struct Asn {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Gln {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

#[derive(Debug, Default)]
pub struct Cys {
    pub χ_1: f64,
}

#[derive(Debug, Default)]
pub struct Sec {
    pub χ_1: f64,
}

#[derive(Debug, Default)]
pub struct Gly {}

#[derive(Debug, Default)]
pub struct Pro {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

#[derive(Debug, Default)]
pub struct Ala {}

#[derive(Debug, Default)]
pub struct Val {
    pub χ_1: f64,
}

#[derive(Debug, Default)]
pub struct Ile {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Leu {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Met {
    pub χ_1: f64,
    pub χ_2: f64,
    pub χ_3: f64,
}

#[derive(Debug, Default)]
pub struct Phe {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Tyr {
    pub χ_1: f64,
    pub χ_2: f64,
}

#[derive(Debug, Default)]
pub struct Trp {
    pub χ_1: f64,
    pub χ_2: f64,
}
