//! This module contains info related to side chains, including their geometry

use crate::{
    chem_definitions::{AminoAcidType, AtomType},
    lin_alg::Vec3,
};

/// An atom, in a sidechain
#[derive(Debug)]
pub struct AtomSidechain {
    /// type of Atom, eg Carbon, Oxygen etc
    pub role: AtomType,
    /// Local position, anchored to Calpha = 0, 0, 0
    /// todo: Orientation of side chain rel to Calpha?
    pub position: Vec3,
}

// the AA-specific structs below specify dihedral angles for each AA instance
// `θ_0` for each is for the bond between the c_alpha, and the first C in the sidechain.

pub struct Arg {
    pub θ_0: f64,
    pub θ_1: f64,
    pub θ_2: f64,
    pub θ_3: f64,
}

pub struct His {
    pub θ_0: f64,
    pub θ_1: f64,
}

pub struct Lys {
    pub θ_0: f64,
    pub θ_1: f64,
    pub θ_2: f64,
    pub θ_3: f64,
}

pub struct Asp {
    pub θ_0: f64,
    pub θ_1: f64,
}

pub struct Glu {
    pub θ_0: f64,
    pub θ_1: f64,
    pub θ_2: f64,
}

pub struct Ser {
    pub θ_0: f64,
}

pub struct Thr {
    pub θ_0: f64,
}

pub struct Asn {
    pub θ_0: f64,
    pub θ_1: f64,
}

pub struct Gln {
    pub θ_0: f64,
    pub θ_1: f64,
    pub θ_2: f64,
}

pub struct Cys {
    pub θ_0: f64,
}

pub struct Sec {
    pub θ_0: f64,
}

pub struct Gly {}

pub struct Pro {}

pub struct Ala {}

pub struct Val {
    pub θ_0: f64,
}

pub struct Ile {
    pub θ_0: f64,
    pub θ_1: f64,
}

pub struct Leu {
    pub θ_0: f64,
    pub θ_1: f64,
}
pub struct Met {
    pub θ_0: f64,
    pub θ_1: f64,
    pub θ_2: f64,
}
pub struct Phe {
    pub θ_0: f64,
    pub θ_1: f64,
}

pub struct Tyr {
    pub θ_0: f64,
    pub θ_1: f64,
}

pub struct Trp {
    pub θ_0: f64,
    pub θ_1: f64,
}

/// Describes the sidechain of a given amino acid
#[derive(Debug)]
pub struct Sidechain {
    atoms: Vec<AtomSidechain>,
}

impl Sidechain {
    pub fn from_aa(aa: AminoAcidType) -> Self {
        let atoms = match aa {
            AminoAcidType::A => vec![AtomSidechain {
                role: AtomType::C,
                position: Vec3::new(0., 0., 0.),
            }],
            AminoAcidType::R => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::N => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::D => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::C => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::S,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::Q => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::E => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::G => vec![], // G has no sidechain. (single H)
            AminoAcidType::H => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::I => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::L => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::K => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::M => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::S,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::F => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::P => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::S => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::T => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::W => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::N,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::Y => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::O,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
            AminoAcidType::V => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
        };

        Self { atoms }
    }
}
