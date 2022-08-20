//! This module contains info related to side chains, including their geometry

use crate::{
    chem_definitions::{AminoAcidType, AtomType},
    coord_gen::{AtomSidechain},
    lin_alg::Vec3,
};

/// Describes the sidechain of a given amino acid
#[derive(Debug)]
pub struct Sidechain {
    atoms: Vec<AtomSidechain>,
}

impl Sidechain {
    pub fn from_aa(aa: AminoAcidType) -> Self {
        let atoms = match aa {
            AminoAcidType::A => vec![
                AtomSidechain {
                    role: AtomType::C,
                    position: Vec3::new(0., 0., 0.),
                },
            ],
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