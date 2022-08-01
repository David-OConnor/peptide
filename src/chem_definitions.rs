//! Names and enumerations of various chemistry definitions.

#![allow(unused)]

#[derive(Clone, Copy, Debug)]
pub enum AtomType {
    C,
    N,
    H,
    O,
    P,
    S,
}

impl AtomType {
    // todo: What is the significance of this? It's a bit nebulous
    pub fn _charge(&self) -> f64 {
        match self {
            Self::C => 4.,
            Self::N => -3.,
            Self::H => 1.,
            Self::O => -2.,
            Self::P => 0., // (5., 3., -3.)
            Self::S => 0., // (-2., 2., 4., 6.)
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum AminoAcidType {
    A,
    R,
    N,
    D,
    C,
    Q,
    E,
    G,
    H,
    I,
    L,
    K,
    M,
    F,
    P,
    S,
    T,
    W,
    Y,
    V,
}

impl AminoAcidType {
    pub fn symbol(&self) -> &str {
        match self {
            Self::A => "Ala",
            Self::R => "Arg",
            Self::N => "Asn",
            Self::D => "Asp",
            Self::C => "Cys",
            Self::Q => "Gln",
            Self::E => "Glu",
            Self::G => "Gly",
            Self::H => "His",
            Self::I => "Ile",
            Self::L => "Leu",
            Self::K => "Lys",
            Self::M => "Met",
            Self::F => "Phe",
            Self::P => "Pro",
            Self::S => "Ser",
            Self::T => "Thr",
            Self::W => "Trp",
            Self::Y => "Tyr",
            Self::V => "Val",
        }
    }
}

#[derive(Clone, Copy, Debug)]
/// The positional role of an atom in an amino acid
pub enum BackboneRole {
    N,
    Cα,
    Cp,
}