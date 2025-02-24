//! Names and enumerations of various chemistry definitions.
#![allow(unused)]

use std::fmt;

#[derive(Clone, Copy, Debug)]
pub enum AtomType {
    C,
    N,
    H,
    O,
    P,
    S,
    Se,
}

impl AtomType {
    // todo: What is the significance of this? It's a bit nebulous
    // pub fn _charge(&self) -> f64 {
    //     match self {
    //         Self::C => 4.,
    //         Self::N => -3.,
    //         Self::H => 1.,
    //         Self::O => -2.,
    //         Self::P => 0., // (5., 3., -3.)
    //         Self::S => 0., // (-2., 2., 4., 6.)
    //     }
    // }

    /// Returns mass, in atomic mass units.
    pub fn mass(&self) -> f64 {
        match self {
            Self::C => 12.011,
            Self::N => 14.007,
            Self::H => 1.0078,
            Self::O => 15.999,
            Self::P => 30.974,
            Self::S => 32.065,
            Self::Se => 78.960,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
/// The positional role of an atom in an amino acid
/// todo: We've overloaded this to include sidechains too! Consider
/// todo how you'd prefer to handle this.
pub enum AtomRole {
    N,
    /// The carbon attom that attaches to the side chain.
    Cα,
    Cp,
    /// Oxygen double-bonded to C'
    O,
    HN,
    HCα,
    CSidechain,
    OSidechain,
    NSidechain,
    SSidechain,
    SeSidechain,
    HSidechain,
}

impl AtomRole {
    pub fn atom_type(&self) -> AtomType {
        match self {
            Self::N => AtomType::N,
            _ => AtomType::C,
        }
    }
}
