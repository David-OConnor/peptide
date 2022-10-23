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
}

#[derive(Clone, Copy, Debug, PartialEq)]
/// Our order follows this Wikipedia example:
/// https://upload.wikimedia.org/wikipedia/commons/4/4f/ProteinogenicAminoAcids.svg
pub enum AminoAcidType {
    Arg,
    His,
    Lys,
    Asp,
    Glu,
    Ser,
    Thr,
    Asn,
    Gln,
    Cys,
    Sec,
    Gly,
    Pro,
    Ala,
    Val,
    Ile,
    Leu,
    Met,
    Phe,
    Tyr,
    Trp,
}

impl fmt::Display for AminoAcidType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let v = match self {
            Self::Arg => "Arg (R)",
            Self::His => "His (H)",
            Self::Lys => "Lys (K)",
            Self::Asp => "Asp (D)",
            Self::Glu => "Glu (E)",
            Self::Ser => "Ser (S)",
            Self::Thr => "Thr (T)",
            Self::Asn => "Asn (N)",
            Self::Gln => "Gln (Q)",
            Self::Cys => "Cys (C)",
            Self::Sec => "Sec (U)",
            Self::Gly => "Gly (G)",
            Self::Pro => "Pro (P)",
            Self::Ala => "Ala (A)",
            Self::Val => "Val (V)",
            Self::Ile => "Ile (I)",
            Self::Leu => "Leu (L)",
            Self::Met => "Met (M)",
            Self::Phe => "Phe (F)",
            Self::Tyr => "Tyr (Y)",
            Self::Trp => "Trp (W)",
        };

        write!(f, "{}", v)
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
