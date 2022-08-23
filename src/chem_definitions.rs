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

impl AminoAcidType {
    pub fn symbol(&self) -> &str {
        match self {
            Self::Ala => "Ala",
            Self::Arg => "Arg",
            Self::Asn => "Asn",
            Self::Asp => "Asp",
            Self::Cys => "Cys",
            Self::Gln => "Gln",
            Self::Glu => "Glu",
            Self::Gly => "Gly",
            Self::His => "His",
            Self::Ile => "Ile",
            Self::Leu => "Leu",
            Self::Lys => "Lys",
            Self::Met => "Met",
            Self::Phe => "Phe",
            Self::Pro => "Pro",
            Self::Ser => "Ser",
            Self::Thr => "Thr",
            Self::Trp => "Trp",
            Self::Tyr => "Tyr",
            Self::Val => "Val",
            Self::Sec => "Sec",
        }
    }
}

#[derive(Clone, Copy, Debug)]
/// The positional role of an atom in an amino acid
/// todo: We've overloaded this to include sidechains too! Consider
/// todo how you'd prefer to handle this.
pub enum BackboneRole {
    N,
    /// The carbon attom that attaches to the side chain.
    Cα,
    Cp,
    /// Oxygen double-bonded to C'
    O,
    CSidechain,
    OSidechain,
    NSidechain,
}

impl BackboneRole {
    pub fn atom_type(&self) -> AtomType {
        match self {
            Self::N => AtomType::N,
            _ => AtomType::C,
        }
    }
}
