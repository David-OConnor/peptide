//! Code for saving to and loading from files. We use a simple binary format.
//!
//! In the future, we may change to a human-readable format.

use std::{
    f32::consts::TAU,
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    str,
};

use crate::{
    kinematics::{ProteinDescription, Residue},
    sidechain::{self, Sidechain},
};

use lin_alg2::f64::Vec3;

// todo: A bit of a sloppy way
const NAME_START_I: usize = 0;
const IDENT_START_I: usize = 64;
const RESIDUE_START_I: usize = 96;

const F64_SIZE: usize = 8;

// Sidechain identifier (0 - 19), then up to 5 dihedral angles.
const SIDECHAIN_SIZE: usize = 1 + F64_SIZE * 5;
const RESIDUE_SIZE: usize = SIDECHAIN_SIZE + 3 * F64_SIZE;

impl Sidechain {
    fn to_bytes(&self) -> [u8; SIDECHAIN_SIZE] {
        [0; SIDECHAIN_SIZE]
    }

    // Convert to a single-byte identifier. Note that we can't use `repr(u8)` due to the enclosed
    // structs.
    fn to_ident(&self) -> u8 {
        match self {
            Self::Arg(_) => 0,
            Self::His(_) => 1,
            Self::Lys(_) => 2,
            Self::Asp(_) => 3,
            Self::Glu(_) => 4,
            Self::Ser(_) => 5,
            Self::Thr(_) => 6,
            Self::Asn(_) => 7,
            Self::Gln(_) => 8,
            Self::Cys(_) => 9,
            Self::Gly(_) => 10,
            Self::Pro(_) => 11,
            Self::Ala(_) => 12,
            Self::Val(_) => 13,
            Self::Ile(_) => 14,
            Self::Leu(_) => 15,
            Self::Met(_) => 16,
            Self::Phe(_) => 17,
            Self::Tyr(_) => 18,
            Self::Trp(_) => 19,
            Self::Sec(_) => 20,
        }
    }
}

impl ProteinDescription {
    /// Save the protein description (metadata, dihedral angles etc) to our custom binary file format.
    pub fn save(&self, filename: &str) {
        let mut file_buf = Vec::new();

        file_buf.resize(RESIDUE_START_I + self.residues.len() * RESIDUE_SIZE, 0);

        let name_buf = self.name.as_bytes();
        let ident_buf = self.pdb_ident.as_bytes();

        // extend_from_slice()

        file_buf[NAME_START_I..NAME_START_I + name_buf.len()].clone_from_slice(&name_buf);
        file_buf[IDENT_START_I..IDENT_START_I + ident_buf.len()].clone_from_slice(&ident_buf);

        let mut i = RESIDUE_START_I;

        for res in &self.residues {
            file_buf[i] = res.sidechain.to_ident();
            i += 1;
            file_buf[i..i + F64_SIZE].clone_from_slice(&res.ω.to_le_bytes());
            i += F64_SIZE;
            file_buf[i..i + F64_SIZE].clone_from_slice(&res.φ.to_le_bytes());
            i += F64_SIZE;
            file_buf[i..i + F64_SIZE].clone_from_slice(&res.ψ.to_le_bytes());
            i += F64_SIZE;
            // file_buf.append(&mut res.sidechain.to_bytes()); // todo
        }

        let f = match File::open(filename) {
            Ok(f_) => f_,
            // todo: This won't handle permission errors etc properly
            Err(_) => File::create(filename).unwrap(),
        };

        let mut writer = BufWriter::new(f);

        writer.write_all(&file_buf).unwrap();
    }

    pub fn load(filename: &str) -> Self {
        let f = File::open(filename).unwrap();
        let mut reader = BufReader::new(f);
        let mut file_buf = Vec::new();

        reader.read_to_end(&mut file_buf).unwrap();

        let name = str::from_utf8(&file_buf[NAME_START_I..IDENT_START_I])
            .unwrap()
            .trim_end_matches(char::from(0))
            .to_owned();
        let pdb_ident = str::from_utf8(&file_buf[IDENT_START_I..RESIDUE_START_I])
            .unwrap()
            .trim_end_matches(char::from(0))
            .to_owned();

        let mut residues = Vec::new();

        let mut i = RESIDUE_START_I;
        while i < file_buf.len() {
            let aa_ident = file_buf[i];
            i += 1;

            let ω = f64::from_le_bytes(file_buf[i..i + F64_SIZE].try_into().unwrap());
            i += F64_SIZE;
            let φ = f64::from_le_bytes(file_buf[i..i + F64_SIZE].try_into().unwrap());
            i += F64_SIZE;
            let ψ = f64::from_le_bytes(file_buf[i..i + F64_SIZE].try_into().unwrap());
            i += F64_SIZE;

            residues.push(Residue {
                ω,
                φ,
                ψ,
                sidechain: Sidechain::His(Default::default()), // todo temp
                dipole: Vec3::new_zero(),
            });
        }

        Self {
            name,
            pdb_ident,
            residues,
        }
    }
}
