//! Code for saving to and loading from files. We use a simple binary format.
//!
//! In the future, we may change to a human-readable format.

use std::{
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    str,
};

use lin_alg2::f64::Vec3;

use crate::{
    sidechain::*,
    types::{ProteinDescription, Residue},
};

// todo: A bit of a sloppy way
const NAME_START_I: usize = 0;
const IDENT_START_I: usize = 64;
const RESIDUE_START_I: usize = 96;

const F64_SIZE: usize = 8;

// Sidechain amino acid identifier, then space for 5 dihedral angles.
const SIDECHAIN_SIZE: usize = 1 + F64_SIZE * 5;
// Sidechain, + the 3 backbone dihedral angles.
const RESIDUE_SIZE: usize = SIDECHAIN_SIZE + 3 * F64_SIZE;

impl Sidechain {
    // Convert to a single-byte identifier. Note that we can't use `repr(u8)` due to the enclosed
    // AA-specific struct in each enum variant.
    fn to_byte_ident(&self) -> u8 {
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

    fn to_bytes(&self) -> [u8; SIDECHAIN_SIZE] {
        let mut result = [0; SIDECHAIN_SIZE];

        result[0] = self.to_byte_ident();

        if let Some(χ) = self.get_χ1() {
            result[1..9].clone_from_slice(&χ.to_le_bytes());
        }
        if let Some(χ) = self.get_χ2() {
            result[9..17].clone_from_slice(&χ.to_le_bytes());
        }
        if let Some(χ) = self.get_χ3() {
            result[17..25].clone_from_slice(&χ.to_le_bytes());
        }
        if let Some(χ) = self.get_χ4() {
            result[25..33].clone_from_slice(&χ.to_le_bytes());
        }
        if let Some(χ) = self.get_χ5() {
            result[33..41].clone_from_slice(&χ.to_le_bytes());
        }

        result
    }

    fn from_bytes(data: &[u8]) -> Self {
        match data[0] {
            // todo: DRY from above ident matching
            0 => {
                Self::Arg(Arg {
                    // todo: Simplifying method like we have for to, so we're not doing this 20 times?
                    χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                    χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
                    χ_3: f64::from_le_bytes(data[17..25].try_into().unwrap()),
                    χ_4: f64::from_le_bytes(data[25..33].try_into().unwrap()),
                    χ_5: f64::from_le_bytes(data[33..41].try_into().unwrap()),
                })
            }
            1 => Self::His(His {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            2 => Self::Lys(Lys {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
                χ_3: f64::from_le_bytes(data[17..25].try_into().unwrap()),
                χ_4: f64::from_le_bytes(data[25..33].try_into().unwrap()),
            }),
            3 => Self::Asp(Asp {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            4 => Self::Glu(Glu {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
                χ_3: f64::from_le_bytes(data[17..25].try_into().unwrap()),
            }),
            5 => Self::Ser(Ser {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
            }),
            6 => Self::Thr(Thr {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
            }),
            7 => Self::Asn(Asn {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            8 => Self::Gln(Gln {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
                χ_3: f64::from_le_bytes(data[17..25].try_into().unwrap()),
            }),
            9 => Self::Cys(Cys {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
            }),
            10 => Self::Gly(Gly {}),
            11 => Self::Pro(Pro {}),
            12 => Self::Ala(Ala {}),
            13 => Self::Val(Val {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
            }),
            14 => Self::Ile(Ile {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            15 => Self::Leu(Leu {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            16 => Self::Met(Met {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
                χ_3: f64::from_le_bytes(data[17..25].try_into().unwrap()),
            }),
            17 => Self::Phe(Phe {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            18 => Self::Tyr(Tyr {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
            }),
            19 => Self::Trp(Trp {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
                χ_2: f64::from_le_bytes(data[9..17].try_into().unwrap()),
                χ_3: f64::from_le_bytes(data[17..25].try_into().unwrap()),
            }),
            20 => Self::Sec(Sec {
                χ_1: f64::from_le_bytes(data[1..9].try_into().unwrap()),
            }),
            // todo temp; fill in the others.
            _ => Self::Gly(Gly {}),
        }
    }
}

impl Residue {
    fn to_bytes(&self) -> [u8; RESIDUE_SIZE] {
        let mut result = [0; RESIDUE_SIZE];

        result[0..F64_SIZE].clone_from_slice(&self.ω.to_le_bytes());
        result[F64_SIZE..2 * F64_SIZE].clone_from_slice(&self.φ.to_le_bytes());
        result[2 * F64_SIZE..3 * F64_SIZE].clone_from_slice(&self.ψ.to_le_bytes());

        result[3 * F64_SIZE..RESIDUE_SIZE].clone_from_slice(&self.sidechain.to_bytes());

        result
    }

    fn from_bytes(data: &[u8]) -> Self {
        Self {
            ω: f64::from_le_bytes(data[0..F64_SIZE].try_into().unwrap()),
            φ: f64::from_le_bytes(data[F64_SIZE..2 * F64_SIZE].try_into().unwrap()),
            ψ: f64::from_le_bytes(data[2 * F64_SIZE..3 * F64_SIZE].try_into().unwrap()),

            sidechain: Sidechain::from_bytes(&data[3 * F64_SIZE..RESIDUE_SIZE]),
            // dipole: Vec3::new_zero(), // not currently saved.
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
            file_buf[i..i + RESIDUE_SIZE].clone_from_slice(&res.to_bytes());

            i += RESIDUE_SIZE;
        }

        // Overwrite teh existing file. Commented-out code below doesn't do anything if
        // the file exists.

        let f = File::create(filename).unwrap();

        // let f = match File::open(filename) {
        //     Ok(f_) => {
        //         println!("YES");
        //         f_
        //     },
        //     // todo: This won't handle permission errors, already exists etc. etc properly
        //     Err(_) => File::create(filename).unwrap(),
        // };

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
            residues.push(Residue::from_bytes(&file_buf[i..i + RESIDUE_SIZE]));

            i += RESIDUE_SIZE
        }

        Self {
            name,
            pdb_ident,
            residues,
        }
    }
}
