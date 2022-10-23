//! This file contains render code that's engine agnostic. Configuration, etc.

use crate::chem_definitions::AtomRole;

use lin_alg2::f64::{Quaternion, Vec3};

// pub const BACKGROUND_COLOR: (f32, f32, f32) = (0.9, 0.9, 0.9);
pub const BACKGROUND_COLOR: (f32, f32, f32) = (0., 0., 0.);

pub const BOND_COLOR_BACKBONE: (f32, f32, f32) = (0.2, 0.2, 0.2);
pub const BOND_COLOR_SIDECHAIN: (f32, f32, f32) = (0.7, 0.65, 0.5);

// Atom colors
pub const CALPHA_COLOR: (f32, f32, f32) = (0.2, 0.8, 0.2);
pub const CP_COLOR: (f32, f32, f32) = (0.4, 0.6, 0.6);
pub const C_SIDECHAIN_COLOR: (f32, f32, f32) = (0.3, 0.5, 0.3);
pub const S_SIDECHAIN_COLOR: (f32, f32, f32) = (0.6, 0.6, 0.);
pub const SE_SIDECHAIN_COLOR: (f32, f32, f32) = (0.0, 0.6, 0.6);

pub const N_COLOR: (f32, f32, f32) = (0., 0., 1.);
pub const O_COLOR: (f32, f32, f32) = (1., 0., 0.);
pub const H_COLOR: (f32, f32, f32) = (0.8, 0.8, 0.8);

// Fraction of C radius.
pub const H_SCALE: f32 = 0.8;

// Note: This active color is deliberately not normalized, so it comes out weaker in saturation,
// but brighter in the avg with the atom color.
pub const ACTIVE_COLOR_ATOM: (f32, f32, f32) = (2., 1., 0.3);

// Shinyness affects specular intensity.
pub const ATOM_SHINYNESS: f32 = 1.;
pub const BOND_SHINYNESS: f32 = 1.;

pub const BOND_RADIUS_BACKBONE: f32 = 0.09;
pub const BOND_RADIUS_SIDECHAIN: f32 = 0.025;
pub const BOND_N_SIDES: usize = 10;

// These sensitivities are in units (position), or radians (orientation) per second.
pub const CAM_MOVE_SENS: f64 = 1.1;
pub const CAM_ROTATE_SENS: f64 = 0.2;
pub const CAM_ROTATE_KEY_SENS: f64 = 0.5;
// Move speed multiplier when the run modifier key is held.
pub const RUN_FACTOR: f64 = 5.;

pub const WINDOW_TITLE: &str = "Peptide";
pub const WINDOW_SIZE_X: f32 = 1_200.0;
pub const WINDOW_SIZE_Y: f32 = 1_000.0;

// Changes the far end of the frustrum; try to have this shortly past the farthest
// view distance you expect. Possible to make this dynamic?
pub const RENDER_DIST: f32 = 200.;

// Render size of an atom, on a [polyhedron] side.
pub const SIDE_LEN: f32 = 0.4;

// // todo: Do we want this in addition to the render lib's cam?
// pub struct Camera {
//     pub position: Vec3,
//     pub orientation: Quaternion,
// }

impl AtomRole {
    pub fn render_color(&self) -> (f32, f32, f32) {
        let cα = CALPHA_COLOR;
        let cp = CP_COLOR;
        let n = N_COLOR;
        let o = O_COLOR;
        let h_n = H_COLOR;
        let h_cα = H_COLOR;
        let cs = C_SIDECHAIN_COLOR;
        let s = S_SIDECHAIN_COLOR;
        let se = SE_SIDECHAIN_COLOR;

        match self {
            Self::Cα => (cα.0, cα.1, cα.2),
            Self::Cp => (cp.0, cp.1, cp.2),
            Self::N => (n.0, n.1, n.2),
            Self::O => (o.0, o.1, o.2),
            Self::HCα => (h_cα.0, h_cα.1, h_cα.2),
            Self::HN => (h_n.0, h_n.1, h_n.2),
            Self::CSidechain => (cs.0, cs.1, cs.2),
            // todo: Consider a diff shade for n and o sidechain colors
            Self::NSidechain => (n.0, n.1, n.2),
            Self::OSidechain => (o.0, o.1, o.2),
            Self::HSidechain => (h_n.0, h_n.1, h_n.2),
            Self::SSidechain => (s.0, s.1, s.2),
            Self::SeSidechain => (se.0, se.1, se.2),
        }
    }
}
