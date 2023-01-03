//! This module contains code for running a simulation over time. Related to the concept of
//! integration.

use core::f64::consts::TAU;

use crate::{
    atom_coords::ProteinCoords,
    bond_vecs::{LEN_O_H, WATER_BOND_H_A, WATER_BOND_H_B, WATER_BOND_M},
    forces::{self, CHARGE_ELECTRON, CHARGE_PROTON, K_C},
    quantum::Electron,
    water::{self, H_MASS, O_H_DIST},
    AminoAcidType, State,
};

use lin_alg2::f64::{Quaternion, Vec3};

// Distance from the origin.
pub const SIM_BOX_DIST: f64 = 40.;

// todo: Util?
/// Get a random value from -0.5 to 0.5
fn rng() -> f64 {
    rand::random::<f64>() - 0.5
}

/// Execute our simulation.
pub fn run(state: &mut State, dt: f32) {
    // Constant term in our noise.
    let dt_modified = state.ui.temperature * state.ui.sim_time_scale * dt as f64;
    // Crude approach, where we add random noise to the angles.

    for res in &mut state.protein.descrip.residues {
        res.ψ += rng() * dt_modified;
        res.φ += rng() * dt_modified;

        crate::clamp_angle(&mut res.φ, res.sidechain.aa_type() == AminoAcidType::Pro);
        crate::clamp_angle(&mut res.ψ, false);

        if let Some(χ) = res.sidechain.get_mut_χ1() {
            *χ += rng() * dt_modified;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ2() {
            *χ += rng() * dt_modified;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ3() {
            *χ += rng() * dt_modified;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ4() {
            *χ += rng() * dt_modified;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ5() {
            *χ += rng() * dt_modified;
            crate::clamp_angle(χ, false);
        }
    }

    state.protein.coords = ProteinCoords::from_descrip(&state.protein.descrip);

    let mut dt_modified = state.ui.sim_time_scale * dt as f64;
    dt_modified *= 1000.; // todo: Temp fudge factor

    // todo: Hack to prevent editing the loop we're itereating through.
    let wm_dup = state.water_env.water_molecules.clone();

    // todo: is thsi right?
    let m_water = water::O_MASS + 2. * water::H_MASS;
    for (i, water) in state.water_env.water_molecules.iter_mut().enumerate() {
        // let v_o = forces::potential(water.o_posit_world, forces::O_CHARGE, &wm_dup);
        // let v_h_a = forces::potential(water.h_a_posit_world, forces::H_CHARGE, &wm_dup);
        // let v_h_b = forces::potential(water.h_b_posit_world, forces::H_CHARGE, &wm_dup);
        // let v_m = forces::potential(water.m_posit_world, forces::M_CHARGE, &wm_dup);

        // let f_o = v_o * forces::CH

        // let potential = v_o + v_h_a + v_h_b + v_m;

        // let f_o = forces::force(water.o_posit, forces::O_CHARGE, &wm_dup, i);
        // let f_h_a = forces::force(water.ha_posit, forces::H_CHARGE, &wm_dup, i);
        // let f_h_b = forces::force(water.hb_posit, forces::H_CHARGE, &wm_dup, i);
        // let f_m = forces::force(water.m_posit, forces::M_CHARGE, &wm_dup, i);

        // let (force, torque) = forces::force_tipt4(water, &wm_dup, i);
        let (force, torque) = forces::water_tipt4(water, &wm_dup, i);

        // todo: Dist of just O?
        // let force = (water.o_posit_world)
        // let a = (f_o + f_h_a + f_h_b + f_m) / m_water;

        // Calculate inertial moments for each hydrogen atom, on the given torque axis.
        // Note that these are already normalized.
        let o_ha_bond_world = water.o_orientation.rotate_vec(WATER_BOND_H_A);
        let o_hb_bond_world = water.o_orientation.rotate_vec(unsafe { WATER_BOND_H_B });

        // Calculate an effective radius; the radius that's going in the direction of rotation.
        // Take the dot product between each hydrogen, and the perpendicular vector of the
        // normalized torque closest to it to find the cosine loss.

        let torque_norm = torque.to_normalized();

        // todo: Can we use this shortcut?
        let r_ha_effective = (1. - (o_ha_bond_world.dot(torque_norm))) * O_H_DIST;
        let r_hb_effective = (1. - (o_hb_bond_world.dot(torque_norm))) * O_H_DIST;

        // // todo: QC this
        // let r_ha_effective = bond_o_ha_world_norm.dot(torque_norm) * LEN_O_H;
        // let r_hb_effective = bond_o_hb_world_norm.dot(torque_norm) * LEN_O_H;

        let l = H_MASS * (r_ha_effective.powi(2) + r_hb_effective.powi(2));

        let a = force / m_water;
        let α = torque / l;

        // todo: Do we want the time sim var in here too?
        water.velocity += a * dt_modified as f64; // todo: Euler integration - not great
        water.angular_velocity += α * dt_modified as f64; // todo: Euler integration - not great

        // todo: Beter way to incorporate temp. Isn't it KE, which is vel sq? so maybe take sqrt of temp?
        // let fudge_factor = 10.;
        // let fudge_factor_rot = 0.1;

        water.o_posit += water.velocity
            // * state.ui.temperature
            * dt_modified as f64;

        let r = Quaternion::from_axis_angle(
            water.angular_velocity.to_normalized(),
            water.angular_velocity.magnitude() * dt_modified,
        );

        water.o_orientation = r * water.o_orientation;

        // todo: Code to re-generate out-of-bond molecules?
    }

    let nuc_dup = state.wavefunction_lab.nuclei.clone();

    for (i, prot_this) in state.wavefunction_lab.nuclei.iter_mut().enumerate() {
        let (force_prot, force_elec) = forces::hydrogen_atoms(
            prot_this,
            &Electron {
                // position: Vec3::new_zero(),
                // velocity: Vec3::new_zero(),
                spin: Default::default(),
            }, // unused
            &nuc_dup,
            &state.wavefunction_lab.electron_posits_dynamic,
            i,
        );

        let a = force_prot / forces::MASS_PROT;
        // todo: dt or dt_modified here, as above?
        prot_this.velocity += a * dt_modified as f64 * 0.001; // todo: Euler integration - not great
    }
    //
    // for (i, elec) in state.wavefunction_lab.electron_posits_dynamic.iter_mut().enumerate() {
    //
    // }

    state.water_env.update_atom_posits();
    state.wavefunction_lab.update_posits(dt_modified);

    // todo: Since we don't have a way to update teh electron centers properly, do this fudge.
    // for (i, elec) in state
    //     .wavefunction_lab
    //     .electrons
    //     .iter_mut()
    //     .enumerate()
    // {
    //     elec.position = state.wavefunction_lab.nuclei[i].position;
    // }
}
