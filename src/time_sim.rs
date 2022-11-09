//! This module contains code for running a simulation over time. Related to the concept of
//! integration.

use crate::{atom_coords::ProteinCoords, forces, AminoAcidType, State};

use lin_alg2::f64::{Quaternion, Vec3};

// todo: Util?
/// Get a random value from -0.5 to 0.5
fn rng() -> f64 {
    rand::random::<f64>() - 0.5
}

/// Execute our simulation.
pub fn run(state: &mut State, dt: f32) {
    // Constant term in our noise.
    let c = state.temperature * state.sim_time_scale * dt as f64;
    // Crude approach, where we add random noise to the angles.

    for res in &mut state.protein_descrip.residues {
        res.ψ += rng() * c;
        res.φ += rng() * c;

        crate::clamp_angle(&mut res.φ, res.sidechain.aa_type() == AminoAcidType::Pro);
        crate::clamp_angle(&mut res.ψ, false);

        if let Some(χ) = res.sidechain.get_mut_χ1() {
            *χ += rng() * c;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ2() {
            *χ += rng() * c;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ3() {
            *χ += rng() * c;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ4() {
            *χ += rng() * c;
            crate::clamp_angle(χ, false);
        }
        if let Some(χ) = res.sidechain.get_mut_χ5() {
            *χ += rng() * c;
            crate::clamp_angle(χ, false);
        }
    }

    state.protein_coords = ProteinCoords::from_descrip(&state.protein_descrip);

    // todo: Hack to prevent editing the loop we're itereating through.
    let wm_dup = state.water_env.water_molecules.clone();

    // todo: is thsi right?
    let m_water = forces::O_MASS + 2. * forces::H_MASS;
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

        let (force, torque) = forces::force(water, &wm_dup, i);

        // todo: Dist of just O?
        // let force = (water.o_posit_world)
        // let a = (f_o + f_h_a + f_h_b + f_m) / m_water;

        let a = force / m_water;
        let ω = torque / m_water;

        water.velocity += a; // todo: Euler integration - not great
        water.angular_velocity += ω; // todo: Euler integration - not great

        // todo: Beter way to incorporate temp. Isn't it KE, which is vel sq? so maybe take sqrt of temp?
        let fudge_factor = 0.01;
        let fudge_factor_rot = 0.001;

        water.o_posit +=
            water.velocity * state.temperature * state.sim_time_scale * dt as f64 * fudge_factor;

        let r = Quaternion::from_axis_angle(
            water.angular_velocity.to_normalized(),
            water.angular_velocity.magnitude() * fudge_factor_rot,
        );

        // water.o_orientation = r * water.o_orientation;

        // todo: Code to re-generate out-of-bond molecules?
    }

    state.water_env.update_atom_posits();
}
