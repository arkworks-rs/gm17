use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{PrimeField, UniformRand};
use ark_poly::GeneralEvaluationDomain;
use ark_std::{cfg_into_iter, vec::Vec};

use crate::{r1cs_to_sap::R1CStoSAP, Proof, ProvingKey};

use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

/// Create a zero-knowledge GrothMaller17 proof.
#[inline]
pub fn create_random_proof<E, C, R>(
    circuit: C,
    params: &ProvingKey<E>,
    rng: &mut R,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    R: Rng,
{
    let d1 = E::Fr::rand(rng);
    let d2 = E::Fr::rand(rng);
    let r = E::Fr::rand(rng);

    create_proof::<E, C>(circuit, params, d1, d2, r)
}

/// Create a GrothMaller17 proof with randomness `d1`, `d2`, `r`.
pub fn create_proof<E, C>(
    circuit: C,
    params: &ProvingKey<E>,
    d1: E::Fr,
    d2: E::Fr,
    r: E::Fr,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    type D<F> = GeneralEvaluationDomain<F>;

    let prover_time = start_timer!(|| "GrothMaller17::Prover");
    let cs = ConstraintSystem::new_ref();

    // Synthesize the circuit.
    let synthesis_time = start_timer!(|| "Constraint synthesis");
    circuit.generate_constraints(cs.clone())?;
    end_timer!(synthesis_time);

    let lc_time = start_timer!(|| "Inlining LCs");
    cs.inline_all_lcs();
    end_timer!(lc_time);

    let witness_map_time = start_timer!(|| "R1CS to SAP witness map");
    let (full_input_assignment, h, _) =
        R1CStoSAP::witness_map::<E::Fr, D<E::Fr>>(cs.clone(), &d1, &d2)?;
    end_timer!(witness_map_time);
    let prover = cs.borrow().unwrap();
    let num_inputs = prover.num_instance_variables;

    let input_assignment = full_input_assignment[1..num_inputs]
        .iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();

    let aux_assignment = cfg_into_iter!(full_input_assignment[num_inputs..])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    drop(full_input_assignment);

    let h_input = h[0..num_inputs]
        .iter()
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    let h_aux = cfg_into_iter!(h[num_inputs..])
        .map(|s| s.into_repr())
        .collect::<Vec<_>>();
    drop(h);

    // Compute A
    let a_acc_time = start_timer!(|| "Compute A");
    let (a_inputs, a_aux) = params.a_query.split_at(num_inputs);
    let a_inputs_acc = VariableBaseMSM::multi_scalar_mul(&a_inputs[1..], &input_assignment);
    let a_aux_acc = VariableBaseMSM::multi_scalar_mul(a_aux, &aux_assignment);

    let r_g = params.g_gamma_z.mul(r);
    let d1_g = params.g_gamma_z.mul(d1);

    let mut g_a = r_g;
    g_a += &params.a_query[0].into_projective();
    g_a += &d1_g;
    g_a += &a_inputs_acc;
    g_a += &a_aux_acc;
    end_timer!(a_acc_time);

    // Compute B
    let b_acc_time = start_timer!(|| "Compute B");

    let (b_inputs, b_aux) = params.b_query.split_at(num_inputs);
    let b_inputs_acc = VariableBaseMSM::multi_scalar_mul(&b_inputs[1..], &input_assignment);
    let b_aux_acc = VariableBaseMSM::multi_scalar_mul(b_aux, &aux_assignment);

    let r_h = params.h_gamma_z.mul(r);
    let d1_h = params.h_gamma_z.mul(d1);

    let mut g_b = r_h;
    g_b += &params.b_query[0].into_projective();
    g_b += &d1_h;
    g_b += &b_inputs_acc;
    g_b += &b_aux_acc;
    end_timer!(b_acc_time);

    // Compute C
    let c_acc_time = start_timer!(|| "Compute C");
    let r_2 = r + &r;
    let r2 = r * &r;
    let d1_r_2 = d1 * &r_2;

    let c1_acc_time = start_timer!(|| "Compute C1");
    let c1_acc = VariableBaseMSM::multi_scalar_mul(&params.c_query_1, &aux_assignment);
    end_timer!(c1_acc_time);

    let c2_acc_time = start_timer!(|| "Compute C2");

    let (c2_inputs, c2_aux) = params.c_query_2.split_at(num_inputs);
    let c2_inputs_acc = VariableBaseMSM::multi_scalar_mul(&c2_inputs[1..], &input_assignment);
    let c2_aux_acc = VariableBaseMSM::multi_scalar_mul(c2_aux, &aux_assignment);

    let c2_acc = c2_inputs_acc + &c2_aux_acc;
    end_timer!(c2_acc_time);

    // Compute G
    let g_acc_time = start_timer!(|| "Compute G");

    let (g_inputs, g_aux) = params.g_gamma2_z_t.split_at(num_inputs);
    let g_inputs_acc = VariableBaseMSM::multi_scalar_mul(g_inputs, &h_input);
    let g_aux_acc = VariableBaseMSM::multi_scalar_mul(g_aux, &h_aux);

    let g_acc = g_inputs_acc + &g_aux_acc;
    end_timer!(g_acc_time);

    let r2_g_gamma2_z2 = params.g_gamma2_z2.mul(r2);
    let r_g_ab_gamma_z = params.g_ab_gamma_z.mul(r);
    let d1_g_ab_gamma_z = params.g_ab_gamma_z.mul(d1);
    let r_c0 = params.c_query_2[0].mul(r);
    let r2_d1_g_gamma2_z2 = params.g_gamma2_z2.mul(d1_r_2);
    let d2_g_gamma2_z_t0 = params.g_gamma2_z_t[0].mul(d2);
    let mut r_c2_exp = c2_acc;
    r_c2_exp *= r;

    let mut g_c = c1_acc;
    g_c += &r2_g_gamma2_z2;
    g_c += &r_g_ab_gamma_z;
    g_c += &d1_g_ab_gamma_z;
    g_c += &r_c0;
    g_c += &r2_d1_g_gamma2_z2;
    g_c += &r_c2_exp;
    g_c += &d2_g_gamma2_z_t0;
    g_c += &g_acc;
    end_timer!(c_acc_time);

    end_timer!(prover_time);

    Ok(Proof {
        a: g_a.into_affine(),
        b: g_b.into_affine(),
        c: g_c.into_affine(),
    })
}
