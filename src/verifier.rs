use ark_ec::{pairing::{Pairing, MillerLoopOutput}, CurveGroup, AffineRepr};
use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};
use core::ops::{Mul, AddAssign};
use num_traits::Zero;
use std::ops::Neg;

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

/// Prepare the verifying key `vk` for use in proof verification.
pub fn prepare_verifying_key<E: Pairing>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    PreparedVerifyingKey {
        vk: vk.clone(),
        g_alpha: vk.g_alpha_g1,
        h_beta: vk.h_beta_g2,
        g_alpha_h_beta_ml: E::miller_loop(vk.g_alpha_g1, vk.h_beta_g2).0,
        g_gamma_pc: vk.g_gamma_g1.into(),
        h_gamma_pc: vk.h_gamma_g2.into(),
        h_pc: vk.h_g2.into(),
        query: vk.query.clone(),
    }
}

/// Verify a GrothMaller17 proof `proof` against the prepared verification key `pvk`,
/// with respect to the instance `public_inputs`.
pub fn verify_proof<E: Pairing>(
    pvk: &PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    public_inputs: &[E::ScalarField],
) -> R1CSResult<bool> {
    if (public_inputs.len() + 1) != pvk.query.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    // e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma}) *
    // e(C, H) where psi = \sum_{i=0}^l input_i pvk.query[i]

    let mut g_psi = pvk.query[0].into_group();
    for (i, b) in public_inputs.iter().zip(pvk.query.iter().skip(1)) {
        g_psi += &b.mul(*i);
    }

    let mut test1_a_g_alpha = proof.a.into_group();
    test1_a_g_alpha.add_assign(&pvk.g_alpha);
    let test1_a_g_alpha = test1_a_g_alpha.into_affine();

    let mut test1_b_h_beta = proof.b.into_group();
    test1_b_h_beta.add_assign(&pvk.h_beta);
    let test1_b_h_beta = test1_b_h_beta.into_affine();

    let test1_r1 = pvk.g_alpha_h_beta_ml;

    let test1_r2_input_g1 = [
        -test1_a_g_alpha.into_group(),
        g_psi.into().into(),
        proof.c.into()
    ];

    let test1_r2_input_g2= [
        test1_b_h_beta.into(),
        pvk.h_gamma_pc.clone(),
        pvk.h_pc.clone()
    ];

    let test1_r2 = E::multi_miller_loop(test1_r2_input_g1.to_owned(), test1_r2_input_g2.to_owned()).0;
    let mlo = MillerLoopOutput(test1_r2 * test1_r1);
    let test1 = E::final_exponentiation(mlo).unwrap();

    // e(A, H^{gamma}) = e(G^{gamma}, B)
    let test2_exp_input_g1= [
        proof.a.into(),
        pvk.g_gamma_pc.clone(),
    ];
    let test2_exp_input_g2= [
        pvk.h_gamma_pc.clone(),
        (-proof.b.into_group()).into().into(),
    ];

    let test2_exp = E::multi_miller_loop(test2_exp_input_g1.to_owned(), test2_exp_input_g2.to_owned()).0;
    let mlo = MillerLoopOutput(test1_r2 * test1_r1);
    let test2 = E::final_exponentiation(mlo).unwrap();

    Ok(test1.is_zero() && test2.is_zero())
}
