use crate::{PreparedVerifyingKey, Proof, VerifyingKey, GM17};
use ark_crypto_primitives::snark::constraints::{CircuitSpecificSetupSNARKGadget, SNARKGadget};
use ark_crypto_primitives::snark::{BooleanInputVar, SNARK};
use ark_ec::{AffineCurve, PairingEngine};
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::CurveVar;
use ark_r1cs_std::{
    alloc::{AllocVar, AllocationMode},
    bits::boolean::Boolean,
    bits::uint8::UInt8,
    eq::EqGadget,
    pairing::PairingVar,
    ToBitsGadget, ToBytesGadget,
};
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::{borrow::Borrow, marker::PhantomData, vec::Vec};

/// The proof variable for the GM17 construction
#[derive(Derivative)]
#[derivative(Clone(bound = "P::G1Var: Clone, P::G2Var: Clone"))]
pub struct ProofVar<E: PairingEngine, P: PairingVar<E>> {
    /// The `A` element in `G1`.
    pub a: P::G1Var,
    /// The `B` element in `G2`.
    pub b: P::G2Var,
    /// The `C` element in `G1`.
    pub c: P::G1Var,
}

/// The verifying key variable for the GM17 construction
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P::G1Var: Clone, P::GTVar: Clone, P::G1PreparedVar: Clone, \
    P::G2PreparedVar: Clone, ")
)]
pub struct VerifyingKeyVar<E: PairingEngine, P: PairingVar<E>> {
    #[doc(hidden)]
    pub h_g2: P::G2Var,
    #[doc(hidden)]
    pub g_alpha_g1: P::G1Var,
    #[doc(hidden)]
    pub h_beta_g2: P::G2Var,
    #[doc(hidden)]
    pub g_gamma_g1: P::G1Var,
    #[doc(hidden)]
    pub h_gamma_g2: P::G2Var,
    #[doc(hidden)]
    pub query: Vec<P::G1Var>,
}

impl<E: PairingEngine, P: PairingVar<E>> VerifyingKeyVar<E, P> {
    /// Prepare the verifying key `vk` for use in proof verification.
    pub fn prepare(&self) -> Result<PreparedVerifyingKeyVar<E, P>, SynthesisError> {
        let g_alpha_pc = P::prepare_g1(&self.g_alpha_g1)?;
        let h_beta_pc = P::prepare_g2(&self.h_beta_g2)?;
        let g_gamma_pc = P::prepare_g1(&self.g_gamma_g1)?;
        let h_gamma_pc = P::prepare_g2(&self.h_gamma_g2)?;
        let h_pc = P::prepare_g2(&self.h_g2)?;
        Ok(PreparedVerifyingKeyVar {
            g_alpha: self.g_alpha_g1.clone(),
            h_beta: self.h_beta_g2.clone(),
            g_alpha_pc,
            h_beta_pc,
            g_gamma_pc,
            h_gamma_pc,
            h_pc,
            query: self.query.clone(),
        })
    }
}

/// Preprocessed verification key parameters variable for the GM17 construction
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P::G1Var: Clone, P::GTVar: Clone, P::G1PreparedVar: Clone, \
    P::G2PreparedVar: Clone, ")
)]
pub struct PreparedVerifyingKeyVar<E: PairingEngine, P: PairingVar<E>> {
    #[doc(hidden)]
    pub g_alpha: P::G1Var,
    #[doc(hidden)]
    pub h_beta: P::G2Var,
    #[doc(hidden)]
    pub g_alpha_pc: P::G1PreparedVar,
    #[doc(hidden)]
    pub h_beta_pc: P::G2PreparedVar,
    #[doc(hidden)]
    pub g_gamma_pc: P::G1PreparedVar,
    #[doc(hidden)]
    pub h_gamma_pc: P::G2PreparedVar,
    #[doc(hidden)]
    pub h_pc: P::G2PreparedVar,
    #[doc(hidden)]
    pub query: Vec<P::G1Var>,
}

/// The SNARK verifier gadget of [[GM17]](https://eprint.iacr.org/2016/260.pdf).
pub struct GM17VerifierGadget<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    _pairing_engine: PhantomData<E>,
    _pairing_gadget: PhantomData<P>,
}

impl<E: PairingEngine, P: PairingVar<E, E::Fq>> SNARKGadget<E::Fr, E::Fq, GM17<E>>
    for GM17VerifierGadget<E, P>
{
    type ProcessedVerifyingKeyVar = PreparedVerifyingKeyVar<E, P>;
    type VerifyingKeyVar = VerifyingKeyVar<E, P>;
    type InputVar = BooleanInputVar<E::Fr, E::Fq>;
    type ProofVar = ProofVar<E, P>;

    type VerifierSize = usize;

    fn verifier_size(circuit_vk: &<GM17<E> as SNARK<E::Fr>>::VerifyingKey) -> Self::VerifierSize {
        circuit_vk.query.len()
    }

    /// Allocates `N::Proof` in `cs` without performing
    /// subgroup checks.
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_proof_unchecked<T: Borrow<Proof<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self::ProofVar, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        f().and_then(|proof| {
            let proof = proof.borrow();
            let a = CurveVar::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "Proof.a"),
                || Ok(proof.a.into_projective()),
                mode,
            )?;
            let b = CurveVar::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "Proof.b"),
                || Ok(proof.b.into_projective()),
                mode,
            )?;
            let c = CurveVar::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "Proof.c"),
                || Ok(proof.c.into_projective()),
                mode,
            )?;
            Ok(ProofVar { a, b, c })
        })
    }

    /// Allocates `N::Proof` in `cs` without performing
    /// subgroup checks.
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_verification_key_unchecked<T: Borrow<VerifyingKey<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self::VerifyingKeyVar, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        f().and_then(|vk| {
            let vk = vk.borrow();
            let g_alpha_g1 = P::G1Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "g_alpha"),
                || Ok(vk.g_alpha_g1.into_projective()),
                mode,
            )?;
            let h_g2 = P::G2Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "h"),
                || Ok(vk.h_g2.into_projective()),
                mode,
            )?;
            let h_beta_g2 = P::G2Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "h_beta"),
                || Ok(vk.h_beta_g2.into_projective()),
                mode,
            )?;
            let g_gamma_g1 = P::G1Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "g_gamma"),
                || Ok(vk.g_gamma_g1.into_projective()),
                mode,
            )?;
            let h_gamma_g2 = P::G2Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "h_gamma"),
                || Ok(vk.h_gamma_g2.into_projective()),
                mode,
            )?;

            let query = vk
                .query
                .iter()
                .map(|g| {
                    P::G1Var::new_variable_omit_prime_order_check(
                        ark_relations::ns!(cs, "g"),
                        || Ok(g.into_projective()),
                        mode,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(VerifyingKeyVar {
                g_alpha_g1,
                h_g2,
                h_beta_g2,
                g_gamma_g1,
                h_gamma_g2,
                query,
            })
        })
    }

    #[tracing::instrument(target = "r1cs", skip(circuit_pvk, x, proof))]
    fn verify_with_processed_vk(
        circuit_pvk: &Self::ProcessedVerifyingKeyVar,
        x: &Self::InputVar,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<E::Fq>, SynthesisError> {
        let circuit_pvk = circuit_pvk.clone();
        // e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma}) *
        // e(C, H) where psi = \sum_{i=0}^l input_i pvk.query[i]
        let g_psi = {
            let mut g_psi = circuit_pvk.query[0].clone();
            let mut input_len = 1;
            let mut input = x.clone().into_iter();
            for (input, b) in input.by_ref().zip(circuit_pvk.query.iter().skip(1)) {
                let input_bits = input.to_bits_le()?;
                g_psi += b.scalar_mul_le(input_bits.iter())?;
                input_len += 1;
            }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == circuit_pvk.query.len() && input.next().is_none());
            g_psi
        };

        let mut test1_a_g_alpha = proof.a.clone();
        test1_a_g_alpha += circuit_pvk.g_alpha.clone();
        let mut test1_b_h_beta = proof.b.clone();
        test1_b_h_beta += circuit_pvk.h_beta.clone();

        let test1_exp = {
            test1_a_g_alpha = test1_a_g_alpha.negate()?;
            let test1_a_g_alpha_prep = P::prepare_g1(&test1_a_g_alpha)?;
            let test1_b_h_beta_prep = P::prepare_g2(&test1_b_h_beta)?;

            let g_psi_prep = P::prepare_g1(&g_psi)?;

            let c_prep = P::prepare_g1(&proof.c)?;

            P::miller_loop(
                &[
                    test1_a_g_alpha_prep,
                    g_psi_prep,
                    c_prep,
                    circuit_pvk.g_alpha_pc.clone(),
                ],
                &[
                    test1_b_h_beta_prep,
                    circuit_pvk.h_gamma_pc.clone(),
                    circuit_pvk.h_pc.clone(),
                    circuit_pvk.h_beta_pc.clone(),
                ],
            )?
        };

        let test1 = P::final_exponentiation(&test1_exp).unwrap();

        // e(A, H^{gamma}) = e(G^{gamma}, B)
        let test2_exp = {
            let a_prep = P::prepare_g1(&proof.a)?;
            // pvk.h_gamma_pc
            //&pvk.g_gamma_pc
            let proof_b = proof.b.negate()?;
            let b_prep = P::prepare_g2(&proof_b)?;
            P::miller_loop(
                &[a_prep, circuit_pvk.g_gamma_pc.clone()],
                &[circuit_pvk.h_gamma_pc, b_prep],
            )?
        };
        let test2 = P::final_exponentiation(&test2_exp)?;

        let one = P::GTVar::one();
        test1.is_eq(&one)?.and(&test2.is_eq(&one)?)
    }

    #[tracing::instrument(target = "r1cs", skip(circuit_vk, x, proof))]
    fn verify(
        circuit_vk: &Self::VerifyingKeyVar,
        x: &Self::InputVar,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<E::Fq>, SynthesisError> {
        let pvk = circuit_vk.prepare()?;
        Self::verify_with_processed_vk(&pvk, x, proof)
    }
}

impl<E, P> CircuitSpecificSetupSNARKGadget<E::Fr, E::Fq, GM17<E>> for GM17VerifierGadget<E, P>
where
    E: PairingEngine,
    P: PairingVar<E, E::Fq>,
{
}

impl<E, P> AllocVar<PreparedVerifyingKey<E>, E::Fq> for PreparedVerifyingKeyVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<PreparedVerifyingKey<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        f().and_then(|pvk| {
            let pvk = pvk.borrow();
            let g_alpha = P::G1Var::new_variable(
                ark_relations::ns!(cs, "g_alpha"),
                || Ok(pvk.g_alpha),
                mode,
            )?;
            let h_beta =
                P::G2Var::new_variable(ark_relations::ns!(cs, "h_beta"), || Ok(pvk.h_beta), mode)?;
            let g_alpha_pc = P::G1PreparedVar::new_variable(
                ark_relations::ns!(cs, "g_alpha_pc"),
                || Ok(pvk.g_alpha.into()),
                mode,
            )?;
            let h_beta_pc = P::G2PreparedVar::new_variable(
                ark_relations::ns!(cs, "h_beta_pc"),
                || Ok(pvk.h_beta.into()),
                mode,
            )?;
            let g_gamma_pc = P::G1PreparedVar::new_variable(
                ark_relations::ns!(cs, "g_gamma_pc"),
                || Ok(&pvk.g_gamma_pc),
                mode,
            )?;
            let h_gamma_pc = P::G2PreparedVar::new_variable(
                ark_relations::ns!(cs, "h_gamma_pc"),
                || Ok(&pvk.h_gamma_pc),
                mode,
            )?;
            let h_pc = P::G2PreparedVar::new_variable(
                ark_relations::ns!(cs, "h_pc"),
                || Ok(&pvk.h_pc),
                mode,
            )?;
            let query = Vec::new_variable(
                ark_relations::ns!(cs, "query"),
                || Ok(pvk.query.clone()),
                mode,
            )?;

            Ok(Self {
                g_alpha,
                h_beta,
                g_alpha_pc,
                h_beta_pc,
                g_gamma_pc,
                h_gamma_pc,
                h_pc,
                query,
            })
        })
    }
}

impl<E, P> AllocVar<VerifyingKey<E>, E::Fq> for VerifyingKeyVar<E, P>
where
    E: PairingEngine,

    P: PairingVar<E>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<VerifyingKey<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        f().and_then(|vk| {
            let vk = vk.borrow();
            let g_alpha_g1 = P::G1Var::new_variable(
                ark_relations::ns!(cs, "g_alpha"),
                || Ok(vk.g_alpha_g1),
                mode,
            )?;
            let h_g2 = P::G2Var::new_variable(ark_relations::ns!(cs, "h"), || Ok(vk.h_g2), mode)?;
            let h_beta_g2 = P::G2Var::new_variable(
                ark_relations::ns!(cs, "h_beta"),
                || Ok(vk.h_beta_g2),
                mode,
            )?;
            let g_gamma_g1 = P::G1Var::new_variable(
                ark_relations::ns!(cs, "g_gamma"),
                || Ok(&vk.g_gamma_g1),
                mode,
            )?;
            let h_gamma_g2 = P::G2Var::new_variable(
                ark_relations::ns!(cs, "h_gamma"),
                || Ok(&vk.h_gamma_g2),
                mode,
            )?;
            let query = Vec::new_variable(
                ark_relations::ns!(cs, "query"),
                || Ok(vk.query.clone()),
                mode,
            )?;
            Ok(Self {
                h_g2,
                g_alpha_g1,
                h_beta_g2,
                g_gamma_g1,
                h_gamma_g2,
                query,
            })
        })
    }
}

impl<E, P> AllocVar<Proof<E>, E::Fq> for ProofVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<Proof<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        f().and_then(|proof| {
            let Proof { a, b, c } = proof.borrow().clone();
            let a = P::G1Var::new_variable(ark_relations::ns!(cs, "a"), || Ok(a), mode)?;
            let b = P::G2Var::new_variable(ark_relations::ns!(cs, "b"), || Ok(b), mode)?;
            let c = P::G1Var::new_variable(ark_relations::ns!(cs, "c"), || Ok(c), mode)?;
            Ok(Self { a, b, c })
        })
    }
}

impl<E, P> ToBytesGadget<E::Fq> for VerifyingKeyVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[inline]
    #[tracing::instrument(target = "r1cs", skip(self))]
    fn to_bytes(&self) -> Result<Vec<UInt8<E::Fq>>, SynthesisError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.h_g2.to_bytes()?);
        bytes.extend_from_slice(&self.g_alpha_g1.to_bytes()?);
        bytes.extend_from_slice(&self.h_beta_g2.to_bytes()?);
        bytes.extend_from_slice(&self.g_gamma_g1.to_bytes()?);
        bytes.extend_from_slice(&self.h_gamma_g2.to_bytes()?);
        for q in &self.query {
            bytes.extend_from_slice(&q.to_bytes()?);
        }
        Ok(bytes)
    }
}

#[cfg(test)]
mod test {
    use crate::{constraints::GM17VerifierGadget, GM17};
    use ark_crypto_primitives::snark::constraints::SNARKGadget;
    use ark_crypto_primitives::snark::{CircuitSpecificSetupSNARK, SNARK};
    use ark_ec::PairingEngine;
    use ark_ff::{test_rng, Field, UniformRand};
    use ark_mnt4_298::{
        constraints::PairingVar as MNT4PairingVar, Fr as MNT4Fr, MNT4_298 as MNT4PairingEngine,
    };
    use ark_mnt6_298::Fr as MNT6Fr;
    use ark_r1cs_std::bits::boolean::Boolean;
    use ark_r1cs_std::{alloc::AllocVar, eq::EqGadget};
    use ark_relations::{
        lc, ns,
        r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError},
    };
    use ark_std::ops::MulAssign;

    #[derive(Copy, Clone)]
    struct Circuit<F: Field> {
        a: Option<F>,
        b: Option<F>,
        num_constraints: usize,
        num_variables: usize,
    }

    impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for Circuit<ConstraintF> {
        fn generate_constraints(
            self,
            cs: ConstraintSystemRef<ConstraintF>,
        ) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
            let c = cs.new_input_variable(|| {
                let mut a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

                a.mul_assign(&b);
                Ok(a)
            })?;

            for _ in 0..(self.num_variables - 3) {
                let _ =
                    cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            }

            for _ in 0..self.num_constraints {
                cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)
                    .unwrap();
            }
            Ok(())
        }
    }

    type TestSNARK = GM17<MNT4PairingEngine>;
    type TestSNARKGadget = GM17VerifierGadget<MNT4PairingEngine, MNT4PairingVar>;

    #[test]
    fn gm17_snark_test() {
        let mut rng = test_rng();
        let a = MNT4Fr::rand(&mut rng);
        let b = MNT4Fr::rand(&mut rng);
        let mut c = a;
        c.mul_assign(&b);

        let circ = Circuit {
            a: Some(a.clone()),
            b: Some(b.clone()),
            num_constraints: 100,
            num_variables: 25,
        };

        let (pk, vk) = TestSNARK::setup(circ, &mut rng).unwrap();

        let proof = TestSNARK::prove(&pk, circ.clone(), &mut rng).unwrap();

        assert!(
            TestSNARK::verify(&vk, &vec![c], &proof).unwrap(),
            "The native verification check fails."
        );

        let cs_sys = ConstraintSystem::<MNT6Fr>::new();
        let cs = ConstraintSystemRef::new(cs_sys);

        let input_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::InputVar::new_input(ns!(cs, "new_input"), || Ok(vec![c]))
        .unwrap();
        let proof_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::ProofVar::new_witness(ns!(cs, "alloc_proof"), || Ok(proof))
        .unwrap();
        let vk_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::VerifyingKeyVar::new_constant(ns!(cs, "alloc_vk"), vk.clone())
        .unwrap();
        <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::verify(&vk_gadget, &input_gadget, &proof_gadget)
        .unwrap()
        .enforce_equal(&Boolean::constant(true))
        .unwrap();

        assert!(
            cs.is_satisfied().unwrap(),
            "Constraints not satisfied: {}",
            cs.which_is_unsatisfied().unwrap().unwrap_or_default()
        );

        let pvk = TestSNARK::process_vk(&vk).unwrap();
        let pvk_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::ProcessedVerifyingKeyVar::new_constant(
            ns!(cs, "alloc_pvk"), pvk.clone()
        )
        .unwrap();
        TestSNARKGadget::verify_with_processed_vk(&pvk_gadget, &input_gadget, &proof_gadget)
            .unwrap()
            .enforce_equal(&Boolean::constant(true))
            .unwrap();

        assert!(
            cs.is_satisfied().unwrap(),
            "Constraints not satisfied: {}",
            cs.which_is_unsatisfied().unwrap().unwrap_or_default()
        );
    }
}
