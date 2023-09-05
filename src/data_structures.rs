use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalSerialize, *};
use ark_std::vec::Vec;

/// A proof in the GM17 SNARK.
#[derive(PartialEq, Eq, Clone, Default, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: Pairing> {
    #[doc(hidden)]
    pub a: E::G1Affine,
    #[doc(hidden)]
    pub b: E::G2Affine,
    #[doc(hidden)]
    pub c: E::G1Affine,
}

/// A verification key in the GM17 SNARK.
#[derive(Eq, PartialEq, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifyingKey<E: Pairing> {
    #[doc(hidden)]
    pub h_g2: E::G2Affine,
    #[doc(hidden)]
    pub g_alpha_g1: E::G1Affine,
    #[doc(hidden)]
    pub h_beta_g2: E::G2Affine,
    #[doc(hidden)]
    pub g_gamma_g1: E::G1Affine,
    #[doc(hidden)]
    pub h_gamma_g2: E::G2Affine,
    #[doc(hidden)]
    pub query: Vec<E::G1Affine>,
}

impl<E: Pairing> Default for VerifyingKey<E> {
    fn default() -> Self {
        Self {
            h_g2: E::G2Affine::default(),
            g_alpha_g1: E::G1Affine::default(),
            h_beta_g2: E::G2Affine::default(),
            g_gamma_g1: E::G1Affine::default(),
            h_gamma_g2: E::G2Affine::default(),
            query: Vec::new(),
        }
    }
}

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(PartialEq, Eq, Clone, CanonicalDeserialize, CanonicalSerialize)]
pub struct PreparedVerifyingKey<E: Pairing> {
    #[doc(hidden)]
    pub vk: VerifyingKey<E>,
    #[doc(hidden)]
    pub g_alpha: E::G1Affine,
    #[doc(hidden)]
    pub h_beta: E::G2Affine,
    #[doc(hidden)]
    pub g_alpha_h_beta_ml: E::TargetField,
    #[doc(hidden)]
    pub g_gamma_pc: E::G1Prepared,
    #[doc(hidden)]
    pub h_gamma_pc: E::G2Prepared,
    #[doc(hidden)]
    pub h_pc: E::G2Prepared,
    #[doc(hidden)]
    pub query: Vec<E::G1Affine>,
}

impl<E: Pairing> Default for PreparedVerifyingKey<E> {
    fn default() -> Self {
        Self {
            vk: VerifyingKey::default(),
            g_alpha: E::G1Affine::default(),
            h_beta: E::G2Affine::default(),
            g_alpha_h_beta_ml: E::TargetField::default(),
            g_gamma_pc: E::G1Prepared::default(),
            h_gamma_pc: E::G2Prepared::default(),
            h_pc: E::G2Prepared::default(),
            query: Vec::new(),
        }
    }
}

impl<E: Pairing> From<PreparedVerifyingKey<E>> for VerifyingKey<E> {
    fn from(other: PreparedVerifyingKey<E>) -> Self {
        other.vk
    }
}

impl<E: Pairing> From<VerifyingKey<E>> for PreparedVerifyingKey<E> {
    fn from(other: VerifyingKey<E>) -> Self {
        crate::prepare_verifying_key(&other)
    }
}

/// Full public (prover and verifier) parameters for the GM17 zkSNARK.
#[derive(PartialEq, Eq, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct ProvingKey<E: Pairing> {
    #[doc(hidden)]
    pub vk: VerifyingKey<E>,
    #[doc(hidden)]
    pub a_query: Vec<E::G1Affine>,
    #[doc(hidden)]
    pub b_query: Vec<E::G2Affine>,
    #[doc(hidden)]
    pub c_query_1: Vec<E::G1Affine>,
    #[doc(hidden)]
    pub c_query_2: Vec<E::G1Affine>,
    #[doc(hidden)]
    pub g_gamma_z: E::G1Affine,
    #[doc(hidden)]
    pub h_gamma_z: E::G2Affine,
    #[doc(hidden)]
    pub g_ab_gamma_z: E::G1Affine,
    #[doc(hidden)]
    pub g_gamma2_z2: E::G1Affine,
    #[doc(hidden)]
    pub g_gamma2_z_t: Vec<E::G1Affine>,
}
