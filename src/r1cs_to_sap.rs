use ark_ff::{Field, PrimeField};
use ark_poly::EvaluationDomain;
use ark_std::{cfg_chunks_mut, cfg_iter, cfg_iter_mut, vec, vec::Vec};

use ark_relations::r1cs::{ConstraintSystemRef, Result as R1CSResult, SynthesisError};
use core::ops::Deref;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub(crate) struct R1CStoSAP;

impl R1CStoSAP {
    #[inline]
    pub(crate) fn instance_map_with_evaluation<F: PrimeField, D: EvaluationDomain<F>>(
        cs: ConstraintSystemRef<F>,
        t: &F,
    ) -> R1CSResult<(Vec<F>, Vec<F>, F, usize, usize)> {
        let matrices = cs.to_matrices().unwrap();
        let num_inputs = cs.num_instance_variables();
        let num_aux = cs.num_witness_variables();
        let num_constraints = cs.num_constraints();

        let domain_size = 2 * num_constraints + 2 * (num_inputs - 1) + 1;
        let domain = D::new(domain_size).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let zt = domain.evaluate_vanishing_polynomial(*t);

        // Evaluate all Lagrange polynomials
        let coefficients_time = start_timer!(|| "Evaluate Lagrange coefficients");
        let u = domain.evaluate_all_lagrange_coefficients(*t);
        end_timer!(coefficients_time);

        let sap_num_variables = 2 * (num_inputs - 1) + num_aux + num_constraints;
        let extra_var_offset = (num_inputs - 1) + num_aux + 1;
        let extra_constr_offset = 2 * num_constraints;
        let extra_var_offset2 = (num_inputs - 1) + num_aux + num_constraints;

        let mut a = vec![F::zero(); sap_num_variables + 1];
        let mut c = vec![F::zero(); sap_num_variables + 1];

        for i in 0..num_constraints {
            let u_2i = u[2 * i];
            let u_2i_plus_1 = u[2 * i + 1];
            let u_add = u_2i + &u_2i_plus_1;
            let u_sub = u_2i - &u_2i_plus_1;

            for &(ref coeff, index) in matrices.a[i].iter() {
                a[index] += &(u_add * coeff);
            }

            for &(ref coeff, index) in matrices.b[i].iter() {
                a[index] += &(u_sub * coeff);
            }

            for &(ref coeff, index) in matrices.c[i].iter() {
                c[index] += &((u_2i * coeff).double().double());
            }
            c[extra_var_offset + i] += &u_add;
        }

        a[0] += &u[extra_constr_offset];
        c[0] += &u[extra_constr_offset];

        for i in 1..num_inputs {
            // First extra constraint

            a[i] += &u[extra_constr_offset + 2 * i - 1];
            a[0] += &u[extra_constr_offset + 2 * i - 1];

            let t_four = u[extra_constr_offset + 2 * i - 1].double().double();

            c[i] += &t_four;
            c[extra_var_offset2 + i] += &u[extra_constr_offset + 2 * i - 1];

            // Second extra constraint

            a[i] += &u[extra_constr_offset + 2 * i];
            a[0] -= &u[extra_constr_offset + 2 * i];
            c[extra_var_offset2 + i] += &u[extra_constr_offset + 2 * i];
        }

        Ok((a, c, zt, sap_num_variables, domain_size))
    }

    #[inline]
    pub(crate) fn witness_map<F: PrimeField, D: EvaluationDomain<F>>(
        prover: ConstraintSystemRef<F>,
        d1: &F,
        d2: &F,
    ) -> R1CSResult<(Vec<F>, Vec<F>, usize)> {
        #[inline]
        fn evaluate_constraint<F: Field>(terms: &[(F, usize)], assignment: &[F]) -> F {
            let mut acc = F::zero();
            for &(coeff, index) in terms {
                let val = assignment[index];
                acc += &(val * &coeff);
            }
            acc
        }

        let zero = F::zero();
        let one = F::one();
        let matrices = prover.to_matrices().unwrap();
        let num_inputs = prover.num_instance_variables();
        let num_aux = prover.num_witness_variables();
        let num_constraints = prover.num_constraints();
        let cs = prover.borrow().unwrap();
        let prover = cs.deref();

        let mut full_input_assignment = prover.instance_assignment.clone();
        full_input_assignment.extend_from_slice(&prover.witness_assignment);

        let temp = cfg_iter!(matrices.a)
            .zip(&matrices.b)
            .map(|(a_i, b_i)| {
                let mut extra_var = evaluate_constraint(&a_i, &full_input_assignment);
                extra_var -= &evaluate_constraint(&b_i, &full_input_assignment);
                extra_var.square_in_place();
                extra_var
            })
            .collect::<Vec<_>>();
        full_input_assignment.extend(temp);

        for i in 1..num_inputs {
            let mut extra_var = full_input_assignment[i];
            extra_var -= &one;
            extra_var.square_in_place();
            full_input_assignment.push(extra_var);
        }

        let domain = D::new(2 * num_constraints + 2 * (num_inputs - 1) + 1)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let coset_domain = domain.get_coset(F::GENERATOR).unwrap();

        let extra_constr_offset = 2 * num_constraints;
        let extra_var_offset = num_inputs + num_aux;
        let extra_var_offset2 = num_inputs + num_aux + num_constraints - 1;

        let mut a = vec![zero; domain_size];
        cfg_chunks_mut!(a[..2 * num_constraints], 2)
            .zip(&matrices.a)
            .zip(&matrices.b)
            .for_each(|((chunk, at_i), bt_i)| {
                chunk[0] = evaluate_constraint(&at_i, &full_input_assignment);
                chunk[0] += &evaluate_constraint(&bt_i, &full_input_assignment);

                chunk[1] = evaluate_constraint(&at_i, &full_input_assignment);
                chunk[1] -= &evaluate_constraint(&bt_i, &full_input_assignment);
            });
        a[extra_constr_offset] = one;
        for i in 1..num_inputs {
            a[extra_constr_offset + 2 * i - 1] = full_input_assignment[i] + &one;
            a[extra_constr_offset + 2 * i] = full_input_assignment[i] - &one;
        }

        domain.ifft_in_place(&mut a);

        let d1_double = d1.double();
        let mut h: Vec<F> = vec![d1_double; domain_size];
        cfg_iter_mut!(h).zip(&a).for_each(|(h_i, a_i)| *h_i *= a_i);
        h[0] -= d2;
        let d1d1 = d1.square();
        h[0] -= &d1d1;
        h.push(d1d1);

        coset_domain.fft_in_place(&mut a);

        let mut aa = domain.mul_polynomials_in_evaluation_domain(&a, &a);
        drop(a);

        let mut c = vec![zero; domain_size];
        cfg_chunks_mut!(c[..2 * num_constraints], 2)
            .enumerate()
            .for_each(|(i, chunk)| {
                let mut tmp = evaluate_constraint(&matrices.c[i], &full_input_assignment);
                tmp.double_in_place();
                tmp.double_in_place();

                let assignment = full_input_assignment[extra_var_offset + i];
                chunk[0] = tmp + &assignment;
                chunk[1] = assignment;
            });
        c[extra_constr_offset] = one;
        for i in 1..num_inputs {
            let mut tmp = full_input_assignment[i];
            tmp.double_in_place();
            tmp.double_in_place();

            let assignment = full_input_assignment[extra_var_offset2 + i];
            c[extra_constr_offset + 2 * i - 1] = tmp + &assignment;
            c[extra_constr_offset + 2 * i] = assignment;
        }

        domain.ifft_in_place(&mut c);
        coset_domain.fft_in_place(&mut c);

        cfg_iter_mut!(aa)
            .zip(c)
            .for_each(|(aa_i, c_i)| *aa_i -= &c_i);

        let i = domain
            .evaluate_vanishing_polynomial(F::GENERATOR)
            .inverse()
            .unwrap();
        ark_std::cfg_iter_mut!(aa).for_each(|eval| *eval *= &i);
        coset_domain.ifft_in_place(&mut aa);

        cfg_iter_mut!(h[..domain_size - 1])
            .enumerate()
            .for_each(|(i, e)| *e += &aa[i]);

        Ok((full_input_assignment, h, domain_size))
    }
}
