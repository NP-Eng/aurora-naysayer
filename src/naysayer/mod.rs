use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::LabeledCommitment;
use ark_poly_naysayer::PCSNaysayer;
use ark_std::rand::RngCore;

use crate::{aurora::{AuroraProof, error::AuroraError}, AuroraR1CS};

pub enum AuroraNaysayerProof {
    /// TODO documentation
    Case1
}

// pub fn naysay<PCS: PolynomialCommitment<F, DensePolynomial<F>>>(
//     &self,
//     vk: &PCS::VerifierKey,
//     instance: Vec<F>,
//     aurora_proof: AuroraProof<F, PCS>,
//     sponge: &mut impl CryptographicSponge,
// ) -> Result<Option<NaysayerProof>, AuroraError<F, PCS>> {
//     assert_padded(&self.r1cs);

//     let matrices = self.r1cs.to_matrices().unwrap();

//     // 0. Initialising sponge with public parameters
//     absorb_public_parameters::<F, PCS>(vk, &matrices, sponge);

//     let ConstraintMatrices {
//         a,
//         b,
//         c,
//         num_constraints: n,
//         num_instance_variables,
//         num_witness_variables,
//         ..
//     } = matrices;

//     let AuroraProof {
//         commitments: com,
//         proof,
//         evals,
//     } = aurora_proof;

//     // Checking instance and witness lengths
//     if instance.len() != self.unpadded_num_instance_variables {
//         return Err(AuroraError::IncorrectInstanceLength {
//             received: instance.len(),
//             expected: self.unpadded_num_instance_variables,
//         });
//     }

//     // Resize the instance to the padded length
//     let mut instance = instance;
//     instance.resize(num_instance_variables, F::ZERO);

//     // Absorb the first 5 commitments
//     sponge.absorb(&com.iter().take(5).collect::<Vec<_>>());

//     let r: F = sponge.squeeze_field_elements(1)[0];

//     // ======================== Verify the proof ========================

//     // Absorb the missing commitments to g1, g2
//     sponge.absorb(&com.iter().skip(5).collect::<Vec<_>>());

//     let a_point: F = sponge.squeeze_field_elements(1)[0];

//     if !PCS::check(vk, &com, &a_point, evals.clone(), &proof, sponge, None).unwrap() {
//         return Ok(false);
//     }

//     // ======================== Zero test ========================

//     // Evaluations of the committed polynomials at a_point
//     let [f_a_a, f_b_a, f_c_a, f_0_a, f_w_a, g_1_a, g_2_a] = evals[..] else {
//         return Err(AuroraError::IncorrectNumberOfEvaluations {
//             received: evals.len(),
//             expected: 7,
//         });
//     };

//     let h = GeneralEvaluationDomain::<F>::new(n).unwrap();

//     let v_h_a = h.evaluate_vanishing_polynomial(a_point);

//     if f_a_a * f_b_a - f_c_a != f_0_a * v_h_a {
//         return Ok(false);
//     }

//     // ======================== Univariate sumcheck test ========================
//     let zero_padded_instance =
//         [instance.clone(), vec![F::ZERO; num_witness_variables]].concat();

//     let lagrange_basis_evals = h.evaluate_all_lagrange_coefficients(a_point);

//     // Returns f(a_point), where f is the unique polynomial of degree < n that
//     // interpolates the given values over h. This requires
//     //  - a one-time evaluation of the Lagrange basis over h at a_point
//     //    (lagrange_basis_evals), which is amortised over all calls
//     //  - a one-time inner product of size n per call.
//     let h_evaluate_interpolator = |evals: &Vec<F>| inner_product(&evals, &lagrange_basis_evals);

//     let v_star_a = h_evaluate_interpolator(&zero_padded_instance);

//     let v_h_in_a: F = h
//         .elements()
//         .take(num_instance_variables)
//         .map(|c| a_point - c)
//         .product();

//     let f_z_a = f_w_a * v_h_in_a + v_star_a;

//     // Computing [1, r, r^2, ... r^(n-1)]
//     let powers_of_r = powers(r, n);
//     let p_r_a = h_evaluate_interpolator(&powers_of_r);

//     let q_ar_a =
//         h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&a, &powers_of_r));
//     let q_br_a =
//         h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&b, &powers_of_r));
//     let q_cr_a =
//         h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&c, &powers_of_r));

//     let r_pow_n = r * powers_of_r[n - 1];

//     let u_a = (p_r_a * f_a_a - q_ar_a * f_z_a)
//         + (p_r_a * f_b_a - q_br_a * f_z_a) * r_pow_n
//         + (p_r_a * f_c_a - q_cr_a * f_z_a) * (r_pow_n * r_pow_n);

//     Ok(u_a == g_1_a * v_h_a + g_2_a * a_point)
// }

/// Naysays the given proof, returning
/// - Err(e) if non-assertion error e was encountered
/// - Ok(None) if no errors were encountered and the proof is valid ("aye")
/// - Ok(Some(naysayer_proof)) if an assertion error was encountered
fn aurora_naysay<F: PrimeField + Absorb, NS: PCSNaysayer<F, DensePolynomial<F>>>(
    aurora_r1cs: &AuroraR1CS<F>,
    vk: &NS::VerifierKey,
    instance: Vec<F>,
    aurora_proof: AuroraProof<F, NS>,
    sponge: &mut impl CryptographicSponge,
) -> Result<Option<AuroraNaysayerProof>, AuroraError<F, NS>>
{
    unimplemented!();
}

/// Verifies the naysayer proof. Returns:
/// - Ok(true) if the original proof is rejected (i.e. the naysayer proof
///   points to a valid issue).
/// - Ok(false) if the original proof is not rejected, i.e. the naysayer
///   proof points to a non-issue
/// - Err if another type of error occurs during verification of the
///   naysayer proof.
fn verify_naysay<'a, F: PrimeField + Absorb, NS: PCSNaysayer<F, DensePolynomial<F>>>(
    vk: &NS::VerifierKey,
    commitments: impl IntoIterator<Item = &'a LabeledCommitment<NS::Commitment>>,
    point: &F,
    values: impl IntoIterator<Item = F>,
    proof_array: &NS::Proof,
    naysayer_proof: &NS::NaysayerProof,
    sponge: &mut impl CryptographicSponge,
    rng: Option<&mut dyn RngCore>,
) -> Result<bool, AuroraError<F, NS>> 
where 
    NS::Commitment: 'a,
{
    unimplemented!();
}
