use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, EvaluationDomain};
use ark_poly_naysayer::PCSNaysayer;
use ark_relations::r1cs::ConstraintMatrices;

use crate::{aurora::{AuroraProof, error::AuroraError, utils::*}, AuroraR1CS};

pub enum AuroraNaysayerProof<F: PrimeField + Absorb, NS: PCSNaysayer<F, DensePolynomial<F>>> {
    // The proof of PCS: open is incorrect, and the NS:NaysayerProof shows it
    PCS(NS::NaysayerProof),
    // The zero check equality does not hold
    ZeroCheck,
    // The univariate sumcheck equality does not hold
    UnivariateSumcheck,
}

/// Naysays the given proof, returning
/// - Err(e) if non-assertion error e was encountered
/// - Ok(None) if no errors were encountered and the proof is valid ("aye")
/// - Ok(Some(naysayer_proof)) if an assertion error was encountered

// TODO It's possible both this and Ligero's naysay will never return Err once
// we have the comprehensive non-assertion-error-handling implementation, in
// which case non-assertion-errors will be one more case of the naysayer proof
// and the return type here will be Option<AuroraNaysayerProof>
fn aurora_naysay<F: PrimeField + Absorb, NS: PCSNaysayer<F, DensePolynomial<F>>>(
    aurora_r1cs: &AuroraR1CS<F>,
    vk: &NS::VerifierKey,
    aurora_proof: AuroraProof<F, NS>,
    instance: Vec<F>,
    sponge: &mut impl CryptographicSponge,
) -> Result<Option<AuroraNaysayerProof<F, NS>>, AuroraError<F, NS>>
{
    
    let r1cs = aurora_r1cs.r1cs();
    
    // TODO Return non-assertion error
    assert_padded(r1cs);

    // TODO Return non-assertion error
    let matrices = r1cs.to_matrices().unwrap();

    // 0. Initialising sponge with public parameters
    absorb_public_parameters::<F, NS>(vk, &matrices, sponge);

    let ConstraintMatrices {
        a,
        b,
        c,
        num_constraints: n,
        num_instance_variables,
        num_witness_variables,
        ..
    } = matrices;

    let AuroraProof {
        commitments: com,
        proof,
        evals,
    } = aurora_proof;

    // Checking instance and witness lengths
    // TODO Return non-assertion error
    if instance.len() != aurora_r1cs.unpadded_num_instance_variables {
        return Err(AuroraError::IncorrectInstanceLength {
            received: instance.len(),
            expected: aurora_r1cs.unpadded_num_instance_variables,
        });
    }

    // Resize the instance to the padded length
    let mut instance = instance;
    instance.resize(num_instance_variables, F::ZERO);

    // Absorb the first 5 commitments
    sponge.absorb(&com.iter().take(5).collect::<Vec<_>>());

    let r: F = sponge.squeeze_field_elements(1)[0];

    // ======================== Verify the proof ========================

    // Absorb the missing commitments to g1, g2
    sponge.absorb(&com.iter().skip(5).collect::<Vec<_>>());

    let a_point: F = sponge.squeeze_field_elements(1)[0];

    // TODO Return non-assertion error
    let pcs_naysay = NS::naysay(vk, &com, &a_point, evals.clone(), &proof, sponge, None).unwrap();

    if let Some(pcs_naysayer_proof) = pcs_naysay {
        return Ok(Some(AuroraNaysayerProof::PCS(pcs_naysayer_proof)));
    }
    
    // ======================== Zero test ========================

    // Evaluations of the committed polynomials at a_point
    // TODO Return non-assertion error
    let [f_a_a, f_b_a, f_c_a, f_0_a, f_w_a, g_1_a, g_2_a] = evals[..] else {
        return Err(AuroraError::IncorrectNumberOfEvaluations {
            received: evals.len(),
            expected: 7,
        });
    };

    // TODO Return non-assertion error
    let h = GeneralEvaluationDomain::<F>::new(n).unwrap();

    let v_h_a = h.evaluate_vanishing_polynomial(a_point);

    if f_a_a * f_b_a - f_c_a != f_0_a * v_h_a {
        return Ok(Some(AuroraNaysayerProof::ZeroCheck));
    }

    // ======================== Univariate sumcheck test ========================
    let zero_padded_instance =
        [instance.clone(), vec![F::ZERO; num_witness_variables]].concat();

    let lagrange_basis_evals = h.evaluate_all_lagrange_coefficients(a_point);

    // Returns f(a_point), where f is the unique polynomial of degree < n that
    // interpolates the given values over h. This requires
    //  - a one-time evaluation of the Lagrange basis over h at a_point
    //    (lagrange_basis_evals), which is amortised over all calls
    //  - a one-time inner product of size n per call.
    let h_evaluate_interpolator = |evals: &Vec<F>| inner_product(&evals, &lagrange_basis_evals);

    let v_star_a = h_evaluate_interpolator(&zero_padded_instance);

    let v_h_in_a: F = h
        .elements()
        .take(num_instance_variables)
        .map(|c| a_point - c)
        .product();

    let f_z_a = f_w_a * v_h_in_a + v_star_a;

    // Computing [1, r, r^2, ... r^(n-1)]
    let powers_of_r = powers(r, n);
    let p_r_a = h_evaluate_interpolator(&powers_of_r);

    let q_ar_a =
        h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&a, &powers_of_r));
    let q_br_a =
        h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&b, &powers_of_r));
    let q_cr_a =
        h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&c, &powers_of_r));

    let r_pow_n = r * powers_of_r[n - 1];

    let u_a = (p_r_a * f_a_a - q_ar_a * f_z_a)
        + (p_r_a * f_b_a - q_br_a * f_z_a) * r_pow_n
        + (p_r_a * f_c_a - q_cr_a * f_z_a) * (r_pow_n * r_pow_n);

    if u_a != g_1_a * v_h_a + g_2_a * a_point {
        return Ok(Some(AuroraNaysayerProof::UnivariateSumcheck));
    }

    Ok(None)
}

/// Verifies the naysayer proof. Returns:
/// - Ok(true) if the original proof is rejected (i.e. the naysayer proof
///   points to a valid issue).
/// - Ok(false) if the original proof is not rejected, i.e. the naysayer
///   proof points to a non-issue
/// - Err if another type of error occurs during verification of the
///   naysayer proof.
fn aurora_verify_naysay<'a, F: PrimeField + Absorb, NS: PCSNaysayer<F, DensePolynomial<F>>>(
    vk: &NS::VerifierKey,
    aurora_r1cs: &AuroraR1CS<F>,
    original_proof: &AuroraProof<F, NS>,
    naysayer_proof: &AuroraNaysayerProof<F, NS>,
    instance: Vec<F>,
    sponge: &mut impl CryptographicSponge,
) -> Result<bool, AuroraError<F, NS>> {

    let r1cs = aurora_r1cs.r1cs();
    let matrices = r1cs.to_matrices().unwrap();

    // 0. Initialising sponge with public parameters
    absorb_public_parameters::<F, NS>(vk, &matrices, sponge);

    let AuroraProof {
        commitments: com,
        proof,
        evals,
    } = original_proof;

    let ConstraintMatrices {
        a,
        b,
        c,
        num_constraints: n,
        num_instance_variables,
        num_witness_variables,
        ..
    } = matrices;

    // Resize the instance to the padded length
    let mut instance = instance;
    instance.resize(num_instance_variables, F::ZERO);

    // Absorb the first 5 commitments
    sponge.absorb(&com.iter().take(5).collect::<Vec<_>>());

    let r: F = sponge.squeeze_field_elements(1)[0];

    // ======================== Verify the PCS proof ========================

    // Absorb the missing commitments to g1, g2
    sponge.absorb(&com.iter().skip(5).collect::<Vec<_>>());

    let a_point: F = sponge.squeeze_field_elements(1)[0];

    // TODO review code structure and decide whether it is better to
    // transform this match into a list of sequential ifs,
    match naysayer_proof {
        AuroraNaysayerProof::PCS(pcs_naysayer_proof) => Ok(NS::verify_naysay(vk, com, &a_point, evals.clone(), proof, pcs_naysayer_proof, sponge, None).unwrap()),
        _ => {

            // ======================== Zero test ========================

            // Evaluations of the committed polynomials at a_point
            // TODO Return non-assertion error
            let [f_a_a, f_b_a, f_c_a, f_0_a, f_w_a, g_1_a, g_2_a] = evals[..] else {
                return Err(AuroraError::IncorrectNumberOfEvaluations {
                    received: evals.len(),
                    expected: 7,
                });
            };

            // TODO Return non-assertion error
            let h = GeneralEvaluationDomain::<F>::new(n).unwrap();

            let v_h_a = h.evaluate_vanishing_polynomial(a_point);

            match naysayer_proof {
                AuroraNaysayerProof::ZeroCheck => Ok(f_a_a * f_b_a - f_c_a != f_0_a * v_h_a),
                AuroraNaysayerProof::UnivariateSumcheck => {
                    let zero_padded_instance =
                        [instance.clone(), vec![F::ZERO; num_witness_variables]].concat();

                    let lagrange_basis_evals = h.evaluate_all_lagrange_coefficients(a_point);

                    let h_evaluate_interpolator = |evals: &Vec<F>| inner_product(&evals, &lagrange_basis_evals);

                    let v_star_a = h_evaluate_interpolator(&zero_padded_instance);

                    let v_h_in_a: F = h
                        .elements()
                        .take(num_instance_variables)
                        .map(|c| a_point - c)
                        .product();

                    let f_z_a = f_w_a * v_h_in_a + v_star_a;

                    // Computing [1, r, r^2, ... r^(n-1)]
                    let powers_of_r = powers(r, n);
                    let p_r_a = h_evaluate_interpolator(&powers_of_r);

                    let q_ar_a =
                        h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&a, &powers_of_r));
                    let q_br_a =
                        h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&b, &powers_of_r));
                    let q_cr_a =
                        h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&c, &powers_of_r));

                    let r_pow_n = r * powers_of_r[n - 1];

                    let u_a = (p_r_a * f_a_a - q_ar_a * f_z_a)
                        + (p_r_a * f_b_a - q_br_a * f_z_a) * r_pow_n
                        + (p_r_a * f_c_a - q_cr_a * f_z_a) * (r_pow_n * r_pow_n);

                    return Ok(u_a != g_1_a * v_h_a + g_2_a * a_point);
                },
                _ => unreachable!(),
            }
        }  
    }
}
