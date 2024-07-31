use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_naysayer::PCSNaysayer;
use ark_relations::r1cs::ConstraintMatrices;
use derivative::Derivative;

use crate::{
    aurora::{error::AuroraError, utils::*, AuroraProof, AuroraVerifierKey},
    AuroraR1CS,
};

#[cfg(any(test, feature = "bench"))]
pub mod tests;

#[derive(Derivative)]
#[derivative(Debug(bound = ""), PartialEq(bound = ""))]
pub enum AuroraNaysayerProof<F: PrimeField + Absorb, NS: PCSNaysayer<F, DensePolynomial<F>>> {
    // The proof of PCS: open is incorrect for one of f_a, f_b, f_c, f_0, f_w and/or g_1
    PCSLarge(NS::NaysayerProof),
    // The proof of PCS: open is incorrect for g_2
    PCSG2(NS::NaysayerProof),
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
pub fn aurora_naysay<F, NS>(
    vk: &AuroraVerifierKey<F, NS>,
    aurora_proof: &AuroraProof<F, NS>,
    instance: &Vec<F>,
    sponge: &mut impl CryptographicSponge,
) -> Result<Option<AuroraNaysayerProof<F, NS>>, AuroraError<F, NS>>
where
    F: PrimeField + Absorb,
    NS: PCSNaysayer<F, DensePolynomial<F>>,
{
    let AuroraVerifierKey {
        r1cs: AuroraR1CS {
            r1cs,
            unpadded_num_instance_variables,
        },
        vk_large,
        vk_small,
    } = vk;

    // TODO Return non-assertion error
    assert_padded(r1cs);

    // TODO Return non-assertion error
    let matrices = r1cs.to_matrices().unwrap();

    // 0. Initialising sponge with public parameters
    absorb_public_parameters::<F, NS>((vk_large, vk_small), &matrices, sponge);

    let ConstraintMatrices {
        a,
        b,
        c,
        num_constraints: n,
        num_instance_variables,
        ..
    } = matrices;

    let AuroraProof {
        large_coms,
        com_g_2,
        large_opening_proof,
        g_2_opening_proof,
        large_evals,
        g_2_a,
    } = aurora_proof;

    // Checking instance and witness lengths
    // TODO Return non-assertion error
    if instance.len() != *unpadded_num_instance_variables {
        return Err(AuroraError::IncorrectInstanceLength {
            received: instance.len(),
            expected: *unpadded_num_instance_variables,
        });
    }

    // ======================== Naysay the proof ========================

    // Absorb the first 5 commitments
    sponge.absorb(&large_coms.iter().take(5).collect::<Vec<_>>());

    let r: F = sponge.squeeze_field_elements(1)[0];

    // Absorb the missing commitments to g1, g2
    sponge.absorb(&large_coms.last().unwrap());
    sponge.absorb(&com_g_2);

    let a_point: F = sponge.squeeze_field_elements(1)[0];

    let pcs_naysay_large = NS::naysay(
        vk_large,
        large_coms,
        &a_point,
        large_evals.clone(),
        &large_opening_proof,
        sponge,
        None,
    )
    .unwrap();

    let g_2_naysay = NS::naysay(
        vk_small,
        [com_g_2],
        &a_point,
        [*g_2_a],
        &g_2_opening_proof,
        sponge,
        None,
    )
    .unwrap();

    if let Some(pcs_naysayer_proof) = pcs_naysay_large {
        return Ok(Some(AuroraNaysayerProof::PCSLarge(pcs_naysayer_proof)));
    }

    if let Some(g_2_naysayer_proof) = g_2_naysay {
        return Ok(Some(AuroraNaysayerProof::PCSG2(g_2_naysayer_proof)));
    }

    // ======================== Zero test ========================

    // Evaluations of the committed polynomials at a_point
    // TODO Return non-assertion error
    let [f_a_a, f_b_a, f_c_a, f_0_a, f_w_a, g_1_a] = large_evals[..] else {
        return Err(AuroraError::IncorrectNumberOfEvaluations {
            received: large_evals.len(),
            expected: 6,
        });
    };

    // TODO Return non-assertion error
    let h = GeneralEvaluationDomain::<F>::new(n).unwrap();

    let v_h_a = h.evaluate_vanishing_polynomial(a_point);

    if f_a_a * f_b_a - f_c_a != f_0_a * v_h_a {
        return Ok(Some(AuroraNaysayerProof::ZeroCheck));
    }

    // ======================== Univariate sumcheck test ========================
    let lagrange_basis_evals = h.evaluate_all_lagrange_coefficients(a_point);

    // Returns f(a_point), where f is the unique polynomial of degree < n that
    // interpolates the given values over h. This requires
    //  - a one-time evaluation of the Lagrange basis over h at a_point
    //    (lagrange_basis_evals), which is amortised over all calls
    //  - a one-time inner product of size n per call.
    let h_evaluate_interpolator =
        |evals: &Vec<F>| inner_product(&evals, &lagrange_basis_evals[0..evals.len()]);

    let v_star_a = h_evaluate_interpolator(&instance);

    let v_h_in_a: F = h
        .elements()
        .take(num_instance_variables)
        .map(|c| a_point - c)
        .product();

    let f_z_a = f_w_a * v_h_in_a + v_star_a;

    // Computing [1, r, r^2, ... r^(n-1)]
    let powers_of_r = powers(r, n);
    let p_r_a = h_evaluate_interpolator(&powers_of_r);

    let q_ar_a = h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&a, &powers_of_r));
    let q_br_a = h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&b, &powers_of_r));
    let q_cr_a = h_evaluate_interpolator(&random_matrix_polynomial_evaluations(&c, &powers_of_r));

    let r_pow_n = r * powers_of_r[n - 1];

    let u_a = (p_r_a * f_a_a - q_ar_a * f_z_a)
        + (p_r_a * f_b_a - q_br_a * f_z_a) * r_pow_n
        + (p_r_a * f_c_a - q_cr_a * f_z_a) * (r_pow_n * r_pow_n);

    if u_a != g_1_a * v_h_a + *g_2_a * a_point {
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
pub fn aurora_verify_naysay<'a, F, NS>(
    vk: &AuroraVerifierKey<F, NS>,
    original_proof: &AuroraProof<F, NS>,
    naysayer_proof: &AuroraNaysayerProof<F, NS>,
    instance: Vec<F>,
    sponge: &mut impl CryptographicSponge,
) -> Result<bool, AuroraError<F, NS>>
where
    F: PrimeField + Absorb,
    NS: PCSNaysayer<F, DensePolynomial<F>>,
{
    let AuroraVerifierKey {
        r1cs: AuroraR1CS { r1cs, .. },
        vk_large,
        vk_small,
    } = vk;

    let matrices = r1cs.to_matrices().unwrap();

    // 0. Initialising sponge with public parameters
    absorb_public_parameters::<F, NS>((&vk_large, &vk_small), &matrices, sponge);

    let AuroraProof {
        large_coms,
        com_g_2,
        large_opening_proof,
        g_2_opening_proof,
        large_evals,
        g_2_a,
    } = original_proof;

    let ConstraintMatrices {
        a,
        b,
        c,
        num_constraints: n,
        num_instance_variables,
        ..
    } = matrices;

    // =================== Naysayer proof for the PCS ===================

    // Absorb the first 5 commitments
    sponge.absorb(&large_coms.iter().take(5).collect::<Vec<_>>());

    let r: F = sponge.squeeze_field_elements(1)[0];

    // Absorb the missing commitments to g1, g2
    sponge.absorb(&large_coms.last().unwrap());
    sponge.absorb(&com_g_2);

    let a_point: F = sponge.squeeze_field_elements(1)[0];

    // TODO review code structure and decide whether it is better to
    // transform this match into a list of sequential ifs,
    match naysayer_proof {
        AuroraNaysayerProof::PCSLarge(pcs_naysayer_proof) => Ok(NS::verify_naysay(
            vk_large,
            large_coms,
            &a_point,
            large_evals.clone(),
            large_opening_proof,
            pcs_naysayer_proof,
            sponge,
            None,
        )
        .unwrap()),
        AuroraNaysayerProof::PCSG2(pcs_naysayer_proof) => Ok(NS::verify_naysay(
            vk_small,
            &[com_g_2.clone()],
            &a_point,
            vec![g_2_a.clone()],
            g_2_opening_proof,
            pcs_naysayer_proof,
            sponge,
            None,
        )
        .unwrap()),
        _ => {
            // ======================== Zero test ========================

            // Evaluations of the committed polynomials at a_point
            // TODO Return non-assertion error
            let [f_a_a, f_b_a, f_c_a, f_0_a, f_w_a, g_1_a] = large_evals[..] else {
                return Err(AuroraError::IncorrectNumberOfEvaluations {
                    received: large_evals.len(),
                    expected: 6,
                });
            };

            // TODO Return non-assertion error
            let h = GeneralEvaluationDomain::<F>::new(n).unwrap();

            let v_h_a = h.evaluate_vanishing_polynomial(a_point);

            match naysayer_proof {
                AuroraNaysayerProof::ZeroCheck => Ok(f_a_a * f_b_a - f_c_a != f_0_a * v_h_a),
                AuroraNaysayerProof::UnivariateSumcheck => {
                    let lagrange_basis_evals = h.evaluate_all_lagrange_coefficients(a_point);

                    let h_evaluate_interpolator = |evals: &Vec<F>| {
                        inner_product(&evals, &lagrange_basis_evals[0..evals.len()])
                    };

                    let v_star_a = h_evaluate_interpolator(&instance);

                    let v_h_in_a: F = h
                        .elements()
                        .take(num_instance_variables)
                        .map(|c| a_point - c)
                        .product();

                    let f_z_a = f_w_a * v_h_in_a + v_star_a;

                    // Computing [1, r, r^2, ... r^(n-1)]
                    let powers_of_r = powers(r, n);
                    let p_r_a = h_evaluate_interpolator(&powers_of_r);

                    let q_ar_a = h_evaluate_interpolator(&random_matrix_polynomial_evaluations(
                        &a,
                        &powers_of_r,
                    ));
                    let q_br_a = h_evaluate_interpolator(&random_matrix_polynomial_evaluations(
                        &b,
                        &powers_of_r,
                    ));
                    let q_cr_a = h_evaluate_interpolator(&random_matrix_polynomial_evaluations(
                        &c,
                        &powers_of_r,
                    ));

                    let r_pow_n = r * powers_of_r[n - 1];

                    let u_a = (p_r_a * f_a_a - q_ar_a * f_z_a)
                        + (p_r_a * f_b_a - q_br_a * f_z_a) * r_pow_n
                        + (p_r_a * f_c_a - q_cr_a * f_z_a) * (r_pow_n * r_pow_n);

                    return Ok(u_a != g_1_a * v_h_a + *g_2_a * a_point);
                }
                _ => unreachable!(),
            }
        }
    }
}
