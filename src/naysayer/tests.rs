
use ark_bn254::Fr;
use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::{Field, PrimeField};
use ark_poly::{univariate::{DensePolynomial, SparsePolynomial}, DenseUVPolynomial};
use ark_poly_commit::{test_sponge, PolynomialCommitment, TestUVLigero};
use ark_std::test_rng;

use crate::{
    aurora::*, naysayer::*, reader::read_constraint_system, TEST_DATA_PATH
};

trait FuzzablePolynomialCommitment<F: PrimeField + Absorb>: PolynomialCommitment<F, DensePolynomial<F>>
{
    fn fuzz_proof(proof: &mut Self::Proof);
}

impl<F: PrimeField + Absorb> FuzzablePolynomialCommitment<F> for TestUVLigero<F> {
    fn fuzz_proof(proof: &mut Self::Proof) {
        proof[0].opening.columns[0][0] += F::ONE;
    }
}

#[derive(PartialEq, Clone, Copy)]
enum AuroraDishonesty {
    // No dishonesty: same algorihm as AuroraR1CS::prove
    None,
    // Tamper with the i-th evaluation sent. Recall the order: 
    // [f_a(a) f_b(a), f_c(a), f_0(a), f_w(a), g_1(a), g_2(a)]
    Evaluation(usize),
    // Fuzz the PCS opening proof
    OpeningProof,
    // Tamper with the witness
    // This will break the zero test in most situations
    Witness,
    // Tamper with the polynomials before committing
    //      These break the zero test
    FA,
    FB,
    FC,
    F0,
    //      These break the unviariate sumcheck
    FW,
    G1,
    G2,
    // Tamper with f_a after absorbing the commitment.
    // This causes the prover and verifier sponges to desynchronise, which in
    // turn breaks the zero test at the squeezed point a
    FAPostAbsorb
}

fn dishonest_aurora_prove<F, PCS>(
    padded_aurora_r1cs: &AuroraR1CS<F>,
    instance: Vec<F>,
    witness: Vec<F>,
    ck: &PCS::CommitterKey,
    vk: &PCS::VerifierKey,
    sponge: &mut impl CryptographicSponge,
    dishonesty: AuroraDishonesty,
) -> AuroraProof<F, PCS>
where
    F: PrimeField + Absorb,
    PCS: FuzzablePolynomialCommitment<F>,
{

    let r1cs = padded_aurora_r1cs.r1cs();
    let matrices = r1cs.to_matrices().unwrap();

    // 0. Initialising sponge with public parameters
    absorb_public_parameters::<F, PCS>(&vk, &matrices, sponge);

    let ConstraintMatrices {
        a,
        b,
        c,
        num_constraints: n,
        num_instance_variables,
        ..
    } = matrices;

    // Resize the instance to the padded length
    let mut instance = instance;
    instance.resize(num_instance_variables, F::ZERO);

    // 1. Constructing committed polynomials
    // Following the notation of the paper
    let h = GeneralEvaluationDomain::new(n).unwrap();
    
    let mut witness = witness;
    
    if dishonesty == AuroraDishonesty::Witness {
        witness[0] += F::ONE;
    }

    let solution = [instance.clone(), witness.clone()].concat();

    // ======================== Computation of f_0 ========================

    // Note we can't compute f_a * f_b using an iFFT
    let mut f_a = matrix_polynomial(&a, &solution, &h);
    let mut f_b = matrix_polynomial(&b, &solution, &h);
    let mut f_c = matrix_polynomial(&c, &solution, &h);

    // Computing f_0 = (f_a * f_b - f_c) / v_h
    let mut f_0 = (&(&f_a * &f_b) - &f_c)
        .divide_by_vanishing_poly(h)
        .unwrap()
        .0;

    // ======================== Computation of f_w ========================

    let zero_padded_witness = [vec![F::ZERO; num_instance_variables], witness.clone()].concat();

    // Return the unique polynomial of degree < n that interpolates the given
    // values over h
    let h_interpolate = |evals: &Vec<F>| DensePolynomial::from_coefficients_vec(h.ifft(&evals));

    // Numerator
    let z_minus_v_star = h_interpolate(&zero_padded_witness);

    // TODO: Is there a more efficient way to compute this?
    // Denominator v_h_in = (x - 1) * (x - zeta^1) * ... * (x - zeta^(k-1))
    let v_h_in = h
        .elements()
        .take(num_instance_variables) // 1, zeta, ..., zeta^(k-1)
        .map(|c| DensePolynomial::from_coefficients_slice(&[-c, F::ONE])) // x - zeta^i
        .reduce(|acc, p| &acc * &p) // multiply together
        .unwrap();

    let mut f_w = &z_minus_v_star / &v_h_in;

    match dishonesty {
        AuroraDishonesty::FA => f_a.coeffs[0] += F::ONE,
        AuroraDishonesty::FB => f_b.coeffs[0] += F::ONE,
        AuroraDishonesty::FC => f_c.coeffs[0] += F::ONE,
        AuroraDishonesty::F0 => f_0.coeffs[0] += F::ONE,
        AuroraDishonesty::FW => f_w.coeffs[0] += F::ONE,
        _ => {}
    };

    // ======================== Computation of f_z ========================

    let f_z = h_interpolate(&solution);

    // ================== Committing to f_a, f_b, f_c, f_0, f_w ==================

    let labeled_polynomials_1 = label_polynomials(&[
        ("f_a", &f_a),
        ("f_b", &f_b),
        ("f_c", &f_c),
        ("f_0", &f_0),
        ("f_w", &f_w),
    ]);

    let (com_1, com_states_1) = PCS::commit(ck, &labeled_polynomials_1, None).unwrap();

    sponge.absorb(&com_1);

    // ======================== Computation of g_1, g_2 ========================

    if dishonesty == AuroraDishonesty::FAPostAbsorb  {
        f_a.coeffs[0] += F::ONE;
    }

    // Randomising polinomial through a squeezed challenge
    let r: F = sponge.squeeze_field_elements(1)[0];

    // Computing [1, r, r^2, ... r^(n-1)]
    let powers_of_r = powers(r, n);
    let p_r = h_interpolate(&powers_of_r);

    let q_ar = h_interpolate(&random_matrix_polynomial_evaluations(&a, &powers_of_r));
    let q_br = h_interpolate(&random_matrix_polynomial_evaluations(&b, &powers_of_r));
    let q_cr = h_interpolate(&random_matrix_polynomial_evaluations(&c, &powers_of_r));

    let r_pow_n = r * powers_of_r[n - 1];

    let u = (&(&p_r * &f_a) - &(&q_ar * &f_z))
        + &(&(&p_r * &f_b) - &(&q_br * &f_z)) * r_pow_n
        + &(&(&p_r * &f_c) - &(&q_cr * &f_z)) * (r_pow_n * r_pow_n);

    // We construct g_1 and g_2 such that u = v_h * g_1 + x * g_2

    // Auxiliary polynomials x and x^n - 1
    let x = DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[(1, F::ONE)]));

    let x_n_minus_1 = DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[(
        n - 1,
        F::ONE,
    )]));

    let dividend = &x_n_minus_1 * &u;

    let (quotient, mut g_2) = dividend.divide_by_vanishing_poly(h).unwrap();

    let mut g_1 = &(&quotient * &x) - &u;

    match dishonesty {
        AuroraDishonesty::G1 => g_1.coeffs[0] += F::ONE,
        AuroraDishonesty::G2 => g_2.coeffs[0] += F::ONE,
        _ => {}
    };

    //==================== Committing to g_1, g_2 ====================

    let labeled_polynomials_2 = label_polynomials(&[("g_1", &g_1), ("g_2", &g_2)]);

    let (com_2, com_states_2) = PCS::commit(ck, &labeled_polynomials_2, None).unwrap();

    sponge.absorb(&com_2);

    let coms = [com_1, com_2].concat();
    let com_states = [com_states_1, com_states_2].concat();
    let labeled_polynomials = [labeled_polynomials_1, labeled_polynomials_2].concat();

    //======================== PCS proof ========================

    let a_point: F = sponge.squeeze_field_elements(1)[0];

    let mut proof = PCS::open(
        ck,
        &labeled_polynomials,
        &coms,
        &a_point,
        sponge,
        &com_states,
        None,
    )
    .unwrap();

    if dishonesty == AuroraDishonesty::OpeningProof {
        PCS::fuzz_proof(&mut proof);
    }

    let mut evals: Vec<F> = labeled_polynomials
        .iter()
        .map(|lp| lp.evaluate(&a_point))
        .collect();

    if let AuroraDishonesty::Evaluation(i) = dishonesty {
        evals[i] += F::ONE;
    }
    
    AuroraProof {
        commitments: coms,
        proof,
        evals,
    }
}

#[test]
fn test_aurora_naysay() {
    let r1cs = read_constraint_system::<Fr>(
        &format!(TEST_DATA_PATH!(), "padding_test.r1cs"),
        &format!(TEST_DATA_PATH!(), "padding_test.wasm"),
    );

    // Instance: (1, a1, a2, b1, b2)
    let sol_c = Fr::from(42) * Fr::from(9 * 289).inverse().unwrap();
    let sol_a2c = Fr::from(9) * sol_c;
    let instance = vec![
        Fr::ONE,
        Fr::from(3),
        Fr::from(9),
        Fr::from(17),
        Fr::from(289),
    ];

    let witness = vec![sol_c, sol_a2c];

    let sponge: PoseidonSponge<Fr> = test_sponge();

    let (aurora_r1cs, ck, vk) =
        AuroraR1CS::setup::<TestUVLigero<Fr>>(r1cs, &mut test_rng()).unwrap();

    let test_aurora_naysay_with = |dishonesty: AuroraDishonesty, expected_naysayer_proof: Option<AuroraNaysayerProof<Fr, TestUVLigero<Fr>>>| {
        let aurora_proof = dishonest_aurora_prove(
            &aurora_r1cs,
            instance.clone(),
            witness.clone(),
            &ck,
            &vk,
            &mut sponge.clone(),
            dishonesty,
        );

        let naysayer_proof = aurora_naysay::<Fr, TestUVLigero<Fr>>(
            &aurora_r1cs,
            &vk,
            aurora_proof.clone(),
            instance.clone(),
            &mut sponge.clone(),
        )
        .unwrap();
            
        assert_eq!(naysayer_proof, expected_naysayer_proof);

        assert_eq!(
            aurora_verify_naysay::<Fr, TestUVLigero<Fr>>(
                &vk,
                &aurora_r1cs,
                &aurora_proof,
                &naysayer_proof.unwrap(),
                instance.clone(),
                &mut sponge.clone(),
            )
            .is_ok(),
            dishonesty == AuroraDishonesty::None
        );
    };
    
    let aurora_proof = aurora_r1cs
        .prove::<TestUVLigero<Fr>>(
            instance.clone(),
            witness.clone(),
            &ck,
            &vk,
            &mut sponge.clone(),
        )
        .unwrap();

    assert!(aurora_r1cs
        .verify::<TestUVLigero<Fr>>(&vk, instance.clone(), aurora_proof.clone(), &mut sponge.clone())
        .unwrap());

    let naysayer_proof = aurora_naysay::<Fr, TestUVLigero<Fr>>(
        &aurora_r1cs,
        &vk,
        aurora_proof.clone(),
        instance.clone(),   
        &mut sponge.clone(),
    ).unwrap();

    assert!(aurora_verify_naysay::<Fr, TestUVLigero<Fr>>(
        &vk,
        &aurora_r1cs,
        &aurora_proof,
        &naysayer_proof.unwrap(),
        instance.clone(),
        &mut sponge.clone(),
    ).unwrap());

    /***************** Case 1 *****************/
    // Honest proof verifies and is not naysaid
    test_aurora_naysay_with(AuroraDishonesty::None, None);
    
    /***************** Case 2 *****************/
    // // (i, x): Set the i-th value in the evaluation vector to x Recall the order
    // // [f_a(a) f_b(a), f_c(a), f_0(a), f_w(a), g_1(a), g_2(a)]
    // Evaluation(usize, F),
    // // Fuzz the PCS opening proof
    // OpeningProof,
    // // Tamper with the witness
    // // This will break the zero test in most situations
    // Witness,
    // // Tamper with the polynomials before committing
    // //      These break the zero test
    // FA,
    // FB,
    // FC,
    // F0,
    // //      These break the unviariate sumcheck
    // FW,
    // G1,
    // G2,
    // // Tamper with f_a after absorbing the commitment.
    // // This causes the prover and verifier sponges to desynchronise, which in
    // // turn breaks the zero test at the squeezed point a
    // FAPostAbsorb
}
