#![cfg(feature = "bench")]
use ark_bn254::Fr;
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ff::Field;
use ark_poly_commit::{linear_codes::LinCodeParametersInfo, test_sponge, TestUVLigero};
use ark_std::test_rng;
use criterion::{criterion_main, BatchSize, Criterion};

use aurora::{
    aurora::{AuroraProof, AuroraVerifierKey},
    naysayer::{
        aurora_naysay, aurora_verify_naysay,
        tests::{dishonest_aurora_prove, AuroraDishonesty},
        AuroraNaysayerProof,
    },
    reader::read_constraint_system,
    AuroraR1CS, TEST_DATA_PATH,
};

fn setup_bench(
    dishonesty: AuroraDishonesty,
) -> (
    AuroraProof<Fr, TestUVLigero<Fr>>,
    Vec<Fr>,
    AuroraVerifierKey<Fr, TestUVLigero<Fr>>,
) {
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

    let (mut pk, mut vk) = AuroraR1CS::setup::<TestUVLigero<Fr>>(r1cs, &mut test_rng()).unwrap();

    pk.ck_large.set_well_formedness(false);
    pk.ck_small.set_well_formedness(false);
    vk.vk_large.set_well_formedness(false);
    vk.vk_small.set_well_formedness(false);

    let proof = dishonest_aurora_prove(
        &pk,
        instance.clone(),
        witness.clone(),
        (&vk.vk_large, &vk.vk_small),
        &mut sponge.clone(),
        dishonesty,
    );

    (proof, instance, vk)
}

fn bench_with_dishonesty(label: &str, dishonesty: AuroraDishonesty) {
    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Verify");

    group.bench_function(label, |b| {
        b.iter_batched(
            || setup_bench(dishonesty),
            |(proof, instance, vk)| {
                AuroraR1CS::verify::<TestUVLigero<Fr>>(
                    &vk,
                    &instance,
                    &proof,
                    &mut test_sponge::<Fr>(),
                )
            },
            BatchSize::SmallInput,
        );
    });

    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Naysay");

    group.bench_function(label, |b| {
        b.iter_batched(
            || setup_bench(dishonesty),
            |(proof, instance, vk)| aurora_naysay(&vk, &proof, &instance, &mut test_sponge::<Fr>()),
            BatchSize::SmallInput,
        );
    });

    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Verify Naysay");

    group.bench_function(label, |b| {
        b.iter_batched(
            || {
                let (proof, instance, vk) = setup_bench(dishonesty);
                let naysayer_proof =
                    aurora_naysay(&vk, &proof, &instance, &mut test_sponge::<Fr>());
                (proof, instance, vk, naysayer_proof.unwrap().unwrap())
            },
            |(proof, instance, vk, naysayer_proof)| {
                aurora_verify_naysay(
                    &vk,
                    &proof,
                    &naysayer_proof,
                    instance,
                    &mut test_sponge::<Fr>(),
                )
            },
            BatchSize::SmallInput,
        );
    });
}

fn bench_aurora() {
    bench_with_dishonesty("Zero Test", AuroraDishonesty::FA);
    bench_with_dishonesty("Univariate Sumcheck Test", AuroraDishonesty::FW);
    bench_with_dishonesty("PCS on f_a", AuroraDishonesty::FAPostAbsorb);
    bench_with_dishonesty("PCS on g_2", AuroraDishonesty::G2PostAbsorb);
}

criterion_main!(bench_aurora);
