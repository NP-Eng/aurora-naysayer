#![cfg(feature = "bench")]
use std::{mem, process::Command};

use ark_bn254::Fr;
use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::Field;
use ark_poly_commit::{linear_codes::LinCodeParametersInfo, TestUVLigero};
use ark_std::test_rng;
use criterion::{criterion_main, BatchSize, Criterion};

use aurora::{
    aurora::{AuroraProof, AuroraVerifierKey},
    naysayer::{
        aurora_naysay, aurora_verify_naysay,
        tests::{dishonest_aurora_prove, AuroraDishonesty},
    },
    reader::read_constraint_system,
    AuroraR1CS, TEST_DATA_PATH,
};

mod sponge;
use sponge::{test_sponge, KeccakSponge};

fn setup_bench(
    dishonesty: AuroraDishonesty,
    num_squarings: usize,
) -> (
    AuroraProof<Fr, TestUVLigero<Fr>>,
    Vec<Fr>,
    AuroraVerifierKey<Fr, TestUVLigero<Fr>>,
) {
    let x = Fr::from(3);
    let mut y = x.clone();
    for _ in 0..num_squarings {
        y.square_in_place();
    }

    // call the python code to generate the circom file.
    // then call `circom` r1cs and wasm files

    // Define the Python script path and arguments
    let python_script = "../scripts/gen_solidity_repeated_squaring.py";
    let circom_directory = "circom";

    // Call the Python script to generate the circom and witness files
    let _ = Command::new("python3")
        .arg(python_script)
        .arg(num_squarings.to_string())
        .output()
        .expect("Failed to execute Python script");

    // Compile the circom file to get the R1CS and WASM files
    let circom_file = format!(
        "{}/repeated_squaring_{}.circom",
        circom_directory, num_squarings
    );

    Command::new("circom")
        .arg(&circom_file)
        .arg("--r1cs")
        .arg("--wasm")
        .arg("-o test-data/")
        .output()
        .expect("Failed to compile circom file");

    let r1cs = read_constraint_system::<Fr>(
        &format!(
            TEST_DATA_PATH!(),
            format!("repeated_squaring_{}.r1cs", num_squarings)
        ),
        &format!(
            TEST_DATA_PATH!(),
            format!(
                "/repeated_squaring_{}_js/repeated_squaring_{}.wasm",
                num_squarings, num_squarings
            )
        ),
        x,
        y,
    );

    let sponge: KeccakSponge = test_sponge();

    let (mut pk, mut vk) =
        AuroraR1CS::setup::<TestUVLigero<Fr>>(r1cs.clone(), &mut test_rng()).unwrap();

    pk.ck_large.set_well_formedness(false);
    pk.ck_small.set_well_formedness(false);
    vk.vk_large.set_well_formedness(false);
    vk.vk_small.set_well_formedness(false);

    let proof = dishonest_aurora_prove(
        &pk,
        r1cs.instance_assignment.clone(),
        r1cs.witness_assignment.clone(),
        (&vk.vk_large, &vk.vk_small),
        &mut sponge.clone(),
        dishonesty,
    );

    (proof, r1cs.instance_assignment.clone(), vk)
}

fn bench_with_dishonesty(label: &str, dishonesty: AuroraDishonesty, n: usize) {
    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Verify");

    group.bench_function(label, |b| {
        b.iter_batched(
            || setup_bench(dishonesty, n),
            |(proof, instance, vk)| {
                AuroraR1CS::verify::<TestUVLigero<Fr>>(
                    &vk,
                    instance.clone(),
                    proof,
                    &mut test_sponge(),
                )
            },
            BatchSize::SmallInput,
        );
    });

    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Naysay");

    group.bench_function(label, |b| {
        b.iter_batched(
            || setup_bench(dishonesty, n),
            |(proof, instance, vk)| aurora_naysay(&vk, proof, instance, &mut test_sponge()),
            BatchSize::SmallInput,
        );
    });

    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Verify Naysay");

    group.bench_function(label, |b| {
        b.iter_batched(
            || {
                let (proof, instance, vk) = setup_bench(dishonesty, n);
                let naysayer_proof =
                    aurora_naysay(&vk, proof.clone(), instance.clone(), &mut test_sponge());
                (proof, instance, vk, naysayer_proof.unwrap().unwrap())
            },
            |(proof, instance, vk, naysayer_proof)| {
                aurora_verify_naysay(&vk, &proof, &naysayer_proof, instance, &mut test_sponge())
            },
            BatchSize::SmallInput,
        );
    });
}

const NUM_SQUARINGS: usize = 10;

fn bench_aurora() {
    bench_with_dishonesty("Zero Test", AuroraDishonesty::FA, NUM_SQUARINGS);
    bench_with_dishonesty(
        "Univariate Sumcheck Test",
        AuroraDishonesty::FW,
        NUM_SQUARINGS,
    );
    bench_with_dishonesty("PCS on f_a", AuroraDishonesty::FAPostAbsorb, NUM_SQUARINGS);
    bench_with_dishonesty("PCS on g_2", AuroraDishonesty::G2PostAbsorb, NUM_SQUARINGS);
}

criterion_main!(bench_aurora);
