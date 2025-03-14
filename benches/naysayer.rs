#![cfg(feature = "bench")]
use std::path::Path;

use ark_bn254::Fr;
use ark_circom::{CircomBuilder, CircomConfig};
use ark_ff::{Field, PrimeField};
use ark_poly_commit::{linear_codes::LinCodeParametersInfo, TestUVLigero};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};
use ark_std::test_rng;
use criterion::{criterion_main, Criterion};

use aurora::{
    aurora::{AuroraProof, AuroraVerifierKey},
    naysayer::{
        aurora_naysay, aurora_verify_naysay,
        tests::{dishonest_aurora_prove, AuroraDishonesty},
    },
    AuroraR1CS, TEST_DATA_PATH,
};

mod sponge;
use num_bigint::{BigInt, BigUint};
use sponge::{test_sponge, KeccakSponge};

pub fn read_constraint_system_and_populate<F: PrimeField>(
    r1cs_file: impl AsRef<Path>,
    wasm_file: impl AsRef<Path>,
    x: F,
    y: F,
) -> ConstraintSystem<F> {
    let cfg = CircomConfig::<F>::new(wasm_file, r1cs_file).unwrap();

    let mut builder = CircomBuilder::new(cfg);
    // this is an ugly hack for now. Ideally, `push_input already accepts `F` elements.
    builder.push_input("x", Into::<BigInt>::into(Into::<BigUint>::into(x)));
    builder.push_input("y", Into::<BigInt>::into(Into::<BigUint>::into(y)));

    let circom = builder.build().unwrap();

    let cs = ConstraintSystem::<F>::new_ref();
    circom.generate_constraints(cs.clone()).unwrap();
    cs.into_inner().unwrap()
}

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

    let r1cs = read_constraint_system_and_populate::<Fr>(
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

    let (proof, instance, vk) = setup_bench(dishonesty, n);

    group.bench_function(label, |b| {
        b.iter(|| {
            AuroraR1CS::verify::<TestUVLigero<Fr>>(&vk, &instance, &proof, &mut test_sponge())
        });
    });

    let mut c = Criterion::default().sample_size(10);
    let mut group = c.benchmark_group("Naysay");

    group.bench_function(label, |b| {
        b.iter(|| aurora_naysay(&vk, &proof, &instance, &mut test_sponge()));
    });

    // shouldn't panic if we're benching an honest case
    match dishonesty {
        AuroraDishonesty::None => return,
        _ => {
            let mut c = Criterion::default().sample_size(10);
            let mut group = c.benchmark_group("Verify Naysay");
            let naysayer_proof = aurora_naysay(&vk, &proof, &instance, &mut test_sponge())
                .unwrap()
                .unwrap();
            group.bench_function(label, |b| {
                b.iter(|| {
                    aurora_verify_naysay(
                        &vk,
                        &proof,
                        &naysayer_proof,
                        &instance,
                        &mut test_sponge(),
                    )
                });
            });
        }
    }
}

const NUM_SQUARINGS: usize = 1 << 10;

fn bench_aurora() {
    bench_with_dishonesty("Honest", AuroraDishonesty::None, NUM_SQUARINGS);
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
