use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use tiny_keccak::{Hasher, Keccak};

#[derive(Clone)]
pub(crate) struct KeccakSponge {
    pub(crate) state: Vec<u8>,
}

impl CryptographicSponge for KeccakSponge {
    type Config = ();

    fn new(_params: &Self::Config) -> Self {
        KeccakSponge { state: vec![] }
    }

    fn absorb(&mut self, input: &impl Absorb) {
        let mut input_bytes = vec![];
        input.to_sponge_bytes(&mut input_bytes);
        self.state.extend_from_slice(&input_bytes);
    }

    fn squeeze_bytes(&mut self, num_bytes: usize) -> Vec<u8> {
        let mut keccak = Keccak::v256();
        let mut output = vec![0u8; num_bytes];
        keccak.update(&self.state);
        keccak.finalize(&mut output);
        self.state = output.clone();
        output
    }

    fn squeeze_bits(&mut self, num_bits: usize) -> Vec<bool> {
        let num_bytes = (num_bits + 7) / 8;
        let tmp = self.squeeze_bytes(num_bytes);
        let dest = tmp
            .iter()
            .flat_map(|byte| (0..8u32).rev().map(move |i| (byte >> i) & 1 == 1))
            .collect::<Vec<_>>();
        dest[..num_bits].to_vec()
    }
}

pub(crate) fn test_sponge() -> KeccakSponge {
    KeccakSponge::new(&())
}
