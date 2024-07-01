pub mod aurora;
pub mod naysayer;
pub mod reader;

#[macro_export]
macro_rules! TEST_DATA_PATH {
    () => {
        concat!(env!("CARGO_MANIFEST_DIR"), "/test-data/{}",)
    };
}

pub use aurora::AuroraR1CS;
