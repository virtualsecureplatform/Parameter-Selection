extern crate concrete_npe;
extern crate concrete_commons;
use concrete_npe as npe;
pub use concrete_commons::dispersion::{DispersionParameter, LogStandardDev, Variance};
pub use concrete_commons::parameters::{
    DecompositionBaseLog, DecompositionLevelCount, GlweDimension, LweDimension, LweSize,
    PlaintextCount, PolynomialSize,
};
pub use concrete_commons::key_kinds::BinaryKeyKind;

fn main(){
	let polynomial_size = PolynomialSize(1024);
	let rlwe_dimension = GlweDimension(1);
	let lwe_dimension = LweDimension(630);
	let level = DecompositionLevelCount(3);
	let base_log = DecompositionBaseLog(7);
	let std = LogStandardDev::from_log_standard_dev(-25.);
	let output_variance = npe::estimate_pbs_noise::<u32, Variance, BinaryKeyKind>(
            lwe_dimension,
            polynomial_size,
            rlwe_dimension,
            base_log,
            level,
            Variance(f64::powi(std.get_standard_dev(), 2)),
        );
	println!("varinace: {}",output_variance.get_variance()/f64::powi(2.0,64));
}