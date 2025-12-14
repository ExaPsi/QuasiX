use quasix_core::linalg::MetricSqrt;

fn assert_send<T: Send>() {}
fn assert_sync<T: Sync>() {}

fn main() {
    assert_send::<MetricSqrt>();
    assert_sync::<MetricSqrt>();
    println!("MetricSqrt is Send + Sync!");
}
