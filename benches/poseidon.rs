#[macro_use]
extern crate criterion;

use criterion::{black_box, Criterion, ParameterizedBenchmark};
use gperftools::profiler::PROFILER;
use poseidon_rs::hash_bytes;
use rand::{thread_rng, Rng};

fn start_profile(stage: &str) {
    PROFILER
        .lock()
        .unwrap()
        .start(format!("./{}.profile", stage))
        .unwrap();
}
fn stop_profile() {
    PROFILER.lock().unwrap().stop().unwrap();
}
fn poseidon_benchmark(c: &mut Criterion) {
    let params = vec![32, 64, 4 * 32, 32 * 10];

    c.bench(
        "hash-poseidon",
        ParameterizedBenchmark::new(
            "non-circuit",
            |b, bytes| {
                let mut rng = thread_rng();
                let data: Vec<u8> = (0..*bytes).map(|_| rng.gen()).collect();

                // start_profile("poseidon");

                b.iter(|| black_box(hash_bytes(&data).unwrap()));
                // stop_profile();
            },
            params,
        ),
    );
}

criterion_group!(benches, poseidon_benchmark);
criterion_main!(benches);
