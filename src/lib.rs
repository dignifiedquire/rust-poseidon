use ff::{Field, PrimeField, PrimeFieldRepr};
use lazy_static::lazy_static;
use paired::bls12_381::{Fr, FrRepr};

lazy_static! {
    static ref CONSTANTS: Constants = generate_constants();
}

const SEED: &str = "poseidon";
const NROUNDSF: usize = 8;
const NROUNDSP: usize = 57;
const T: usize = 6;

#[derive(Debug)]
pub struct Constants {
    c: Vec<Fr>,
    m: Vec<Vec<Fr>>,
}

/// Hash bytes using blake2b and convert the result into an `Fr` field element.
fn blake2b(input: &[u8]) -> Fr {
    let mut h = blake2b_simd::Params::new();
    h.hash_length(32);
    let mut res = h.hash(input).as_ref().to_vec();

    // strip last two bits, to ensure result is in Fr.
    // TODO: proper mod r
    res[31] &= 0b0011_1111;
    assert_eq!(res.len(), 32);

    let mut repr = FrRepr::default();
    repr.read_le(std::io::Cursor::new(res)).unwrap();
    Fr::from_repr(repr).unwrap()
}

fn generate_constants() -> Constants {
    let c = get_pseudo_random(format!("{}{}", SEED, "_constants"), NROUNDSF + NROUNDSP);
    let m = get_mds();

    Constants { c, m }
}

fn get_pseudo_random(seed: String, n: usize) -> Vec<Fr> {
    let mut h = blake2b(seed.as_bytes());

    let mut res = Vec::with_capacity(n);
    for _ in 0..n {
        let mut bytes = vec![0u8; 32];
        h.into_repr().write_le(&mut bytes).unwrap();
        res.push(h);

        h = blake2b(&bytes);
    }
    res
}

fn nonce_to_string(n: usize) -> String {
    let mut r = format!("{}", n);
    while r.len() < 4 {
        r = format!("0{}", r);
    }
    r
}

fn get_mds() -> Vec<Vec<Fr>> {
    let mut nonce = 0;
    let mut cauchy_matrix =
        get_pseudo_random(format!("{}_matrix_{}", SEED, nonce_to_string(nonce)), T * 2);
    while !are_all_different(&cauchy_matrix) {
        nonce = nonce + 1;
        cauchy_matrix =
            get_pseudo_random(format!("{}_matrix_{}", SEED, nonce_to_string(nonce)), T * 2);
    }
    let mut m = Vec::with_capacity(T);
    for i in 0..T {
        let mut mi = Vec::with_capacity(T);
        for j in 0..T {
            let mut a = cauchy_matrix[i].clone();
            a.sub_assign(&cauchy_matrix[T + j]);
            let res = a.inverse().unwrap();

            mi.push(res);
        }
        m.push(mi);
    }
    m
}

fn are_all_different(v: &[Fr]) -> bool {
    for (i, vi) in v.iter().enumerate() {
        if vi.is_zero() {
            return false;
        }
        for j in i + 1..v.len() {
            if vi == &v[j] {
                return false;
            }
        }
    }
    true
}

#[inline]
fn ark(state: &mut [Fr], c: &Fr) {
    for s in state.iter_mut() {
        s.add_assign(c);
    }
}

#[inline]
fn cubic(a: &Fr) -> Fr {
    let mut res = a.clone();
    res.square(); // ^2
    res.square(); // ^4
    res.mul_assign(a); // ^5

    res
}

fn sbox(state: &mut [Fr], i: usize) {
    if i < NROUNDSF / 2 || i >= NROUNDSF / 2 + NROUNDSP {
        for j in 0..T {
            state[j] = cubic(&state[j]);
        }
    } else {
        state[0] = cubic(&state[0]);
    }
}

fn mix(state: &mut [Fr], old_state: &[Fr], m: &[Vec<Fr>]) {
    for i in 0..state.len() {
        state[i] = Fr::zero();
        for j in 0..state.len() {
            let mut v = m[i][j].clone();
            v.mul_assign(&old_state[j]);

            state[i].add_assign(&v);
        }
    }
}

fn poseidon_hash(state: &mut [Fr]) -> Result<Fr, String> {
    assert_eq!(state.len(), T, "Wrong inputs length");

    let mut old_state = vec![Fr::zero(); state.len()];
    for i in 0..(NROUNDSF + NROUNDSP) {
        ark(state, &CONSTANTS.c[i]);
        sbox(state, i);
        old_state.copy_from_slice(state);
        mix(state, &old_state, &CONSTANTS.m);
    }

    Ok(state[0].clone())
}

pub fn hash(inp: Vec<Fr>) -> Result<Fr, String> {
    let mut r = Fr::zero();
    let mut five_elems = Vec::with_capacity(T);

    for i in (0..inp.len()).step_by(5) {
        for j in 0..5 {
            if i + j < inp.len() {
                five_elems.push(inp[i + j].clone());
            } else {
                five_elems.push(Fr::zero());
            }
        }
        five_elems.resize(T, Fr::zero());
        let ph = poseidon_hash(&mut five_elems);
        match ph {
            Result::Err(err) => return Err(err.to_string()),
            Result::Ok(res) => r.add_assign(&res),
        }
        five_elems.clear();
    }
    Ok(r)
}

pub fn hash_bytes(b: &[u8]) -> Result<Fr, String> {
    let n = 31;

    let ints = b
        .chunks(n)
        .map(|chunk| {
            let mut chunk = chunk.to_vec();
            chunk.resize(32, 0);

            let mut repr = FrRepr::default();
            repr.read_le(std::io::Cursor::new(chunk)).unwrap();
            Fr::from_repr(repr).unwrap()
        })
        .collect();

    hash(ints)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blake2_version() {
        let h = blake2b(b"poseidon_constants");
        assert_eq!(
            h,
            Fr::from_str("e57ba154fb2c47811dc1a2369b27e25a44915b4e4ece4eb8ec74850cb78e01b1")
                .unwrap()
        );
    }

    #[test]
    fn test_poseidon_hash() {
        let b1 = Fr::from_str("1").unwrap();
        let b2 = Fr::from_str("2").unwrap();
        let mut big_arr: Vec<Fr> = Vec::new();
        big_arr.push(b1.clone());
        big_arr.push(b2.clone());
        big_arr.resize(T, Fr::zero());
        let h = poseidon_hash(&mut big_arr).unwrap();
        assert_eq!(
            h,
            Fr::from_str(
                "12242166908188651009877250812424843524687801523336557272219921456462821518061"
            )
            .unwrap()
        );

        let b3 = Fr::from_str("3").unwrap();
        let b4 = Fr::from_str("4").unwrap();
        let mut big_arr34: Vec<Fr> = Vec::new();
        big_arr34.push(b3.clone());
        big_arr34.push(b4.clone());
        big_arr34.resize(T, Fr::zero());
        let h34 = poseidon_hash(&mut big_arr34).unwrap();
        assert_eq!(
            h34,
            Fr::from_str(
                "17185195740979599334254027721507328033796809509313949281114643312710535000993"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_hash_bytes() {
        let msg = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.";
        let h = hash_bytes(msg.as_bytes()).unwrap();
        assert_eq!(
            h,
            Fr::from_str(
                "11821124228916291136371255062457365369197326845706357273715164664419275913793"
            )
            .unwrap()
        );

        let msg2 = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum. Lorem ipsum dolor sit amet.";
        let h2 = hash_bytes(msg2.as_bytes()).unwrap();
        assert_eq!(
            h2,
            Fr::from_str(
                "10747013384255785702102976082726575658403084163954725275481577373644732938016"
            )
            .unwrap()
        );
    }
}
