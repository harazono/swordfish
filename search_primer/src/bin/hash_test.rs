use std::collections::HashSet;

fn main() {
    let mut hash_from_u128: HashSet<u128> = HashSet::new();
    for i in 0..u128::MAX{
        hash_from_u128.insert(i);
        eprintln!("hash_count: {}", i);
    }
}