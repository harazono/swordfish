use std::collections::HashSet;
use std::sync::Mutex;

fn main(){
    eprintln!("start allocation");
    let h_cbf_h_oyadama: Mutex<HashSet<u128>> = Mutex::new(HashSet::with_capacity((u32::MAX >> 1) as usize));
    eprintln!("finish allocation");

}