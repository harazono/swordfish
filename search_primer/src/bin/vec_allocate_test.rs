fn main() {
    pub const BLOOMFILTER_TABLE_SIZE: usize = u32::MAX as usize + 1;
    eprintln!("Allocating Vec<[u32; BLOOMFILTER_TABLE_SIZE]> where BLOOMFILTER_TABLE_SIZE = {}", BLOOMFILTER_TABLE_SIZE);
    let mut ret_array: Vec<u32> = Vec::with_capacity(BLOOMFILTER_TABLE_SIZE);
    eprintln!("finish allocating");
}