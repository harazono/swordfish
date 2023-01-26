use search_primer_and_probe::counting_bloomfilter_util::{L_LEN, M_LEN, R_LEN, HASHSET_SIZE};
use search_primer_and_probe::sequence_encoder_util::{LmrTuple};
use getopts::Options;
use std::env;
use std::fs::File;
use std::io::{Write, BufWriter, Read, BufReader};
use std::collections::HashSet;

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optmulti("i", "inputfiles", "set input file names", "NAME");
    opts.optopt("o", "outputfile", "set output file name", "NAME");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };

    if matches.opt_present("h") {
        println!("{}", opts.usage(&format!("Usage: {} [options]", args[0])));
        return;
    }

    let mut lmr_set: HashSet<LmrTuple> = HashSet::with_capacity(HASHSET_SIZE);


    let files = matches.opt_strs("i");
    for file in files {
        let f: File = File::open(&file).unwrap();
        let mut reader = BufReader::new(f);
        let mut buf: [u8; 24] = [0; 24];
        loop {
            match reader.read(&mut buf).unwrap() {
                0 => break,
                _n => {
                    let buf_l = &buf[0..8];
                    let buf_m = &buf[8..16];
                    let buf_r = &buf[16..24];
                    let mut l: u64 = 0;
                    let mut m: u64 = 0;
                    let mut r: u64 = 0;
                    for i in 0..8{
                        l <<= 4;
                        m <<= 4;
                        r <<= 4;
                        l += buf_l[7 - i] as u64;
                        m += buf_m[7 - i] as u64;
                        r += buf_r[7 - i] as u64;
                    }
                    let tmp_lmr_tuple = LmrTuple::new(l, m, r);
                    lmr_set.insert(tmp_lmr_tuple);
                }
            }
        }
    }
    let mut lmr_vec: Vec<LmrTuple> = lmr_set.into_iter().collect();
    lmr_vec.sort();


    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    }else{
        "lmrtuple_binary_decoder_out.bin".to_string()
    };

    let mut w = BufWriter::new(File::create(&output_file).unwrap());

    for each_lmr in &lmr_vec{
        w.write(&each_lmr.decode_as_single_vec()).unwrap();
        w.write(b"\n").unwrap();

    }
}
