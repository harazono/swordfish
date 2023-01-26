use search_primer_and_probe::counting_bloomfilter_util::{HASHSET_SIZE};
use search_primer_and_probe::sequence_encoder_util::{LmrTuple};
use getopts::Options;
use std::env;
use std::fs::File;
use std::io::{Write, BufWriter, Read, BufReader};
use std::collections::HashSet;
use std::mem;

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
                    let l_buf = &buf[0..8];
                    let m_buf = &buf[8..16];
                    let r_buf = &buf[16..24];
                    let mut l_u64: u64 = l_buf[7] as u64;
                    let mut m_u64: u64 = m_buf[7] as u64;
                    let mut r_u64: u64 = r_buf[7] as u64;
                    for i in 0..7{
                        l_u64 <<= 8;
                        m_u64 <<= 8;
                        r_u64 <<= 8;
                        l_u64 += l_buf[6 - i] as u64;
                        m_u64 += m_buf[6 - i] as u64;
                        r_u64 += r_buf[6 - i] as u64;
                    }

                    //let u64s: [u64; 3] = unsafe {mem::transmute(buf)};
                    let tmp_lmr_tuple = LmrTuple::new(l_u64, m_u64, r_u64);
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
