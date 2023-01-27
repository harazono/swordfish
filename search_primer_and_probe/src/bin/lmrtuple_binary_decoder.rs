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

    let mut lmr_set: HashSet<LmrTuple> = HashSet::with_capacity(10000);
    //let mut lmr_vec: Vec<LmrTuple> = Vec::with_capacity(10000);

    let files = matches.opt_strs("i");
    for file in files {
        let f: File = File::open(&file).unwrap();
        let mut reader = BufReader::new(f);
        let mut buffer = [0u8; 24];

        loop {
            let result = reader.by_ref().take(24).read_exact(&mut buffer);
            match result {
                Ok(val) => {},
                Err(err) => break,
            }
            /*
            if bytes_read == 0 {
                break;
            }
            if bytes_read != 24 {
                eprintln!("never reached!");
            }
            */
            let mut l: u64 = buffer[0] as u64;
            let mut m: u64 = buffer[8] as u64;
            let mut r: u64 = buffer[16] as u64;
            for i in 1..8{
                l <<= 8;
                m <<= 8;
                r <<= 8;
                l += buffer[i + 0] as u64;
                m += buffer[i + 8] as u64;
                r += buffer[i + 16] as u64;
            }
            //let u64s: [u64; 3] = unsafe {mem::transmute(buffer)};
            let tmp_lmr_tuple = LmrTuple::new(l, m, r);
            // /*
            eprint!("{:?} ", String::from_utf8(tmp_lmr_tuple.decode_as_single_vec()).unwrap());
            eprint!("{:08b} ", buffer[0]);
            eprint!("{:064b} ", l);
            eprint!("{:064b} ", m);
            eprint!("{:064b} ", r);
            eprintln!();
            // */
            lmr_set.insert(tmp_lmr_tuple);
            //lmr_vec.push(tmp_lmr_tuple);
        }
    }
    //Vecに突っ込むと4こ要素が増える謎

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
