use search_primer_and_probe::counting_bloomfilter_util::{L_LEN, M_LEN, R_LEN};
use search_primer_and_probe::sequence_encoder_util::{LmrTuple};
extern crate getopts;
use std::{env, process};
use std::io::{Write, BufWriter};
use std::io::{Read,  BufReader};
use std::fs;
use std::fs::File;
use getopts::Options;


fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn blast_formatter(sequence: &LmrTuple) -> String{
    let (l_vec, m_vec, r_vec) = sequence.decode_as_triple_vec();
    let name = sequence.id();
    let n_str: &str = std::str::from_utf8(&name).unwrap();
    let l_str: &str = std::str::from_utf8(&l_vec).unwrap();
    let m_str: &str = std::str::from_utf8(&m_vec).unwrap();
    let r_str: &str = std::str::from_utf8(&r_vec).unwrap();
    let fasta_fmt = format!(">{}-L\n{}>{}-M\n{}\n>{}-R\n{}\n", n_str, l_str, n_str, m_str, n_str, r_str);
    return fasta_fmt;
}



fn main(){
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("a", "output1", "sequence.fasta", "FILENAME");
    opts.optopt("b", "output2", "namelist.txt", "FILENAME");

    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };

    let output_file_1 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap()
    }else{
        "out_1.fasta".to_string()
    };
    let output_file_2 = if matches.opt_present("b") {
        matches.opt_str("b").unwrap()
    }else{
        "out_2.txt".to_string()
    };

    let mut w1 = BufWriter::new(fs::File::create(&output_file_1).unwrap());
    let mut w2 = BufWriter::new(fs::File::create(&output_file_2).unwrap());


    let f: File = File::open(&input_file).unwrap();
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
                for i in 0..4{
                    l <<= 2;
                    m <<= 2;
                    r <<= 2;
                    l += buf_l[7 - i] as u64;
                    m += buf_m[7 - i] as u64;
                    r += buf_r[7 - i] as u64;
                }
                let tmp_lmr_tuple = LmrTuple::new(l, m, r);
            writeln!(&mut w1, "{}", blast_formatter(&tmp_lmr_tuple)).unwrap();
                writeln!(&mut w2, "{:?}", tmp_lmr_tuple.id()).unwrap();
            }
        }
    }
}