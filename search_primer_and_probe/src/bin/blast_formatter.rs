use lmr_tuple_count::counting_bloomfilter_util::L_LEN;
use lmr_tuple_count::counting_bloomfilter_util::M_LEN;
use lmr_tuple_count::counting_bloomfilter_util::R_LEN;
use lmr_tuple_count::sequence_encoder_util::LmrTuple;
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
    let (l_vec, m_vec, r_vec) = sequence.decode();
    let name = sequence.as_vec();
    let l_str: &str = std::str::from_utf8(&l_vec).unwrap();
    let m_str: &str = std::str::from_utf8(&m_vec).unwrap();
    let r_str: &str = std::str::from_utf8(&r_vec).unwrap();
    let fasta_fmt = format!(">{:X?}-L\n{}>{:X?}-M\n{}\n>{:X?}-R\n{}\n", name, l_str, name, m_str, name, r_str);
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
    let mut buf: [u8; L_LEN + M_LEN + R_LEN] = [0; L_LEN + M_LEN + R_LEN];
    let mut tmp_seq_as_u128: u128 = 0;


    let mut buf_l: [u8; L_LEN] = [0; L_LEN];
    let mut buf_m: [u8; M_LEN] = [0; M_LEN];
    let mut buf_r: [u8; R_LEN] = [0; R_LEN];



    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                //let buf = &buf[..n];
                let buf_l = &buf[0..L_LEN];
                let buf_m = &buf[L_LEN..M_LEN];
                let buf_r = &buf[(L_LEN + M_LEN)..R_LEN];
/*             
                for i in 0..16 {
                    tmp_seq_as_u128 <<= 8;
                    tmp_seq_as_u128 += u128::from(buf[i]);
                }
 */
                let tmp_lmr_tuple = LmrTuple::new_from_bytes(buf_l, buf_m, buf_r);
                writeln!(&mut w1, "{}", blast_formatter(&tmp_lmr_tuple)).unwrap();
                writeln!(&mut w2, "{:X?}", tmp_lmr_tuple.as_vec()).unwrap();
            }
        }
    }
}