extern crate kmer_count;
extern crate getopts;
use std::{env, process};
use std::fs::File;
use std::io::{Read, BufReader};
use std::io::prelude::*;
use std::process::{Command, Stdio};
use std::thread;
use std::sync::Arc;
use std::sync::Mutex;
use getopts::Options;
use kmer_count::sequence_encoder_util::LmrTuple;

fn sequence_dump(sequence: LmrTuple) -> String{
    let (l_vec, m_vec, r_vec) = LmrTuple.decode();
    let l_str: &str = std::str::from_utf8(&l_vec).unwrap();
    let m_str: &str = std::str::from_utf8(&m_vec).unwrap();
    let r_str: &str = std::str::from_utf8(&r_vec).unwrap();
    let ret_fmt = format!("{} {} {}\n", l_str, m_str, r_str);
    return ret_fmt;
}



fn sequences_dump(sequences: &Vec<u128>) -> Vec<String>{
    let mut str_vec: Vec<String> = Vec::new();
    for each_seq in sequences {
        str_vec.push(sequence_dump(each_seq));
    }
    return str_vec;
}




fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE ", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main(){
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("t", "thread", "number of thread to use for radix sort. default value is 8.", "THREAD");


    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let thread_number: usize = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    }else{
        4
    };

    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };
    let f: File = File::open(&input_file).unwrap();
    let mut reader = BufReader::new(f);
    let mut buf: [u8; 16] = [0; 16];
    let mut tmp_seq_as_u128: u128 = 0;
    let mut candidates: Vec<u128> = Vec::new();

    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                let buf = &buf[..n];
                for i in 0..16 {
                    tmp_seq_as_u128 <<= 8;
                    tmp_seq_as_u128 += u128::from(buf[i]);
                }
                println!("{:?}", sequence_dump(&tmp_seq_as_u128));
            }
        }
    }
/* 
    let sequences_dump_string: Vec<String> = sequences_dump(&candidates);
    for candidate in sequences_dump_string{
        println!("{}", candidate);
    }
*/
}