use search_primer_and_probe::counting_bloomfilter_util::{L_LEN, M_LEN, R_LEN};
use search_primer_and_probe::sequence_encoder_util::{LmrTuple};
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



fn primer3_core_input_sequence(sequences: &Vec<LmrTuple>) -> Vec<String>{
    let mut str_vec: Vec<String> = Vec::new();
    let many_n = "N".to_string().repeat(50);
    for each_seq in sequences {
        let (l_vec, m_vec, r_vec) = each_seq.decode_as_triple_vec();
        let name = each_seq.id();
        let n_str: &str = std::str::from_utf8(&name).unwrap();
        let l_str: &str = std::str::from_utf8(&l_vec).unwrap();
        let m_str: &str = std::str::from_utf8(&m_vec).unwrap();
        let r_str: &str = std::str::from_utf8(&r_vec).unwrap();
    
        let sequence_with_internal_n = format!("{}{}{}{}{}", l_str, many_n, m_str, many_n, r_str);
        let primer3_fmt_str = format!("SEQUENCE_ID={}
SEQUENCE_TEMPLATE={}
PRIMER_TASK=pick_pcr_primers
PRIMER_OPT_SIZE=27
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=31
PRIMER_PRODUCT_SIZE_RANGE=101-200 201-301
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_OPT_TM=65.0
PRIMER_MAX_TM=70.0
=", n_str, sequence_with_internal_n);
        str_vec.push(primer3_fmt_str);
    }
    return str_vec;
}


fn execute_primer3(formatted_string: String) -> String{
    let process = match Command::new("primer3_core")
    .stdin(Stdio::piped())
    .stdout(Stdio::piped())
    .spawn() {
        Err(why) => panic!("couldn't spawn primer3: {}", why),
        Ok(process) => process,
    };

    match process.stdin.as_ref().unwrap().write_all(formatted_string.as_bytes()) {
        Err(why) => panic!("couldn't write to primer3_core stdin: {}", why),
        Ok(_) => eprintln!("sent pangram to primer3_core"),
    }

    let output = process.wait_with_output().expect("Failed to wait on child");
    let result = String::from_utf8(output.stdout).unwrap();
    return result;
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
    let mut buf: [u8; L_LEN + M_LEN + R_LEN] = [0; L_LEN + M_LEN + R_LEN];
    let mut candidates: Vec<LmrTuple> = Vec::new();

    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            _n => {
                let buf_l = &buf[0..L_LEN];
                let buf_m = &buf[L_LEN..M_LEN];
                let buf_r = &buf[(L_LEN + M_LEN)..(L_LEN + M_LEN + R_LEN)];

                let tmp_lmr_tuple = LmrTuple::new_from_bytes(buf_l, buf_m, buf_r);
                candidates.push(tmp_lmr_tuple);
            }
        }
    }

    let primer3_fmt_string: Vec<String> = primer3_core_input_sequence(&candidates);

    let mut chunks_of_input: Vec<String> = Vec::new();
    for _i in 0..thread_number{
        chunks_of_input.push(String::new());
    }
    for (index, string) in primer3_fmt_string.iter().enumerate(){
        chunks_of_input[index % thread_number] += string;
        chunks_of_input[index % thread_number] += "\n";
    }

    let arc_chunks_of_input: Arc<Vec<String>> = Arc::new(chunks_of_input);
    let final_result: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let mut children = Vec::new();
    for i in 0..thread_number{
        let chunks_of_input  = Arc::clone(&arc_chunks_of_input);
        let arc_final_result = Arc::clone(&final_result);
        children.push(
            thread::spawn(move|| {
                let primer3_results: String = execute_primer3((*chunks_of_input[i]).to_string());
                arc_final_result.lock().unwrap().push(primer3_results);
        })
        );
    }
    for child in children{
        let _ = child.join();
    }
    for i in final_result.lock().unwrap().iter(){
        println!("{}", i);
    }
}