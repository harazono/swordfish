extern crate search_primer;
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
use search_primer::sequence_encoder_util::{decode_u128_l, decode_u128_r};
use std::time::Instant;


fn primer3_core_input_sequence(sequences: &Vec<u128>, library_file_name: &Option<String>) -> Vec<String>{
    let mut str_vec: Vec<String> = Vec::new();
    let many_n = "N".to_string().repeat(50);
    eprintln!("primer3_core_input_sequence: sequense length...{}", sequences.len());

    for each_seq in sequences {
        let l_u8_array = decode_u128_l(&each_seq);
        let r_u8_array = decode_u128_r(&each_seq);
        let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
        let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
        let sequence_with_internal_n = format!("{}{}{}", l_str, many_n, r_str);
        let mut primer3_fmt_str = format!("SEQUENCE_ID={:0x}
SEQUENCE_TEMPLATE={}
PRIMER_TASK=pick_pcr_primers
PRIMER_OPT_SIZE=27
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=31
PRIMER_PRODUCT_SIZE_RANGE=101-200 201-301
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_OPT_TM=65.0
PRIMER_MAX_TM=70.0", each_seq, sequence_with_internal_n);
        // Check if library_file_name is Some or None
        match library_file_name.as_ref() {
            Some(file_name) => {
                // If Some, append the file name to the string
                primer3_fmt_str.push_str(&format!("\nPRIMER_MISPRIMING_LIBRARY={}\n=", file_name));
            }
            None => {
                primer3_fmt_str.push_str("\n=");
            }
        }
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
    //println!("{}", result);
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
    opts.optopt("l", "library", "library file name which will be used for PRIMER_MISPRIMING_LIBRARY", "LIBRARY");


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

    let library_file_name: Option<String> = matches.opt_str("l");



    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };
    eprintln!("start  loading {:?}", &input_file);
    let f: File = File::open(&input_file).unwrap();
    eprintln!("finish loading {:?}", &input_file);
    let mut reader = BufReader::new(f);
    let mut buffer: [u8; 16] = [0; 16];
    let mut tmp_seq_as_u128: u128;
    let mut candidates: Vec<u128> = Vec::new();
    loop {
        let result = reader.by_ref().take(16).read_exact(&mut buffer);
        match result {
            Ok(_val) => {},
            Err(_err) => break,
        }
        tmp_seq_as_u128 = u128::from_le_bytes(buffer);
        candidates.push(tmp_seq_as_u128);
    }

    eprintln!("start formatting string");
    let primer3_fmt_string: Vec<String> = primer3_core_input_sequence(&candidates, &library_file_name);

    let mut chunks_of_input: Vec<Vec<String>> = Vec::new();
    for _i in 0..thread_number{
        chunks_of_input.push(Vec::new());
    }
    for (index, string) in primer3_fmt_string.iter().enumerate(){
        chunks_of_input[index % thread_number].push(string.clone());
    }

    let arc_chunks_of_input: Arc<Vec<Vec<Vec<String>>>> = Arc::new(chunks_of_input.into_iter().map(|v| v.chunks(5000).map(|c| c.to_vec()).collect()).collect());
    let final_result: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let mut children = Vec::new();

    for i in 0..thread_number {
        let arc_final_result = Arc::clone(&final_result);
        let chunk = arc_chunks_of_input[i].clone(); // Clone the chunk for the thread
        let thread_id = i;
        let total_chunks = arc_chunks_of_input.len();
        children.push(
            thread::spawn(move|| {
                let start_time = Instant::now();
                eprintln!("start calling primer3_core in iteration {} (Thread ID: {}, Chunk: {} of {})", i + 1, thread_id, i + 1, total_chunks);
                let mut primer3_results = Vec::new();
                // Process each element in the chunk
                for chunk in chunk {
                    for (j, input) in chunk.iter().enumerate() {
                        primer3_results.push(execute_primer3(input.clone()));
                        let progress = ((j + 1) as f64 / chunk.len() as f64) * 100.0;
                        eprintln!("Thread ID: {}, Progress: {:.2}%", thread_id, progress);
                    }
                }
                {
                    let mut final_result = arc_final_result.lock().unwrap();
                    for result in primer3_results {
                        final_result.push(result);
                    }
                }
                let duration = start_time.elapsed();
                eprintln!("finish calling primer3_core in iteration {} (Thread ID: {}). Time elapsed: {:.2} seconds", i + 1, thread_id, duration.as_secs_f64());
            })
        );
    }

    for child in children {
        // Wait for the thread to finish. Returns a result.
        let _ = child.join();
    }

    for i in final_result.lock().unwrap().iter(){
        println!("{}", i);
    }
}