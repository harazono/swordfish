extern crate bio;
extern crate rdxsort;
extern crate getopts;
use getopts::Options;
use std::{env, process};
use std::fs;
use std::io::{BufWriter, Write};
use std::collections::HashSet;
use std::thread;
use std::sync::{Mutex, Arc};
use std::iter::zip;
use voracious_radix_sort::{RadixSort};
use search_primer::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use search_primer::counting_bloomfilter_util::{L_LEN, R_LEN, HASHSET_SIZE};
use search_primer::counting_bloomfilter_util::{build_counting_bloom_filter, number_of_high_occurence_lr_tuple};
use search_primer::sequence_encoder_util::{decode_u128_2_dna_seq};
use search_primer::sequence_encoder_util::DnaSequence;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use std::fs::File;
use crate::bio::io::fasta::FastaRead;



fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "output", "set output file name", "NAME");
    opts.optopt("t", "thread", "number of threads to use for radix sort. default value is 8.", "THREAD");
    opts.optopt("a", "threshold", "threshold for hyper log counter. default value is 8.", "THRESHOLD");
    opts.optflag("b", "binary", "outputs binary file");
    opts.optflag("r", "only-num", "outputs only total number of k-mer");
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

    let threads: usize = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    }else{
        8
    };

    let threshold:u32 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u32>().unwrap()
    }else{
        8
    };

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    }else{
        format!("{:?}_threshold{}_threads{}.out", input_file, threshold, threads)
    };
    eprintln!("input  file: {:?}",  input_file);
    eprintln!("loading {:?} done", input_file);


    let file = File::open(&input_file).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", input_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }

/*
ここにマルチスレッド処理を書く
*/

    let chunk_size: usize = sequences.len() / (threads - 1);
    let sequences_ref = &sequences;
    let mut cbf_oyadama: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];

    thread::scope(|scope|{
        let mut children_1 = Vec::new();
        for i in 1..threads {
            children_1.push(
                scope.spawn(move || 
                    {
                        let start_idx: usize = (i - 1) * chunk_size;
                        let end_idx: usize;
                        if i != threads - 1{
                            end_idx = i * chunk_size;
                        }else{
                            end_idx = sequences_ref.len() - 1;
                        }
                        eprintln!("start calling build_counting_bloom_filter[{}]", i);
                        let cbf: Vec<u32> = build_counting_bloom_filter(sequences_ref, start_idx, end_idx, i);
                        eprintln!("finish calling build_counting_bloom_filter[{}]", i);
                        cbf
                    }
                )
            )
        }
        for child in children_1{
            let cbf = child.join().unwrap();
            zip(cbf_oyadama.iter_mut(), cbf).for_each(|(x, y)| *x = x.checked_add(y).unwrap_or(u32::MAX));
        }
    });
    let h_cbf_h_oyadama: Arc<Mutex<HashSet<u128>>> = Arc::new(Mutex::new(HashSet::with_capacity(HASHSET_SIZE)));
    let cbf_oyadama_ref = &cbf_oyadama;
    let h_cbf_h_oyadama_ref = &h_cbf_h_oyadama;


    thread::scope(|scope|{
        let mut children_2 = Vec::new();
        for i in 1..threads {
            children_2.push(
                scope.spawn(move ||
                    {
                        let start_idx: usize = (i - 1) * chunk_size;
                        let end_idx: usize;
                        if i != threads - 1{
                            end_idx = i * chunk_size;
                        }else{
                            end_idx = sequences_ref.len() - 1;
                        }
                        eprintln!("thread [{}]: start calling number_of_high_occurence_lr_tuple", i);
                        let h_cbf_h: HashSet<u128> = number_of_high_occurence_lr_tuple(cbf_oyadama_ref, sequences_ref, start_idx, end_idx, threshold, i);
                        h_cbf_h_oyadama_ref.lock().unwrap().extend(&h_cbf_h);
                        eprintln!("thread [{}]: finish calling number_of_high_occurence_lr_tuple", i);
                    }
                )
            )
        }
        for child in children_2{
            let _ = child.join();
        }
    });


    let mut high_occurence_lr_tuple: Vec<u128> = Vec::from_iter(h_cbf_h_oyadama.lock().unwrap().clone());
    high_occurence_lr_tuple.sort();
 /*
ここまで
*/
    let mut w = BufWriter::new(fs::File::create(&output_file).unwrap());
    let mut previous_lr_tuple: u128 = 0;
    let mut cnt = 0;
    let mut buf_array: [u8; 16] = [0; 16];
    let mut buf_num: u128;

    if matches.opt_present("r") {
        eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
        for each_lr_tuple in &high_occurence_lr_tuple{
            if previous_lr_tuple != *each_lr_tuple{
                cnt += 1;
            }
            previous_lr_tuple = *each_lr_tuple;
        }
        writeln!(&mut w, "k-mer count: {}\tthreshold: {}\tinput file {:?}", cnt, threshold, &input_file).unwrap();
    }
    if !matches.opt_present("r") && matches.opt_present("b"){
        eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
        for each_lr_tuple in &high_occurence_lr_tuple{
            if previous_lr_tuple != *each_lr_tuple{
                cnt += 1;
                buf_num = *each_lr_tuple;
                for i in 0..16{
                    buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
                    buf_num >>= 8;
                }
                w.write(&buf_array).unwrap();
            }
            previous_lr_tuple = *each_lr_tuple;
        }
    }
    if !matches.opt_present("r") && !matches.opt_present("b"){
        eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
        for each_lr_tuple in &high_occurence_lr_tuple{
            if previous_lr_tuple != *each_lr_tuple{
                cnt += 1;
                writeln!(&mut w, "{:?}", String::from_utf8(decode_u128_2_dna_seq(&each_lr_tuple, 54)).unwrap()).unwrap();
            }
            previous_lr_tuple = *each_lr_tuple;
        }
    }

    eprintln!("finish writing to output file: {:?}", &output_file);
    eprint!("L: {}\tR: {}\tthreshold:{}({}x63)\tcardinarity: {}\t", L_LEN, R_LEN, threshold, threshold / 63, cnt);
    eprintln!("threads: {}\tinput file {:?}", threads, &input_file);

}