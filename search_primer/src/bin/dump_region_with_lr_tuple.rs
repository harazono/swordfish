extern crate bio;
extern crate rdxsort;
extern crate getopts;
use getopts::Options;
use std::{env, process};
use std::thread;
use std::fs::File;
use std::io::{Write, Read};
use std::sync::Arc;
use search_primer::counting_bloomfilter_util::aggregate_length_between_lr_tuple;
use search_primer::sequence_encoder_util::DnaSequence;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use crate::bio::io::fasta::FastaRead;
use std::io::BufReader;
use search_primer::sequence_encoder_util::{decode_u128_l, decode_u128_r};

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("i", "input", "set input file name (u128 binary)", "NAME");
    opts.optopt("o", "output", "set output file name", "NAME");
    opts.optopt("t", "thread", "number of threads to use for radix sort. default value is 8.", "THREAD");
    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let lr_tuple_filename = if matches.opt_present("i") {
        matches.opt_str("i").unwrap()
    }else{
        print_usage(&program, &opts);
        return;
    };


    let ngsread_input_file = if !matches.free.is_empty() {
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


    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    }else{
        "region.out".to_string()
    };



    eprintln!("Input  file name: {:?}", ngsread_input_file);
    eprintln!("Output file name: {:?}", output_file);
    eprintln!("Number of threads: {:?}", threads);


    let mut lr_tuple: Vec<(Vec<u8>, DnaSequence, DnaSequence)> = Vec::new();

    let f: File = File::open(&lr_tuple_filename).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = [0u8; 16];
    loop {
        let result = reader.by_ref().take(16).read_exact(&mut buffer);
        match result {
            Ok(_val) => {},
            Err(_err) => break,
        }
        let tmp_val: u128 = u128::from_be_bytes(buffer);
        let l            = &decode_u128_l(&tmp_val).to_vec();
        let r            = &decode_u128_r(&tmp_val).to_vec();
        let left_primer  = DnaSequence::new(l);
        let right_primer = DnaSequence::new(r);
        let hex_string   = format!("0x{:016x}", tmp_val);
        lr_tuple.push((hex_string.into(), left_primer, right_primer));
    }

    eprintln!("Number of lr_tuple: {:?}", &lr_tuple.len());

    for each_lr_tuple in &lr_tuple {
        eprintln!("lr_tuple: {}", String::from_utf8(each_lr_tuple.1.decode(0, each_lr_tuple.1.len())).unwrap());
        eprintln!("lr_tuple: {}", String::from_utf8(each_lr_tuple.2.decode(0, each_lr_tuple.2.len())).unwrap());
    }



    let ngsread_file = File::open(&ngsread_input_file).expect("Error during opening the file");
    let mut reader = faReader::new(ngsread_file);
    let mut record = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", ngsread_input_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }
    let chunk_size: usize = sequences.len() / (threads - 1);
    let sequences_ref = Arc::new(sequences);
    let lr_tuple_ref    = Arc::new(lr_tuple);
    //let mut cbf_oyadama: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];

    //let mut product_size_hashmap = HashMap::<u32, usize>::new();
    thread::scope(|scope|{
        let mut children_1 = Vec::new();
        for i in 1..threads {
            let sequences_ref = Arc::clone(&sequences_ref);
            let lr_tuple_ref    = Arc::clone(&lr_tuple_ref);
            children_1.push(
                scope.spawn(move || 
                    {   
                        let start_idx: usize = (i - 1) * chunk_size;
                        let end_idx  : usize;
                        if i != threads - 1{
                            end_idx = i * chunk_size;
                        }else{
                            end_idx = sequences_ref.len() - 1;
                        }
                        let lr_tuple_ref_mine = (*lr_tuple_ref).clone();
                        let slice_sequences = Vec::from(sequences_ref[start_idx..end_idx].to_vec());

                        eprintln!("start calling aggregate_length_between_lr_tuple[{}], # of sequence: {}", i, &slice_sequences.len());
                        let cbf: Vec<u8> = aggregate_length_between_lr_tuple(&slice_sequences, i, &lr_tuple_ref_mine, usize::MAX);
                        eprintln!("finish calling aggregate_length_between_lr_tuple[{}]", i);
                        return cbf
                    }
                )
            )
        }
        eprintln!("start  writing to output file: {:?}", &output_file);
        let mut file = File::create(&output_file).unwrap();
        for child in children_1 {
            let result = child.join(); // Resultを返さず、直接戻り値を取得
            file.write_all(&result.unwrap()).unwrap();
        }
        eprintln!("finish writing to output file: {:?}", &output_file);
    });
}
