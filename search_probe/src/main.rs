extern crate bio;
extern crate rdxsort;
extern crate getopts;
use getopts::Options;
use std::{env, process};
use std::fs;
use std::thread;
use std::fs::File;
use std::iter::zip;
use std::io::{BufWriter, Write};
use std::collections::HashSet;
use std::collections::HashMap;
use std::sync::{Mutex, Arc};
use search_probe::find_taqman_probe::BLOOMFILTER_TABLE_SIZE;
use search_probe::find_taqman_probe::{PROBE_LEN, HASHSET_SIZE};
use search_probe::find_taqman_probe::{build_counting_bloom_filter, number_of_high_occurence_kmer, aggregate_length_between_primer};
use search_probe::sequence_encoder_util::{decode_u128_2_dna_seq};
use search_probe::sequence_encoder_util::DnaSequence;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use crate::bio::io::fasta::FastaRead;
use std::io::BufReader;
use std::io::BufRead;

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
    opts.optopt("s", "triming_size", "each primer will be trimmed to this size. 3' side will be remain.", "TRIMSIZE");
    opts.optopt("l", "max_product_size", "max product size of PCR", "PRODUCT SIZE");
    opts.optopt("p", "primer", "input primers (TSV file).", "TSV FILE");
    opts.optflag("e", "extract", "extract genomic region where primer is located");
    opts.optflag("b", "binary", "outputs binary file");
    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

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

    let triming_size: usize = if matches.opt_present("s") {
        matches.opt_str("s").unwrap().parse::<usize>().unwrap()
    }else{
        15
    };

    let threshold:u32 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u32>().unwrap()
    }else{
        8
    };
    
    let primer_filename = if matches.opt_present("p") {
        matches.opt_str("p").unwrap()
    }else{
        print_usage(&program, &opts);
        return;
    };

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    }else{
        format!("{:?}_threshold{}_threads{}.out", ngsread_input_file, threshold, threads)
    };

    let max_product_size = if matches.opt_present("l") {
        matches.opt_str("l").unwrap().parse::<usize>().unwrap()
    }else{
        usize::MAX
    };


    eprintln!("Input  file name: {:?}", ngsread_input_file);
    eprintln!("Output file name: {:?}", output_file);
    eprintln!("Number of threads: {:?}", threads);
    eprintln!("Triming size: {:?}", triming_size);
    eprintln!("Input primers file: {:?}", primer_filename);
    eprintln!("Outputs binary file: {:?}", matches.opt_present("b"));
	eprintln!("Extract: {:?}", matches.opt_present("e"));

    /*
    primer id       left primer     right primer    primer left Tm  primer right Tm primer pair product Tm
    2baf2cca8e913286d099b61c1665bc2_1       AGGTGGTTAGTATAGGGATGGCAC        GCGGTTAGTCGACGCGCTTGAC  61.490  66.591  74.2
    bcaba327c863b6bb4266d8705996f0a_2       TAGGGTGGATAGCTTAGACGATGTCGG     CTTAACGGCGCGGTTAGTCGAC  65.125  63.999  75.4
    bcaba327c863b6b9099b61c1665bc2b_0       TAGGGTGGATAGCTTAGACGATGTCGG     CCTTAACGGCGCGGTTAGTCGAC 65.125  65.907  75.7
    2cca8e913286f2adc1665bc2b29e9a9_4       TGGCACATAGGACGTTAGGGT   GCCCGCCAGCCTACCTTAAC    60.898  63.217  75.4
    9cd45706d73b94fa15d233d86fdb60_3        TATCCACCCTAACGTCCTATGTGCCA      GAAACGTCGAATATCTGAGGGTCCA       64.991  62.400  72.6
    aebcb32a3a44ca19b61c1665bc2b29e_1       GGTGGTTAGTATAGGGATGGCAC CCTACCTTAACGGCGCGGTTAGTC        60.244  65.268  75.4
    baf2cca8e913286f61c1665bc2b29e9_2       GTGGTTAGTATAGGGATGGCAC  CTACCTTAACGGCGCGGTTAGTCG        57.996  65.492  74.0
    ae8c9f218edaed0a6d8705996f0aca7_3       TGGATAGCTTAGACGATGTCGGTGTCA     CCTACCTTAACGGCGCGGTTAGTC        65.224  65.268  75.1
    baf2cca8e913286ed099b61c1665bc2_2       GTGGTTAGTATAGGGATGGCAC  GCGGTTAGTCGACGCGCTTGA   57.996  65.776  74.2
    gc016 check_crossing_reaction 23/04/07 23:33:13$ 
    */

    let primer_file = File::open(&primer_filename).expect("Error during opening the file");
    //let primer_reader = BufReader::new(primer_file);

    // 各行を読み取り、タブで分割する
    let mut primer: Vec<(Vec<u8>, DnaSequence, DnaSequence)> = Vec::new();

    for result in BufReader::new(File::open(&primer_filename).unwrap()).lines() {
        let line = result.unwrap();
        if line.starts_with('p') {continue}
        let fields: Vec<&str> = line.split('\t').collect();
        // 1, 2, 3カラムのデータをu8型に変換してベクターに格納する
        //pub fn subsequence(&self, ranges: Vec<[usize; 2]>) -> DnaSequence{
        /*
        L primerは、添字が若い方が5'
        つまり前半数塩基をトリミングする
        R primerも、添字が若い方が5'
        つまり前半数塩基をトリミングする
        */
        let primer_id            = Vec::from(fields[0].as_bytes());
        let left_primer_seq      = &fields[1][fields[1].len() - triming_size..]; // 後ろから15文字を取得
        let right_primer_seq     = &fields[2][fields[2].len() - triming_size..]; // 後ろから15文字を取得
        let left_primer          = DnaSequence::new(&left_primer_seq.into());
        let right_primer         = DnaSequence::new(&right_primer_seq.into());

        //let left_primer          = DnaSequence::new(&Vec::from(fields[1].as_bytes()));
        //let right_primer         = DnaSequence::new(&Vec::from(fields[2].as_bytes()));
        let left_primer_revcomp  = DnaSequence::new(&Vec::from(fields[1].as_bytes())).reverse_complement();
        let right_primer_revcomp = DnaSequence::new(&Vec::from(fields[2].as_bytes())).reverse_complement();
        let id_1: Vec<u8> = [&primer_id, &b"LeftPrimerForward_LeftPrimerRevcomp"[..]].concat();
        let id_2: Vec<u8> = [&primer_id, &b"LeftPrimerForward_RightPrimerRevcomp"[..]].concat();
        let id_4: Vec<u8> = [&primer_id, &b"RightPrimerRevcomp_RightPrimerForward"[..]].concat();

        primer.push((id_1, left_primer.clone(), left_primer_revcomp.clone()));
        primer.push((id_2, left_primer.clone(), right_primer_revcomp.clone()));
        //primer.push((id_3, right_primer.clone(), left_primer_revcomp.clone()));
        primer.push((id_4, right_primer_revcomp.clone(), right_primer.clone()));

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
    let primer_ref    = Arc::new(primer);
    let mut cbf_oyadama: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];

    if matches.opt_present("e") {
        //let mut product_size_hashmap = HashMap::<u32, usize>::new();
        thread::scope(|scope|{
            let mut children_1 = Vec::new();
            for i in 1..threads {
                let sequences_ref = Arc::clone(&sequences_ref);
                let primer_ref    = Arc::clone(&primer_ref);
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
                            let primer_ref_mine = (*primer_ref).clone();
                            let slice_sequences = Vec::from(sequences_ref[start_idx..end_idx].to_vec());

                            eprintln!("start calling aggregate_length_between_primer[{}], # of sequence: {}", i, &slice_sequences.len());
                            let cbf: Vec<u8> = aggregate_length_between_primer(&slice_sequences, i, &primer_ref_mine, max_product_size);
                            eprintln!("finish calling aggregate_length_between_primer[{}]", i);
                            return cbf
                        }
                    )
                )
            }
            let mut results: Vec<u8> = Vec::new();
            for child in children_1 {
                match child.join() {
                    Ok(result) => {
                        results.extend(result);
                    }
                    Err(e) => {
                        eprintln!("Error: {:?}", e);
                    }
                }
            }
    
            let mut file = File::create(&output_file).unwrap();
            file.write_all(&results).unwrap();

            eprintln!("finish writing to output file: {:?}", &output_file);

        }
    );
    }else{
        thread::scope(|scope|{
            let mut children_1 = Vec::new();
            for i in 1..threads {
                let sequences_ref = Arc::clone(&sequences_ref);
                let primer_ref    = Arc::clone(&primer_ref);

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
                            let cbf: Vec<u32> = build_counting_bloom_filter(&sequences_ref, start_idx, end_idx, i, &primer_ref);
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
                let sequences_ref = Arc::clone(&sequences_ref);
                let primer_ref    = Arc::clone(&primer_ref);

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
                            eprintln!("thread [{}]: start calling number_of_high_occurence_kmer", i);
                            let h_cbf_h: HashSet<u128> = number_of_high_occurence_kmer(cbf_oyadama_ref, &sequences_ref, start_idx, end_idx, threshold, i, &primer_ref);
                            //h_cbf_h_oyadama = h_cbf_h_oyadama_ref.lock().unwrap().union(&h_cbf_h);
                            h_cbf_h_oyadama_ref.lock().unwrap().extend(&h_cbf_h);
                            eprintln!("thread [{}]: finish calling number_of_high_occurence_kmer", i);
                        }
                    )
                )
            }
            for child in children_2{
                let _ = child.join();
            }
        });


        let high_occurence_kmer: Vec<u128> = Vec::from_iter(h_cbf_h_oyadama.lock().unwrap().clone());
        let mut w = BufWriter::new(fs::File::create(&output_file).unwrap());

        let mut previous_kmer: u128 = 0;
        let mut cnt = 0;
        let mut buf_array: [u8; 16] = [0; 16];
        let mut buf_num: u128;

        if matches.opt_present("r") {
            eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
            for each_kmer in &high_occurence_kmer{
                if previous_kmer != *each_kmer{
                    cnt += 1;
                }
                previous_kmer = *each_kmer;
            }
            writeln!(&mut w, "k-mer count: {}\tthreshold: {}\tinput file {:?}", cnt, threshold, &ngsread_input_file).unwrap();
        }
        if !matches.opt_present("r") && matches.opt_present("b"){
            eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
            for each_kmer in &high_occurence_kmer{
                if previous_kmer != *each_kmer{
                    cnt += 1;
                    buf_num = *each_kmer;
                    for i in 0..16{
                        buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
                        buf_num >>= 8;
                    }
                    w.write(&buf_array).unwrap();
                }
                previous_kmer = *each_kmer;
            }
        }
        if !matches.opt_present("r") && !matches.opt_present("b"){
            eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
            for each_kmer in &high_occurence_kmer{
                if previous_kmer != *each_kmer{
                    cnt += 1;
                    writeln!(&mut w, "{:?}", String::from_utf8(decode_u128_2_dna_seq(&each_kmer, PROBE_LEN)).unwrap()).unwrap();
                }
                previous_kmer = *each_kmer;
            }
        }


        eprintln!("finish writing to output file: {:?}", &output_file);
        eprintln!("threshold: {}({}x63)\tcardinarity: {}", threshold, threshold / 63, cnt);
        eprintln!("total cardinarity: {}", cnt);
        eprintln!("threads: {}\tthreshold: {}\tinput file {:?}", threads, threshold, &ngsread_input_file);
    }
}