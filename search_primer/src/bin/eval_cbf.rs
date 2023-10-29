use search_primer::counting_bloomfilter_util::{hash_from_u128, BLOOMFILTER_TABLE_SIZE};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::BufWriter;
use std::io::Write;
use std::iter::zip;
use std::sync::{Arc, Mutex};
use std::thread;

extern crate bio;
extern crate getopts;
extern crate rdxsort;
use getopts::Options;
use std::{env, process};

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn calc_index(thread_id: usize, threads: usize, limit: usize) -> (usize, usize) {
    // 各スレッドが担当する基本の範囲を計算
    let base_chunk_size = limit / threads;
    let remainder = limit % threads;

    // thread_id=0が余りの部分を担当する
    if thread_id == 0 {
        return (0, base_chunk_size + remainder);
    }

    // それ以外のスレッドは基本の範囲を担当
    let start = thread_id * base_chunk_size + remainder;
    let end = start + base_chunk_size;

    (start, end)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optopt("o", "output", "set output file name prefix.", "NAME");
    opts.optopt(
        "t",
        "thread",
        "number of threads to use for radix sort. default value is 8.",
        "THREAD",
    );

    opts.optopt(
        "l",
        "limit",
        "0..l is used for hash function. default value is 10_000.",
        "LIMIT",
    );

    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f.to_string())
        }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let threads: usize = if matches.opt_present("t") {
        std::cmp::max(matches.opt_str("t").unwrap().parse::<usize>().unwrap(), 2)
    } else {
        8
    };

    let limit: usize = if matches.opt_present("l") {
        matches.opt_str("l").unwrap().parse::<usize>().unwrap()
    } else {
        10_000
    };

    let output_file: String = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        format!("eval_cbf",)
    };

    let mut cbf_oyadama: Vec<u16> = vec![0u16; BLOOMFILTER_TABLE_SIZE];

    thread::scope(|scope: &thread::Scope<'_, '_>| {
        let mut children_1: Vec<thread::ScopedJoinHandle<'_, Vec<u16>>> = Vec::new();
        for thread_idx in 0..(threads - 1) {
            children_1.push(scope.spawn(move || {
                let mut index_table_in_a_thread: Vec<u16> = vec![0u16; BLOOMFILTER_TABLE_SIZE];
                let (start, end) = calc_index(thread_idx, threads, limit);
                eprintln!(
                    "thread {:02} start: {:04} end: {:04}",
                    thread_idx, start, end
                );
                for hash_src in start..=end {
                    if (hash_src - start) % 100000 == 0 {
                        eprintln!(
                            "thread {:02} hash_src: {:08} in {:08}",
                            thread_idx,
                            hash_src,
                            end - start
                        );
                    }
                    let hash_values = hash_from_u128(hash_src as u128, BLOOMFILTER_TABLE_SIZE);
                    for i in 0..8 {
                        *index_table_in_a_thread
                            .get_mut(hash_values[i] as usize)
                            .unwrap() += 1;
                    }
                }
                index_table_in_a_thread
            }))
        }
        for child in children_1 {
            let cbf: Vec<u16> = child.join().unwrap();
            zip(cbf_oyadama.iter_mut(), cbf)
                .for_each(|(x, y)| *x = x.checked_add(y).unwrap_or(u16::MAX));
        }
    });

    let mut index_counter: HashMap<u32, usize> = HashMap::new();
    for i in 0..BLOOMFILTER_TABLE_SIZE {
        let accumurated_val: u32 = cbf_oyadama[i] as u32;
        index_counter.insert(
            accumurated_val,
            index_counter.get(&accumurated_val).unwrap_or(&0) + 1,
        );
    }
    let output_file1: String = String::from(&output_file) + "_hist.csv";
    let mut w1: BufWriter<File> = BufWriter::new(fs::File::create(&output_file1).unwrap());
    let mut sorted_vec: Vec<(&u32, &usize)> = index_counter.iter().collect();
    sorted_vec.sort_by_key(|k| k.0);
    for each_cbf_hist in &index_counter {
        writeln!(&mut w1, "{},{}", each_cbf_hist.0, each_cbf_hist.1).unwrap();
    }
    let oyadama_file: String = String::from(&output_file) + "_cbf.txt";
    let mut w2: BufWriter<File> = BufWriter::new(fs::File::create(&oyadama_file).unwrap());

    for value in &cbf_oyadama {
        writeln!(&mut w2, "{}", value).unwrap();
    }
}
