extern crate getopts;

use getopts::Options;
use std::env;
use std::fs::File;
use std::io::{Read, Write, BufWriter};
use std::collections::BinaryHeap;
use std::cmp::Ordering;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    opts.optopt("a", "file1", "first sorted file", "FILE");
    opts.optopt("b", "file2", "second sorted file", "FILE");
    opts.optopt("o", "output", "output file", "OUTPUT");
    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("{}", f.to_string()),
    };

    if matches.opt_present("h") {
        print_usage(&opts);
        return;
    }

    let file1_path = matches.opt_str("a").expect("Please provide the first file using -a or --file1");
    let file2_path = matches.opt_str("b").expect("Please provide the second file using -b or --file2");
    let output_path = matches.opt_str("o").unwrap_or("merged_output.txt".to_string());

    merge_and_uniq(&file1_path, &file2_path, &output_path);
}

fn print_usage(opts: &Options) {
    let brief = "Usage: merge_and_uniq -a FILE1 -b FILE2 [-o OUTPUT]";
    print!("{}", opts.usage(&brief));
}

fn merge_and_uniq(file1_path: &str, file2_path: &str, output_path: &str) {
    let mut file1 = File::open(file1_path).expect("Failed to open the first file");
    let mut file2 = File::open(file2_path).expect("Failed to open the second file");

    let mut numbers = BinaryHeap::new();

    loop {
        let mut buffer = [0u8; 16];
        match file1.read_exact(&mut buffer) {
            Ok(_) => {
                let number = u128::from_be_bytes(buffer);
                numbers.push(Reverse(number));
            }
            Err(_) => break,
        }
    }

    loop {
        let mut buffer = [0u8; 16];
        match file2.read_exact(&mut buffer) {
            Ok(_) => {
                let number = u128::from_be_bytes(buffer);
                numbers.push(Reverse(number));
            }
            Err(_) => break,
        }
    }

    let output_file = File::create(output_path).expect("Failed to create the output file");
    let mut writer = BufWriter::new(output_file);

    let mut last_written: Option<u128> = None;

    while let Some(Reverse(number)) = numbers.pop() {
        if last_written.as_ref() != Some(&number) {
            writer.write_all(&number.to_le_bytes()).expect("Failed to write to the output file");
            last_written = Some(number);
        }
    }
}

#[derive(Eq, PartialEq, Ord, PartialOrd)]
struct Reverse(u128);
