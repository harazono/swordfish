use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write, BufRead};
use getopts::Options;
use std::env;

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} -i INPUT -o OUTPUT -p PREFIX [-n NUMFILES | -s FILESIZE]", program);
    print!("{}", opts.usage(&brief));
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("i", "input", "set input file name", "NAME");
    opts.optopt("o", "output", "set output directory", "DIR");
    opts.optopt("p", "prefix", "set output file prefix", "PREFIX");
    opts.optopt("n", "numfiles", "set number of output files", "NUM");
    opts.optopt("s", "filesize", "set max size of each output file", "SIZE");
    opts.optflag("h", "help", "print this help menu");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            eprintln!("{}", f.to_string());
            print_usage(&program, opts);
            return Ok(());
        }
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return Ok(());
    }
    let input_file_name = matches.opt_str("i").unwrap();
    let output_dir      = matches.opt_str("o").unwrap();
    let output_prefix   = matches.opt_str("p").unwrap();
    let num_files       = matches.opt_str("n").map(|n| n.parse::<usize>().unwrap());
    let file_size       = matches.opt_str("s").map(|s| s.parse::<usize>().unwrap());

    if num_files.is_some() && file_size.is_some() {
        eprintln!("Error: Cannot specify both number of files and file size");
        print_usage(&program, opts);
        return Ok(());
    }

    
    // Read all files
    let file = File::open(input_file_name)?;
    let reader = BufReader::new(file);
    let mut intermediate_files = Vec::new();
    let mut set: HashSet<u128> = HashSet::new();

    for (index, line) in reader.lines().enumerate() {
        let path = line?;
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);

        loop {
            let mut buffer = [0; 16]; // u128 is 16 bytes
            match reader.read_exact(&mut buffer) {
                Ok(_) => {
                    let value = u128::from_le_bytes(buffer);
                    set.insert(value);
                }
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
        }
        eprintln!("Finished reading file {}: {}", index + 1, path);

        // If we've read 20 files, write the set to an intermediate file
        if (index + 1) % 20 == 0 {
            let intermediate_file_name = format!("intermediate_{}.bin", (index + 1) / 20);
            let mut writer = BufWriter::new(File::create(&intermediate_file_name)?);
            for value in &set {
                writer.write_all(&value.to_le_bytes())?;
            }
            intermediate_files.push(intermediate_file_name);
            set.clear(); // Clear the set for the next batch of files
        }
    }

    // Don't forget the last batch if it's less than 20 files
    if !set.is_empty() {
        let intermediate_file_name = format!("intermediate_{}.bin", intermediate_files.len() + 1);
        let mut writer = BufWriter::new(File::create(&intermediate_file_name)?);
        for value in &set {
            writer.write_all(&value.to_le_bytes())?;
        }
        intermediate_files.push(intermediate_file_name);
    }

    // Now merge the intermediate files
    let mut set = HashSet::new();
    for (index, file_name) in intermediate_files.iter().enumerate() {
        let file = File::open(file_name)?;
        let mut reader = BufReader::new(file);

        loop {
            let mut buffer = [0; 16]; // u128 is 16 bytes
            match reader.read_exact(&mut buffer) {
                Ok(_) => {
                    let value = u128::from_le_bytes(buffer);
                    set.insert(value);
                }
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
        }
        eprintln!("Finished merging file {}: {}", index + 1, file_name);
    }

    // Write to output files
    let values_per_file = match (num_files, file_size) {
        (Some(n), None) => (set.len() + n - 1) / n, // ceil division
        (None, Some(s)) => s / 16, // each value is 16 bytes
        (None, None) => (set.len() + 4999) / 5000, // default to 5000 files
        _ => unreachable!(), // This should never happen
    };
    let mut file_index = 0;
    let mut values_in_file = 0;
    let mut writer = BufWriter::new(File::create(format!("{}/{}{}.bin", output_dir, output_prefix, file_index))?);

    for (index, value) in set.into_iter().enumerate() {
        if values_in_file == values_per_file {
            file_index += 1;
            values_in_file = 0;
            writer = BufWriter::new(File::create(format!("{}/{}{}.bin", output_dir, output_prefix, file_index))?);
        }
        writer.write_all(&value.to_le_bytes())?;
        values_in_file += 1;

        if index % 1000 == 0 {
            eprintln!("Written {} values...", index);
        }
    }

    Ok(())
}
