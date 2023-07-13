use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write, BufRead};
use getopts::Options;
use std::env;

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    opts.optopt("i", "input", "set input file name", "NAME");
    opts.optopt("o", "output", "set output directory", "DIR");
    opts.optopt("p", "prefix", "set output file prefix", "PREFIX");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("{}", f.to_string()),
    };
    let input_file_name = matches.opt_str("i").unwrap();
    let output_dir = matches.opt_str("o").unwrap();
    let output_prefix = matches.opt_str("p").unwrap();

    let mut set = HashSet::new();

    // Read all files
    let file = File::open(input_file_name)?;
    let reader = BufReader::new(file);
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
    }

    // Write to output files
    let values_per_file = (set.len() + 4999) / 5000; // ceil division
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
