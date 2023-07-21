use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

fn merge_files(file1: &Path, file2: &Path, output: &Path) -> std::io::Result<()> {
    let mut reader1 = BufReader::new(File::open(file1)?);
    let mut reader2 = BufReader::new(File::open(file2)?);
    let mut writer = BufWriter::new(File::create(output)?);

    let mut set = HashSet::new();

    let mut buffer = [0; 16]; // u128 is 16 bytes
    while let Ok(_) = reader1.read_exact(&mut buffer) {
        let value = u128::from_be_bytes(buffer);
        set.insert(value);
    }

    while let Ok(_) = reader2.read_exact(&mut buffer) {
        let value = u128::from_be_bytes(buffer);
        set.insert(value);
    }

    let mut values: Vec<_> = set.into_iter().collect();
    values.sort_unstable();

    for value in values {
        writer.write_all(&value.to_be_bytes())?;
    }

    Ok(())
}

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <file1> <file2> <output>", args[0]);
        return;
    }

    let file1 = Path::new(&args[1]);
    let file2 = Path::new(&args[2]);
    let output = Path::new(&args[3]);

    if let Err(e) = merge_files(file1, file2, output) {
        eprintln!("Error: {}", e);
    }
}
