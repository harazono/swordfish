use binary_heap_plus::BinaryHeap;
use std::cmp::Ordering;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

const U128_SIZE: usize = 16; // u128 is 16 bytes

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let mut input_files: Vec<String> = Vec::new();
    let mut output_file: String = String::new();

    // Parse command line arguments
    let mut args_iter: std::iter::Skip<std::slice::Iter<'_, String>> = args.iter().skip(1);
    while let Some(arg) = args_iter.next() {
        match arg.as_str() {
            "-o" => output_file = args_iter.next().unwrap().clone(),
            "-i" => input_files.push(args_iter.next().unwrap().clone()),
            _ => {}
        }
    }

    if input_files.is_empty() || output_file.is_empty() {
        eprintln!("Usage: merge_u128 -i <input_files> -o <output_file>");
        return Ok(());
    }

    // Open all input files and create BufReaders
    let mut readers: Vec<_> = input_files
        .iter()
        .map(|file: &String| BufReader::new(File::open(file).expect("Unable to open input file")))
        .collect();

    // Create a BufWriter for the output file
    let mut writer: BufWriter<File> = BufWriter::new(File::create(&output_file)?);

    // Create a binary heap for merging
    let mut heap = BinaryHeap::new_by(|a: &(u128, usize), b: &(u128, usize)| {
        a.0.cmp(&b.0).reverse() // We want a min-heap
    });

    // Initialize the heap with the first element from each reader
    for (index, reader) in readers.iter_mut().enumerate() {
        let mut buffer: [u8; 16] = [0u8; U128_SIZE];
        if reader.read_exact(&mut buffer).is_ok() {
            heap.push((u128::from_le_bytes(buffer), index));
        }
    }

    // Track the last number written to avoid duplicates
    let mut last_written: Option<u128> = None;

    // Iterate over the heap and write unique values to the output file
    while let Some((number, index)) = heap.pop() {
        if last_written != Some(number) {
            writer.write_all(&number.to_le_bytes())?;
            last_written = Some(number);
        }
        // Read the next number from the reader that provided the last number
        let mut buffer: [u8; 16] = [0u8; U128_SIZE];
        if readers[index].read_exact(&mut buffer).is_ok() {
            heap.push((u128::from_le_bytes(buffer), index));
        }
    }

    // Flush the writer to ensure all data is written to the output file
    writer.flush()?;

    // Output results to stderr in markdown format
    eprintln!("| Input File | Number of u128 integers |");
    eprintln!("|------------|------------------------|");
    for file in &input_files {
        eprintln!("| {} | {} |", file, count_u128_integers(file)?);
    }
    eprintln!("| Output File | {} |", output_file);
    eprintln!(
        "| Number of u128 integers | {} |",
        last_written.map_or(0, |_| heap.len())
    );

    Ok(())
}

fn count_u128_integers(file_name: &str) -> io::Result<usize> {
    let file: File = File::open(file_name)?;
    let metadata: std::fs::Metadata = file.metadata()?;
    Ok(metadata.len() as usize / U128_SIZE)
}
