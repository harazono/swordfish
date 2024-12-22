use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::thread;

use bio::io::fasta;
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::distributions::Alphanumeric;
use rand::{thread_rng, Rng};

fn get_db_taxids(db: &str) -> HashSet<String> {
    let output = Command::new("blastdbcmd")
        .args(&["-db", db, "-entry", "all", "-outfmt", "%T"])
        .output()
        .expect("failed to execute process");
    let stdout = String::from_utf8_lossy(&output.stdout);
    stdout.lines().map(|s| s.to_string()).collect()
}

fn main() {
    // 引数の解析
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("Usage: {} <tgt_db_list_filename> <taxon_id_to_be_ignored> <fasta_file> <slots> <output_dir>", args[0]);
        std::process::exit(1);
    }
    let tgt_db_list_filename = &args[1];
    let taxon_id_to_be_ignored = &args[2];
    let fasta_file = &args[3];
    let slots: usize = args[4].parse().expect("Invalid slots value");
    let output_dir = &args[5];

    // 除外する taxon ID の読み込み
    let negative_taxon_ids: HashSet<String> =
        BufReader::new(File::open(taxon_id_to_be_ignored).expect("Failed to open taxon ID file"))
            .lines()
            .map(|l| l.expect("Failed to read line").trim().to_string())
            .collect();

    // 出力ディレクトリの作成
    std::fs::create_dir_all(output_dir).expect("Failed to create output directory");
    let log_dir = Path::new(output_dir).join("logs");
    std::fs::create_dir_all(&log_dir).expect("Failed to create log directory");

    // 出力ファイルの設定
    let output_file_path = Path::new(output_dir).join("blastn_results.gz");

    // データベースリストの読み込み
    let db_list: Vec<String> = BufReader::new(
        File::open(tgt_db_list_filename).expect("Failed to open database list file"),
    )
    .lines()
    .map(|l| l.expect("Failed to read line").trim().to_string())
    .collect();

    // FASTA ファイルの読み込みと重複除去
    let reader = fasta::Reader::from_file(fasta_file).expect("Failed to open FASTA file");
    let mut seen_sequences = HashSet::new();
    let mut seq_name_mapping = Vec::new();
    let unique_sequences: Vec<bio::io::fasta::Record> = reader
        .records()
        .filter_map(|record| {
            let record = record.expect("Failed to read record");
            let seq = record.seq().to_vec();
            if !seen_sequences.contains(&seq) {
                seen_sequences.insert(seq);
                let new_seq_name: String = thread_rng()
                    .sample_iter(&Alphanumeric)
                    .take(32)
                    .map(char::from)
                    .collect();
                seq_name_mapping.push((new_seq_name.clone(), record.id().to_string()));
                Some(bio::io::fasta::Record::with_attrs(
                    &new_seq_name,
                    None,
                    &record.seq(),
                ))
            } else {
                None
            }
        })
        .collect();

    // 一時 FASTA ファイルの作成
    let temp_fasta_file = Path::new(output_dir).join(format!(
        "unique_sequences_{}.fasta",
        thread_rng()
            .sample_iter(&Alphanumeric)
            .take(8)
            .map(char::from)
            .collect::<String>()
    ));
    let mut writer =
        fasta::Writer::to_file(&temp_fasta_file).expect("Failed to create temporary FASTA file");
    for record in &unique_sequences {
        writer.write_record(record).expect("Failed to write record");
    }

    // BLAST 検索の実行
    thread::scope(|scope| {
        let temp_files: Vec<_> = db_list
            .into_iter()
            .map(|db| {
                let log_file = log_dir.join(format!("blastn_{}.log", db));
                let exec_cmd_logfile = log_dir.join(format!("blastn_{}_exec_cmd.txt", db));
                let temp_output = Path::new(output_dir).join(format!("{}.gz", db));

                let taxon_ids_in_db = get_db_taxids(&db);
                let negative_taxon_ids_for_this_db: Vec<&str> = negative_taxon_ids
                    .intersection(&taxon_ids_in_db)
                    .map(|s| s.as_str())
                    .collect();
                let negative_taxids_arg = if !negative_taxon_ids_for_this_db.is_empty() {
                    format!("-negative_taxids {}", negative_taxon_ids_for_this_db.join(","))
                } else {
                    String::new()
                };

                let execute_cmd = format!(
                    "export BLASTDB=/usr/local/db/blast/ncbi/v5/; \
                    echo '{}' > {:?}; \
                    blastn \
                    -task blastn-short \
                    -query {:?} \
                    -db {} \
                    -outfmt '6 qseqid sseqid sacc qlen qstart qend slen sstart send qseq sseq evalue length staxid staxids ssciname scomname' \
                    -num_threads {} \
                    -dust no \
                    -soft_masking false \
                    -word_size 15 \
                    -best_hit_overhang 0.01 \
                    -best_hit_score_edge 0.49 \
                    -perc_identity 100 \
                    {} \
                    2>> {:?}",
                    negative_taxids_arg,
                    log_file,
                    temp_fasta_file,
                    db,
                    slots,
                    negative_taxids_arg,
                    log_file
                );

                // スレッドを作成して BLAST 検索を実行
                scope.spawn(move || {
                    // コマンドの実行
                    let output = Command::new("sh")
                        .arg("-c")
                        .arg(execute_cmd)
                        .output()
                        .expect("failed to execute process");

                    // ログの書き出し
                    let mut log_writer =
                        BufWriter::new(File::create(&exec_cmd_logfile).expect("Failed to create log file"));
                    log_writer
                        .write_all(&output.stderr)
                        .expect("Failed to write to log file");

                    temp_output
                })
                .join()
                .unwrap()
            })
            .collect();

        // BLAST 結果の結合と出力
        let f_out = File::create(output_file_path).expect("Failed to create output file");
        let mut encoder = GzEncoder::new(f_out, Compression::default());
        for temp_file in temp_files {
            let f_in = File::open(temp_file).expect("Failed to open temporary file");
            let reader = BufReader::new(f_in);
            for line in reader.lines() {
                let line = line.expect("Failed to read line");
                let mut fields: Vec<&str> = line.split('\t').collect();
                let new_seq_name = fields[0];
                let original_seq_name = seq_name_mapping
                    .iter()
                    .find(|(k, _)| k == &new_seq_name)
                    .map(|(_, v)| v)
                    .expect("Failed to find original sequence name");
                fields[0] = original_seq_name;
                let new_line = fields.join("\t");
                writeln!(encoder, "{}", new_line).expect("Failed to write to output file");
            }
        }
    });
}
