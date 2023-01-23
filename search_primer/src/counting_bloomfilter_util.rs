pub const L_LEN: usize = 31;
pub const R_LEN: usize = 31;
pub const HASHSET_SIZE: usize = (u32::MAX >> 4) as usize;
const dupplications: u32 = 1;
use crate::sequence_encoder_util::{DnaSequence, decode_u128_2_dna_seq};
use std::fs::File;
use sha2::Sha256;
use sha2::Digest;

use std::time::{Instant};
use std::cmp;
use std::collections::HashSet;

//use bio::io::fastq::Reader as fqReader;
//use bio::io::fastq::Record as fqRecord;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
//use bio::io::fastq::FastqRead;
use bio::io::fasta::FastaRead;

//pub const BLOOMFILTER_TABLE_SIZE: usize = u32::MAX as usize;
pub const BLOOMFILTER_TABLE_SIZE: usize = (u32::MAX >> 1) as usize;

//全てのL, Rと、hash値を出力する
//部分配列のdecoderを書き、テストする
pub fn build_counting_bloom_filter(sequences: &Vec<DnaSequence>, start_idx: usize, end_idx: usize, thread_id: usize) -> Vec<u32>{
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let chunk_max: usize = 500;

    let mut loop_cnt:usize = 0;
    eprintln!("Allocating Vec<u32> where BLOOMFILTER_TABLE_SIZE = {}", BLOOMFILTER_TABLE_SIZE);
    //let mut ret_array: Vec<u32> = Vec::with_capacity(BLOOMFILTER_TABLE_SIZE);
    let mut ret_array: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];
    eprintln!("Filling Vec<u32; {}> with 0", BLOOMFILTER_TABLE_SIZE);
    eprintln!("finish allocating");

    let start_time = Instant::now();
    let mut previous_time = start_time.elapsed();

    'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
        let mut add_bloom_filter_cnt: usize = 0;
        let mut l_window_cnt: usize         = 0;
        loop_cnt += 1;
        l_window_start = 0;
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len() + 1{
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let (l_has_poly_base,     l_offset_1) = current_sequence.has_poly_base(    l_window_start, l_window_end);
            let (l_has_simple_repeat, l_offset_2) = current_sequence.has_simple_repeat(l_window_start, l_window_end);
            let (l_has_2base_repeat,  l_offset_3) = current_sequence.has_2base_repeat( l_window_start, l_window_end);
            //eprintln!("l_has_poly_base: {}, {}, l_window_start: {}, l_window_end: {}",l_has_poly_base, String::from_utf8(decode_u128_2_dna_seq(&current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end]]), 19)).unwrap(), l_window_start, l_window_end);
            if l_has_poly_base|l_has_simple_repeat|l_has_2base_repeat {
                l_window_start += cmp::max(cmp::max(l_offset_1, l_offset_2), l_offset_3) + 1;
                //l_window_start += 1;
                continue 'each_l_window;
            }
            r_window_start = l_window_end;
            'each_r_window: loop{
                r_window_end = r_window_start + R_LEN;
                if r_window_end >= current_sequence.len() + 1{
                    let end = start_time.elapsed();
                    eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
                    previous_time = end;
                    continue 'each_read;
                }
                if r_window_end - l_window_start > chunk_max - R_LEN{
                    break 'each_r_window;
                }
                let (m_has_poly_base,     m_offset_1) = current_sequence.has_poly_base(    r_window_start, r_window_end);
                let (m_has_simple_repeat, m_offset_2) = current_sequence.has_simple_repeat(r_window_start, r_window_end);
                let (m_has_2base_repeat,  m_offset_3) = current_sequence.has_2base_repeat( r_window_start, r_window_end);
                if m_has_poly_base|m_has_simple_repeat|m_has_2base_repeat {
                    r_window_start += cmp::max(cmp::max(m_offset_1, m_offset_2), m_offset_3) + 1;
                    //r_window_start += 1;
                    continue 'each_r_window;
                }

                add_bloom_filter_cnt += 1;
                //assert!(l_has_poly_base == false, "assertion failed");
                //assert!(m_has_poly_base == false, "assertion failed");
                //assert!(r_has_poly_base == false, "assertion failed");
                let lmr_string: u128 = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end], [r_window_start, r_window_end]]);
                let table_indice:[u32;8] = hash_from_u128(lmr_string);//u128を受けてhashを返す関数
                for i in 0..8{
                    let idx: usize = table_indice[i] as usize;
                    if ret_array[idx] == u32::MAX{
                        eprintln!("index {} reaches u32::MAX", idx);
                    }else{
                        ret_array[idx] += 1;
                    }
                }




                r_window_start += 1;
            }
            l_window_start += 1;
        }
        let end = start_time.elapsed();
        eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
        previous_time = end;
    }
    return ret_array;
}

//BLOOMFILTER_TABLE_SIZEの範囲内で柔軟にhash値を返すようにする。
fn hash_from_u128(source: u128) -> [u32; 8]{
    let mut ret_val: [u32;8] = [0;8];
    let mut hasher = Sha256::new();
    let mut u8_array: [u8; 16] = [0; 16];
    let mut src_copy: u128 = source;
    for i in 0..16{
        u8_array[i] = (src_copy & 255).try_into().unwrap();
        src_copy >>= 8;
    }
    hasher.update(u8_array);
    let result = hasher.finalize();
    let sha256_bit_array = result.as_slice();//&[u8;32]
    for i in 0..8{
        for j in 0..4{
            ret_val[i] += (sha256_bit_array[i * 4 + j] as u32);
            ret_val[i] <<= 8;
        }
        ret_val[i] %= (BLOOMFILTER_TABLE_SIZE as u32);
    }
    return ret_val;
}

fn count_occurence_from_counting_bloomfilter_table(counting_bloomfilter_table: &Vec<u32>, indice: [u32; 8]) -> u32{
    let mut retval: u32 = u32::MAX;
    for index in indice{
        if counting_bloomfilter_table[index as usize] < retval{
            retval = counting_bloomfilter_table[index as usize];
        }
    }
    return retval;
}


pub fn number_of_high_occurence_kmer(source_table: &Vec<u32>, sequences: &Vec<DnaSequence>, start_idx: usize, end_idx: usize, threshold: u32, thread_id: usize) -> HashSet<u128>{
    let mut ret_table: HashSet<u128> = HashSet::with_capacity(HASHSET_SIZE);
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let chunk_max: usize = 141;
    let mut ho_lmr: usize = 0;


    let start = Instant::now();
    let mut previous_time = start.elapsed();
    let mut loop_cnt:usize = 0;
    'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
        let mut add_bloom_filter_cnt: usize = 0;
        let mut l_window_cnt: usize         = 0;
        loop_cnt += 1;
        l_window_start = 0;
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len() + 1{
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let (l_has_poly_base, l_offset_1)     = current_sequence.has_poly_base(l_window_start, l_window_end);
            let (l_has_simple_repeat, l_offset_2) = current_sequence.has_simple_repeat(l_window_start, l_window_end);
            let (l_has_2base_repeat, l_offset_3)  = current_sequence.has_2base_repeat(l_window_start, l_window_end);
            if l_has_poly_base|l_has_simple_repeat|l_has_2base_repeat {
                l_window_start += cmp::max(cmp::max(l_offset_1, l_offset_2), l_offset_3) + 1;
                //l_window_start += 1;
                continue 'each_l_window;
            }
            r_window_start = l_window_end;
            'each_r_window: loop{
                r_window_end = r_window_start + R_LEN;
                if r_window_end >= current_sequence.len() + 1{
                    let end = start.elapsed();
                    eprintln!("2nd loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}\tho_lmr: {}", thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt, ho_lmr);
                    previous_time = end;
                    continue 'each_read;
                }
                if r_window_end - l_window_start > chunk_max - R_LEN{
                    break 'each_r_window;
                }
                let (m_has_poly_base, m_offset_1)     = current_sequence.has_poly_base(r_window_start, r_window_end);
                let (m_has_simple_repeat, m_offset_2) = current_sequence.has_simple_repeat(r_window_start, r_window_end);
                let (m_has_2base_repeat, m_offset_3)  = current_sequence.has_2base_repeat(r_window_start, r_window_end);
                if m_has_poly_base|m_has_simple_repeat|m_has_2base_repeat {
                    r_window_start += cmp::max(cmp::max(m_offset_1, m_offset_2), m_offset_3) + 1;
                    //r_window_start += 1;
                    continue 'each_r_window;
                }
                add_bloom_filter_cnt += 1;
                let lmr_string:u128 = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end], [r_window_start, r_window_end]]);
                let table_indice:[u32;8] = hash_from_u128(lmr_string);//u128を受けてhashを返す関数
                let occurence: u32 = count_occurence_from_counting_bloomfilter_table(source_table, table_indice);
                if occurence >= threshold * dupplications{
                    ret_table.insert(lmr_string);
                    ho_lmr += 1;
                }

                r_window_start += 1;
            }
            l_window_start += 1;
        }
        let end = start.elapsed();
        eprintln!("2nd loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}\tho_lmr: {}", thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt, ho_lmr);
        previous_time = end;
    }
    return ret_table;
}
