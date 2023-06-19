pub const PROBE_LEN: usize = 30;
pub const HASHSET_SIZE: usize = (u32::MAX >> 4) as usize;
const DUPPLICATION: u32 = 1;
use crate::sequence_encoder_util::{DnaSequence, decode_u128_2_dna_seq};
use sha2::Sha256;
use sha2::Digest;

use std::time::{Instant};
use std::collections::HashSet;
pub const BLOOMFILTER_TABLE_SIZE: usize = (u32::MAX >> 1) as usize;

pub fn build_counting_bloom_filter(sequences: &Vec<DnaSequence>, start_idx: usize, end_idx: usize, thread_id: usize, primer: &Vec<(Vec<u8>, DnaSequence, DnaSequence)>) -> Vec<u32>{
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut m_window_start: usize;
    let mut m_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let mut primer_l_seq:   u128;
    let mut primer_l_size:  usize;
    let mut primer_r_seq:   u128;
    let mut primer_r_size:  usize;
    let mut mask_l:         u128;
    let mut mask_r:         u128;
    let chunk_max = 200;
    let mut loop_cnt:usize = 0;
    eprintln!("Allocating Vec<u32> where BLOOMFILTER_TABLE_SIZE = {}", BLOOMFILTER_TABLE_SIZE);
    let mut ret_array: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];
    eprintln!("Filling Vec<u32; {}> with 0", BLOOMFILTER_TABLE_SIZE);
    eprintln!("finish allocating");

    let start_time = Instant::now();
    let mut previous_time = start_time.elapsed();
    'each_primer: for current_primer in primer.iter() {
        primer_l_size = current_primer.1.len();
        primer_l_seq  = current_primer.1.subsequence_as_u128(vec!([0, primer_l_size]));
        primer_r_size = current_primer.2.len();
        primer_r_seq  = current_primer.2.subsequence_as_u128(vec!([0, primer_r_size]));
        mask_l        = u128::MAX >> (64 - primer_l_size) * 2;
        mask_r        = u128::MAX >> (64 - primer_r_size) * 2;
    

        'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
            let mut add_bloom_filter_cnt: usize = 0;
            let mut l_window_cnt: usize         = 0;
            loop_cnt += 1;
            l_window_start = 0;
            if current_sequence.len() < primer_l_size || current_sequence.len() <  primer_r_size ||current_sequence.len() < PROBE_LEN{
                continue 'each_read;
            }
            'each_l_window: loop{
                l_window_end = l_window_start + primer_l_size;
                if l_window_end >= current_sequence.len() + 1{
                    break 'each_l_window;
                }
                l_window_cnt += 1;
                let l_window_as_u128: u128 = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end]]);
                if (l_window_as_u128 & mask_l) != primer_l_seq{
                    l_window_start += 1;
                    continue 'each_l_window;
                }
                r_window_start = l_window_end + PROBE_LEN;
                'each_r_window: loop{
                    r_window_end = r_window_start + primer_r_size;
                    if r_window_end >= current_sequence.len() + 1{
                        let end = start_time.elapsed();
                        eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
                        previous_time = end;
                        continue 'each_read;
                    }
                    if r_window_end - l_window_start > chunk_max{
                        break 'each_l_window;
                    }
                    let r_window_as_u128: u128 = current_sequence.subsequence_as_u128(vec![[r_window_start, r_window_end]]);
                    if (r_window_as_u128 & mask_r) != primer_r_seq {
                        r_window_start += 1;
                        continue 'each_r_window;
                    }
                    //ここまでで、LとRが一致してる。この下のループでProbeを探索する。
                    m_window_start = l_window_end;
                    'each_m_window: loop{
                        m_window_end = m_window_start + PROBE_LEN;
                        if m_window_end >= r_window_start{
                            break 'each_m_window;
                        }
                        //print!("poly base check...");
                        let (m_has_repeat_bool, m_has_repeat_offset) = current_sequence.has_repeat(m_window_start, m_window_end);
                        if m_has_repeat_bool {
                            m_window_start += m_has_repeat_offset + 1;
                            continue 'each_m_window;
                        }
                        //println!("OK");
                        //ここからcounting bloom filterに追加していく。
                        add_bloom_filter_cnt += 1;
                        let probe_candidate: u128 = current_sequence.subsequence_as_u128(vec![[m_window_start, m_window_end]]);
                        let table_indice:  [u32;8] = hash_from_u128(probe_candidate);//u128を受けてhashを返す関数
                        println!("{:?}", String::from_utf8(decode_u128_2_dna_seq(&probe_candidate, PROBE_LEN)).unwrap());
                        for i in 0..8{
                            let idx: usize = table_indice[i] as usize;
                            if ret_array[idx] == u32::MAX{
                                eprintln!("index {} reaches u32::MAX", idx);
                            }else{
                                ret_array[idx] += 1;
                            }
                        }
                        m_window_start += 1;
                    }
                    r_window_start += 1;
                }
                l_window_start += 1;
            }
            let end = start_time.elapsed();
            eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
            previous_time = end;
        }
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
            ret_val[i] += sha256_bit_array[i * 4 + j] as u32;
            ret_val[i] <<= 8;
        }
        ret_val[i] %= BLOOMFILTER_TABLE_SIZE as u32;
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

pub fn number_of_high_occurence_kmer(source_table: &Vec<u32>, sequences: &Vec<DnaSequence>, start_idx: usize, end_idx: usize, threshold: u32, thread_id: usize, primer: &Vec<(Vec<u8>, DnaSequence, DnaSequence)>) -> HashSet<u128>{
    let mut ret_table: HashSet<u128> = HashSet::with_capacity(HASHSET_SIZE);
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut m_window_start: usize;
    let mut m_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let mut primer_l_seq:   u128;
    let mut primer_l_size:  usize;
    let mut primer_r_seq:   u128;
    let mut primer_r_size:  usize;
    let mut mask_l:         u128;
    let mut mask_r:         u128;
    let chunk_max = 200;


    let start_time = Instant::now();
    let mut previous_time = start_time.elapsed();

    let mut loop_cnt:usize = 0;
    'each_primer: for current_primer in primer.iter() {
        primer_l_size = current_primer.1.len();
        primer_l_seq  = current_primer.1.subsequence_as_u128(vec!([0, primer_l_size]));
        primer_r_size = current_primer.2.len();
        primer_r_seq  = current_primer.2.subsequence_as_u128(vec!([0, primer_r_size]));
        mask_l        = u128::MAX >> (64 - primer_l_size) * 2;
        mask_r        = u128::MAX >> (64 - primer_r_size) * 2;
    

        'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
            let mut add_bloom_filter_cnt: usize = 0;
            let mut l_window_cnt: usize         = 0;
            loop_cnt += 1;
            l_window_start = 0;
            if current_sequence.len() < primer_l_size || current_sequence.len() <  primer_r_size ||current_sequence.len() < PROBE_LEN{
                continue 'each_read;
            }

            'each_l_window: loop{
                l_window_end = l_window_start + primer_l_size;
                if l_window_end >= current_sequence.len() + 1{
                    break 'each_l_window;
                }
                l_window_cnt += 1;
                let l_window_as_u128: u128 = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end]]);
                if l_window_as_u128 & mask_l != primer_l_seq{
                    l_window_start += 1;
                    continue 'each_l_window;
                }
                r_window_start = l_window_end + PROBE_LEN;
                'each_r_window: loop{
                    r_window_end = r_window_start + primer_r_size;
                    if r_window_end >= current_sequence.len() + 1{
                        let end = start_time.elapsed();
                        eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
                        previous_time = end;
                        continue 'each_read;
                    }
                    if r_window_end - l_window_start > chunk_max{
                        break 'each_l_window;
                    }
                    let r_window_as_u128: u128 = current_sequence.subsequence_as_u128(vec![[r_window_start, r_window_end]]);

                    if r_window_as_u128 & mask_r != primer_r_seq {
                        r_window_start += 1;
                        continue 'each_r_window;
                    }
                    //ここまでで、LとRが一致してる。この下のループでPROBE Probeを探索する。
                    m_window_start = l_window_end;
                    'each_m_window: loop{
                        m_window_end = m_window_start + PROBE_LEN;
                        if m_window_end >= r_window_start{
                            break 'each_m_window;
                        }
                        let (m_has_repeat_bool, m_has_repeat_offset) = current_sequence.has_repeat(m_window_start, m_window_end);
                        if m_has_repeat_bool {
                            m_window_start += m_has_repeat_offset + 1;
                            continue 'each_m_window;
                        }
                        //ここからcounting bloom filterに追加していく。
                        add_bloom_filter_cnt += 1;
                        let probe_candidate: u128    = current_sequence.subsequence_as_u128(vec![[m_window_start, m_window_end]]);
                        let table_indice:     [u32;8] = hash_from_u128(probe_candidate);//u128を受けてhashを返す関数
                        let occurence:        u32     = count_occurence_from_counting_bloomfilter_table(source_table, table_indice);
                        if occurence >= threshold * DUPPLICATION{
                            ret_table.insert(probe_candidate);
                        }
                        m_window_start += 1;
                    }
                    r_window_start += 1;
                }
                l_window_start += 1;
            }
            let end = start_time.elapsed();
            eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
            previous_time = end;
        }
    }
    return ret_table;
}


pub fn aggregate_length_between_primer(sequences: &Vec<DnaSequence>, thread_id: usize, primer: &Vec<(Vec<u8>, DnaSequence, DnaSequence)>, product_size_max: usize) -> Vec<u8>{
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let mut primer_l_seq:   u128;
    let mut primer_l_size:  usize;
    let mut primer_r_seq:   u128;
    let mut primer_r_size:  usize;
    let mut mask_l:         u128;
    let mut mask_r:         u128;
    let mut primer_id:      Vec<u8>;
    let mut loop_cnt:usize = 0;
    let mut ret_array: Vec<u8> = Vec::new();
    let mut lr_hit_counter:usize = 0;
    let mut l_hit_counter:usize  = 0;
    eprintln!("[{}]primer pairs: {}", thread_id, primer.len());

    let start_time = Instant::now();
    let mut previous_time = start_time.elapsed();
    'each_primer: for current_primer in primer.iter() {
        primer_l_size = current_primer.1.len();
        primer_l_seq  = current_primer.1.subsequence_as_u128(vec!([0, primer_l_size]));
        primer_r_size = current_primer.2.len();
        primer_r_seq  = current_primer.2.subsequence_as_u128(vec!([0, primer_r_size]));
        mask_l        = u128::MAX >> (64 - primer_l_size) * 2;
        mask_r        = u128::MAX >> (64 - primer_r_size) * 2;
        loop_cnt += 1;
        

        'each_read: for current_sequence in sequences.iter() {
            l_window_start = 0;
            'each_l_window: loop{
                l_window_end = l_window_start + primer_l_size;
                if l_window_end >= current_sequence.len() + 1{
                    break 'each_l_window;
                }
                //l_window_cnt += 1;
                let l_window_as_u128: u128 = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end]]);
                if (l_window_as_u128 & mask_l) != primer_l_seq{
                    l_window_start += 1;
                    continue 'each_l_window;
                }
                r_window_start = l_window_end + PROBE_LEN;
                l_hit_counter += 1;
                'each_r_window: loop{
                    r_window_end = r_window_start + primer_r_size;
                    if r_window_end >= current_sequence.len() + 1{
                        let end = start_time.elapsed();
                        //eprintln!("loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos());
                        previous_time = end;
                        continue 'each_read;
                    }
                    if r_window_end - l_window_start > product_size_max{
                        break 'each_l_window;
                    }
                    let r_window_as_u128: u128 = current_sequence.subsequence_as_u128(vec![[r_window_start, r_window_end]]);
                    if (r_window_as_u128 & mask_r) != primer_r_seq {
                        r_window_start += 1;
                        continue 'each_r_window;
                    }
                    //ここまでで、LとRが一致してる
                    let length: u32 = (r_window_end - l_window_start) as u32;
                    //ret_array.push(length);
                    let primer_id          = String::from_utf8(current_primer.0.clone()).unwrap();
                    let primer_id_str      = format!("{}", primer_id);
                    let sequence_slice     = String::from_utf8(current_sequence.decode(l_window_start, r_window_end)).unwrap();
                    let sequence_slice_str = format!("{}", sequence_slice);
                    let result_str = format!(">{}_{}\n{}\n", primer_id_str, r_window_end - l_window_start, sequence_slice_str);
                    let result_bytes = result_str.into_bytes();
                    // Add the bytes to ret_array
                    ret_array.extend(result_bytes);
                    lr_hit_counter += 1;
                    r_window_start += 1;
                }
                l_window_start += 1;
            }
        }
        let end = start_time.elapsed();
        
        eprintln!("loop[{:02?}]: {:06?}\t{:09?}\t{}\t{}\tsec: {}.{:03}",thread_id, primer.len(), loop_cnt, lr_hit_counter, l_hit_counter, end.as_secs() - previous_time.as_secs(), end.subsec_nanos() - previous_time.subsec_nanos());
        previous_time = end;
        lr_hit_counter = 0;
        l_hit_counter = 0;
    }
    return ret_array;
}
