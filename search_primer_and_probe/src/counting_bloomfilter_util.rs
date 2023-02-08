pub const L_LEN: usize = 32;
pub const M_LEN: usize = 32;
pub const R_LEN: usize = 32;
pub const HASHSET_SIZE: usize = (u32::MAX >> 4) as usize;
pub const BLOOMFILTER_TABLE_SIZE: usize = (u32::MAX >> 1) as usize;
//const length: usize = 141;
const DUPPLICATION: u32 = 1;
use crate::sequence_encoder_util::{DnaSequence, LmrTuple};

use std::time::{Instant};
use std::collections::HashSet;


//全てのL, M, Rと、hash値を出力する
//部分配列のdecoderを書き、テストする
pub fn build_counting_bloom_filter(sequences: &Vec<DnaSequence>, start_idx: usize, end_idx: usize, length: usize, thread_id: usize) -> Vec<u32>{
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut m_window_start: usize;
    let mut m_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;

    let mut loop_cnt:usize = 0;
    eprintln!("Allocating Vec<u32> where BLOOMFILTER_TABLE_SIZE = {}", BLOOMFILTER_TABLE_SIZE);
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
        if current_sequence.len() < L_LEN + M_LEN + R_LEN{
            continue 'each_read;
        }
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len() + 1{
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let (l_has_repeat_bool, l_has_repeat_offset) = current_sequence.has_repeat(l_window_start, l_window_end);
            if l_has_repeat_bool {
                l_window_start += l_has_repeat_offset + 1;
                continue 'each_l_window;
            }
            m_window_start = l_window_end;
            'each_m_window: loop{
                m_window_end = m_window_start + M_LEN;
                if m_window_end >= current_sequence.len() + 1{
                    break 'each_m_window;
                }
                if m_window_end - l_window_start >= length - R_LEN{
                    break 'each_m_window;
                }
                let (m_has_repeat_bool, m_has_repeat_offset) = current_sequence.has_repeat(m_window_start, m_window_end);
                if m_has_repeat_bool {
                    m_window_start += m_has_repeat_offset + 1;
                    continue 'each_m_window;
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
                    if r_window_end - l_window_start >= length{
                        break 'each_r_window;
                    }
                    let (r_has_repeat_bool, r_has_repeat_offset) = current_sequence.has_repeat(r_window_start, r_window_end);
                    if r_has_repeat_bool {
                        r_window_start += r_has_repeat_offset + 1;
                        continue 'each_r_window;
                    }
                    add_bloom_filter_cnt += 1;
                    let lmr_string: LmrTuple = current_sequence.subsequence_as_lmrtuple([[l_window_start, l_window_end], [r_window_start, r_window_end], [r_window_start, r_window_end]]);
                    let table_indice:[u32;8] = lmr_string.hash();
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
                m_window_start += 1;
            }
            l_window_start += 1;
        }
        let end = start_time.elapsed();
        eprintln!("1st loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
        previous_time = end;
    }
    return ret_array;
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

pub fn number_of_high_occurence_lmr_tuple(source_table: &Vec<u32>, sequences: &Vec<DnaSequence>, start_idx: usize, end_idx: usize, threshold: u32, length: usize, thread_id: usize) -> HashSet<LmrTuple>{
    let mut ret_table: HashSet<LmrTuple> = HashSet::with_capacity(HASHSET_SIZE);
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut m_window_start: usize;
    let mut m_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let mut ho_lmr: usize = 0;

    let start = Instant::now();
    let mut previous_time = start.elapsed();
    let mut loop_cnt:usize = 0;
    'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
        let mut add_bloom_filter_cnt: usize = 0;
        let mut l_window_cnt: usize         = 0;
        loop_cnt += 1;
        l_window_start = 0;
        if current_sequence.len() < L_LEN + M_LEN + R_LEN{
            continue 'each_read;
        }
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len() + 1{
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let (l_has_repeat_bool, l_has_repeat_offset) = current_sequence.has_repeat(l_window_start, l_window_end);
            if l_has_repeat_bool {
                l_window_start += l_has_repeat_offset + 1;
                continue 'each_l_window;
            }
            m_window_start = l_window_end;
            'each_m_window: loop{
                m_window_end = m_window_start + M_LEN;
                if m_window_end >= current_sequence.len() + 1{
                    let end = start.elapsed();
                    eprintln!("2nd loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}\tho_lmr: {}", thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt, ho_lmr);
                    previous_time = end;
                    continue 'each_read;
                }
                if m_window_end - l_window_start >= length - R_LEN{
                    break 'each_m_window;
                }
                let (m_has_repeat_bool, m_has_repeat_offset) = current_sequence.has_repeat(m_window_start, m_window_end);
                if m_has_repeat_bool {
                    m_window_start += m_has_repeat_offset + 1;
                    continue 'each_m_window;
                }
                r_window_start = m_window_end;
                'each_r_window: loop{
                    r_window_end = r_window_start + R_LEN;
                    if r_window_end >= current_sequence.len() + 1{
                        let end = start.elapsed();
                        eprintln!("2nd loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}\tho_lmr: {}", thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt, ho_lmr);
                        previous_time = end;
                        continue 'each_read;
                    }
                    if r_window_end - l_window_start > length{
                        break 'each_r_window;
                    }
                    let (r_has_repeat_bool, r_has_repeat_offset) = current_sequence.has_repeat(r_window_start, r_window_end);
                    if r_has_repeat_bool {
                        r_window_start += r_has_repeat_offset + 1;
                        continue 'each_r_window;
                    }
                    add_bloom_filter_cnt += 1;
                    let lmr_string: LmrTuple = current_sequence.subsequence_as_lmrtuple([[l_window_start, l_window_end], [r_window_start, r_window_end], [r_window_start, r_window_end]]);
                    let table_indice:[u32;8] = lmr_string.hash();
                    let occurence: u32 = count_occurence_from_counting_bloomfilter_table(source_table, table_indice);
                    if occurence >= threshold * DUPPLICATION{
                        if ret_table.len() < HASHSET_SIZE{
                            ret_table.insert(lmr_string);
                            ho_lmr += 1;
                        }else{
                            panic!("reached to the maximum size of has table");
                        }
                    }
                    r_window_start += 1;
                }
                m_window_start += 1;
            }
            l_window_start += 1;
        }
        let end = start.elapsed();
        eprintln!("2nd loop[{:02?}]({}-{}, length is {}): {:09?}\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}\tho_lmr: {}", thread_id, start_idx, end_idx, end_idx - start_idx, loop_cnt, current_sequence.len(), end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt, ho_lmr);
        previous_time = end;
    }
    return ret_table;
}
