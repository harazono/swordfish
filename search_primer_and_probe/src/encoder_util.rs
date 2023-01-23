

pub fn encode_dna_seq_2_u64(sequence: &[u8]) -> u64{
    let mut result: u64 = 0;
    for each_base in sequence.iter(){
        match each_base{
            b'A' => {result |= 0;}
            b'C' => {result |= 1;}
            b'G' => {result |= 2;}
            b'T' => {result |= 3;}
            _   => {panic!("Unexpected character: {}", each_base);}
        }
        result = result << 2;
    }
    result = result >> 2;
    return result;
}

pub fn encode_dna_seq_2_u128(sequence: &[u8]) -> u128{
    let mut result: u128 = 0;
    for each_base in sequence.iter(){
        match each_base{
            b'A' => {result |= 0;}
            b'C' => {result |= 1;}
            b'G' => {result |= 2;}
            b'T' => {result |= 3;}
            _   => {panic!("Unexpected character: {}", each_base);}
        }
        result = result << 2;
    }
    result = result >> 2;
    return result;
}

pub fn encode_dna_seq_2_u8_vec(sequence: &[u8]) -> Vec<u8>{
    let mut result: Vec<u8> = Vec::new();
    for base in 0..(sequence.len() as usize){
        match base{
            0 => {result.push('A' as u8);}
            1 => {result.push('C' as u8);}
            2 => {result.push('G' as u8);}
            3 => {result.push('T' as u8);}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;

}



pub fn decode_u64_2_dna_seq(source:u64, index: usize, length: usize) ->u8{
    let mut result: u8 = 0;
    let mut tmp: u64 = source << (64 - length * 2);
    tmp = tmp << index * 2;
    tmp = tmp >> 62;
    match tmp{
        0 => {result = b'A';}
        1 => {result = b'C';}
        2 => {result = b'G';}
        3 => {result = b'T';}
        _ => {panic!("Never reached!!!tmp: {}", tmp);}
    }
    return result;
}

pub fn decode_u128_2_dna_seq(source:&u128) -> [u8; L_LEN + R_LEN]{
    let mut result: [u8; L_LEN + R_LEN] = [b'X'; L_LEN + R_LEN];
    let mut tmp: u128 = source.clone() << 20;
    let mut base: u128;
    for i in 0..L_LEN + R_LEN{
        base = (tmp & 0xC000_0000_0000_0000_0000_0000_0000_0000) >> 126;
        //println!("{}, {:#130b}", base, tmp);
        tmp = tmp << 2;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!tmp: {}", tmp);}
        }
    }
    return result;
}
pub fn decode_u128_2_occurence(source: &u128) -> u32{
    return TryFrom::try_from(source >> (L_LEN + R_LEN) * 2).unwrap();
}
