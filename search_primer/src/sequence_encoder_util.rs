use crate::counting_bloomfilter_util::L_LEN;
use crate::counting_bloomfilter_util::R_LEN;
//use crate::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use std::cmp;
//use std::hash::Hash;



pub fn decode_u128_2_dna_seq(source:&u128, char_size: usize) -> Vec<u8>{
    let mut result: Vec<u8> = Vec::new();
    let mut base;
    for i in 0..char_size{
        base = source >> 2 * (char_size - 1 - i) & 3;
        match base{
            0 => {result.push(b'A');}
            1 => {result.push(b'C');}
            2 => {result.push(b'G');}
            3 => {result.push(b'T');}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}

pub fn decode_u128_l(source: &u128) -> [u8; L_LEN]{
    let mut result: [u8; L_LEN] = [b'X'; L_LEN];
    let mut base;
    for i in 0..L_LEN{
        base = source >> (((L_LEN + R_LEN) - i - 1) * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}


pub fn decode_u128_r(source: &u128) -> [u8; R_LEN]{
    let mut result: [u8; R_LEN] = [b'X'; R_LEN];
    let mut base;
    for i in 0..R_LEN{
        base = source >> ((R_LEN - i - 1) * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}




pub struct DnaSequence{
    length:   usize,
    sequence: Vec<u64>
}

impl DnaSequence{
    pub fn new(source: &Vec<u8>) -> DnaSequence{
        let mut retval:Vec<u64> = Vec::new();
        let mut buf: u64 = 0;
        let mut cnt: u64 = 0;
        for each_base in source.iter(){
            match each_base{
                b'A' => {buf |= 0;}
                b'C' => {buf |= 1;}
                b'G' => {buf |= 2;}
                b'T' => {buf |= 3;}
                _   => {panic!("Unexpected character: {}", each_base);}
            }
            //buf <<= 2;//?????????????????????????????????????????????????????????
            if cnt == 31{//cnt == 31????????????buf???????????????????????????
                retval.push(buf);
                buf = 0;
                cnt = 0;
            }else{
                buf <<= 2;
                cnt += 1;
            }
        }

        if cnt != 32{
            buf <<= 2 * (31 - cnt);//???????????????bit????????????????????????
            retval.push(buf);//32???????????????????????????????????????????????????push?????????
        }
/*
        for each_buf in retval.iter(){
            println!("{:064b}", each_buf);
        }
*/
        return DnaSequence{sequence: retval, length: source.len()};
    }

    pub fn encode(source: &Vec<u8>) -> DnaSequence{
        return DnaSequence::new(source);
    }

    pub fn len(&self) -> usize{
        return self.length;
    }

    pub fn decode(&self, start: usize, end: usize) -> Vec<u8>{
        assert!(start < end, "DnaSequence::decode assertion failed: {} !< {}", start, end);
        assert!(end <= self.length, "DnaSequence::decode assertion failed: {} !< {}", end, self.length);
        let mut retval = Vec::new();
        let mut buf: u8;
        for i in start..end{
            buf = ((self.sequence[i / 32] >> (2 * (31 - i % 32))) & 3).try_into().unwrap();
            match buf{
                0 => {retval.push('A' as u8);}
                1 => {retval.push('C' as u8);}
                2 => {retval.push('G' as u8);}
                3 => {retval.push('T' as u8);}
                _ => {panic!("Never reached!!!buf: {}", buf);}
            }
        }
        return retval;
    }

    pub fn subsequence(&self, ranges: Vec<[usize; 2]>) -> DnaSequence{
        let mut buf: u64 = 0;
        let mut cnt: u64 = 0;
        let mut retval: Vec<u64> = Vec::new();
        let mut length: usize = 0;
        for each_range in ranges.iter(){
            let start: usize = each_range[0];
            let end:   usize = each_range[1];
            assert!(start < end, "DnaSequence::subsequence assertion failed: {} !< {}", start, end);
            //assert!(start >= 0, "DnaSequence::subsequence assertion failed: {} >= 0", start);
            assert!(end < self.length, "DnaSequence::subsequence assertion failed: {} < {}", end, self.length);
            length += end - start;
            for i in start..end{
                buf += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
                buf <<= 2;
                cnt += 1;
                if cnt == 31{
                    retval.push(buf);
                    buf = 0;
                    cnt = 0;
                }
            }
        }
        if cnt != 31{
            buf <<= 2 * (31 - cnt);//???????????????bit????????????????????????
            retval.push(buf);//32???????????????????????????????????????????????????push?????????
        }
        return DnaSequence{sequence: retval, length: length}
    }

//subsequence_as_u128????????????
    pub fn subsequence_as_u128(&self, ranges: Vec<[usize; 2]>) -> u128{
        let mut buf: u128 = 0;
        let mut cnt: usize = 0;
        for each_range in ranges.iter(){
            let start: usize = each_range[0];
            let end:   usize = each_range[1];
            assert!(start < end, "DnaSequence::subsequence_as_u128 assertion failed: {} !< {}", start, end);
            assert!(end <= self.length, "DnaSequence::subsequence_as_u128 assertion failed: {} < {}", end, self.length);
            for i in start..end{
                buf <<= 2;
                cnt += 1;
                buf += ((self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3) as u128;
            }
            assert!(cnt <= 64, "DnaSequence::subsequence_as_u128 assertion failed: too many DNA size. cnt reaches {}", cnt);
        }
        return buf
    }

    pub fn has_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        let has_one_base_repeat: (bool, usize)   = self.has_one_base_repeat(start, end);
        let has_two_base_repeat: (bool, usize)   = self.has_two_base_repeat(start, end);
        let has_three_base_repeat: (bool, usize) = self.has_three_base_repeat(start, end);
        let retval_bool: bool  = has_one_base_repeat.0 | has_two_base_repeat.0 | has_three_base_repeat.0;
        let retval_base: usize = cmp::max(has_one_base_repeat.1, cmp::max(has_two_base_repeat.1, has_three_base_repeat.1));
        return (retval_bool, retval_base)
    }

    pub fn has_one_base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(start < end, "DnaSequence::has_one_base_repeat assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_one_base_repeat assertion failed: {} - {} < 4", end, start);
        assert!(end - start <= 32, "DnaSequence::has_one_base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_one_base_repeat assertion failed: end coordinate must be smaller than length of the sequence. start: {}, end: {}, self.lngth: {}", start, end, self.length);
        let zero_ichi: u64  = 0x5555555555555554;

        let upper_mask:u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let mut left_bits: u128 = self.sequence[start/32] as u128;
        left_bits <<= 64 - 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 64 - 2 * (start % 32);
        let mut original = u64::try_from(left_bits).unwrap();
        original <<= 2 * (end % 32);
        original += self.sequence[end / 32] >> (64 - 2 * (end % 32));
        let val1  = original;
        let val2  = original << 2;
        let val3  = val1 ^ val2;
        let val4  = val3 >> 1;
        let val5  = val3 | val4;
        let val6  = !val5;
        let val7  = val6 & zero_ichi;//??????1bit??????0????????????
        let val8  = val7 << 2;
        let val9  = val7 << 4;
        let val10 = val7 & val8 & val9;
        let val11 = (val10 << (2 * (32 + start - end))).leading_zeros() / 2;

        #[cfg(test)]{
            println!("has_one_base_repeat");
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start:    {}", start);
            println!("end:      {}", end);
            println!(" 0101: {:064b}", zero_ichi);
            println!(" val1: {:064b}", val1);
            println!(" val2: {:064b}", val2);
            println!(" val3: {:064b}", val3);
            println!(" val4: {:064b}", val4);
            println!(" val5: {:064b}", val5);
            println!(" val6: {:064b}", val6);
            println!(" val7: {:064b}", val7);
            println!(" val8: {:064b}", val8);
            println!(" val9: {:064b}", val9);
            println!("val10: {:064b}", val10);
            println!("val11: {}", val11);
        }
        if val11 == 32{
            return (false, 0)
        }else{
            return (true, val11.try_into().unwrap())
        }
        //shift???????????????????????????????????????????????????
    }


    pub fn has_two_base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(start < end, "DnaSequence::has_two_base_repeat assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_two_base_repeat assertion failed: {} - {} < 4", end, start);
        assert!(end - start <= 32, "DnaSequence::has_two_base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_two_base_repeat assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let zero_ichi: u64 = 0x5555_5555_5555_5555 & !63;
        let upper_mask:u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let mut left_bits: u128 = self.sequence[start/32] as u128;
        left_bits <<= 64 - 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 64 - 2 * (start % 32);
        let mut original = u64::try_from(left_bits).unwrap();
        original <<= 2 * (end % 32);
        original += self.sequence[end / 32] >> (64 - 2 * (end % 32));
        //??????????????????original???????????????????????????????????????????????????
        let val1 = original;
        let val2 = original << 4;
        let val3 = val1 ^ val2 & !15;
        let val4 = val3 >> 1;
        let val5 = val3 | val4;
        let val6 = !val5;
        let val7 = val6 & zero_ichi;
        let last  = val7 & 
                    val7 << 2 &
                    val7 << 4 &
                    val7 << 6 &
                    val7 << 8 &
                    val7 << 10 &
                    val7 << 12 &
                    val7 << 14 ;
        let leading0 = (last << (2 * (32 + start - end))).leading_zeros() / 2;
        #[cfg(test)]{
            println!("has_two_base_repeat");
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start: {}", start);
            println!("end:   {}", end);
            println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
            println!("val4:  {:064b}", val4);
            println!("val5:  {:064b}", val5);
            println!("val6:  {:064b}", val6);
            println!("val7:  {:064b}", val7);
            println!("last:  {:064b}", last);
            println!("leading0: {}", leading0);
        }
        if leading0 == 32{
            return (false, 0)
        }else{
            return (true, leading0.try_into().unwrap())
        }
    }

    pub fn has_three_base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(start < end, "DnaSequence::has_three_base_repeat assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_three_base_repeat assertion failed: {} - {} < 4", end, start);
        assert!(end - start <= 32, "DnaSequence::has_three_base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_three_base_repeat assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let zero_ichi: u64  = 0x5555555555555555 & !63;
        let upper_mask:u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let mut left_bits: u128 = self.sequence[start/32] as u128;
        left_bits <<= 64 - 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 64 - 2 * (start % 32);
        let mut original = u64::try_from(left_bits).unwrap();
        original <<= 2 * (end % 32);
        original += self.sequence[end / 32] >> (64 - 2 * (end % 32));
        
        //??????????????????original???????????????????????????????????????????????????
        let val1 = original;
        let val2 = original << 6;
        let val3 = val1 ^ val2;
        let val4 = val3 >> 1;
        let val5 = val3 | val4;
        let val6 = !val5;
        let val7 = val6 & zero_ichi;
        let val8 = val7 << 2;
        let val9 = val7 << 4;
        let val10 = val7 << 6;
        let val11 = val7 << 8;
        let val12 = val7 << 10;
        let last  = val7 & val8 & val9 & val10 & val11 & val12;
        let leading0 = (last << (2 * (32 + start - end))).leading_zeros() / 2;

        #[cfg(test)]{
            println!("has_three_base_repeat");
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start: {}", start);
            println!("end:   {}", end);
            println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
            println!("val4:  {:064b}", val4);
            println!("val5:  {:064b}", val5);
            println!("val6:  {:064b}", val6);
            println!("val7:  {:064b}", val7);
            println!("val8:  {:064b}", val8);
            println!("val9:  {:064b}", val9);
            println!("val10: {:064b}", val10);
            println!("val11: {:064b}", val11);
            println!("val12: {:064b}", val12);
            println!("last:  {:064b}", last);
        }
        if leading0 == 32{
            return (false, 0)
        }else{
            return (true, leading0.try_into().unwrap())
        }
    }
}


#[cfg(test)]
mod tests{
    use crate::sequence_encoder_util::DnaSequence;
    use ::function_name::named;
/*
*
*Encode Test
*
*/
    #[test]
    #[named]
    fn encode_test_4A(){
        let source: Vec<u8> = vec![b'A', b'A', b'A', b'A'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x0000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_4C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x5500000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_4G(){
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0xAA00000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_4T(){
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0xFF00000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_ACGT(){
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x1b00000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_16C(){
        let source: String = "CCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555500000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_31C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555554, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_32C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_33C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x4000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_34C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_63C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5555555555555554, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_64C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5555555555555555, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_65C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[2] == 0x4000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_32N(){
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0b1000000110000111101111110111000111001100000011010111110101111101, "{} failed", function_name!());
    }//                               G A A C G A C T G T T T C T A C T A T A A A T C C T T C C T T C
    #[test]
    #[named]
    fn encode_test_120N(){
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x8187bf71cc0d7d7d, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x725cd3f7a2d7eb81, "{} failed", function_name!());
        assert!(obj.sequence[2] == 0xeca09de04446f57e, "{} failed", function_name!());
        assert!(obj.sequence[3] == 0x8f6c5ce0c75b0000, "{} failed", function_name!());
    }

/*
*
*has_poly_base_test
*4?????????true
*
*/

#[test]
    #[named]
    fn has_poly_base_test_8C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (true, 0),  "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_8G(){
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (true, 0),  "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_8T(){
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (true, 0),  "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_8N(){
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_120N(){
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 9)   == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 19)  == (true, 9),  "{} failed", function_name!());
        assert!(obj.has_poly_base(9, 20)  == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(28, 39) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_1(){
        let source: String = "ATTCATACTTAATACTGTATCAGTTGA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_2(){
        let source: String = "TTCATACTTAATACTGTATCAGTTGAG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_3(){
        let source: String = "TCATACTTAATACTGTATCAGTTGAGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_4(){
        let source: String = "TTCGAAAATCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (true, 4), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_5(){
        let source: String = "TTCGAAATTCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_6(){
        let source: String = "TTCGTAATTCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_7(){
        let source: String = "AAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (true, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_8(){
        let source: String = "AAAACGTCGTCGTCGCATACGATCGAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (true, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_9(){
        let source: String = "AAACGTCGTCGTCGCATACGAATCGAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_10(){
        let source: String = "GTCGTCGTCGCATACGAATCGATAAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (true, 23), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_11(){
        let source: String = "GTCGTCGTCGCATACGAATCGATTAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 23) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_12(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 30) == (true, 23), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 19) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_13(){
        let source: String = "TAATGTATTCACTAACAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 19) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_14(){
        let source: String = "TAATGTATTCACTAACCAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 19) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_27N_15(){
        let source: String = "CTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 19) == (true, 10), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_16(){
        let source: String = "GCTTGTATACAGGGGATTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 19) == (true, 11), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_27N_17(){
        let source: String = "CTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 19) == (true, 10), "{} failed", function_name!());
    }


    #[test]
    #[named]
    fn has_poly_base_test_27N_18(){
        let source: String = "AACGTACGTACGTACGTACGTACGTACGTACGTCTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(30, 52) == (true, 13), "{} failed", function_name!());
    }


/*
*
*has_simple_repeat
*
*/
    #[test]
    #[named]
    fn has_simple_repeat_27N_1(){
        let source: String = "TCATATGCAACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 27) == (true, 7), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_simple_repeat_27N_1_a(){
        let source: String = "ACGTACGTCAACAAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 15) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_simple_repeat_27N_2(){
        let source: String = "TCATATGCTACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_simple_repeat_27N_3(){
        let source: String = "CAACAACTGC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 10) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_simple_repeat_27N_4(){
        let source: String = "ACGTACGTCAACAAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 15) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_simple_repeat_27N_5(){
        let source: String = "CGTACGTCAACAACA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 15) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_simple_repeat_27N_6(){
        let source: String = "GTACGTCAACAACAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 15) == (true, 6), "{} failed", function_name!());
    }
/*
*
*has_2base_repeat
*
*/
    #[test]
    #[named]
    fn has_2base_repeat_27N_1(){
        let source: String = "TCATATGCTACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2(){
        let source: String = "CGTACGCTTATATATATATATACCGCA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 8), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2x6_1(){
        let source: String = "TATATATATATAGCCCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2x6_2(){
        let source: String = "GCCCGCACGTACGCTATATATATATAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 14), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2x5_1(){
        let source: String = "TATATATATACTGCCCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2x5_2(){
        let source: String = "CTGCCCGCACGTACGCTATATATATAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 16), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2x4_1(){
        let source: String = "TATATATACTGCCCGAACACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_2x4_2(){
        let source: String = "CTGCCCGAACACGTACGCTATATATAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_3(){
        let source: String = "GAATTCGTTACGTAACGACGCGCGCGC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_5(){
        let source: String = "ACGTTATACCTGTACCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_6(){
        let source: String = "GCGGTATATACGTACCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_2base_repeat_27N_7(){
        let source: String = "CGGTATATACGTACCGCACGCACACAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }
/*

has_poly_base_or_simple_repeat test
*/
    #[test]
    #[named]
    fn has_poly_base_or_simple_repeat_1(){
        let source: String = "CTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base_or_simple_repeat(0, 19) == true, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_or_simple_repeat_2(){
        let source: String = "CATCACCAATTATTGGTCCTAATGTA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base_or_simple_repeat(0, 26) == false, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_or_simple_repeat_3(){
        let source: String = "ATCCTCAGCTGCTTGTATA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base_or_simple_repeat(0, 19) == false, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_or_simple_repeat_4(){
        for string in ["ATCCTCAGCTGCTTGTATA"]{
            let source: String = string.to_string();
            let v: Vec<u8> = source.into_bytes();
            let obj = DnaSequence::new(&v);
            assert!(obj.has_poly_base_or_simple_repeat(0, 19) == false, "{} failed", function_name!());
        }
    }

/*
*
*Decode test
*
*/
    #[test]
    #[named]
    fn decode_test_8C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.decode(0, 4) == vec![67, 67, 67, 67            ], "{} failed", function_name!());
        assert!(obj.decode(0, 5) == vec![67, 67, 67, 67, 67        ], "{} failed", function_name!());
        assert!(obj.decode(0, 6) == vec![67, 67, 67, 67, 67, 67    ], "{} failed", function_name!());
        assert!(obj.decode(0, 7) == vec![67, 67, 67, 67, 67, 67, 67], "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn decode_test_120N(){
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.decode(0, 120) == v, "{} failed", function_name!());
    }





/*
*
*Subsequence Test
*
*/
    #[test]
    #[named]
    fn subsequence_as_u128_test_64G(){
        let source: String = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 64 as usize]]);
        let expectedvalue: u128 = 0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_10G(){
        let source: String = "GGGGGGGGGG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 10 as usize]]);
        let expectedvalue: u128 = 0b10101010101010101010;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_12N(){
        let source: String = "ACGTAACCGGTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 12 as usize]]);
        let expectedvalue: u128 = 0b000110110000010110101111;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }


    #[test]
    #[named]
    fn subsequence_as_u128_test_1(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }
    

    #[test]
    #[named]
    fn subsequence_as_u128_test_2(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 19], [19, 45], [45, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_3(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 27], [64, 84]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGAAAAAAAAAAAAAAAAAAAA".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 47]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }
/*
    #[test]
    #[named]
    fn subsequence_as_u128_test_4(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn subsequence_as_u128_test_5(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn subsequence_as_u128_test_6(){
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }
 */
}

