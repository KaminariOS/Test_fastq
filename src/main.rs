use bio::io::fastq;
use bio::utils::TextSlice;
use std::fs::File;
use std::slice::Iter;
use std::collections::HashMap;

mod cli;
use cli::Args;
use clap::Parser;

const DECODE_MAP: [char; 4] = ['A', 'C', 'G', 'T'];

#[derive(Clone, Copy)]
enum Code {
    R = -1,
    I = -2,
    O = -3,
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

use Code::*;

const ENCODE_MAP: [Code; 256] = [
    O, O, O, O, O, O, O, O, O, O, I, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, R, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O,
    O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
    O, A, R, C, R, O, O, G, R, O, O, R, O, R, R, O,
    O, O, R, R, T, O, R, R, R, R, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
    O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O
];

struct KmerIterator<'s, const K: usize> {
    buffer: u64,
    seq: Iter<'s, u8>,
    filled: usize
}

impl<'s, const K: usize> KmerIterator<'s, K> {
   pub fn new(seq: TextSlice<'s>) -> Self {
       Self {
           buffer: 0,
           seq: seq.iter(),
           filled: 0
       }
   }

   fn refill_buffer(&mut self) -> bool {
    while self.filled < K {
        let code;
        if let Some(&c) = self.seq.next() {
            code = c;
        } else {
            return false;
        }
        let encode = ENCODE_MAP[code as usize] as i8;
        if encode < 0 {
            self.filled = 0;
        } else {
            self.buffer <<= 2;
            self.buffer |= encode as u64; 
            self.filled += 1;
        }
    }
    // self.buffer &= (0 << (2 * K));
    true
   }
}

type Kmer = u64;

impl<const K: usize> Iterator for KmerIterator<'_, K> {
    type Item = Kmer;

    fn next(&mut self) -> Option<Self::Item> {
        let filled = self.refill_buffer();
        if filled {
            let val = self.buffer;
            self.filled -= 1;
            return Some(val);
        }
        None
    }
}

fn main() {
    let args = Args::parse();
    // let file = File::open("1_control_psbA3_2019_minq7.fastq").unwrap();
    let file = File::open(args.filename).unwrap();
    let reader = fastq::Reader::new(file);

    let mut nb_reads = 0;
    let mut nb_bases = 0;
    let mut strip = 0;
    let mut total = 0usize;
    let mut kmer_count = 0;
    // let mut map = HashMap::with_capacity(1 << 20);
    for result in reader.records() {
        let record = result.expect("Error during fastq record parsing");
        let seq = record.seq();
        let kmer_iterator = 
            KmerIterator::<32>::new(seq);
        for c in seq {
            if ENCODE_MAP[*c as usize] as i8 >= 0 {
                strip += 1;
            } else {
                if strip >= 32 {
                    total += strip - 31;
                }
                strip = 0;
            }
        }
        for kmer in kmer_iterator {
            kmer_count += 1;
            // *map.entry(kmer).or_insert(0) += 1usize;
        }
        nb_reads += 1;
        nb_bases += seq.len();
        if strip >= 32 {
            total += strip - 31;
        }
        strip = 0;
    }
    assert_eq!(kmer_count, total);
    let total_size = total * 8;
    println!("Number of total 32 mer: {}; size: {} M", total, total_size / (1024 * 1024));
    println!("Number of reads: {}", nb_reads);
    println!("Number of bases: {}", nb_bases);

    let file = File::create("map.json").unwrap();
    // serde_json::to_writer(file, &map).unwrap();
}
