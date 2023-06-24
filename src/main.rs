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

fn main() {
    use bio::io::fastq;
    use std::fs::File;

    // let file = File::open("1_control_psbA3_2019_minq7.fastq").unwrap();
    let file = File::open("../../kmer_dataset/SRR1513870.fastq").unwrap();
    let reader = fastq::Reader::new(file);

    let mut nb_reads = 0;
    let mut nb_bases = 0;
    let mut strip = 0;
    let mut total = 0usize;
    for result in reader.records() {
        let record = result.expect("Error during fastq record parsing");
        for c in record.seq() {
            if ENCODE_MAP[*c as usize] as i8 >= 0 {
                strip += 1;
            } else {
                if strip >= 32 {
                    total += strip - 31;
                }
                strip = 0;
            }
        }
        nb_reads += 1;
        nb_bases += record.seq().len();
        if strip >= 32 {
            total += strip - 32;
        }
        strip = 0;
    }

    let total_size = total * 8;
    println!("Number of total 32 mer: {}; size: {} M", total, total_size / (1024 * 1024));
    println!("Number of reads: {}", nb_reads);
    println!("Number of bases: {}", nb_bases);
}
