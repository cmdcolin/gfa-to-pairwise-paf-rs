use std::env;
use std::process;

use crate::nodes::parse_int;

pub struct Args {
    pub gfa: Option<String>,
    pub reference: Vec<u8>,
    pub queries: Option<Vec<Vec<u8>>>,
    pub max_gap: i64,
    pub min_block: i64,
    pub chrom_sizes_dir: Option<String>,
    pub pair_x: bool,
    pub hold_queries: bool,
    pub contig_lengths: Option<String>,
}

const USAGE: &str =
    "usage: gfa-to-pairwise-paf [-h] [--version] --reference REFERENCE [--queries QUERIES]
                           [--max-gap MAX_GAP] [--min-block MIN_BLOCK]
                           [--chrom-sizes-dir CHROM_SIZES_DIR] [--no-x]
                           [--hold-queries] [--contig-lengths CONTIG_LENGTHS]
                           [gfa]";

const HELP: &str = "
Unpack pairwise alignments out of a pangenome graph's GFA: one PAF record per
chain of graph nodes a query haplotype shares with the reference path, with a
cg:Z: CIGAR over =/X/I/D. An X says the graph put different sequence between
two anchors; no base is compared. https://github.com/cmdcolin/gfa-to-pairwise-paf-rs

positional arguments:
  gfa                   GFA file, gz accepted; default stdin (pigz -dc file.gfa.gz | ...
                        is faster for big files)

options:
  -h, --help            show this help message and exit
  --version             show the version and exit
  --reference REFERENCE
                        reference sample#hap, e.g. GRCh38#0 or K12#0 (a bare sample
                        means #0)
  --queries QUERIES     comma-separated query sample#hap list; default every other
                        sample#hap seen
  --max-gap MAX_GAP     a chain may skip up to this many private bp on either side
                        between two anchors (default 10000)
  --min-block MIN_BLOCK
                        drop records spanning fewer reference bp than this
  --chrom-sizes-dir CHROM_SIZES_DIR
                        write <sample>.<hap>.chrom.sizes per query here
  --no-x                write private runs as I then D instead of pairing them as X
  --hold-queries        align every query walk after the whole file is read, for a
                        file whose reference walks come after query walks that share
                        their nodes
  --contig-lengths CONTIG_LENGTHS
                        chrom.sizes or .fai giving exact query contig lengths, keyed
                        by contig or sample#hap#contig";

const OPTIONS: [&str; 10] = [
    "--help",
    "--version",
    "--reference",
    "--queries",
    "--max-gap",
    "--min-block",
    "--chrom-sizes-dir",
    "--no-x",
    "--hold-queries",
    "--contig-lengths",
];

fn fail(message: &str) -> ! {
    eprintln!("{USAGE}\ngfa-to-pairwise-paf: error: {message}");
    process::exit(2);
}

// argparse's reading of sample names: `HG01109.1` means HG01109#1, and a bare
// sample means haplotype 0.
pub fn pansn_of_prefix(prefix: &str) -> Vec<u8> {
    if prefix.contains('#') {
        return prefix.as_bytes().to_vec();
    }
    match prefix.rsplit_once('.') {
        Some((sample, hap)) if !hap.is_empty() && hap.bytes().all(|b| b.is_ascii_digit()) => {
            format!("{sample}#{hap}").into_bytes()
        }
        _ => format!("{prefix}#0").into_bytes(),
    }
}

fn full_name(arg: &str) -> &'static str {
    if let Some(&exact) = OPTIONS.iter().find(|&&o| o == arg) {
        return exact;
    }
    let matches: Vec<_> = OPTIONS.iter().filter(|o| o.starts_with(arg)).collect();
    match matches.as_slice() {
        [one] => one,
        [] => fail(&format!("unrecognized arguments: {arg}")),
        many => fail(&format!(
            "ambiguous option: {arg} could match {}",
            many.iter().map(|o| **o).collect::<Vec<_>>().join(", ")
        )),
    }
}

fn looks_negative(arg: &str) -> bool {
    arg.strip_prefix('-')
        .is_some_and(|rest| !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()))
}

fn integer(option: &str, value: &str) -> i64 {
    parse_int(value.as_bytes())
        .unwrap_or_else(|_| fail(&format!("argument {option}: invalid int value: '{value}'")))
}

pub fn parse() -> Args {
    let mut args = Args {
        gfa: None,
        reference: Vec::new(),
        queries: None,
        max_gap: 10000,
        min_block: 0,
        chrom_sizes_dir: None,
        pair_x: true,
        hold_queries: false,
        contig_lengths: None,
    };
    let mut reference = None;
    let mut iter = env::args().skip(1);
    let mut positional_only = false;
    while let Some(arg) = iter.next() {
        if positional_only || !arg.starts_with('-') || arg == "-" || looks_negative(&arg) {
            if args.gfa.is_some() {
                fail(&format!("unrecognized arguments: {arg}"));
            }
            args.gfa = Some(arg);
            continue;
        }
        if arg == "--" {
            positional_only = true;
            continue;
        }
        if arg == "-h" {
            println!("{USAGE}\n{HELP}");
            process::exit(0);
        }
        let (name, inline) = match arg.split_once('=') {
            Some((name, value)) if name.starts_with("--") => (name, Some(value.to_string())),
            _ => (arg.as_str(), None),
        };
        if !name.starts_with("--") {
            fail(&format!("unrecognized arguments: {arg}"));
        }
        let option = full_name(name);
        let flag = matches!(option, "--help" | "--version" | "--no-x" | "--hold-queries");
        if let Some(extra) = inline.as_ref().filter(|_| flag) {
            fail(&format!(
                "argument {option}: ignored explicit argument '{extra}'"
            ));
        }
        let mut value = || match inline.clone() {
            Some(value) => value,
            None => match iter.next() {
                Some(next) if !next.starts_with('-') || looks_negative(&next) => next,
                _ => fail(&format!("argument {option}: expected one argument")),
            },
        };
        match option {
            "--help" => {
                println!("{USAGE}\n{HELP}");
                process::exit(0);
            }
            "--version" => {
                println!("gfa-to-pairwise-paf {}", env!("CARGO_PKG_VERSION"));
                process::exit(0);
            }
            "--reference" => reference = Some(value()),
            "--queries" => {
                args.queries = Some(
                    value()
                        .split(',')
                        .map(str::trim)
                        .filter(|p| !p.is_empty())
                        .map(pansn_of_prefix)
                        .collect(),
                )
            }
            "--max-gap" => args.max_gap = integer(option, &value()),
            "--min-block" => args.min_block = integer(option, &value()),
            "--chrom-sizes-dir" => args.chrom_sizes_dir = Some(value()),
            "--no-x" => args.pair_x = false,
            "--hold-queries" => args.hold_queries = true,
            "--contig-lengths" => args.contig_lengths = Some(value()),
            _ => unreachable!(),
        }
    }
    let Some(reference) = reference else {
        fail("the following arguments are required: --reference");
    };
    args.reference = pansn_of_prefix(&reference);
    args
}
