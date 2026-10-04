// Every row worked out by hand over a graph small enough to check on paper,
// ported assertion for assertion from the Python tool's tests. Nothing
// downstream of the converter would catch a flipped chain's coordinates and
// CIGAR direction, or where a chain breaks, so these pin each one.
//
// GRCh38 chrA walks 1 2 4 5 6 7 8 9 8 10 (59 bp; node 8 twice, at 33 and 49).
// HG01109 ctgA arrives as two W pieces: the first takes the SNP allele 3 for 2,
// skips 5, inserts 11, and reaches 8 the first time; the second starts at
// offset 50 on 9 and reaches 8 the second time. HG01123 ctgB traverses 7 and 5
// backwards, skipping 6 between them. HG00097 ctgC is one shared node.

use std::fs;
use std::io::Write;
use std::process::{Command, Output, Stdio};

const WITH_WALKS: &str = include_str!("data/walks.gfa");
const WITH_PATHS: &str = include_str!("data/paths.gfa");
const QUERIES_FIRST: &str = include_str!("data/queries-first.gfa");
const LATE_REFERENCE: &str = include_str!("data/late-reference.gfa");

type Row = Vec<String>;

macro_rules! row {
    ($($field:expr),* $(,)?) => { vec![$($field.to_string()),*] };
}

// qname qlen qstart qend strand tname tlen tstart tend matches blocklen mapq cg
fn ctga_first_piece() -> Row {
    row![
        "HG01109#1#ctgA",
        70,
        0,
        50,
        "+",
        "GRCh38#0#chrA",
        59,
        0,
        49,
        45,
        53,
        255,
        "cg:Z:10=1X5=3D8=4I22="
    ]
}
fn ctga_second_piece() -> Row {
    row![
        "HG01109#1#ctgA",
        70,
        50,
        70,
        "+",
        "GRCh38#0#chrA",
        59,
        39,
        59,
        20,
        20,
        255,
        "cg:Z:20="
    ]
}
fn ctgb_forward() -> Row {
    row![
        "HG01123#1#ctgB",
        23,
        0,
        10,
        "+",
        "GRCh38#0#chrA",
        59,
        0,
        10,
        10,
        10,
        255,
        "cg:Z:10="
    ]
}
// 7 then 5 backwards: reference 16..33 with the 8 bp of node 6 deleted, and
// the CIGAR reads along the reference, so node 5 comes first
fn ctgb_inverted() -> Row {
    row![
        "HG01123#1#ctgB",
        23,
        10,
        19,
        "-",
        "GRCh38#0#chrA",
        59,
        16,
        33,
        9,
        17,
        255,
        "cg:Z:3=8D6="
    ]
}
fn ctgb_tail() -> Row {
    row![
        "HG01123#1#ctgB",
        23,
        19,
        23,
        "+",
        "GRCh38#0#chrA",
        59,
        55,
        59,
        4,
        4,
        255,
        "cg:Z:4="
    ]
}
fn ctgc() -> Row {
    row![
        "HG00097#1#ctgC",
        10,
        0,
        10,
        "+",
        "GRCh38#0#chrA",
        59,
        0,
        10,
        10,
        10,
        255,
        "cg:Z:10="
    ]
}
fn all_rows() -> Vec<Row> {
    vec![
        ctgc(),
        ctga_first_piece(),
        ctga_second_piece(),
        ctgb_forward(),
        ctgb_inverted(),
        ctgb_tail(),
    ]
}

fn identical(q: &str, qlen: i64, qs: i64, qe: i64, strand: &str, ts: i64, te: i64) -> Row {
    let n = qe - qs;
    row![
        q,
        qlen,
        qs,
        qe,
        strand,
        "GRCh38#0#chrA",
        59,
        ts,
        te,
        n,
        n,
        255,
        format!("cg:Z:{n}=")
    ]
}

fn with_qlen(mut fields: Row, qlen: i64) -> Row {
    fields[1] = qlen.to_string();
    fields
}

struct Run {
    code: i32,
    stdout: String,
    stderr: String,
}

fn convert(gfa: &str, args: &[&str]) -> Run {
    let mut child = Command::new(env!("CARGO_BIN_EXE_gfa-to-pairwise-paf"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(gfa.as_bytes())
        .unwrap();
    let Output {
        status,
        stdout,
        stderr,
    } = child.wait_with_output().unwrap();
    Run {
        code: status.code().unwrap(),
        stdout: String::from_utf8(stdout).unwrap(),
        stderr: String::from_utf8(stderr).unwrap(),
    }
}

fn paf_rows(stdout: &str) -> Vec<Row> {
    let mut rows: Vec<Row> = stdout
        .split('\n')
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    rows.sort_by_key(|r| (r[0].clone(), r[2].parse::<i64>().unwrap()));
    rows
}

fn succeeded(run: &Run) {
    assert_eq!(run.code, 0, "{}", run.stderr);
}

#[test]
fn every_row_by_hand() {
    let run = convert(WITH_WALKS, &["--reference", "GRCh38#0"]);
    succeeded(&run);
    assert_eq!(paf_rows(&run.stdout), all_rows());
    assert!(
        run.stderr
            .contains("HG01109#1: 2 walks, 9 anchors -> 2 chains, 65 bp =, 73 columns")
    );
    assert!(
        run.stderr
            .contains("11 nodes, 10 GRCh38#0 steps on 1 walks")
    );
}

#[test]
fn p_lines_and_bare_reference() {
    let run = convert(WITH_PATHS, &["--reference", "GRCh38"]);
    succeeded(&run);
    assert_eq!(paf_rows(&run.stdout), all_rows());
}

#[test]
fn max_gap_breaks_a_chain() {
    let run = convert(
        WITH_WALKS,
        &[
            "--reference",
            "GRCh38",
            "--queries",
            "HG01109#1,HG01123#1",
            "--max-gap",
            "2",
        ],
    );
    succeeded(&run);
    assert_eq!(
        paf_rows(&run.stdout),
        vec![
            row![
                "HG01109#1#ctgA",
                70,
                0,
                16,
                "+",
                "GRCh38#0#chrA",
                59,
                0,
                16,
                15,
                16,
                255,
                "cg:Z:10=1X5="
            ],
            identical("HG01109#1#ctgA", 70, 16, 24, "+", 19, 27),
            identical("HG01109#1#ctgA", 70, 28, 50, "+", 27, 49),
            ctga_second_piece(),
            ctgb_forward(),
            identical("HG01123#1#ctgB", 23, 10, 16, "-", 27, 33),
            identical("HG01123#1#ctgB", 23, 16, 19, "-", 16, 19),
            ctgb_tail(),
        ]
    );
}

#[test]
fn no_x_writes_insertion_then_deletion() {
    let run = convert(
        WITH_WALKS,
        &["--reference", "GRCh38", "--queries", "HG01109#1", "--no-x"],
    );
    succeeded(&run);
    assert_eq!(
        paf_rows(&run.stdout),
        vec![
            row![
                "HG01109#1#ctgA",
                70,
                0,
                50,
                "+",
                "GRCh38#0#chrA",
                59,
                0,
                49,
                45,
                54,
                255,
                "cg:Z:10=1I1D5=3D8=4I22="
            ],
            ctga_second_piece(),
        ]
    );
}

#[test]
fn min_block_and_queries() {
    let run = convert(
        WITH_WALKS,
        &[
            "--reference",
            "GRCh38",
            "--queries",
            "HG01123#1",
            "--min-block",
            "15",
        ],
    );
    succeeded(&run);
    assert_eq!(paf_rows(&run.stdout), vec![ctgb_inverted()]);
}

#[test]
fn chrom_sizes_per_query() {
    let tmp = std::env::temp_dir().join(format!("gfa-to-pairwise-paf-{}", std::process::id()));
    fs::create_dir_all(&tmp).unwrap();
    let lengths = tmp.join("lengths.txt");
    fs::write(&lengths, "ctgA\t100\n").unwrap();
    let sizes = tmp.join("sizes");
    let run = convert(
        WITH_WALKS,
        &[
            "--reference",
            "GRCh38",
            "--chrom-sizes-dir",
            sizes.to_str().unwrap(),
            "--contig-lengths",
            lengths.to_str().unwrap(),
        ],
    );
    succeeded(&run);
    for (name, expected) in [
        ("HG01109.1.chrom.sizes", "ctgA\t100\n"),
        ("HG01123.1.chrom.sizes", "ctgB\t23\n"),
        ("HG00097.1.chrom.sizes", "ctgC\t10\n"),
    ] {
        assert_eq!(fs::read_to_string(sizes.join(name)).unwrap(), expected);
    }
    fs::remove_dir_all(&tmp).unwrap();
    assert_eq!(
        paf_rows(&run.stdout)
            .into_iter()
            .filter(|r| r[0] == "HG01109#1#ctgA")
            .collect::<Vec<_>>(),
        vec![
            with_qlen(ctga_first_piece(), 100),
            with_qlen(ctga_second_piece(), 100)
        ]
    );
}

#[test]
fn query_walks_ahead_of_the_reference() {
    let held = convert(QUERIES_FIRST, &["--reference", "GRCh38"]);
    succeeded(&held);
    assert_eq!(paf_rows(&held.stdout), all_rows());

    // minigraph-cactus writes one chromosome's S, L and W lines after
    // another's, so a second reference contig can follow the query walks;
    // node 11 was private to ctgA when it was aligned, which the guard catches
    let refused = convert(LATE_REFERENCE, &["--reference", "GRCh38"]);
    assert_ne!(refused.code, 0);
    assert!(
        refused
            .stderr
            .contains("GRCh38#0 chrB:0 arrived after a query walk that visits 1 of its nodes")
    );

    let held_all = convert(LATE_REFERENCE, &["--reference", "GRCh38", "--hold-queries"]);
    succeeded(&held_all);
    assert_eq!(
        paf_rows(&held_all.stdout)
            .into_iter()
            .filter(|r| r[0] == "HG01109#1#ctgA")
            .collect::<Vec<_>>(),
        vec![
            row![
                "HG01109#1#ctgA",
                70,
                0,
                24,
                "+",
                "GRCh38#0#chrA",
                59,
                0,
                27,
                23,
                27,
                255,
                "cg:Z:10=1X5=3D8="
            ],
            row![
                "HG01109#1#ctgA",
                70,
                24,
                28,
                "+",
                "GRCh38#0#chrB",
                4,
                0,
                4,
                4,
                4,
                255,
                "cg:Z:4="
            ],
            row![
                "HG01109#1#ctgA",
                70,
                28,
                50,
                "+",
                "GRCh38#0#chrA",
                59,
                27,
                49,
                22,
                22,
                255,
                "cg:Z:22="
            ],
            ctga_second_piece(),
        ]
    );
}

#[test]
fn wrong_reference_fails() {
    let run = convert(WITH_WALKS, &["--reference", "CHM13"]);
    assert_ne!(run.code, 0);
    assert!(run.stderr.contains("no CHM13#0 walk in the input"));
}
