mod align;
mod cli;
mod nodes;
mod reference;

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process;
use std::time::Instant;

use flate2::read::MultiGzDecoder;

use align::{Converter, QueryWalk, Settings};
use cli::Args;
use nodes::{Nodes, parse_int};
use reference::Reference;

const PAF_MAPQ: u32 = 255;

struct Held {
    name: Vec<u8>,
    contig: Vec<u8>,
    start: i64,
    end: i64,
    steps: Vec<i32>,
}

fn open_input(path: Option<&str>) -> Result<Box<dyn BufRead>, String> {
    let mut raw: Box<dyn Read> = match path {
        None | Some("-") => Box::new(io::stdin()),
        Some(path) => Box::new(File::open(path).map_err(|e| format!("{path}: {e}"))?),
    };
    let mut magic = Vec::with_capacity(2);
    (&mut raw)
        .take(2)
        .read_to_end(&mut magic)
        .map_err(|e| e.to_string())?;
    let gzipped = magic == [0x1f, 0x8b];
    let stream = io::Cursor::new(magic).chain(raw);
    Ok(if gzipped {
        Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(stream),
        ))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, stream))
    })
}

fn find(line: &[u8], byte: u8, from: usize) -> Option<usize> {
    line.get(from..)?
        .iter()
        .position(|&b| b == byte)
        .map(|n| from + n)
}

fn field_tab(line: &[u8], from: usize) -> Result<usize, String> {
    find(line, b'\t', from).ok_or_else(|| {
        format!(
            "too few fields: {}",
            String::from_utf8_lossy(&line[..line.len().min(80)])
        )
    })
}

fn bound(text: &[u8]) -> Result<i64, String> {
    if text == b"*" { Ok(0) } else { parse_int(text) }
}

fn segment_length(line: &[u8], tab1: usize) -> Result<u32, String> {
    let end = find(line, b'\t', tab1 + 1).unwrap_or(line.len() - 1) as i64;
    let mut length = end - tab1 as i64 - 1;
    if length == 1 && line[tab1 + 1] == b'*' {
        length = line_length_tag(line);
    }
    u32::try_from(length).map_err(|_| {
        format!(
            "segment {} has length {length}",
            String::from_utf8_lossy(&line[2..tab1])
        )
    })
}

fn line_length_tag(line: &[u8]) -> i64 {
    let mut from = 0;
    while let Some(at) = line[from..]
        .windows(6)
        .position(|w| w == b"\tLN:i:")
        .map(|n| from + n + 6)
    {
        let digits = line[at..].iter().take_while(|b| b.is_ascii_digit()).count();
        if digits > 0 {
            return parse_int(&line[at..at + digits]).unwrap_or(i64::MAX);
        }
        from = at;
    }
    0
}

fn path_name_parts(name: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let fields: Vec<_> = name.split(|&b| b == b'#').collect();
    if fields.len() >= 3 {
        (
            [fields[0], b"#", fields[1]].concat(),
            fields[2..].join(&b'#'),
        )
    } else {
        ([name, b"#0"].concat(), name.to_vec())
    }
}

fn read_contig_lengths(path: Option<&str>) -> Result<HashMap<Vec<u8>, i64>, String> {
    let mut lengths = HashMap::new();
    if let Some(path) = path {
        let text = fs::read(path).map_err(|e| format!("{path}: {e}"))?;
        for line in text.split_inclusive(|&b| b == b'\n') {
            let mut fields = line
                .split(|b| b" \t\n\r\x0b\x0c".contains(b))
                .filter(|f| !f.is_empty());
            if let (Some(name), Some(length)) = (fields.next(), fields.next()) {
                lengths.insert(name.to_vec(), parse_int(length)?);
            }
        }
    }
    Ok(lengths)
}

struct Run {
    args: Args,
    started: Instant,
    nodes: Nodes,
    reference: Reference,
    converter: Converter,
    wanted: Option<HashSet<Vec<u8>>>,
    held: Vec<Held>,
    steps: Vec<i32>,
}

impl Run {
    fn walk(&mut self, name: Vec<u8>, contig: &[u8], start: i64, end: i64) -> Result<(), String> {
        let lengths = &self.nodes.lengths;
        if name == self.args.reference {
            let late = self
                .reference
                .index_walk(contig, start, &self.steps, lengths)?;
            let walk = self.reference.walks.last().unwrap();
            eprintln!(
                "{}#{}: {} steps, {} bp; {:.0}s",
                String::from_utf8_lossy(&self.args.reference),
                String::from_utf8_lossy(contig),
                walk.hi - walk.lo,
                self.reference.lengths[walk.contig],
                self.started.elapsed().as_secs_f64()
            );
            if late > 0 {
                return Err(format!(
                    "{} {}:{start} arrived after a query walk that visits {late} of its nodes had \
                     already been aligned; re-run with --hold-queries",
                    String::from_utf8_lossy(&self.args.reference),
                    String::from_utf8_lossy(contig),
                ));
            }
        } else if self.args.hold_queries || self.reference.walks.is_empty() {
            self.held.push(Held {
                name,
                contig: contig.to_vec(),
                start,
                end,
                steps: self.steps.clone(),
            });
        } else {
            let walk = QueryWalk {
                name: &name,
                contig,
                start,
                end,
                steps: &self.steps,
            };
            self.converter.align(&mut self.reference, lengths, walk)?;
        }
        Ok(())
    }

    fn wants(&self, name: &[u8]) -> bool {
        name == self.args.reference.as_slice()
            || self.wanted.as_ref().is_none_or(|w| w.contains(name))
    }

    fn line(&mut self, line: &[u8]) -> Result<(), String> {
        match line[0] {
            b'S' => {
                let tab1 = field_tab(line, 2)?;
                let length = segment_length(line, tab1)?;
                self.nodes.add(&line[2..tab1], length)?;
            }
            b'W' => {
                let tab1 = field_tab(line, 2)?;
                let tab2 = field_tab(line, tab1 + 1)?;
                let name = [&line[2..tab1], b"#", &line[tab1 + 1..tab2]].concat();
                if self.wants(&name) {
                    let rest = &line[tab2 + 1..];
                    let rest = &rest
                        [..rest.len() - rest.iter().rev().take_while(|&&b| b == b'\n').count()];
                    let fields: Vec<_> = rest.splitn(5, |&b| b == b'\t').collect();
                    if fields.len() < 4 {
                        return Err(format!(
                            "W line for {} has too few fields",
                            String::from_utf8_lossy(&name)
                        ));
                    }
                    let (start, end) = (bound(fields[1])?, bound(fields[2])?);
                    self.nodes.walk_steps(fields[3], &mut self.steps)?;
                    self.walk(name, fields[0], start, end)?;
                }
            }
            b'P' => {
                let tab1 = field_tab(line, 2)?;
                let (name, contig) = path_name_parts(&line[2..tab1]);
                if self.wants(&name) {
                    let end = find(line, b'\t', tab1 + 1).unwrap_or(line.len() - 1);
                    let path = line.get(tab1 + 1..end).unwrap_or_default();
                    self.nodes.path_steps(path, &mut self.steps)?;
                    self.walk(name, &contig, 0, 0)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn write_paf(
        &self,
        out: &mut impl Write,
        contig_lengths: &HashMap<Vec<u8>, i64>,
    ) -> io::Result<()> {
        let reference = &self.args.reference;
        let lookup = |full: &[u8], contig: &[u8], fallback: i64| {
            contig_lengths
                .get(full)
                .or_else(|| contig_lengths.get(contig))
                .copied()
                .unwrap_or(fallback)
        };
        for query in &self.converter.queries {
            for contig in &query.contigs {
                let qname = [&query.name[..], b"#", &contig.name].concat();
                let qlen = lookup(&qname, &contig.name, contig.length);
                for row in &contig.rows {
                    let ref_contig = &self.reference.contigs[row.contig];
                    let tname = [&reference[..], b"#", ref_contig].concat();
                    let tlen = lookup(&tname, ref_contig, self.reference.lengths[row.contig]);
                    out.write_all(&qname)?;
                    write!(
                        out,
                        "\t{qlen}\t{}\t{}\t{}\t",
                        row.qstart,
                        row.qend,
                        if row.flipped { '-' } else { '+' }
                    )?;
                    out.write_all(&tname)?;
                    write!(
                        out,
                        "\t{tlen}\t{}\t{}\t{}\t{}\t{PAF_MAPQ}\tcg:Z:",
                        row.tstart, row.tend, row.matches, row.columns
                    )?;
                    out.write_all(&row.cigar)?;
                    out.write_all(b"\n")?;
                }
            }
        }
        out.flush()
    }

    fn write_chrom_sizes(
        &self,
        directory: &str,
        contig_lengths: &HashMap<Vec<u8>, i64>,
    ) -> io::Result<()> {
        fs::create_dir_all(directory)?;
        for query in &self.converter.queries {
            let split = query.name.iter().position(|&b| b == b'#').unwrap();
            let (sample, hap) = (&query.name[..split], &query.name[split + 1..]);
            let mut sizes: Vec<(&[u8], i64)> = query
                .contigs
                .iter()
                .map(|contig| {
                    let full = [&query.name[..], b"#", &contig.name].concat();
                    let size = contig_lengths
                        .get(&full)
                        .or_else(|| contig_lengths.get(&contig.name))
                        .copied()
                        .unwrap_or(contig.length);
                    (contig.name.as_slice(), size)
                })
                .collect();
            sizes.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
            let file = [sample, b".", hap, b".chrom.sizes"].concat();
            let mut out = BufWriter::new(File::create(
                Path::new(directory).join(String::from_utf8_lossy(&file).as_ref()),
            )?);
            for (contig, size) in sizes {
                out.write_all(contig)?;
                writeln!(out, "\t{size}")?;
            }
            out.flush()?;
        }
        Ok(())
    }
}

fn run(args: Args) -> Result<(), String> {
    let started = Instant::now();
    let mut input = open_input(args.gfa.as_deref())?;
    let contig_lengths = read_contig_lengths(args.contig_lengths.as_deref())?;
    let settings = Settings {
        max_gap: args.max_gap,
        min_block: args.min_block,
        pair_x: args.pair_x,
    };
    let wanted = args.queries.as_ref().map(|q| q.iter().cloned().collect());
    let mut run = Run {
        args,
        started,
        nodes: Nodes::new(),
        reference: Reference::new(),
        converter: Converter::new(settings),
        wanted,
        held: Vec::new(),
        steps: Vec::new(),
    };
    let mut consumed: u64 = 0;
    let mut line = Vec::with_capacity(1 << 20);
    loop {
        line.clear();
        let n = input
            .read_until(b'\n', &mut line)
            .map_err(|e| format!("reading the GFA: {e}"))?;
        if n == 0 {
            break;
        }
        consumed += n as u64;
        run.line(&line)?;
    }
    let reference_name = String::from_utf8_lossy(&run.args.reference).into_owned();
    if run.reference.walks.is_empty() {
        return Err(format!(
            "no {reference_name} walk in the input; is --reference right? (a bare sample means haplotype 0)"
        ));
    }
    for held in std::mem::take(&mut run.held) {
        let walk = QueryWalk {
            name: &held.name,
            contig: &held.contig,
            start: held.start,
            end: held.end,
            steps: &held.steps,
        };
        run.converter
            .align(&mut run.reference, &run.nodes.lengths, walk)?;
    }
    let stdout = io::stdout();
    run.write_paf(
        &mut BufWriter::with_capacity(1 << 20, stdout.lock()),
        &contig_lengths,
    )
    .map_err(|e| format!("writing PAF: {e}"))?;
    if let Some(directory) = run.args.chrom_sizes_dir.as_deref() {
        run.write_chrom_sizes(directory, &contig_lengths)
            .map_err(|e| format!("{directory}: {e}"))?;
    }
    let elapsed = started.elapsed().as_secs_f64();
    for query in &run.converter.queries {
        eprintln!(
            "{}: {} walks, {} anchors -> {} chains, {} bp =, {} columns",
            String::from_utf8_lossy(&query.name),
            query.walks,
            query.anchors,
            query.chains,
            query.matches,
            query.columns
        );
    }
    let unseen: Vec<_> = run
        .args
        .queries
        .iter()
        .flatten()
        .filter(|q| !run.converter.has_query(q))
        .map(|q| String::from_utf8_lossy(q).into_owned())
        .collect();
    let megabytes = consumed as f64 / 1e6;
    eprintln!(
        "{} nodes, {} {reference_name} steps on {} walks; {megabytes:.0} MB in {elapsed:.1}s ({:.0} MB/s){}",
        run.nodes.count,
        run.reference.offsets.len(),
        run.reference.walks.len(),
        megabytes / elapsed.max(1e-9),
        if unseen.is_empty() {
            String::new()
        } else {
            format!("; no walks for {}", unseen.join(","))
        }
    );
    Ok(())
}

fn main() {
    if let Err(message) = run(cli::parse()) {
        eprintln!("{message}");
        process::exit(1);
    }
}
