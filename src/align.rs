use std::collections::HashMap;

use crate::reference::{NODE_EXTRA, Reference, node_of};

const MATCH: u8 = b'=';
const MISMATCH: u8 = b'X';
const INSERTION: u8 = b'I';
const DELETION: u8 = b'D';

pub struct QueryWalk<'a> {
    pub name: &'a [u8],
    pub contig: &'a [u8],
    pub start: i64,
    pub end: i64,
    pub steps: &'a [i32],
}

pub struct Settings {
    pub max_gap: i64,
    pub min_block: i64,
    pub pair_x: bool,
}

pub struct Row {
    pub qstart: i64,
    pub qend: i64,
    pub flipped: bool,
    pub contig: usize,
    pub tstart: i64,
    pub tend: i64,
    pub matches: i64,
    pub columns: i64,
    pub cigar: Vec<u8>,
}

pub struct QueryContig {
    pub name: Vec<u8>,
    pub rows: Vec<Row>,
    pub length: i64,
}

pub struct Query {
    pub name: Vec<u8>,
    pub walks: u64,
    pub anchors: u64,
    pub chains: u64,
    pub matches: i64,
    pub columns: i64,
    pub contigs: Vec<QueryContig>,
    contig_ids: HashMap<Vec<u8>, usize>,
}

struct Chain {
    lo: usize,
    hi: usize,
    contig: usize,
    flipped: bool,
    qstart: i64,
    tfixed: i64,
    runs: Vec<(i64, u8)>,
}

fn add_run(runs: &mut Vec<(i64, u8)>, length: i64, op: u8) {
    match runs.last_mut() {
        Some(last) if last.1 == op => last.0 += length,
        _ => runs.push((length, op)),
    }
}

impl Query {
    fn new(name: Vec<u8>) -> Self {
        Query {
            name,
            walks: 0,
            anchors: 0,
            chains: 0,
            matches: 0,
            columns: 0,
            contigs: Vec::new(),
            contig_ids: HashMap::new(),
        }
    }

    fn contig(&mut self, name: &[u8]) -> usize {
        if let Some(&id) = self.contig_ids.get(name) {
            return id;
        }
        self.contig_ids.insert(name.to_vec(), self.contigs.len());
        self.contigs.push(QueryContig {
            name: name.to_vec(),
            rows: Vec::new(),
            length: 0,
        });
        self.contigs.len() - 1
    }

    fn emit(&mut self, contig: usize, chain: Chain, qend: i64, tmoving: i64, min_block: i64) {
        let (tstart, tend) = if chain.flipped {
            (tmoving, chain.tfixed)
        } else {
            (chain.tfixed, tmoving)
        };
        if tend - tstart < min_block.max(1) {
            return;
        }
        let mut runs = chain.runs;
        if chain.flipped {
            runs.reverse();
        }
        let matches = runs.iter().filter(|r| r.1 == MATCH).map(|r| r.0).sum();
        let columns = runs.iter().map(|r| r.0).sum();
        let mut cigar = Vec::with_capacity(runs.len() * 4);
        for (length, op) in runs {
            cigar.extend_from_slice(length.to_string().as_bytes());
            cigar.push(op);
        }
        self.contigs[contig].rows.push(Row {
            qstart: chain.qstart,
            qend,
            flipped: chain.flipped,
            contig: chain.contig,
            tstart,
            tend,
            matches,
            columns,
            cigar,
        });
        self.chains += 1;
        self.matches += matches;
        self.columns += columns;
    }
}

pub struct Converter {
    pub settings: Settings,
    pub queries: Vec<Query>,
    query_ids: HashMap<Vec<u8>, usize>,
}

impl Converter {
    pub fn new(settings: Settings) -> Self {
        Converter {
            settings,
            queries: Vec::new(),
            query_ids: HashMap::new(),
        }
    }

    pub fn has_query(&self, name: &[u8]) -> bool {
        self.query_ids.contains_key(name)
    }

    pub fn align(
        &mut self,
        reference: &mut Reference,
        lengths: &[u32],
        walk: QueryWalk,
    ) -> Result<(), String> {
        let id = match self.query_ids.get(walk.name) {
            Some(&id) => id,
            None => {
                self.query_ids
                    .insert(walk.name.to_vec(), self.queries.len());
                self.queries.push(Query::new(walk.name.to_vec()));
                self.queries.len() - 1
            }
        };
        let query = &mut self.queries[id];
        let contig = query.contig(walk.contig);
        let walked = align_walk(
            &self.settings,
            reference,
            lengths,
            query,
            contig,
            walk.start,
            walk.steps,
        )?;
        let length = &mut query.contigs[contig].length;
        *length = (*length).max(walked).max(walk.end);
        Ok(())
    }
}

// Follows one walk. The chain's hot state lives in locals, and the common
// step, the next reference node with no private bp before it, only lengthens
// the open `=` run.
fn align_walk(
    settings: &Settings,
    reference: &mut Reference,
    lengths: &[u32],
    query: &mut Query,
    contig: usize,
    start: i64,
    steps: &[i32],
) -> Result<i64, String> {
    reference.grow(lengths.len());
    let Settings {
        max_gap,
        min_block,
        pair_x,
    } = *settings;
    let mut anchors = 0;
    let mut q = start;
    let mut chain: Option<Chain> = None;
    let mut flipped = false;
    let mut last: i64 = 0;
    let mut qend: i64 = 0;
    let mut tmoving: i64 = 0;
    let mut expect: i64 = 0;
    let mut step4: i64 = 4;
    let mut edge: i64 = 0;
    for &step in steps {
        let node = node_of(step, lengths.len())?;
        let reversed = step < 0;
        let length = i64::from(lengths[node]);
        let p = reference.pos[node];
        if p == 0 {
            reference.visited[node] = 1;
            q += length;
            continue;
        }
        anchors += 1;
        let orientation = i64::from(reversed ^ flipped);
        if q == qend && i64::from(p) == expect | orientation {
            if let Some(chain) = chain.as_mut() {
                chain.runs.last_mut().unwrap().0 += length;
            }
            qend = q + length;
            last += step4.signum();
            expect = if last == edge { -1 } else { expect + step4 };
            tmoving += if flipped { -length } else { length };
            q += length;
            continue;
        }
        let mut chosen = None;
        if let Some(open) = chain.as_mut() {
            for c in reference.occurrences(node, p) {
                let rank = i64::from(c >> 2) - 1;
                let ahead = if flipped { rank < last } else { rank > last };
                let nearer =
                    chosen.is_none_or(|best| if flipped { rank > best } else { rank < best });
                if ((c & 1 != 0) ^ reversed) == flipped
                    && (open.lo as i64) <= rank
                    && rank < open.hi as i64
                    && ahead
                    && nearer
                {
                    chosen = Some(rank);
                }
            }
            if let Some(rank) = chosen {
                let t = i64::from(reference.offsets[rank as usize]);
                let qgap = q - qend;
                let rgap = if flipped {
                    tmoving - (t + length)
                } else {
                    t - tmoving
                };
                if qgap <= max_gap && rgap <= max_gap {
                    let runs = &mut open.runs;
                    if pair_x {
                        let x = qgap.min(rgap);
                        if x != 0 {
                            add_run(runs, x, MISMATCH);
                        }
                        if qgap > x {
                            add_run(runs, qgap - x, INSERTION);
                        }
                        if rgap > x {
                            add_run(runs, rgap - x, DELETION);
                        }
                    } else {
                        if qgap != 0 {
                            add_run(runs, qgap, INSERTION);
                        }
                        if rgap != 0 {
                            add_run(runs, rgap, DELETION);
                        }
                    }
                    add_run(runs, length, MATCH);
                    last = rank;
                    qend = q + length;
                    tmoving = if flipped { t } else { t + length };
                    expect = if last == edge {
                        -1
                    } else {
                        ((last + 1) << 2) + step4
                    };
                } else {
                    chosen = None;
                }
            }
        }
        if chosen.is_none() {
            if let Some(done) = chain.take() {
                query.emit(contig, done, qend, tmoving, min_block);
            }
            let c = p & !NODE_EXTRA;
            last = i64::from(c >> 2) - 1;
            flipped = (c & 1 != 0) ^ reversed;
            let walk = reference.walk_of(last as usize);
            let t = i64::from(reference.offsets[last as usize]);
            chain = Some(Chain {
                lo: walk.lo,
                hi: walk.hi,
                contig: walk.contig,
                flipped,
                qstart: q,
                tfixed: if flipped { t + length } else { t },
                runs: vec![(length, MATCH)],
            });
            qend = q + length;
            tmoving = if flipped { t } else { t + length };
            step4 = if flipped { -4 } else { 4 };
            edge = if flipped {
                walk.lo as i64
            } else {
                walk.hi as i64 - 1
            };
            expect = if last == edge {
                -1
            } else {
                ((last + 1) << 2) + step4
            };
        }
        q += length;
    }
    if let Some(done) = chain {
        query.emit(contig, done, qend, tmoving, min_block);
    }
    query.walks += 1;
    query.anchors += anchors;
    Ok(q)
}
