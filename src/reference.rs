use std::collections::HashMap;

pub const NODE_EXTRA: u32 = 2;
const MAX_RANK: usize = (u32::MAX >> 2) as usize - 1;

pub struct Walk {
    pub lo: usize,
    pub hi: usize,
    pub contig: usize,
}

// Every reference walk laid end to end as one ranked step list. A node's
// `pos` packs its first visit's rank + 1 above two flag bits: reversed, and
// NODE_EXTRA when later visits wait in `extra`.
pub struct Reference {
    pub pos: Vec<u32>,
    pub visited: Vec<u8>,
    pub offsets: Vec<u32>,
    extra: HashMap<usize, Vec<u32>>,
    pub walks: Vec<Walk>,
    pub contigs: Vec<Vec<u8>>,
    contig_ids: HashMap<Vec<u8>, usize>,
    pub lengths: Vec<i64>,
}

impl Reference {
    pub fn new() -> Self {
        Reference {
            pos: Vec::new(),
            visited: Vec::new(),
            offsets: Vec::new(),
            extra: HashMap::new(),
            walks: Vec::new(),
            contigs: Vec::new(),
            contig_ids: HashMap::new(),
            lengths: Vec::new(),
        }
    }

    pub fn grow(&mut self, nodes: usize) {
        if self.pos.len() < nodes {
            self.pos.resize(nodes, 0);
            self.visited.resize(nodes, 0);
        }
    }

    // Ranks the walk's steps after every earlier reference walk, and returns
    // how many of its steps land on nodes an aligned query walked as private.
    pub fn index_walk(
        &mut self,
        contig: &[u8],
        start: i64,
        steps: &[i32],
        lengths: &[u32],
    ) -> Result<u64, String> {
        self.grow(lengths.len());
        let lo = self.offsets.len();
        let mut rank = lo;
        let mut offset = start;
        let mut late = 0;
        for &step in steps {
            let node = node_of(step, lengths.len())?;
            if rank >= MAX_RANK {
                return Err("the reference has more than 2^30 steps".to_string());
            }
            let packed = (((rank + 1) << 2) as u32) | u32::from(step < 0);
            if self.pos[node] != 0 {
                self.extra.entry(node).or_default().push(packed);
                self.pos[node] |= NODE_EXTRA;
            } else {
                self.pos[node] = packed;
            }
            late += u64::from(self.visited[node]);
            self.offsets.push(u32::try_from(offset).map_err(|_| {
                format!(
                    "reference offset {offset} on {} does not fit 32 bits",
                    String::from_utf8_lossy(contig)
                )
            })?);
            offset += i64::from(lengths[node]);
            rank += 1;
        }
        let id = match self.contig_ids.get(contig) {
            Some(&id) => id,
            None => {
                self.contig_ids.insert(contig.to_vec(), self.contigs.len());
                self.contigs.push(contig.to_vec());
                self.lengths.push(0);
                self.contigs.len() - 1
            }
        };
        self.walks.push(Walk {
            lo,
            hi: rank,
            contig: id,
        });
        self.lengths[id] = self.lengths[id].max(offset);
        Ok(late)
    }

    pub fn occurrences(&self, node: usize, first: u32) -> impl Iterator<Item = u32> + '_ {
        let extra = if first & NODE_EXTRA != 0 {
            self.extra[&node].as_slice()
        } else {
            &[]
        };
        std::iter::once(first & !NODE_EXTRA).chain(extra.iter().copied())
    }

    pub fn walk_of(&self, rank: usize) -> &Walk {
        &self.walks[self.walks.partition_point(|walk| walk.lo <= rank) - 1]
    }
}

pub fn node_of(step: i32, nodes: usize) -> Result<usize, String> {
    let node = step.unsigned_abs() as usize;
    if node < nodes {
        Ok(node)
    } else {
        Err(format!(
            "a walk visits segment {node}, which no S line defines"
        ))
    }
}
