use std::collections::HashMap;

const MAX_NODES: usize = 1 << 31;

pub struct Nodes {
    pub lengths: Vec<u32>,
    index: Option<HashMap<Vec<u8>, u32>>,
    pub count: u64,
}

impl Nodes {
    pub fn new() -> Self {
        Nodes {
            lengths: vec![0],
            index: None,
            count: 0,
        }
    }

    pub fn add(&mut self, name: &[u8], length: u32) -> Result<(), String> {
        self.count += 1;
        if self.index.is_none() && is_dense_id(name) {
            let node = parse_int(name)? as usize;
            if node >= MAX_NODES {
                return Err(format!(
                    "segment id {} is past 2^31",
                    String::from_utf8_lossy(name)
                ));
            }
            if node >= self.lengths.len() {
                self.lengths.resize(node + 1, 0);
            }
            self.lengths[node] = length;
            return Ok(());
        }
        let lengths = &mut self.lengths;
        let index = self.index.get_or_insert_with(|| {
            (1..lengths.len())
                .map(|i| (i.to_string().into_bytes(), i as u32))
                .collect()
        });
        match index.get(name) {
            Some(&node) => lengths[node as usize] = length,
            None => {
                if lengths.len() >= MAX_NODES {
                    return Err("more than 2^31 segments".to_string());
                }
                index.insert(name.to_vec(), lengths.len() as u32);
                lengths.push(length);
            }
        }
        Ok(())
    }

    pub fn walk_steps(&self, walk: &[u8], steps: &mut Vec<i32>) -> Result<(), String> {
        steps.clear();
        let is_sep = |b: &u8| *b == b'>' || *b == b'<';
        let Some(mut at) = walk.iter().position(is_sep) else {
            return Ok(());
        };
        while at < walk.len() {
            let reversed = walk[at] == b'<';
            let start = at + 1;
            let end = walk[start..]
                .iter()
                .position(is_sep)
                .map_or(walk.len(), |n| start + n);
            let token = &walk[start..end];
            match &self.index {
                None => steps.push(dense_step(token, reversed)?),
                Some(index) if !token.is_empty() => {
                    let node = lookup(index, token)?;
                    steps.push(if reversed { -node } else { node });
                }
                Some(_) => {}
            }
            at = end;
        }
        Ok(())
    }

    pub fn path_steps(&self, path: &[u8], steps: &mut Vec<i32>) -> Result<(), String> {
        steps.clear();
        for token in path.split(|&b| b == b',') {
            let forward = token.last() == Some(&b'+');
            let name = &token[..token.len().saturating_sub(1)];
            let step = match &self.index {
                None => {
                    let value = parse_int(name)?;
                    to_step(if forward { value } else { -value })?
                }
                Some(index) => {
                    let node = lookup(index, name)?;
                    if forward { node } else { -node }
                }
            };
            steps.push(step);
        }
        Ok(())
    }
}

fn is_dense_id(name: &[u8]) -> bool {
    !name.is_empty() && name[0] != b'0' && name.iter().all(u8::is_ascii_digit)
}

fn lookup(index: &HashMap<Vec<u8>, u32>, name: &[u8]) -> Result<i32, String> {
    index.get(name).map(|&node| node as i32).ok_or_else(|| {
        format!(
            "a walk visits segment {}, which no S line before it defines",
            String::from_utf8_lossy(name)
        )
    })
}

fn to_step(value: i64) -> Result<i32, String> {
    i32::try_from(value).map_err(|_| format!("segment id {value} does not fit 32 bits"))
}

fn dense_step(token: &[u8], reversed: bool) -> Result<i32, String> {
    if !token.is_empty() && token.len() < 10 && token.iter().all(u8::is_ascii_digit) {
        let value = token
            .iter()
            .fold(0i32, |n, &b| n * 10 + i32::from(b - b'0'));
        return Ok(if reversed { -value } else { value });
    }
    let value = if reversed {
        parse_int(&[b"-", token].concat())?
    } else {
        parse_int(token)?
    };
    to_step(value)
}

// Python's int() over bytes: surrounding ASCII whitespace, a sign, and
// underscores between digits.
pub fn parse_int(text: &[u8]) -> Result<i64, String> {
    let invalid = || {
        format!(
            "invalid literal for int() with base 10: {:?}",
            String::from_utf8_lossy(text)
        )
    };
    let space = |b: &u8| b" \t\n\r\x0b\x0c".contains(b);
    let start = text.iter().position(|b| !space(b)).unwrap_or(text.len());
    let end = text
        .iter()
        .rposition(|b| !space(b))
        .map_or(start, |end| end + 1);
    let trimmed = &text[start..end];
    let (negative, digits) = match trimmed.first() {
        Some(b'-') => (true, &trimmed[1..]),
        Some(b'+') => (false, &trimmed[1..]),
        _ => (false, trimmed),
    };
    if digits.first().is_none_or(|b| !b.is_ascii_digit()) || digits.last() == Some(&b'_') {
        return Err(invalid());
    }
    let mut value: i64 = 0;
    let mut underscore = false;
    for &b in digits {
        if b == b'_' {
            if underscore {
                return Err(invalid());
            }
            underscore = true;
            continue;
        }
        if !b.is_ascii_digit() {
            return Err(invalid());
        }
        underscore = false;
        value = value
            .checked_mul(10)
            .and_then(|v| v.checked_add(i64::from(b - b'0')))
            .ok_or_else(|| format!("{} does not fit 64 bits", String::from_utf8_lossy(text)))?;
    }
    Ok(if negative { -value } else { value })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_int_follows_python() {
        assert_eq!(parse_int(b"12"), Ok(12));
        assert_eq!(parse_int(b" -1_000\r\n"), Ok(-1000));
        assert_eq!(parse_int(b"+7"), Ok(7));
        assert_eq!(parse_int(b"007"), Ok(7));
        assert_eq!(parse_int(b"\x0b5\x0c"), Ok(5));
        for bad in [
            &b""[..],
            b"-",
            b"1__0",
            b"_1",
            b"1_",
            b"1a",
            b"--1",
            b"\x0b",
        ] {
            assert!(parse_int(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn walk_steps_dense_and_indexed() {
        let mut nodes = Nodes::new();
        for (name, length) in [(&b"1"[..], 3), (b"2", 4), (b"4", 1)] {
            nodes.add(name, length).unwrap();
        }
        let mut steps = Vec::new();
        nodes.walk_steps(b">1<2>-4", &mut steps).unwrap();
        assert_eq!(steps, [1, -2, -4]);
        nodes.add(b"s9", 2).unwrap();
        nodes.walk_steps(b"x>s9<4>>1", &mut steps).unwrap();
        assert_eq!(steps, [5, -4, 1]);
        assert_eq!(nodes.lengths, [0, 3, 4, 0, 1, 2]);
        assert!(nodes.walk_steps(b">s10", &mut steps).is_err());
    }
}
