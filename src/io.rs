/// TODO: This whole module to be removed and `fgoxide` crate imported
/// TODO: instead once it is published to crates.io
use anyhow::Result;
use csv::{QuoteStyle, ReaderBuilder, WriterBuilder};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde::{de::DeserializeOwned, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// The set of file extensions to treat as GZIPPED
const GZIP_EXTENSIONS: [&str; 2] = ["gz", "bgz"];

pub struct Io {
    compression: Compression,
}

/** Returns a Default implementation that will compress to gzip level 5. */
impl Default for Io {
    fn default() -> Self {
        Io::new(5)
    }
}

impl Io {
    /// Creates a new Io instance with the given compression level.
    fn new(compression: u32) -> Io {
        Io { compression: flate2::Compression::new(compression) }
    }

    /// Returns true if the path ends with a recognized GZIP file extension
    fn is_gzip_path<P: AsRef<Path>>(p: &P) -> bool {
        if let Some(ext) = p.as_ref().extension() {
            match ext.to_str() {
                Some(x) => GZIP_EXTENSIONS.contains(&x),
                None => false,
            }
        } else {
            false
        }
    }

    /// Opens a file for reading.  Transparently handles reading gzipped files based
    /// extension.
    pub fn new_reader<P>(&self, p: &P) -> Result<BufReader<Box<dyn Read>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(p)?;
        let read: Box<dyn Read> =
            if Io::is_gzip_path(p) { Box::new(GzDecoder::new(file)) } else { Box::new(file) };

        Ok(BufReader::new(read))
    }

    /// Opens a file for writing. Transparently handles writing GZIP'd data if the file
    /// ends with a recognized GZIP extension.
    pub fn new_writer<P>(&self, p: &P) -> Result<BufWriter<Box<dyn Write>>>
    where
        P: AsRef<Path>,
    {
        let file = File::create(p)?;
        let write: Box<dyn Write> = if Io::is_gzip_path(p) {
            Box::new(GzEncoder::new(file, self.compression))
        } else {
            Box::new(file)
        };

        Ok(BufWriter::new(write))
    }

    /// Reads lines from a file into a Vec
    pub fn read_lines<P>(&self, p: &P) -> Result<Vec<String>>
    where
        P: AsRef<Path>,
    {
        let r = self.new_reader(p)?;
        let mut v = Vec::new();
        for result in r.lines() {
            v.push(result?);
        }

        Ok(v)
    }

    /// Writes all the lines from an iterable of string-like values to a file, separated by new lines.
    pub fn write_lines<P, S>(&self, p: &P, lines: impl IntoIterator<Item = S>) -> Result<()>
    where
        P: AsRef<Path>,
        S: AsRef<str>,
    {
        let mut out = self.new_writer(p)?;
        for line in lines {
            out.write_all(line.as_ref().as_bytes())?;
            out.write_all(&[b'\n'])?;
        }

        out.flush()?;
        Ok(())
    }
}

/// Unit-struct that contains associated functions for reading and writing Structs to/from
/// delimited files.  Structs should use serde's Serialize/Deserialize derive macros in
/// order to be used with these functions.
pub struct DelimFile {
    io: Io,
}

/// Generates a default implementation that uses the default Io instance
impl Default for DelimFile {
    fn default() -> Self {
        DelimFile { io: Io::default() }
    }
}

impl DelimFile {
    /// Writes a series of one or more structs to a delimited file.  If `quote` is true then fields
    /// will be quoted as necessary, otherwise they will never be quoted.
    pub fn write<S, P>(
        &self,
        path: &P,
        recs: impl IntoIterator<Item = S>,
        delimiter: u8,
        quote: bool,
    ) -> Result<()>
    where
        S: Serialize,
        P: AsRef<Path>,
    {
        let write = self.io.new_writer(path)?;

        let mut writer = WriterBuilder::new()
            .delimiter(delimiter)
            .has_headers(true)
            .quote_style(if quote { QuoteStyle::Necessary } else { QuoteStyle::Never })
            .from_writer(write);

        for rec in recs {
            writer.serialize(rec)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Writes structs implementing `[Serialize]` to a file with tab separators between fields.
    pub fn write_tsv<S, P>(&self, path: &P, recs: impl IntoIterator<Item = S>) -> Result<()>
    where
        S: Serialize,
        P: AsRef<Path>,
    {
        self.write(path, recs, b'\t', true)
    }

    /// Writes structs implementing `[Serialize]` to a file with comma separators between fields.
    pub fn write_csv<S, P>(&self, path: &P, recs: impl IntoIterator<Item = S>) -> Result<()>
    where
        S: Serialize,
        P: AsRef<Path>,
    {
        self.write(path, recs, b',', true)
    }

    /// Writes a series of one or more structs to a delimited file.  If `quote` is true then fields
    /// will be quoted as necessary, otherwise they will never be quoted.
    pub fn read<D, P>(&self, path: &P, delimiter: u8, quote: bool) -> Result<Vec<D>>
    where
        D: DeserializeOwned,
        P: AsRef<Path>,
    {
        let read = self.io.new_reader(path)?;

        let mut reader = ReaderBuilder::new()
            .delimiter(delimiter)
            .has_headers(true)
            .quoting(quote)
            .from_reader(read);

        let mut results = vec![];

        for result in reader.deserialize::<D>() {
            let rec = result?;
            results.push(rec);
        }

        Ok(results)
    }

    /// Reads structs implementing `[Deserialize]` from a file with tab separators between fields.
    pub fn read_tsv<D, P>(&self, path: &P) -> Result<Vec<D>>
    where
        D: DeserializeOwned,
        P: AsRef<Path>,
    {
        self.read(path, b'\t', true)
    }

    /// Reads structs implementing `[Deserialize]` from a file with tab separators between fields.
    pub fn read_csv<D, P>(&self, path: &P) -> Result<Vec<D>>
    where
        D: DeserializeOwned,
        P: AsRef<Path>,
    {
        self.read(path, b',', true)
    }
}

#[cfg(test)]
mod tests {
    use crate::io::{DelimFile, Io};
    use serde::{Deserialize, Serialize};
    use tempfile::TempDir;

    /// Record type used in testing DelimFile
    #[derive(Debug, Serialize, Deserialize, PartialEq)]
    struct Rec {
        s: String,
        i: usize,
        b: bool,
        o: Option<f64>,
    }

    #[test]
    fn test_reading_and_writing_lines_to_file() {
        let lines = vec!["foo", "bar,splat,whee", "baz\twhoopsie"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("strs.txt");
        let f2 = tempdir.path().join("Strings.txt");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let strings: Vec<String> = lines.iter().map(|l| l.to_string()).collect();
        io.write_lines(&f2, &strings).unwrap();

        let r1 = io.read_lines(&f1).unwrap();
        let r2 = io.read_lines(&f2).unwrap();

        assert_eq!(r1, lines);
        assert_eq!(r2, lines);
    }

    #[test]
    fn test_reading_and_writing_gzip_files() {
        let lines = vec!["foo", "bar", "baz"];
        let tempdir = TempDir::new().unwrap();
        let text = tempdir.path().join("text.txt");
        let gzipped = tempdir.path().join("gzipped.txt.gz");

        let io = Io::default();
        io.write_lines(&text, &mut lines.iter()).unwrap();
        io.write_lines(&gzipped, &mut lines.iter()).unwrap();

        let r1 = io.read_lines(&text).unwrap();
        let r2 = io.read_lines(&gzipped).unwrap();

        assert_eq!(r1, lines);
        assert_eq!(r2, lines);

        // Also check that we actually wrote gzipped data to the gzip file!
        assert_ne!(text.metadata().unwrap().len(), gzipped.metadata().unwrap().len());
    }

    #[test]
    fn test_reading_and_writing_empty_delim_file() {
        let recs: Vec<Rec> = vec![];
        let tmp = TempDir::new().unwrap();
        let csv = tmp.path().join("recs.csv");
        let tsv = tmp.path().join("recs.tsv.gz");

        let df = DelimFile::default();
        df.write_csv(&csv, &recs).unwrap();
        df.write_tsv(&tsv, &recs).unwrap();
        let from_csv: Vec<Rec> = df.read_csv(&csv).unwrap();
        let from_tsv: Vec<Rec> = df.read_tsv(&tsv).unwrap();

        assert_eq!(from_csv, recs);
        assert_eq!(from_tsv, recs);
    }

    #[test]
    fn test_reading_and_writing_delim_file() {
        let recs: Vec<Rec> = vec![
            Rec { s: "Hello".to_string(), i: 123, b: true, o: None },
            Rec { s: "A,B,C".to_string(), i: 456, b: false, o: Some(123.45) },
        ];
        let tmp = TempDir::new().unwrap();
        let csv = tmp.path().join("recs.csv");
        let tsv = tmp.path().join("recs.tsv.gz");

        let df = DelimFile::default();
        df.write_csv(&csv, &recs).unwrap();
        df.write_tsv(&tsv, &recs).unwrap();
        let from_csv: Vec<Rec> = df.read_csv(&csv).unwrap();
        let from_tsv: Vec<Rec> = df.read_tsv(&tsv).unwrap();

        assert_eq!(from_csv, recs);
        assert_eq!(from_tsv, recs);
    }
}
