use crate::command::Command;
use crate::guide::*;
use ahash::AHashMap;
use anyhow::{Context, Result};
use clap::Parser;
use fastq::{parse_path, Record};
use fgoxide::io::{DelimFile, Io};
use itertools::Itertools;
use log::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

type GuideMap<'a> = AHashMap<Vec<u8>, &'a Guide>;

/// Counts the guides observed in a CRISPR screen, starting from one or more FASTQs.  FASTQs are
/// one per sample and currently only single-end FASTQ inputs are supported.
///
/// A set of sample IDs may be provided using `--samples id1 id2 ..`.  If provided it must have the
/// same number of values as input FASTQs.  If not provided the FASTQ names are used minus any
/// fastq/fq/gz suffixes.
///
/// Automatically determines the range of valid offsets within the sequencing reads where the
/// guide sequences are located, independently for each FASTQ input.  The first `offset-sample-size`
/// reads from each FASTQ are examined to determine the offsets at which guides are found. When
/// processing the full FASTQ, checks only those offsets that accounted for at least
/// `offset-min-fraction` of the first `offset-sample-size` reads.
///
/// Matching by default allows for one mismatch (and no indels) between the read sub-sequence
/// and the expected guide sequences.  Exact matching may be enabled by specifying the
/// `--exact-match` option.
///
/// Optionally lists may be provided of essential genes, nonessential genes and control guide ids,
/// as well as a regular expression to be used to identify control guides.  Using this information
/// guides are classified as either Essential, Nonessential, Control, or Other.
///
/// Three output files are generated.  The first is named `{output}.counts.txt` and contains columns
/// for the guide id, the gene targeted by the guide and one count column per input FASTQ with
/// raw/un-normalized counts.  The second, `{output}.extended-counts.txt` is identical to the first
/// except for having a `guide_type` column inserted as the third column.  Finally
/// `{output}.stats.txt` contains basic QC statistics per input FASTQ on the matching process.
#[derive(Parser, Debug)]
pub(crate) struct Count {
    /// Input fastq file(s)
    #[clap(long, short = 'i', required = true, multiple_values = true)]
    input: Vec<PathBuf>,

    /// Sample names corresponding to the input fastqs. If provided must be the same length as
    /// input.  Otherwise will be inferred from input file names.
    #[clap(long, short = 's', multiple_values = true)]
    samples: Vec<String>,

    /// Path to the guide library metadata.  May be a tab- or comma-separated file.  Must have
    /// a header line, and the first three fields must be (in order): i) the ID of the guide,
    /// ii) the base sequence of the guide, iii) the gene the guide targets.
    #[clap(long, short = 'l')]
    library: PathBuf,

    /// Optional path to file with list of essential genes.  Gene names should appear one
    /// per line and are case sensitive. If the file has multiple tab-separated columns, the first
    /// column is used.
    #[clap(long, short = 'e')]
    essential_genes: Option<PathBuf>,

    /// Optional path to file with list of nonessential genes.  Gene names should appear one
    /// per line and are case sensitive.  If the file has multiple tab-separated columns, the first
    /// column is used.
    #[clap(long, short = 'n')]
    nonessential_genes: Option<PathBuf>,

    /// Optional path to file with list control guide IDs.  IDs should appear one
    /// per line and are case sensitive.  If the file has multiple tab-separated columns, the first
    /// column is used.
    #[clap(long, short = 'c')]
    control_guides: Option<PathBuf>,

    /// Optional regular expression used to ID control guides. Pattern is matched, case
    /// insensitive, to guide IDs and Gene names.
    #[clap(long, short = 'C')]
    control_pattern: Option<String>,

    /// Perform exact matching only, don't allow mismatches between reads and guides.
    #[clap(long, short = 'x')]
    exact_match: bool,

    /// The number of reads to be examined when determining the offsets at which guides may
    /// be found in the input reads.
    #[clap(long, short = 'N', default_value = "100000")]
    offset_sample_size: u64,

    /// After sampling the first `offset_sample_size` reads, use offsets that
    #[clap(long, short = 'f', default_value = "0.0025")]
    offset_min_fraction: f64,

    /// Path prefix to use for all output files
    #[clap(long, short = 'o')]
    output: String,
}

/// Simple Command impl that just receives params and delegates off to other functions
impl Command for Count {
    /// execute function that is called from the command line parser
    fn execute(&self) -> Result<()> {
        // Auto-fill the sample names if not given
        let sample_ids = if self.samples.is_empty() {
            self.input
                .iter()
                .enumerate()
                .map(|(idx, fq)| Count::sample_name(fq, idx + 1))
                .collect_vec()
        } else {
            assert_eq!(
                self.samples.len(),
                self.input.len(),
                "Different numbers of --samples and --input."
            );
            self.samples.clone()
        };

        // Load up the library and guide lookup
        let library = GuideLibrary::from_files(
            &self.library,
            &self.essential_genes,
            &self.nonessential_genes,
            &self.control_guides,
            &self.control_pattern,
        )?;
        let lookup = Count::build_lookup(&library, !self.exact_match);

        // Generate the counts per sample
        let results = self
            .input
            .iter()
            .zip(sample_ids)
            .map(|(fq, sample)| {
                let prefix_info = Count::determine_prefixes(
                    fq,
                    sample.as_str(),
                    &library,
                    &lookup,
                    self.offset_sample_size,
                    self.offset_min_fraction,
                )
                .expect("Failed to determine offsets.");

                Count::count_reads(fq, sample.as_str(), &library, &lookup, &prefix_info)
                    .expect("Failed to count guide.")
            })
            .collect_vec();

        // Write the outputs
        let counts_file = PathBuf::from(format!("{}.counts.txt", self.output));
        let ext_counts_file = PathBuf::from(format!("{}.extended-counts.txt", self.output));
        let stats_file = PathBuf::from(format!("{}.stats.txt", self.output));

        Count::write_counts(&counts_file, &library, &results, false)?;
        Count::write_counts(&ext_counts_file, &library, &results, true)?;
        Count::write_stats(&stats_file, &library, &results)?;
        Ok(())
    }
}

/// Implementation of the Count command and related functions.
impl Count {
    /// Returns a sample name given a fastq file. Strips off any .gz and fastq-like
    /// suffixes.  If the file doesn't have a valid filename, will return a name
    /// based on the index passed in.
    fn sample_name(p: &Path, idx: usize) -> String {
        if let Some(os_name) = p.file_name() {
            if let Some(name) = os_name.to_str() {
                return name
                    .trim_end_matches(".gz")
                    .trim_end_matches(".fastq")
                    .trim_end_matches(".fq")
                    .to_string();
            }
        }

        format!("s{}", idx)
    }

    /// Builds a lookup from a Vec<u8> of bases to Guides.  The resulting HashMap will contain
    /// keys for every exact guide sequence (in upper case).  If `allow_mismatch` is true,
    /// the map will also contain keys for every one-mismatch version of every guide with the
    /// exception of any sequences that match equally well to multiple guides.
    fn build_lookup(library: &GuideLibrary, allow_mismatch: bool) -> GuideMap {
        info!("Building lookup.");
        let mut lookup = GuideMap::default();
        let mut dupes = HashSet::new();

        lookup.reserve(library.len());

        if allow_mismatch {
            lookup.reserve(library.len() + library.guide_length * 3);

            for guide in library.guides.iter() {
                let bases = &guide.bases;

                for i in 0..bases.len() {
                    for b in [b'A', b'C', b'G', b'T'] {
                        if bases[i] != b {
                            let mut modded = bases.clone();
                            modded[i] = b;

                            let prev = lookup.insert(modded, guide);
                            if prev.is_some() {
                                let mut dupe = bases.clone();
                                dupe[i] = b;
                                dupes.insert(dupe);
                            }
                        }
                    }
                }
            }

            // Make sure no duplicated sequences remain in the lookup
            for dupe in dupes.into_iter() {
                lookup.remove(&dupe);
            }
        }

        // Insert all the exact matches last so they're always present
        for guide in library.guides.iter() {
            lookup.insert(guide.bases.clone(), guide);
        }

        info!("Lookup built with {} entries.", lookup.len());
        lookup
    }

    /// Goes through an input fastq file and determines the set of prefix-lengths that
    /// occur before the guide sequence is observed.  Samples the first `sample_size` reads
    /// from the FASTQ and checks all possible prefixes.  Returns the set of prefixes where
    /// each prefix individually accounts for >= `min_fraction` of the reads that matched
    /// to a guide.
    fn determine_prefixes(
        fastq: &Path,
        sample: &str,
        library: &GuideLibrary,
        lookup: &GuideMap,
        sample_size: u64,
        min_fraction: f64,
    ) -> Result<PrefixInfo> {
        let guide_length = library.guide_length;
        let mut prefix_lengths = vec![0u64; 500];
        let mut count = 0u64;

        // Parse the first `sample_size` records to find exact match guides and
        // extract the sequence that precedes the guide
        parse_path(Some(fastq), |parser| {
            parser
                .each(|rec| {
                    let read_bases = rec.seq();
                    let read_length = read_bases.len();

                    if read_length >= guide_length {
                        for trim in 0..=(read_length - guide_length) {
                            let bases = &read_bases[trim..trim + guide_length];

                            if lookup.contains_key(bases) {
                                prefix_lengths[trim] += 1;
                            }
                        }
                    }

                    count += 1;
                    count < sample_size
                })
                .expect("Failed to parse.");
        })
        .context(format!("Failed to read {:?}", fastq))?;

        let total_matched: u64 = prefix_lengths.iter().sum();
        let fraction_matched = total_matched as f64 / count as f64;
        info!(
            "In {:?} examined {} reads for guide start position and matched {} ({:.4}).",
            fastq, count, total_matched, fraction_matched
        );

        // Tuple of offset -> count where count is > 0
        let non_zeros =
            prefix_lengths.iter().copied().enumerate().filter(|(_idx, n)| *n > 0).collect_vec();

        info!(
            "{} read offsets: {}",
            sample,
            non_zeros.iter().map(|(o, n)| format!("{}->{}", o, n)).join(", ")
        );

        // Filter to just those trim lengths that have at least min_fraction of the data each
        let trims_to_return: Vec<usize> = non_zeros
            .into_iter()
            .filter(|(_idx, n)| *n as f64 / total_matched as f64 >= min_fraction)
            .map(|(idx, _n)| idx)
            .collect();

        let info = PrefixInfo { lengths: trims_to_return };
        Ok(info)
    }

    /// Generates a set of guide counts for a single input FASTQ given a guide lookup and a set of
    /// read offsets/prefixes to check.
    ///
    /// Returns a CountResult which contains a count of the total number of reads in the FASTQ
    /// and a Map of Guide to count of tht guide.  The map will contain an entry for every guide
    /// including those with zero counts.
    fn count_reads<'a, P>(
        fastq: &P,
        sample: &str,
        library: &'a GuideLibrary,
        lookup: &GuideMap,
        prefix_info: &PrefixInfo,
    ) -> Result<CountResult<'a>>
    where
        P: AsRef<Path>,
    {
        let mut count: u64 = 0;
        let mut counts: Vec<u64> = vec![0; library.len()];
        let fastq_path = fastq.as_ref();

        // TODO: remove empty file checking once this is moved over to seq_io instead of fastq
        let fastq_size = std::fs::metadata(fastq).map(|m| m.len()).unwrap_or(0);
        let empty_fastq = fastq_path.is_file() && fastq_path.exists() && fastq_size == 0;

        if !empty_fastq {
            parse_path(Some(fastq), |parser| {
                parser
                    .each(|rec| {
                        let read_bases = rec.seq();
                        let read_length = read_bases.len();
                        let guide_length = library.guide_length;

                        for trim in prefix_info.lengths.iter() {
                            if trim + guide_length <= read_length {
                                let bases = &read_bases[*trim..(*trim + guide_length)];
                                if let Some(guide) = lookup.get(bases) {
                                    counts[guide.index] += 1;
                                }
                            }
                        }

                        count += 1;
                        if count % 10_000_000 == 0 {
                            info!("Processed {}m reads from {:?}.", count / 1_000_000, fastq_path);
                        }

                        true
                    })
                    .expect("Failed to parse.");
            })?;
        }

        let count_map: HashMap<&Guide, u64> =
            library.guides.iter().map(|g| (g, counts[g.index])).collect();

        let count_result = CountResult {
            source: fastq_path.to_str().unwrap_or("").to_string(),
            sample: sample.to_string(),
            counts: count_map,
            total_reads: count,
        };

        info!(
            "Processed {} reads and matched {} ({:.4}) from {:?}.",
            count,
            count_result.mapped_reads(),
            count_result.mapped_frac(),
            fastq_path
        );

        Ok(count_result)
    }

    /// Writes out the counts matrix given one or more CountResults. If extended is false, the
    /// columns produced are "guide", "gene" and then one column per sample with counts for each
    /// sample.  If extended is true, an additional "guide_type" is inserted after "gene".
    fn write_counts(
        path: &Path,
        library: &GuideLibrary,
        counts: &[CountResult],
        extended: bool,
    ) -> Result<()> {
        let mut writer = Io::default().new_writer(&path)?;
        let sep = "\t".as_bytes();
        let newline = "\n".as_bytes();

        // Output the header
        let mut header_fields = vec!["guide", "gene"];
        if extended {
            header_fields.extend_from_slice(&["guide_type"]);
        }

        header_fields.extend(counts.iter().map(|c| c.sample.as_str()));
        writer.write_all(header_fields.join("\t").as_bytes())?;
        writer.write_all(newline)?;

        // Output the body
        for guide in library.guides.iter() {
            writer.write_all(guide.id.as_bytes())?;
            writer.write_all(sep)?;
            writer.write_all(guide.gene.as_bytes())?;

            if extended {
                writer.write_all(sep)?;
                writer.write_all(guide.kind.to_string().as_bytes())?;
            }

            for sample in counts {
                let n = sample.counts.get(&guide).copied().unwrap_or(0);
                writer.write_all(sep)?;
                writer.write_all(n.to_string().as_bytes())?;
            }

            writer.write_all(newline)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Generates and writes out some simple per-sample QC statistics
    fn write_stats(
        stats_file: &Path,
        library: &GuideLibrary,
        results: &[CountResult],
    ) -> Result<()> {
        let recs = results
            .iter()
            .map(|r| CountStats {
                file: r.source.clone(),
                label: r.sample.to_string(),
                total_guides: library.len() as u64,
                total_reads: r.total_reads,
                mapped_reads: r.mapped_reads(),
                frac_mapped: Count::round(r.mapped_frac(), 4),
                mean_reads_per_guide: Count::round(
                    r.mapped_reads() as f64 / library.len() as f64,
                    2,
                ),
                zero_read_guides: r.counts.values().filter(|n| **n == 0).count() as u64,
                mean_reads_essential: Count::compute_mean_cov(r, GuideType::Essential),
                mean_reads_nonessential: Count::compute_mean_cov(r, GuideType::Nonessential),
                mean_reads_other: Count::compute_mean_cov(r, GuideType::Other),
                mean_reads_control: Count::compute_mean_cov(r, GuideType::Control),
            })
            .collect_vec();

        DelimFile::default().write_tsv(&stats_file, recs)?;
        Ok(())
    }

    /// Simple method to round f64s to a maximum number of decimal places
    fn round(f: f64, dp: i32) -> f64 {
        let factor = 10f64.powi(dp);
        (f * factor).round() / factor
    }

    /// Computes the mean coverage of a subset of guides based on guide type. If no guides of the
    /// type exist, returns 0.
    fn compute_mean_cov(counts: &CountResult, kind: GuideType) -> f64 {
        let subset = counts
            .counts
            .iter()
            .filter_map(|(g, n)| if g.kind == kind { Some(*n as f64) } else { None })
            .collect_vec();

        if subset.is_empty() {
            0.0
        } else {
            let total: f64 = subset.iter().sum::<f64>();
            Count::round(total / subset.len() as f64, 2)
        }
    }
}

/// Struct to hold information about the length and bases of the sequence before the guides
struct PrefixInfo {
    pub lengths: Vec<usize>,
}

/// Struct to hold the results of counting a single fastq
struct CountResult<'a> {
    source: String,
    sample: String,
    total_reads: u64,
    counts: HashMap<&'a Guide, u64>,
}

/// Struct to output stats on each sample mapped
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
struct CountStats {
    file: String,
    label: String,
    total_guides: u64,
    total_reads: u64,
    mapped_reads: u64,
    frac_mapped: f64,
    mean_reads_per_guide: f64,
    mean_reads_essential: f64,
    mean_reads_nonessential: f64,
    mean_reads_control: f64,
    mean_reads_other: f64,
    zero_read_guides: u64,
}

impl<'a> CountResult<'a> {
    /// Returns the total number of reads that mapped to a guide.
    pub fn mapped_reads(&self) -> u64 {
        self.counts.values().sum()
    }

    /// Returns the fraction of reads that mapped to a guide.
    pub fn mapped_frac(&self) -> f64 {
        self.mapped_reads() as f64 / self.total_reads as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgoxide::io::{DelimFile, Io};
    use tempfile::TempDir;

    #[test]
    fn test_sample_name() {
        assert_eq!(Count::sample_name(PathBuf::from("/foo/splat.fq").as_path(), 1), "splat");
        assert_eq!(Count::sample_name(PathBuf::from("/foo/splat.fq.gz").as_path(), 1), "splat");
        assert_eq!(Count::sample_name(PathBuf::from("/foo/splat.fastq").as_path(), 1), "splat");
        assert_eq!(Count::sample_name(PathBuf::from("/foo/splat.fastq.gz").as_path(), 1), "splat");
        assert_eq!(Count::sample_name(PathBuf::new().as_path(), 1), "s1");
    }

    #[test]
    fn test_build_lookup() {
        let g1 = Guide::new(0, "g0", "AAAAAAAAAA", "gene-A", GuideType::Other);
        let g2 = Guide::new(1, "g1", "GGGGGGGGGG", "gene-G", GuideType::Other);
        let g3 = Guide::new(2, "g2", "AGGGGGGGGG", "gene-AG", GuideType::Other);

        // Build with one guide and no mismatches
        let library = GuideLibrary::new(vec![g1.clone()]).unwrap();
        let lookup = Count::build_lookup(&library, false);
        assert_eq!(lookup.len(), 1);
        assert_eq!(lookup[&g1.bases], &g1);

        // One guide with mismatches
        let lookup = Count::build_lookup(&library, true);
        assert_eq!(lookup.len(), 31); // original plus three mismatches x ten positions
        assert_eq!(lookup[&g1.bases], &g1);
        assert_eq!(lookup["AAAACAAAAA".as_bytes()], &g1);

        // Two guides without mismatches
        let library = GuideLibrary::new(vec![g1.clone(), g2.clone()]).unwrap();
        let lookup = Count::build_lookup(&library, false);
        assert_eq!(lookup.len(), 2);
        assert_eq!(lookup[&g1.bases], &g1);
        assert_eq!(lookup[&g2.bases], &g2);

        // Two guides with mismatches and no collisions
        let lookup = Count::build_lookup(&library, true);
        assert_eq!(lookup.len(), 62);

        // Two guides with mismatches and collisions!
        let library = GuideLibrary::new(vec![g2.clone(), g3.clone()]).unwrap();
        let lookup = Count::build_lookup(&library, true);
        assert_eq!(lookup.len(), 56);
        assert_eq!(lookup[&g2.bases], &g2); // collision shouldn't override perfect match
        assert_eq!(lookup[&g3.bases], &g3); // collision shouldn't override perfect match
        assert!(!lookup.contains_key("CGGGGGGGGG".as_bytes())); // ambiguous
        assert!(!lookup.contains_key("TGGGGGGGGG".as_bytes())); // ambiguous
    }

    /// Helper function to generate guide library determine_prefixes and counting tests
    fn test_library() -> GuideLibrary {
        GuideLibrary::new(vec![
            Guide::new(0, "g1", "ACGTACGT", "AAA", GuideType::Other),
            Guide::new(1, "g1", "AACCGGTT", "AAA", GuideType::Other),
            Guide::new(2, "g1", "AAACAGAT", "AAA", GuideType::Other),
        ])
        .unwrap()
    }

    /// Helper function to write a bunch of reads to a FASTQ file
    fn write_fastq(reads: &[String], path: PathBuf) -> PathBuf {
        let mut file = Io::default().new_writer(&path).unwrap();

        for (idx, read) in reads.iter().enumerate() {
            let quals = vec![b'#'; read.len()];

            file.write_all(format!("@q{}\n", idx).as_bytes()).unwrap();
            file.write_all(format!("{}\n", read).as_bytes()).unwrap();
            file.write_all("+\n".as_bytes()).unwrap();
            file.write_all(quals.as_slice()).unwrap();
            file.write_all("\n".as_bytes()).unwrap();
        }

        file.flush().unwrap();
        path
    }

    #[test]
    fn test_determine_prefixes_finds_offset_zero() {
        let library = test_library();
        let lookup = Count::build_lookup(&library, false);
        let reads = vec![
            format!("{}tttttttttt", library.guides[0].bases_str),
            format!("{}tttttttttt", library.guides[1].bases_str),
            format!("{}tttttttttt", library.guides[2].bases_str),
        ];
        let tempdir = TempDir::new().unwrap();
        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));
        let prefixes = Count::determine_prefixes(&fastq, "s", &library, &lookup, 100, 0.0).unwrap();
        assert_eq!(prefixes.lengths, vec![0]);
    }

    #[test]
    fn test_determine_prefixes_finds_last_offset() {
        let library = test_library();
        let lookup = Count::build_lookup(&library, false);
        let reads = vec![
            format!("tttttttttt{}", library.guides[0].bases_str),
            format!("tttttttttt{}", library.guides[1].bases_str),
            format!("tttttttttt{}", library.guides[2].bases_str),
        ];
        let tempdir = TempDir::new().unwrap();
        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));
        let prefixes = Count::determine_prefixes(&fastq, "s", &library, &lookup, 100, 0.0).unwrap();
        assert_eq!(prefixes.lengths, vec![10]);
    }

    #[test]
    fn test_determine_prefixes_actually_subsamples() {
        let library = test_library();
        let lookup = Count::build_lookup(&library, false);

        let mut reads = vec![];
        for _ in 0..100 {
            reads.push(format!("ttttt{}ttttt", library.guides[0].bases_str))
        }
        for _ in 0..100 {
            reads.push(format!("tttttt{}tttt", library.guides[0].bases_str))
        }

        let tempdir = TempDir::new().unwrap();
        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));
        let p1 = Count::determine_prefixes(&fastq, "s", &library, &lookup, 200, 0.01).unwrap();
        let p2 = Count::determine_prefixes(&fastq, "s", &library, &lookup, 100, 0.01).unwrap();
        assert_eq!(p1.lengths, vec![5, 6]);
        assert_eq!(p2.lengths, vec![5]);
    }

    #[test]
    fn test_determine_prefixes_min_fraction() {
        let library = test_library();
        let lookup = Count::build_lookup(&library, false);

        let mut reads = vec![];
        for _ in 0..50 {
            // offset = 5
            reads.push(format!("ttttt{}ttttt", library.guides[0].bases_str))
        }
        for _ in 0..30 {
            // offset = 6
            reads.push(format!("tttttt{}tttt", library.guides[0].bases_str))
        }
        for _ in 0..20 {
            // offset = 7
            reads.push(format!("ttttttt{}ttt", library.guides[0].bases_str))
        }
        for _ in 0..100 {
            // doesn't match to guides
            reads.push("tttttttttttttttttttt".to_string());
        }

        let tempdir = TempDir::new().unwrap();
        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));
        let p1 = Count::determine_prefixes(&fastq, "s", &library, &lookup, 200, 0.01).unwrap();
        let p2 = Count::determine_prefixes(&fastq, "s", &library, &lookup, 200, 0.25).unwrap();
        let p3 = Count::determine_prefixes(&fastq, "s", &library, &lookup, 200, 0.50).unwrap();

        assert_eq!(p1.lengths, vec![5, 6, 7]);
        assert_eq!(p2.lengths, vec![5, 6]);
        assert_eq!(p3.lengths, vec![5]);
    }

    #[test]
    fn test_count_reads_handles_empty_fastq() {
        let library = test_library();
        let lookup = Count::build_lookup(&library, false);
        let prefixes = PrefixInfo { lengths: vec![0, 1, 2] };
        let tempdir = TempDir::new().unwrap();
        let fastq = write_fastq(&[], tempdir.path().join("in.fastq"));

        let counts = Count::count_reads(&fastq, "s1", &library, &lookup, &prefixes).unwrap();
        assert_eq!(counts.sample, "s1");
        assert_eq!(counts.total_reads, 0);
        for guide in library.guides.iter() {
            assert_eq!(counts.counts[guide], 0);
        }
    }

    #[test]
    fn test_count_reads() {
        let library = test_library();
        let lookup = Count::build_lookup(&library, false);
        let prefixes = PrefixInfo { lengths: vec![4, 5, 6] };
        let tempdir = TempDir::new().unwrap();
        let mut reads = vec![];

        // Create a bunch of reads, with each guide getting:
        //   100 * 3 * 1_based_guide_index
        for prefix in ["tttt", "ttttt", "tttttt"] {
            for _ in 0..100 {
                for guide in library.guides.iter() {
                    for _ in 0..(guide.index + 1) {
                        reads.push(format!("{}{}ggggg", prefix, guide.bases_str));
                    }
                }
            }
        }

        // Add a few reads at different offsets that won't count
        reads.push(format!("{}ttttttttt", library.guides[0].bases_str));
        reads.push(format!("t{}tttttttt", library.guides[0].bases_str));
        reads.push(format!("tt{}ttttttt", library.guides[0].bases_str));
        reads.push(format!("ttt{}tttttt", library.guides[0].bases_str));

        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));
        let counts = Count::count_reads(&fastq, "s1", &library, &lookup, &prefixes).unwrap();
        assert_eq!(counts.sample, "s1");
        assert_eq!(counts.total_reads, reads.len() as u64);
        for guide in library.guides.iter() {
            let expected = 100 * 3 * (guide.index + 1) as u64;
            assert_eq!(counts.counts[guide], expected);
        }
    }

    #[test]
    fn test_end_to_end() {
        let tempdir = TempDir::new().unwrap();

        // Write the library to disk
        let library = test_library();
        let library_path = tempdir.path().join("library.txt");
        let mut lib_lines = vec!["guide\tbases\tgene".to_string()];
        for guide in library.guides.iter() {
            lib_lines.push(format!("{}\t{}\t{}", guide.id, guide.bases_str, guide.gene));
        }
        Io::default().write_lines(&library_path, &lib_lines).unwrap();

        // Generate a fastq of reads to count
        let mut reads = vec![];

        // Create a bunch of reads, with each guide getting:
        //   100 * 3 * 1_based_guide_index
        for prefix in ["tttt", "ttttt", "tttttt"] {
            for _ in 0..100 {
                for guide in library.guides.iter() {
                    for _ in 0..(guide.index + 1) {
                        reads.push(format!("{}{}ggggg", prefix, guide.bases_str));
                    }
                }
            }
        }

        // Add a few reads at different offsets that won't count
        reads.push(format!("{}ttttttttt", library.guides[0].bases_str));
        reads.push(format!("t{}tttttttt", library.guides[0].bases_str));
        reads.push(format!("tt{}ttttttt", library.guides[0].bases_str));
        reads.push(format!("ttt{}tttttt", library.guides[0].bases_str));

        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));

        // Run the count command
        let prefix = tempdir.path().join("out").to_str().unwrap().to_string();
        let counts = tempdir.path().join("out.counts.txt");
        let stats = tempdir.path().join("out.stats.txt");

        let cmd = Count {
            library: library_path,
            input: vec![fastq],
            essential_genes: None,
            nonessential_genes: None,
            control_guides: None,
            control_pattern: None,
            samples: vec!["sample1".to_string()],
            output: prefix,
            exact_match: false,
            offset_min_fraction: 0.005,
            offset_sample_size: 100000,
        };

        cmd.execute().unwrap();

        assert!(counts.exists());
        assert!(stats.exists());

        // Read the stats back in
        let stat_records: Vec<CountStats> = DelimFile::default().read_tsv(&stats).unwrap();

        assert_eq!(stat_records.len(), 1);
        assert_eq!(stat_records[0].label, "sample1");
        assert_eq!(stat_records[0].total_guides, 3);
        assert_eq!(stat_records[0].total_reads, 1804);
        assert_eq!(stat_records[0].mapped_reads, 1800);
        assert!((stat_records[0].mean_reads_per_guide - 600.0).abs() <= 0.01);
        assert!((stat_records[0].frac_mapped - 1800f64 / 1804f64).abs() <= 0.01);
        assert_eq!(stat_records[0].zero_read_guides, 0);
    }

    #[test]
    fn test_reads_shorter_than_guides_ok() {
        let tempdir = TempDir::new().unwrap();

        // Write the library to disk
        let library = test_library();
        let library_path = tempdir.path().join("library.txt");
        let mut lib_lines = vec!["guide\tbases\tgene".to_string()];
        for guide in library.guides.iter() {
            lib_lines.push(format!("{}\t{}\t{}", guide.id, guide.bases_str, guide.gene));
        }
        Io::default().write_lines(&library_path, &lib_lines).unwrap();

        // Generate a fastq of reads to count
        let mut reads = vec![];
        reads.push("A".to_string());
        reads.push("AC".to_string());
        reads.push("ACG".to_string());
        reads.push("ACGT".to_string());
        reads.push("ACGTA".to_string());
        reads.push("ACGTAC".to_string());
        reads.push("ACGTACG".to_string());
        reads.push("ACGTACGT".to_string());
        reads.push("ACGTACGTA".to_string());
        reads.push("ACGTACGTAC".to_string());

        // Create a bunch of reads all with the guide at offset=5, with each guide getting:
        for _ in 0..100 {
            for guide in library.guides.iter() {
                for _ in 0..=guide.index {
                    reads.push(format!("ttttt{}ggggg", guide.bases_str));
                }
            }
        }

        // Lastly create some reads that are longer than a guide length but shorter than
        // guide-length + max prefix
        reads.push("tttttACGTACGT".to_string());
        reads.push("tttttACGTACGTA".to_string());
        reads.push("tttttACGTACGTAC".to_string());

        let fastq = write_fastq(&reads, tempdir.path().join("in.fastq"));

        // Run the count command
        let prefix = tempdir.path().join("out").to_str().unwrap().to_string();
        let counts = tempdir.path().join("out.counts.txt");
        let stats = tempdir.path().join("out.stats.txt");

        let cmd = Count {
            library: library_path,
            input: vec![fastq],
            essential_genes: None,
            nonessential_genes: None,
            control_guides: None,
            control_pattern: None,
            samples: vec!["sample1".to_string()],
            output: prefix,
            exact_match: false,
            offset_min_fraction: 0.005,
            offset_sample_size: 100000,
        };

        cmd.execute().unwrap();
        assert!(counts.exists());
        assert!(stats.exists());
    }
}
