use anyhow::{anyhow, bail, Result};
use fgoxide::io::Io;
use itertools::Itertools;
use log::*;
use regex::RegexBuilder;
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::path::Path;


/// Guides can either target essential genes, non-essential genes, control sequences
/// or other genes.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum GuideType {
    Essential,
    Nonessential,
    Control,
    Other,
}

/// Implement Disaply for GuideType
impl Display for GuideType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// A struct to represent a CRISPR guide
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Guide {
    pub index: usize,
    pub id: String,
    pub kind: GuideType,
    pub bases: Vec<u8>,
    pub bases_str: String,
    pub gene: String,
}

impl Guide {
    /// Generates a new guide that stores the bases in upper-case and ensures that the
    /// bases and bases_str are in sync.
    pub fn new<S: Into<String>>(index: usize, id: S, bases: S, gene: S, kind: GuideType) -> Guide {
        let bases_upper = bases.into().to_uppercase();
        Guide {
            index,
            id: id.into(),
            kind,
            bases: bases_upper.as_bytes().to_vec(),
            bases_str: bases_upper,
            gene: gene.into(),
        }
    }

    /// Returns the length of the guide sequence
    fn len(&self) -> usize {
        self.bases.len()
    }
}

/// Struct representing a guide library used in a crispr screen
pub struct GuideLibrary {
    pub guides: Vec<Guide>,
    pub guide_length: usize,
}

impl GuideLibrary {
    /// Constructs a new guide library from a set of guides.  Will return an error if:
    ///   - the guides have a mix of lengths
    ///   - the guide sequences are not all A/C/G/T
    ///   - there are non-unique guide sequences
    pub fn new(guides: Vec<Guide>) -> Result<GuideLibrary> {
        let lengths: HashSet<usize> = guides.iter().map(|g| g.len()).collect();
        let unique: HashSet<&Vec<u8>> = guides.iter().map(|g| &g.bases).collect();
        let genes: HashSet<&str> = guides.iter().map(|g| g.gene.as_str()).collect();
        let bad = guides.iter().filter(|g| !GuideLibrary::is_acgt(&g.bases)).collect_vec();

        if guides.is_empty() {
            Ok(GuideLibrary { guides, guide_length: 0 })
        } else if lengths.len() != 1 {
            Err(anyhow!("More than one guide length found: {}.", lengths.iter().join(", ")))
        } else if !bad.is_empty() {
            Err(anyhow!("{} guides had non-ACGT bases in their sequence.", bad.len()))
        } else if unique.len() < guides.len() {
            Err(anyhow!(
                "Guide library had {} guides but only {} unique sequences.",
                guides.len(),
                unique.len()
            ))
        } else {
            info!(
                "Loaded library with {} guides for {} genes; {}=essential, {}=nonessential, {}=control, {}=other.",
                guides.len(),
                genes.len(),
                guides.iter().filter(|g| g.kind == GuideType::Essential).count(),
                guides.iter().filter(|g| g.kind == GuideType::Nonessential).count(),
                guides.iter().filter(|g| g.kind == GuideType::Control).count(),
                guides.iter().filter(|g| g.kind == GuideType::Other).count(),
            );

            Ok(GuideLibrary {
                guides,
                guide_length: *lengths.iter().next().expect("Where'd it go?"),
            })
        }
    }

    /// Reads a guide library from a file.  The file:
    ///   - May be either tab- or comma-delimited
    ///   - Must have a header row
    ///   - Must have at least three columns with the first three columns in order being:
    ///     - a unique ID for the guide
    ///     - the sequence of the guide
    ///     - the gene, or other target, of the guide
    pub fn from_file<P>(path: &P) -> Result<GuideLibrary>
    where
        P: AsRef<Path>,
    {
        let no_pattern: Option<&str> = None;
        GuideLibrary::from_files(path, &None, &None, &None, &no_pattern)
    }

    pub fn from_files<P, S>(
        lib_path: &P,
        essential_gene_path: &Option<P>,
        non_essential_gene_path: &Option<P>,
        control_guide_list_path: &Option<P>,
        control_pattern: &Option<S>,
    ) -> Result<GuideLibrary>
    where
        P: AsRef<Path>,
        S: AsRef<str>,
    {
        let essentials = GuideLibrary::read_to_set(essential_gene_path)?;
        let non_essentials = GuideLibrary::read_to_set(non_essential_gene_path)?;
        let control_guides = GuideLibrary::read_to_set(control_guide_list_path)?;
        let control_regex = if let Some(p) = control_pattern {
            Some(RegexBuilder::new(p.as_ref()).case_insensitive(true).build()?)
        } else {
            None
        };

        let lines = Io::default().read_lines(lib_path)?;

        if lines.len() < 2 {
            GuideLibrary::new(vec![])
        } else {
            let delim: char = if lines[0].chars().filter(|ch| *ch == '\t').count() >= 2 {
                '\t'
            } else if lines[0].chars().filter(|ch| *ch == ',').count() >= 2 {
                ','
            } else {
                bail!("couldn't detect delimiter from first line of {:?}", lib_path.as_ref());
            };

            // Read in the guides
            let mut guides = Vec::with_capacity(1024);
            let mut idx: usize = 0;
            for line in lines.iter().skip(1) {
                let trimmed = line.trim();

                if !trimmed.is_empty() {
                    let fields = trimmed.split(delim).collect_vec();
                    if fields.len() < 3 {
                        bail!("Too few fields in line: '{}'", line);
                    }

                    let guide_id = fields[0];
                    let bases = fields[1];
                    let gene = fields[2];

                    let kind = if essentials.contains(gene) {
                        GuideType::Essential
                    } else if non_essentials.contains(gene) {
                        GuideType::Nonessential
                    } else if control_guides.contains(guide_id)
                        || control_regex
                            .as_ref()
                            .filter(|re| re.is_match(guide_id) || re.is_match(gene))
                            .is_some()
                    {
                        GuideType::Control
                    } else {
                        GuideType::Other
                    };

                    let guide = Guide::new(idx, guide_id, bases, gene, kind);
                    guides.push(guide);
                    idx += 1;
                }
            }

            GuideLibrary::new(guides)
        }
    }

    /// Reads all lines from a file, trims them, extracts the first tab-separated field and then
    /// returns the unique set of values as a set.
    fn read_to_set<P: AsRef<Path>>(path: &Option<P>) -> Result<HashSet<String>> {
        let items: HashSet<String> = match path {
            None => HashSet::new(),
            Some(p) => Io::default()
                .read_lines(p)?
                .into_iter()
                .map(|line| line.trim().to_string())
                .filter(|line| !line.is_empty())
                .map(|line| line.split('\t').next().unwrap().to_string())
                .collect(),
        };

        Ok(items)
    }

    /// Returns true if the sequence is all upper-case ACGT bases
    fn is_acgt(bases: &[u8]) -> bool {
        bases.iter().copied().all(|b| b == b'A' || b == b'C' || b == b'G' || b == b'T')
    }

    /// Returns the number of guides in the library
    pub fn len(&self) -> usize {
        self.guides.len()
    }

    /// True if there are no guides in the library, false otherwise
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use tempfile::TempDir;

    const CSV_LIBRARY: &str = "\
id,bases,gene
a1,CGATCGCTTAAGCTAGCA,FOO
a2,ATGCTAGATCGCGCTATT,FOO
a3,GGCTTCTAGATCGCTATA,Control
";

    const TSV_LIBRARY: &str = "\
id\tbases\tgene
b1\tCGATCGCTTAAGCTAGCA\tFOO
b2\tATGCTAGATCGCGCTATT\tFOO
b3\tGGCTTCTAGATCGCTATA\tControl
";

    #[test]
    fn test_guide_uppercases_sequence() {
        let g1 = Guide::new(0, "foo-1", "AAAAACCCCCGGGGGTTTTT", "FOO", GuideType::Other);
        assert_eq!(g1.index, 0);
        assert_eq!(g1.id, "foo-1");
        assert_eq!(g1.bases_str, "AAAAACCCCCGGGGGTTTTT");
        assert_eq!(g1.bases, "AAAAACCCCCGGGGGTTTTT".as_bytes());
        assert_eq!(g1.gene, "FOO");
        assert_eq!(g1.len(), 20);

        let g2 = Guide::new(0, "foo-2", "aAaAcCcCgGgGtTtTacgt", "FOO", GuideType::Other);
        assert_eq!(g2.index, 0);
        assert_eq!(g2.id, "foo-2");
        assert_eq!(g2.bases_str, "AAAACCCCGGGGTTTTACGT");
        assert_eq!(g2.bases, "AAAACCCCGGGGTTTTACGT".as_bytes());
        assert_eq!(g2.gene, "FOO");
        assert_eq!(g2.len(), 20);
    }

    #[test]
    fn test_is_all_acgt() {
        // True cases
        assert!(GuideLibrary::is_acgt("".as_bytes()));
        assert!(GuideLibrary::is_acgt("AACGCTGACTGA".as_bytes()));

        // False cases
        assert!(!GuideLibrary::is_acgt("N".as_bytes()));
        assert!(!GuideLibrary::is_acgt("AC GT".as_bytes()));
        assert!(!GuideLibrary::is_acgt("AC-GT".as_bytes()));
        assert!(!GuideLibrary::is_acgt("acgt".as_bytes()));
    }

    #[test]
    fn test_guide_library_positive() {
        let g1 = Guide::new(0, "foo-1", "ACGTCAGCATGCATGACGTT", "FOO", GuideType::Other);
        let g2 = Guide::new(1, "foo-2", "GCTAGACTGGACTCTAATGC", "FOO", GuideType::Other);

        let l1 = GuideLibrary::new(vec![]).unwrap();
        assert_eq!(l1.len(), 0);
        assert!(l1.is_empty());

        let l2 = GuideLibrary::new(vec![g1.clone(), g2.clone()]).unwrap();
        assert_eq!(l2.len(), 2);
        assert_eq!(l2.guides, vec![g1, g2]);
        assert_eq!(l2.guide_length, 20);
    }

    #[test]
    fn test_guide_library_rejects_mixed_length() {
        let g1 = Guide::new(0, "foo-1", "ACGTCAGCATGCATGACGTT", "FOO", GuideType::Other);
        let g2 = Guide::new(1, "foo-2", "GCTAGACTGGACTCTAATGCC", "FOO", GuideType::Other);

        let result = GuideLibrary::new(vec![g1, g2]);
        assert!(result.err().unwrap().to_string().contains("More than one guide length found"));
    }

    #[test]
    fn test_guide_library_rejects_invalid_sequences() {
        let g1 = Guide::new(0, "foo-1", "ACGTCAGCANNNATGACGTT", "FOO", GuideType::Other);
        let g2 = Guide::new(1, "foo-2", "hello!", "FOO", GuideType::Other);

        assert!(GuideLibrary::new(vec![g1]).err().unwrap().to_string().contains("non-ACGT"));
        assert!(GuideLibrary::new(vec![g2]).err().unwrap().to_string().contains("non-ACGT"));
    }

    #[test]
    fn test_guide_library_rejects_duplicate_sequences() {
        let g1 = Guide::new(0, "foo-1", "ACGTCAGCATGCATGACGTT", "FOO", GuideType::Other);
        let g2 = Guide::new(1, "foo-2", "ACGTCAGCATGCATGACGTT", "FOO", GuideType::Other);
        assert!(GuideLibrary::new(vec![g1, g2]).err().unwrap().to_string().contains("unique"));
    }

    #[test]
    fn test_reading_guide_library_from_csv_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("lib.csv");
        Io::default().write_lines(&path, vec![CSV_LIBRARY]).unwrap();
        let lib = GuideLibrary::from_file(&path).unwrap();

        assert_eq!(lib.len(), 3);
        assert_eq!(lib.guides.iter().map(|g| g.id.as_str()).collect_vec(), vec!["a1", "a2", "a3"]);
    }

    #[test]
    fn test_reading_guide_library_from_tsv_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("lib.tsv");
        Io::default().write_lines(&path, vec![TSV_LIBRARY]).unwrap();
        let lib = GuideLibrary::from_file(&path).unwrap();

        assert_eq!(lib.len(), 3);
        assert_eq!(lib.guides.iter().map(|g| g.id.as_str()).collect_vec(), vec!["b1", "b2", "b3"]);
    }

    #[test]
    fn test_load_guide_library_from_files() {
        let tmp = TempDir::new().unwrap();
        let lib_path = tmp.path().join("library.tsv.gz");
        let ess_path = tmp.path().join("essential.txt");
        let non_path = tmp.path().join("non-essential.txt");
        let ctl_path = tmp.path().join("control-guides.txt");

        let io = Io::default();
        io.write_lines(
            &lib_path,
            vec![
                "guide\tbases\tgene",
                "g1.1\tAAAAAAAAAA\tG1",
                "g1.2\tCCCCCCCCCC\tG1",
                "g2.1\tGGGGGGGGGG\tG2",
                "g2.2\tTTTTTTTTTT\tG2",
                "g3.1\tACACACACAC\tG3",
                "g3.2\tAGAGAGAGAG\tG3",
                "g4.1\tATATATATAT\tG4",
                "g4.2\tCACACACACA\tG4",
                "c1\tGGAGGAGGAG\tnon-target1",
                "c2\tGGTGGTGGTG\tnon-target2",
                "c3\tAAGTAAGTCC\tcontrol",
            ],
        )
        .unwrap();

        io.write_lines(&ess_path, vec!["G1", "G3"]).unwrap();
        io.write_lines(&non_path, vec!["G4"]).unwrap();
        io.write_lines(&ctl_path, vec!["c1", "c2"]).unwrap();

        let lib = GuideLibrary::from_files(
            &lib_path,
            &Some(ess_path),
            &Some(non_path),
            &Some(ctl_path),
            &Some("Control"),
        )
        .unwrap();

        let map: HashMap<&str, &GuideType> =
            lib.guides.iter().map(|g| (g.id.as_str(), &g.kind)).collect();
        assert_eq!(map["g1.1"], &GuideType::Essential);
        assert_eq!(map["g1.2"], &GuideType::Essential);
        assert_eq!(map["g2.1"], &GuideType::Other);
        assert_eq!(map["g2.2"], &GuideType::Other);
        assert_eq!(map["g3.1"], &GuideType::Essential);
        assert_eq!(map["g3.2"], &GuideType::Essential);
        assert_eq!(map["g4.1"], &GuideType::Nonessential);
        assert_eq!(map["g4.2"], &GuideType::Nonessential);
        assert_eq!(map["c1"], &GuideType::Control);
        assert_eq!(map["c2"], &GuideType::Control);
        assert_eq!(map["c3"], &GuideType::Control);
    }
}
