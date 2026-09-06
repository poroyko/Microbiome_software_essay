# Python and R Packages for Microbiome Analysis — Chapter 8: Outbreak Phylogenetics and Long-Read Assembly Refinement

## Introduction

Chapter 7 established the isolate WGS toolchain through assembly, annotation, resistance detection, serotyping, and pangenome comparison. This chapter extends that toolchain in two directions that arose naturally from Chapter 7's own findings. **Outbreak phylogenetics** picks up directly where Roary and Panaroo left off: once a set of related isolates has been compared at the gene-content level, the next question in an outbreak or transmission investigation is usually finer-grained — which isolates are most closely related by single-nucleotide differences, and how should recombination (which can mimic close relatedness without reflecting recent shared ancestry) be handled before building a tree. This is covered here through **Snippy** (rapid SNP calling), **Gubbins** (recombination detection and removal), and **IQ-TREE** (maximum-likelihood tree building). **Long-read assembly refinement** picks up the thread left open explicitly by Chapter 7's own Unicycler entry, where Unicycler's author was noted to now recommend two newer tools for "long-read-first" assembly given improved Nanopore sequencing quality: **Trycycler** and **Polypolish** are covered here in full. As in Chapter 7, this chapter's research surfaced a striking, self-reinforcing pattern: nearly every tool here was authored by one of two small, closely linked groups (Torsten Seemann's bacterial genomics toolkit, and Ryan Wick/Kathryn Holt's long-read assembly lineage), each of which has been shown across this series to practice unusually explicit, public tool succession.

## Outbreak Phylogenetics: SNP Calling, Recombination, and Tree Building

### Snippy

Snippy rapidly calls variants (SNPs and small insertions/deletions) between a haploid reference genome and a set of sequencing reads or assembled contigs, and — critically for outbreak investigation — can combine the results from many isolates called against the same reference into a single "core SNP" alignment suitable for building a high-resolution phylogeny directly.

**Prerequisites.** Snippy is a Perl-based command-line tool distributed via GitHub and Conda (from the same author, Torsten Seemann, as Prokka in Chapter 5), wrapping BWA-MEM for read alignment, FreeBayes for variant calling, and, if given an annotated reference, SnpEff to determine each variant's predicted effect on genes and other features.

**Popularity.** Snippy has no dedicated peer-reviewed publication — like vegan, the `microbiome` R package, and scikit-bio elsewhere in this series, it is cited as software directly via its GitHub repository — yet it is one of the most widely taught and used tools specifically for bacterial outbreak SNP analysis, appearing in numerous public health bioinformatics training curricula (including Galaxy Training Network materials) as the standard first step before phylogenetic tree building.

**Efficiency.** Snippy's own documentation emphasizes speed as its core design goal (its name is explicitly derived from being "fast" and "SNAPPY"), reporting typical bacterial genome variant-calling runtimes on the order of a few minutes using multiple CPU cores — a design priority directly suited to outbreak investigations, where results are often needed on an urgent, near-real-time basis rather than as part of a leisurely research timeline.

- **Citation:** Cited as software (no dedicated peer-reviewed publication); see [github.com/tseemann/snippy](https://github.com/tseemann/snippy)
- **GitHub:** [github.com/tseemann/snippy](https://github.com/tseemann/snippy)

### Gubbins

Gubbins identifies and removes the effects of recombination from a bacterial whole-genome alignment before phylogenetic analysis, addressing a specific problem that plain SNP-distance approaches (like Snippy's core-SNP alignment above) do not solve on their own: a single recombination event can introduce many SNPs at once, which — if left uncorrected — would be misread by a phylogenetic method as evidence of many independent point mutations, distorting both tree topology and branch lengths.

**Prerequisites.** Gubbins is implemented in Python and C, distributed via GitHub and Conda, and internally depends on a maximum-likelihood tree-building tool for its iterative ancestral reconstruction step — historically RAxML by default, with IQ-TREE (below) supported as an alternative since Gubbins version 3.

**Popularity.** Gubbins has been applied to reconstruct the recent evolutionary history of numerous major bacterial pathogens since its 2015 publication, including *Streptococcus pneumoniae*, *Vibrio cholerae*, *Mycobacterium abscessus*, multiple *Salmonella* serovars, *Escherichia coli*, and *Neisseria gonorrhoeae* — a breadth of application the original publication itself cites as evidence of the tool's applicability across very different bacterial recombination rates and mechanisms.

**Efficiency.** The original publication's simulations report a positive predictive value greater than 99.5% for correctly attributing base substitutions to recombination rather than point mutation, with a false negative rate peaking at just over 3% — figures that directly translate into the tool's stated purpose of significantly improving the accuracy of phylogenetic topology and branch-length estimation relative to methods that ignore recombination entirely.

- **Citation:** Croucher NJ, Page AJ, Connor TR, Delaney AJ, Keane JA, Bentley SD, Parkhill J, Harris SR. Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Res.* 2015;43(3):e15. PMID: [25414349](https://pubmed.ncbi.nlm.nih.gov/25414349/)
- **GitHub:** [github.com/nickjcroucher/gubbins](https://github.com/nickjcroucher/gubbins)

### IQ-TREE

IQ-TREE builds maximum-likelihood phylogenetic trees, integrating fast automatic model selection (ModelFinder), an efficient stochastic tree-search algorithm, and an ultrafast bootstrap approximation for assessing branch support — the general-purpose tree-building engine that both stands alone for outbreak phylogenies (e.g., from a Gubbins-corrected or Snippy-derived core-genome alignment) and, as noted above, can serve as Gubbins' internal reconstruction engine.

**Prerequisites.** IQ-TREE is implemented in C++ and distributed via GitHub, Conda, and a web server interface, requiring no external dependencies for its core maximum-likelihood inference, though its AliSim module for alignment simulation and various downstream analyses can draw on additional data.

**Popularity.** IQ-TREE is, per its own 2020 publication, user-friendly and widely used for phylogenetic inference broadly (not limited to bacterial outbreak contexts), and has continued rapid, active development well beyond the version most commonly cited: a third major version, IQ-TREE 3, was published in 2026, alongside a growing body of methods work on phylogenetic accuracy at "pandemic scale" from an overlapping set of authors.

**Efficiency.** An independent study cited in IQ-TREE's own publication found that its tree-search algorithm achieves good performance in both computing time and likelihood maximization compared with other popular maximum-likelihood phylogenetics software, including RAxML and PhyML, while its ultrafast bootstrap approximation is reported to run 10 to 40 times faster than RAxML's rapid bootstrap while producing less biased support values — two distinct, independently framed speed claims (tree search, and bootstrap support) rather than one blanket performance figure.

- **Citation:** Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol.* 2020;37(5):1530–1534. PMID: [32011700](https://pubmed.ncbi.nlm.nih.gov/32011700/)
- **GitHub:** [github.com/iqtree/iqtree2](https://github.com/iqtree/iqtree2)

## Long-Read Assembly Refinement

### Trycycler

Trycycler produces a consensus long-read bacterial genome assembly by combining multiple independent input assemblies of the same genome (e.g., generated with different assemblers or different read subsets), clustering their contigs and generating a consensus sequence for each cluster — exploiting the fact that while individual long-read assemblies almost always contain some errors, different assemblies of the same genome typically make different errors, so a consensus across several is more accurate than any single one.

**Prerequisites.** Trycycler is distributed via GitHub and Conda, and — unusually for a tool in this series — explicitly requires human judgment and manual intervention at several pipeline steps (clustering review and contig reconciliation), meaning its output is not fully deterministic and different users can produce slightly different assemblies from identical input data.

**Popularity.** Trycycler was developed by the same author (Ryan Wick) and lab as Unicycler (Chapter 7), and its publication explicitly frames it as the recommended long-read-first alternative for the higher-depth, higher-accuracy long-read sequencing now widely available — precisely the succession Chapter 7's Unicycler entry anticipated; a 2025 successor tool, Autocycler, has already been published specifically to automate the manual-intervention steps Trycycler requires, extending consensus assembly to the larger, less hands-on-friendly datasets now common in bacterial genomics.

**Efficiency.** The original publication's benchmarking, using both simulated and real sequencing reads, found that Trycycler consensus assemblies contained fewer errors than assemblies constructed with any single long-read assembler tested, and that Trycycler assemblies combined with post-assembly polishing (using Medaka and Pilon) were the most accurate genomes produced in the study overall — with the explicit caveat, stated candidly by the authors themselves, that because Trycycler requires manual intervention, different users converge on similarly accurate but not identical assemblies from the same data.

- **Citation:** Wick RR, Judd LM, Cerdeira LT, Hawkey J, Méric G, Vezina B, Wyres KL, Holt KE. Trycycler: consensus long-read assemblies for bacterial genomes. *Genome Biol.* 2021;22(1):266. PMID: [34521459](https://pubmed.ncbi.nlm.nih.gov/34521459/)
- **GitHub:** [github.com/rrwick/Trycycler](https://github.com/rrwick/Trycycler)

### Polypolish

Polypolish corrects residual errors in a long-read-only genome assembly using short-read (Illumina) data, specifically targeting errors in repetitive genomic regions that other short-read polishing tools struggle with, by using all valid alignment positions for each short read (rather than only its single best alignment) to disambiguate which copy of a repeated sequence a given read actually supports.

**Prerequisites.** Polypolish is distributed via GitHub and Conda, requiring only a long-read assembly and a set of short reads aligned to it (as a standard SAM file) as input, with no additional reference database.

**Popularity.** Polypolish is explicitly positioned as one component of a broader, published, author-recommended pipeline for assembling bacterial genomes "to perfection" — Trycycler long-read assembly, Medaka long-read polishing, Polypolish short-read polishing, followed by additional short-read polishing tools and manual curation — rather than as a complete standalone solution, and has itself already received a follow-up publication (2024) from an overlapping author group describing a further-updated version (v0.6.0).

**Efficiency.** The original publication reports that Polypolish performed well in benchmarking using both simulated and real reads and was very unlikely to introduce new errors during polishing, while also candidly noting that Polypolish can fix errors other polishers cannot, and vice versa — meaning the authors' own stated best practice is to use Polypolish in combination with other polishing tools rather than as a sole correction step, a direct parallel to Gubbins-plus-IQ-TREE and Snippy-plus-Gubbins above, where this chapter's tools are consistently designed to be combined rather than chosen between.

- **Citation:** Wick RR, Holt KE. Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLoS Comput Biol.* 2022;18(1):e1009802. PMID: [35073327](https://pubmed.ncbi.nlm.nih.gov/35073327/)
- **GitHub:** [github.com/rrwick/Polypolish](https://github.com/rrwick/Polypolish)

## Comparative Assessment

**Two small author groups produced nearly every tool in this chapter, and both groups practice the same unusually explicit succession pattern documented repeatedly across this series.** Snippy and Gubbins both trace to the Sanger Institute bacterial genomics tradition (Snippy by Torsten Seemann, also responsible for Prokka in Chapter 5; Gubbins co-authored by several of the same researchers behind Roary in Chapter 7); Trycycler and Polypolish both trace to Ryan Wick and Kathryn Holt's group, also responsible for Unicycler in Chapter 7. Wick's own tool lineage is now explicit across three chapters of this series: Unicycler (2017, short-read-first) → Trycycler (2021, long-read-first consensus, requiring manual steps) → Autocycler (2025, automating those manual steps) → Hybracter (an overlapping-author tool for fast, automatic long-read-first assembly at scale, referenced directly in Polypolish's own citation network). Few individual research groups in this entire eight-chapter survey have authored as continuously evolving and explicitly self-documented a tool succession as this one.

**This chapter makes the "combine, don't choose" pattern more explicit than any prior chapter.** Snippy's core-SNP output is designed to feed into a tree-building step; Gubbins is designed to sit between that alignment and a final phylogeny specifically to correct for recombination; IQ-TREE can be either that final tree-building step or Gubbins' own internal engine; Trycycler's consensus output is explicitly meant to be polished further (with Medaka, then Polypolish, then still more tools) rather than used directly. Every tool in this chapter is, by its own authors' explicit statements, one stage in a multi-tool pipeline rather than a complete answer on its own — a pattern present throughout this series (DAS Tool, DRAM, VirSorter2→CheckV→DRAM-v) but never before stated this consistently, this many times, by the tools' own original publications.

**Snippy's lack of a peer-reviewed publication is worth setting against Gubbins' and IQ-TREE's rigorous, heavily benchmarked ones.** This is the same citation-practice distinction noted for vegan and scikit-bio in earlier chapters, but here it sits directly beside two tools (Gubbins, IQ-TREE) that are among the most methodologically detailed and independently benchmarked in this entire series. A reader assembling an outbreak phylogenetics pipeline should not read the absence of a Snippy paper as an absence of rigor in the underlying variant-calling methods it wraps (BWA-MEM, FreeBayes) — but should recognize that Snippy's own specific implementation choices have not been independently peer-reviewed in the way its two most common downstream partners' have.

## Conclusion

This chapter closes two threads left open by Chapter 7: the pangenome-level comparison provided by Roary and Panaroo naturally leads to finer SNP-level outbreak phylogenetics (Snippy, Gubbins, IQ-TREE) when the question shifts from "what genes differ" to "how closely related are these specific isolates," and Unicycler's own author's redirection toward long-read-first assembly is followed through directly via Trycycler and Polypolish. Both threads reveal the same two small, prolific author groups behind an unusually large fraction of the modern bacterial genomics toolkit, both practicing the same pattern of explicit, public tool succession documented repeatedly across this eight-chapter series — from Prokka to Bakta, from ANCOM to ANCOM-BC, from Unicycler to Trycycler to Autocycler. Combined with the seven chapters preceding it, this chapter extends the series' central finding into its final domain: whether the subject is a microbial community or a single outbreak isolate, the modern computational toolkit for studying it is best understood not as a static catalog of individual programs, but as a continuously evolving lineage of tools in explicit, often self-documented conversation with one another.

## References

1. Croucher NJ, Page AJ, Connor TR, Delaney AJ, Keane JA, Bentley SD, Parkhill J, Harris SR. Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Res.* 2015;43(3):e15. PMID: 25414349. doi:10.1093/nar/gku1196
2. Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol.* 2020;37(5):1530–1534. PMID: 32011700. doi:10.1093/molbev/msaa015
3. Wick RR, Judd LM, Cerdeira LT, Hawkey J, Méric G, Vezina B, Wyres KL, Holt KE. Trycycler: consensus long-read assemblies for bacterial genomes. *Genome Biol.* 2021;22(1):266. PMID: 34521459. doi:10.1186/s13059-021-02483-z
4. Wick RR, Holt KE. Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLoS Comput Biol.* 2022;18(1):e1009802. PMID: 35073327. doi:10.1371/journal.pcbi.1009802

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| Snippy | Perl | https://github.com/tseemann/snippy |
| Gubbins | Python/C | https://github.com/nickjcroucher/gubbins |
| IQ-TREE | C++ | https://github.com/iqtree/iqtree2 |
| Trycycler | Python | https://github.com/rrwick/Trycycler |
| Polypolish | Rust | https://github.com/rrwick/Polypolish |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above. Snippy has no dedicated peer-reviewed publication and is cited as software directly, consistent with this series' practice of stating this explicitly rather than fabricating a citation.*
