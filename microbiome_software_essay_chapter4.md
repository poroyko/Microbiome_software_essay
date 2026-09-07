# Python and R Packages for Microbiome Analysis — Chapter 4: Genome-Resolved Metagenomics (MAG Assembly, Binning, and Quality)

## Introduction

Every tool in Chapters 1 through 3 operated on a common currency: a feature table of taxa or gene counts, derived from marker-gene (16S rRNA) sequencing or from mapping shotgun reads against existing reference databases. Genome-resolved metagenomics asks a fundamentally different question — rather than measuring *how much* of a known organism is present, it attempts to reconstruct the actual draft genomes of the organisms living in a sample directly from shotgun sequencing reads, including organisms with no prior reference genome at all. This produces metagenome-assembled genomes (MAGs), and requires an entirely different software pipeline: assembling reads into contigs, grouping ("binning") those contigs into putative genomes, refining and evaluating the quality of those bins, removing redundant genomes recovered from related samples, and finally assigning each surviving genome a taxonomic identity. This chapter covers ten tools spanning that complete pipeline: assembly (**MEGAHIT**, **metaSPAdes**), binning (**MetaBAT2**, **MaxBin2**, **CONCOCT**), bin refinement (**DAS Tool**), quality assessment (**CheckM**, **CheckM2**), dereplication (**dRep**), and taxonomic classification (**GTDB-Tk**) — of which MetaBAT2, CheckM2, and GTDB-Tk were specifically requested, with the remaining seven added to complete the workflow they sit inside, including every tool this chapter's own text names as a direct comparator or dependency (MaxBin2 and CONCOCT were named as DAS Tool's example inputs; CheckM was named as dRep's own dependency). Unlike the R-heavy statistical ecosystem of Chapters 1–3, this is a field dominated by C++ and Python, for reasons this chapter's Comparative Assessment addresses directly.

## Assembly: MEGAHIT and metaSPAdes

MEGAHIT assembles short sequencing reads into longer contiguous sequences (contigs) using a succinct de Bruijn graph representation, a data structure chosen specifically because it uses substantially less memory than conventional de Bruijn graph implementations — the critical bottleneck for assembling the very large, complex datasets typical of environmental metagenomes.

**Prerequisites.** MEGAHIT is implemented in C++ (with an optional GPU-accelerated code path for its graph-construction step) and is distributed as source code or precompiled binaries via GitHub and Conda; it requires no companion database, unlike several downstream tools in this chapter.

**Popularity.** MEGAHIT is one of the two most widely used metagenome assemblers in current practice (alongside metaSPAdes), and is frequently the first tool invoked in a genome-resolved metagenomics pipeline specifically because of its ability to handle very large datasets on a single compute node without requiring a computing cluster.

**Efficiency.** The original publication reports assembling a 252-billion-base-pair soil metagenome — at the time, an unusually large and complex dataset — in 44.1 hours with a GPU or 99.6 hours without one, using at most 260–345 GB of memory on a single server; the authors specifically contrast this against competing assemblers of the era, which could not complete this assembly at all within comparable resource budgets, establishing memory efficiency, not raw speed alone, as MEGAHIT's central contribution.

- **Citation:** Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics.* 2015;31(10):1674–1676. PMID: [25609793](https://pubmed.ncbi.nlm.nih.gov/25609793/)
- **GitHub:** [github.com/voutcn/megahit](https://github.com/voutcn/megahit)

metaSPAdes extends the SPAdes genome assembler (Chapter 7's standard tool for single-isolate assembly) to metagenomic data, applying computational strategies originally developed for single-cell and diploid-genome assembly — settings that, like metagenomes, involve highly non-uniform coverage — to the more complex problem of assembling reads from an unknown number of related but distinct organisms at once.

**Prerequisites.** metaSPAdes ships as a mode within the same SPAdes toolkit and codebase (C++/Python) described in Chapter 7, requiring no separate installation; selecting metagenomic assembly is a matter of passing the `--meta` flag rather than installing different software.

**Popularity.** metaSPAdes and MEGAHIT are, per independent benchmarking literature, considered the two leading metagenome assemblers in current use, with metaSPAdes generally favored when assembly contiguity is the priority and MEGAHIT favored when memory efficiency and very large sample sizes are the priority.

**Efficiency.** The original publication benchmarks metaSPAdes against MEGAHIT and other contemporaneous metagenome assemblers (IDBA-UD, Ray-Meta) across diverse real datasets, reporting that metaSPAdes produces high-quality assemblies competitive with or exceeding these alternatives; independent practice has since established that this quality advantage typically comes at the cost of substantially higher memory and runtime requirements than MEGAHIT, making the choice between the two a genuine trade-off rather than one tool strictly dominating the other.

- **Citation:** Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. *Genome Res.* 2017;27(5):824–834. PMID: [28298430](https://pubmed.ncbi.nlm.nih.gov/28298430/)
- **GitHub:** [github.com/ablab/spades](https://github.com/ablab/spades) (same repository as SPAdes; invoked via the `--meta` flag)

## Binning: MetaBAT2, MaxBin2, and CONCOCT

MetaBAT2 groups assembled contigs into putative genome bins using two complementary signals: tetranucleotide frequency (a compositional signature that tends to be consistent within, but differ between, genomes) and per-sample contig coverage depth (contigs from the same genome tend to be present at similar abundance across samples) — the same two-feature strategy underlying most modern binning tools, refined here through an adaptive algorithm that removes the need for manual parameter tuning present in the original MetaBAT.

**Prerequisites.** MetaBAT2 is implemented in C++; notably, its canonical source repository is hosted on Bitbucket rather than GitHub (`bitbucket.org/berkeleylab/metabat`), a detail worth stating plainly for anyone searching GitHub specifically and coming up empty, echoing the ANCOM situation from Chapter 2 where the "obvious" hosting platform was not where the authoritative source actually lived.

**Popularity.** MetaBAT2's publication has accumulated over 2,500 citations by independent tracking, and its predecessor, MetaBAT, is described in the tool's own follow-up publication as having become one of the most popular binning tools available, largely credited to its computational efficiency and ease of use on large, multi-sample binning experiments.

**Efficiency.** The original publication reports that MetaBAT2, benchmarked against alternative binning tools (including MaxBin2, CONCOCT, and its own predecessor) across more than 100 real-world metagenome assemblies, achieved superior accuracy and computing speed, consistently recovering more bins meeting 95%, 70%, and 50% completeness thresholds than competing tools, while completing a typical metagenome assembly's binning in only a few minutes on a single commodity workstation.

- **Citation:** Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ.* 2019;7:e7359. PMID: [31388474](https://pubmed.ncbi.nlm.nih.gov/31388474/)
- **Repository:** No canonical GitHub repository; source hosted at [bitbucket.org/berkeleylab/metabat](https://bitbucket.org/berkeleylab/metabat)

MaxBin2 bins contigs using a similar feature set to MetaBAT2 (tetranucleotide frequency and coverage) but via a distinct algorithmic strategy: an expectation-maximization procedure that also incorporates single-copy marker genes directly into the binning process itself, both to guide clustering and to estimate the number of genomes present, rather than requiring that number as a separate input.

**Prerequisites.** MaxBin2 is distributed via SourceForge rather than GitHub — a second binning tool in this chapter, alongside MetaBAT2, whose canonical source lives outside GitHub — and is installable via Conda, requiring Perl and a standard bioinformatics toolchain (Bowtie2, FragGeneScan) for its coverage-estimation and marker-gene-detection steps.

**Popularity.** MaxBin2 is, alongside MetaBAT2 and CONCOCT, one of the three most established and commonly used k-mer-composition-and-abundance binners per independent methods surveys, and is one of the two binners most frequently paired with MetaBAT2 as input to consensus tools like DAS Tool (below).

**Efficiency.** The original publication reports that MaxBin 2.0 is highly accurate in recovering individual genomes from both simulated datasets and real environmental metagenomes, and specifically demonstrates that binning multiple related samples together (co-assembly binning) recovers more genomes than binning any single sample alone — a direct multi-sample efficiency gain the expectation-maximization approach was extended specifically to capture.

- **Citation:** Wu YW, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics.* 2016;32(4):605–607. PMID: [26515820](https://pubmed.ncbi.nlm.nih.gov/26515820/)
- **Repository:** No canonical GitHub repository; source hosted at [sourceforge.net/projects/maxbin2](https://sourceforge.net/projects/maxbin2)

CONCOCT bins contigs using unsupervised Gaussian mixture modeling over the same two feature types (composition and multi-sample coverage), with a particular emphasis on scaling gracefully to datasets spanning many samples at once — the paper's own validation case study reconstructs genomes from a real *E. coli* O104:H4 outbreak investigation, tying its binning method directly to a real epidemiological application.

**Prerequisites.** CONCOCT is a Python package distributed via GitHub and Conda, depending on standard scientific Python libraries (scikit-learn, for its Gaussian mixture model implementation) alongside a read aligner for generating the multi-sample coverage tables it requires as input.

**Popularity.** CONCOCT is the third member of the "big three" composition-and-coverage binners discussed throughout this chapter, and — unlike MetaBAT2 and MaxBin2 — has a canonical GitHub repository, making it the only one of the three binners compared here for which "search GitHub" is a reliable way to find the authoritative source.

**Efficiency.** The original publication demonstrates high recall and precision on both artificial datasets and real human gut metagenome data, and its outbreak case study specifically shows that CONCOCT can resolve closely related genomic variants (distinguishing the outbreak strain from background *E. coli* diversity) purely from coverage and composition patterns across samples, without prior reference genomes for the outbreak strain itself.

- **Citation:** Alneberg J, Bjarnason BS, de Bruijn I, Schirmer M, Quick J, Ijaz UZ, Lahti L, Loman NJ, Andersson AF, Quince C. Binning metagenomic contigs by coverage and composition. *Nat Methods.* 2014;11(11):1144–1146. PMID: [25218180](https://pubmed.ncbi.nlm.nih.gov/25218180/)
- **GitHub:** [github.com/BinPro/CONCOCT](https://github.com/BinPro/CONCOCT)

## Bin Refinement: DAS Tool

Different binning algorithms make different mistakes on the same data, and no single binner consistently outperforms the others across every sample type. DAS Tool addresses this directly by taking the output of multiple independent binning runs (e.g., MetaBAT2, MaxBin2, and CONCOCT all applied to the same assembly) and using a dereplication, aggregation, and scoring strategy — based on the presence of single-copy marker genes — to select the best non-redundant set of bins across all of them, rather than trusting any one binner's output alone.

**Prerequisites.** DAS Tool's core scoring algorithm is implemented in R, invoked via a command-line wrapper, and depends on a chosen external search engine (DIAMOND by default, with BLASTP or USEARCH as alternatives) to identify the single-copy marker genes its scoring relies on.

**Popularity.** DAS Tool is one of two commonly used bin-refinement strategies, the other being the refinement module bundled with the metaWRAP pipeline; both are frequently used as points of comparison for newer refinement tools, including a 2022 successor, MAGScoT, that specifically benchmarks itself against DAS Tool's iterative scoring approach.

**Efficiency.** The original publication reports that DAS Tool applied to a constructed (ground-truth-known) community generated more accurate bins than any single automated binning method, and when applied to real environmental and host-associated samples of varying complexity, recovered substantially more near-complete genomes — including previously unreported lineages — than any individual binner alone, directly demonstrating the value of combining, rather than choosing between, competing binning algorithms.

- **Citation:** Sieber CMK, Probst AJ, Sharrar A, Thomas BC, Hess M, Tringe SG, Banfield JF. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nat Microbiol.* 2018;3(7):836–843. PMID: [29807988](https://pubmed.ncbi.nlm.nih.gov/29807988/)
- **GitHub:** [github.com/cmks/DAS_Tool](https://github.com/cmks/DAS_Tool)

## Quality Assessment: CheckM and CheckM2

Before any recovered bin can be trusted as a genuine draft genome, its completeness (what fraction of the expected genome was recovered) and contamination (how much of the bin's content actually belongs to other organisms) must be estimated. CheckM addresses this using collocated sets of single-copy marker genes specific to a bin's estimated phylogenetic lineage, placing each bin within a reference genome tree to select an appropriate, lineage-specific marker set rather than relying on one universal set applied indiscriminately across all bacteria and archaea.

**Prerequisites.** CheckM is a Python package (requiring Python 3, as of v1.1.x) distributed via GitHub and Conda, depending on Prodigal (Chapter 5) for gene prediction and pplacer for its phylogenetic placement step — a direct, documented dependency on the same foundational gene-caller that underlies Prokka, Bakta, and CheckV elsewhere in this series.

**Popularity.** CheckM has been directly depended upon by other tools discussed in this chapter — dRep uses it by default for genome quality filtering during its scoring step — and its benchmarking on 2,281 isolate genomes from major reference-genome projects (GEBA, HMP among them) made it a standard quality-assessment step across genome-resolved metagenomics for the better part of a decade before CheckM2's release.

**Efficiency.** The original publication reports that, using lineage-specific marker sets, 96% of a 2,281-genome real-world benchmark set were correctly estimated at ≥95% complete with ≤5% contamination; CheckM's own GitHub repository states plainly that CheckM v1 is no longer under active feature development (only compatibility fixes), and directly recommends CheckM2 as the generally more accurate successor — while also noting that running both in parallel remains informative, since the two use genuinely independent methodologies (lineage-specific marker sets versus universally trained machine learning models) that can each catch errors the other misses.

- **Citation:** Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Res.* 2015;25(7):1043–1055. PMID: [25977477](https://pubmed.ncbi.nlm.nih.gov/25977477/)
- **GitHub:** [github.com/Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM)

CheckM2 replaces the lineage-specific marker-gene approach of the original CheckM with universally trained machine learning models, applied regardless of taxonomic lineage — a change specifically intended to improve accuracy on novel or poorly represented lineages that lack enough reference genomes to build a reliable lineage-specific marker set.

**Prerequisites.** CheckM2 is a Python package distributed via GitHub and Conda, and depends on a DIAMOND-searchable reference database that must be downloaded separately before first use.

**Popularity.** CheckM2 has rapidly become a standard component of genome-resolved metagenomics pipelines since its 2023 publication, and is frequently invoked as the quality-assessment step inside newer bin-refinement tools (such as Binette, mentioned in community binning tutorials) that use CheckM2's completeness and contamination scores to choose among candidate hybrid bins.

**Efficiency.** The original publication demonstrates, using both synthetic and experimental data, that CheckM2 outperforms existing genome-quality tools in both accuracy and computational speed, and specifically shows that CheckM2 remains accurate on organisms with reduced genome size and unusual biology — such as Patescibacteria and the DPANN archaeal superphylum — where the original CheckM's lineage-specific marker approach is more prone to error due to the sparse reference genome sets available for these groups.

- **Citation:** Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods.* 2023;20(8):1203–1212. PMID: [37500759](https://pubmed.ncbi.nlm.nih.gov/37500759/)
- **GitHub:** [github.com/chklovski/CheckM2](https://github.com/chklovski/CheckM2)

## Dereplication: dRep

Genome-resolved metagenomic studies spanning multiple samples routinely recover the same organism's genome, in slightly different draft-quality versions, from several samples independently. dRep identifies these redundant genome sets using a two-step comparison strategy — a fast, approximate genome-distance estimate followed by a slower, precise average nucleotide identity (ANI) calculation only within candidate clusters — and selects the single best representative genome from each set, producing a non-redundant collection for downstream analysis.

**Prerequisites.** dRep is a Python package installable via pip or Conda, and depends on CheckM for genome quality filtering during its scoring step (an optional dependency that can be bypassed if genome quality metrics, e.g. from CheckM2, are supplied directly).

**Popularity.** dRep's documentation describes it as part of a now-standard genome-resolved metagenomic workflow — bin, dereplicate to species-representative genomes at a 95% ANI threshold, then map all samples back against the dereplicated set — a pattern that has become common enough that dRep's own maintainer references a companion perspective piece specifically titled "To Dereplicate or Not To Dereplicate?" addressing when this step is and is not appropriate.

**Efficiency.** The original publication reports that dRep achieves a 28-fold increase in speed with perfect recall and precision when benchmarked against previously developed comparison algorithms, and demonstrates that using dRep for genome recovery from time-series data recovers significantly more and higher-quality genomes than assembling all time points together (co-assembly); more recently, dRep version 4 replaced its default distance-estimation tools with a newer method (skani), reported to be substantially faster still, with memory usage that scales roughly linearly rather than quadratically with the number of genomes compared — a meaningful practical improvement as genome collections have grown into the tens of thousands.

- **Citation:** Olm MR, Brown CT, Brooks B, Banfield JF. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. *ISME J.* 2017;11(12):2864–2868. PMID: [28742071](https://pubmed.ncbi.nlm.nih.gov/28742071/)
- **GitHub:** [github.com/MrOlm/drep](https://github.com/MrOlm/drep)

## Taxonomic Classification: GTDB-Tk

Once a set of high-quality, dereplicated MAGs has been produced, GTDB-Tk assigns each one an objective taxonomic classification based on the Genome Taxonomy Database (GTDB) — a genome-phylogeny-based taxonomy that has substantially revised, and in many cases replaced, traditional NCBI taxonomy for bacterial and archaeal classification, precisely because it is built to be applied consistently to MAGs that may represent entirely novel, previously unnamed lineages.

**Prerequisites.** GTDB-Tk is implemented in Python and depends on a substantial reference package (the GTDB reference tree and associated marker sets) that must be downloaded separately; the original version's memory requirement — approximately 320 GB, driven by the size of the bacterial reference tree — was significant enough to warrant a dedicated follow-up release, GTDB-Tk v2, specifically to address it.

**Popularity.** GTDB and GTDB-Tk have been widely adopted by the microbiology community for taxonomic classification of both isolate genomes and MAGs recovered from environmental and human-associated samples, and are now the standard reference framework cited by large-scale genome catalog efforts, including a unified catalog of over 200,000 reference genomes from the human gut microbiome.

**Efficiency.** GTDB-Tk v2's central contribution is a divide-and-conquer strategy that first places a query genome into a bacterial reference tree containing only family-level representatives, then places it more precisely into an appropriate class-level subtree containing full species representatives — the original publication demonstrates that this two-stage approach produces classifications nearly equivalent to the original full-tree method while providing a substantial reduction in required memory, making the tool practical on hardware that could not run the original version at all.

- **Citations:** Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics.* 2020;36(6):1925–1927. PMID: [31730192](https://pubmed.ncbi.nlm.nih.gov/31730192/). Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics.* 2022;38(23):5315–5316. PMID: [36218463](https://pubmed.ncbi.nlm.nih.gov/36218463/)
- **GitHub:** [github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)

## Comparative Assessment

**This is the first chapter in this series where R plays a minor, rather than central, role — and the reason is structural, not incidental.** Of the ten tools covered here, six (MEGAHIT, metaSPAdes, MetaBAT2, MaxBin2, CheckM, CheckM2 — plus GTDB-Tk and dRep, both Python) are implemented in C++, Perl, or Python, and only one (DAS Tool) is R, and even DAS Tool's R script functions as a scoring engine invoked from a command-line wrapper rather than as an interactive analysis environment. This contrasts sharply with Chapters 1–3, where R's statistical and community-ecology heritage made it the natural home for downstream analysis of already-summarized feature tables. Genome-resolved metagenomics instead begins from raw sequencing reads and must perform graph construction, sequence comparison, and machine-learning classification at a scale (assemblies with hundreds of gigabases, reference trees with hundreds of thousands of nodes) where C++'s raw performance and Python's machine-learning ecosystem are better structural fits than R's vectorized-statistics design — a genuine example of task shaping language choice, rather than either language being categorically better.

**A third of this chapter's tools have no canonical presence on GitHub at all, which is a higher proportion than any other chapter in this series.** MetaBAT2 and CheckV-adjacent tools' Bitbucket hosting was already notable in earlier chapters, but this chapter adds a new pattern: MaxBin2's canonical source lives on SourceForge, a third hosting platform beyond GitHub and Bitbucket entirely. Of the three "big three" binners compared directly in this chapter, only CONCOCT has a straightforward canonical GitHub repository — a genuinely useful, concrete finding for anyone assembling a binning pipeline and expecting all three tools to be found the same way.

**Memory, not runtime, is the recurring efficiency bottleneck across this entire pipeline.** MEGAHIT's central innovation (succinct de Bruijn graphs) exists specifically to reduce assembly memory footprint, a design point sharpened by direct comparison with metaSPAdes, which achieves competitive or superior assembly quality at a real, acknowledged cost in memory and runtime; GTDB-Tk v2 exists specifically to reduce its predecessor's ~320 GB requirement; dRep v4's adoption of skani is framed explicitly around better memory scaling with genome count. This is a different pattern from the amplicon-sequencing tools of Chapters 1–3, where the more common efficiency story was runtime or statistical power (Kraken 2's 85% memory reduction from Chapter 1 is the one clear exception, and notably, Kraken 2 is also implemented in C++ and also targets shotgun sequencing data at genome-comparison scale) — suggesting memory pressure specifically accompanies genome- and read-scale operations, rather than being a universal concern across all microbiome bioinformatics.

**Each tool in this pipeline exists because the previous step is imperfect, not merely to add a feature.** DAS Tool exists because no single binner is reliably best across sample types, a claim now directly visible in this chapter's own text by comparing MetaBAT2, MaxBin2, and CONCOCT's genuinely distinct algorithmic strategies (adaptive distance-based clustering, expectation-maximization, Gaussian mixture modeling) rather than taking DAS Tool's premise on faith; CheckM2 exists because the original CheckM's lineage-specific markers fail on novel lineages, with CheckM's own GitHub repository now explicitly recommending its successor while still endorsing running both in parallel; GTDB-Tk v2 exists because GTDB-Tk v1's memory requirements limited its adoption. This chapter's pipeline is therefore best understood not as ten independent tools that happen to be used in sequence, but as a series of direct responses to specific, documented failure modes of the step before — a pattern consistent with, but more tightly sequential than, the tool-succession pattern (DEICODE→Gemelli, ANCOM→ANCOM-BC, SourceTracker→SourceTracker2/FEAST) observed throughout Chapters 1–3.

## Conclusion

Genome-resolved metagenomics answers a question the marker-gene and reference-mapping approaches of Chapters 1–3 cannot: what are the actual genomes of the organisms in a sample, including those with no existing reference at all? The ten-tool pipeline in this chapter — assemble (MEGAHIT, metaSPAdes), bin (MetaBAT2, MaxBin2, CONCOCT), refine (DAS Tool), assess quality (CheckM, CheckM2), dereplicate (dRep), and classify (GTDB-Tk) — represents a mature, widely adopted workflow for answering it, but one built almost entirely outside the R ecosystem that dominated the rest of this series, for reasons rooted in the raw scale of the underlying computation rather than any inherent limitation of R itself. Combined with the amplicon-sequencing, differential-abundance, diversity-estimation, source-tracking, and ordination-visualization tools surveyed in Chapters 1 through 3, this chapter completes a genuinely broad map of the modern microbiome bioinformatics landscape — and, consistent with every prior chapter's closing observation, the tools here too are best understood in relation to one another rather than as standalone choices.

## References

1. Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics.* 2015;31(10):1674–1676. PMID: 25609793. doi:10.1093/bioinformatics/btv033
2. Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. *Genome Res.* 2017;27(5):824–834. PMID: 28298430. doi:10.1101/gr.213959.116
3. Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ.* 2019;7:e7359. PMID: 31388474. doi:10.7717/peerj.7359
4. Wu YW, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics.* 2016;32(4):605–607. PMID: 26515820. doi:10.1093/bioinformatics/btv638
5. Alneberg J, Bjarnason BS, de Bruijn I, Schirmer M, Quick J, Ijaz UZ, Lahti L, Loman NJ, Andersson AF, Quince C. Binning metagenomic contigs by coverage and composition. *Nat Methods.* 2014;11(11):1144–1146. PMID: 25218180. doi:10.1038/nmeth.3103
6. Sieber CMK, Probst AJ, Sharrar A, Thomas BC, Hess M, Tringe SG, Banfield JF. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nat Microbiol.* 2018;3(7):836–843. PMID: 29807988. doi:10.1038/s41564-018-0171-1
7. Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Res.* 2015;25(7):1043–1055. PMID: 25977477. doi:10.1101/gr.186072.114
8. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods.* 2023;20(8):1203–1212. PMID: 37500759. doi:10.1038/s41592-023-01940-w
9. Olm MR, Brown CT, Brooks B, Banfield JF. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. *ISME J.* 2017;11(12):2864–2868. PMID: 28742071. doi:10.1038/ismej.2017.126
10. Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics.* 2020;36(6):1925–1927. PMID: 31730192. doi:10.1093/bioinformatics/btz848
11. Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics.* 2022;38(23):5315–5316. PMID: 36218463. doi:10.1093/bioinformatics/btac672

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| MEGAHIT | C++ | https://github.com/voutcn/megahit |
| metaSPAdes | C++/Python | https://github.com/ablab/spades |
| MetaBAT2 | C++ | No canonical GitHub repo; source at [bitbucket.org/berkeleylab/metabat](https://bitbucket.org/berkeleylab/metabat) |
| MaxBin2 | Perl | No canonical GitHub repo; source at [sourceforge.net/projects/maxbin2](https://sourceforge.net/projects/maxbin2) |
| CONCOCT | Python | https://github.com/BinPro/CONCOCT |
| DAS Tool | R (CLI wrapper) | https://github.com/cmks/DAS_Tool |
| CheckM | Python | https://github.com/Ecogenomics/CheckM |
| CheckM2 | Python | https://github.com/chklovski/CheckM2 |
| dRep | Python | https://github.com/MrOlm/drep |
| GTDB-Tk | Python | https://github.com/Ecogenomics/GTDBTk |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above.*
