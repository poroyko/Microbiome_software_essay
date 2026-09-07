# Python and R Packages for Microbiome Analysis — Chapter 6: Beyond Species-Level MAGs (Viral Metagenomics and Strain Tracking)

## Introduction

Every tool in Chapters 4 and 5 operated at the level of the bacterial or archaeal species-representative genome: bin contigs into a MAG, assess its quality, classify it taxonomically, annotate its functional potential. This chapter examines two directions of departure from that species-level unit, each pushing resolution in a different way. **Viral/phage metagenomics** goes sideways: bacterial and archaeal genome-recovery pipelines are not built to detect or classify viral sequences at all, requiring an entirely separate identification, quality-assessment, and taxonomy toolchain, covered here through **VirSorter2**, **CheckV**, **geNomad**, and **vConTACT2** — the third of these added specifically because CheckV's own documentation names it, by URL, as the tool to use for the one task (virus prediction) CheckV itself does not perform, a self-identified gap in this chapter's original scope rather than an externally imposed one. **Strain-level tracking** goes deeper: even a single well-classified bacterial species in a MAG collection is not one genetically uniform population, and two tools, **inStrain** and **StrainPhlAn**, recover the sub-species genetic variation and individual strain identities that species-level binning and taxonomy necessarily collapse away. A third theme emerges directly from the research for this chapter rather than being planned in advance: an unusually large share of these tools' canonical source code lives on Bitbucket rather than GitHub, concentrated specifically within two closely linked research groups (the Sullivan Lab/MAVERICLab for viral ecology tools, and the Banfield Lab for genome-resolved and strain-level metagenomics) — worth noting explicitly here since it directly affects where a reader should actually look for the software.

## Viral/Phage Metagenomics

### VirSorter2

VirSorter2 identifies viral sequences (DNA and RNA, including several viral groups beyond the well-studied Caudovirales bacteriophages) within assembled metagenomic contigs, using a multi-classifier, machine-learning approach that estimates "viralness" from genomic features including structural, functional, and taxonomic annotation alongside known viral hallmark genes, rather than relying on similarity to a single reference database alone.

**Prerequisites.** VirSorter2 is a Python package built as a wrapper around Snakemake for reproducible, cluster-friendly pipeline execution, installed via Conda (mamba is explicitly recommended by its own documentation) and requiring a downloadable reference database.

**Popularity.** VirSorter2 is a standard first step in modern viral metagenomics workflows and is explicitly designed to hand off directly to CheckV and DRAM-v (Chapter 5) in a single published standard operating procedure from its own developing lab, making it less a standalone tool than the entry point to a documented three-tool pipeline.

**Efficiency.** The original publication reports that, when benchmarked against genomes from both isolated and uncultivated viruses, VirSorter2 uniquely maintained high accuracy (F1-score > 0.8) consistently across major viral groups, while competing tools under-detected viruses outside Caudovirales specifically because those tools' reference databases were themselves biased toward that better-studied group — a finding about database representativeness as much as about algorithmic performance.

- **Citation:** Guo J, Bolduc B, Zayed AA, Varsani A, Dominguez-Huerta G, Delmont TO, Pratama AA, Gazitúa MC, Vik D, Sullivan MB, Roux S. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome.* 2021;9(1):37. PMID: [33522966](https://pubmed.ncbi.nlm.nih.gov/33522966/)
- **Repository:** [github.com/jiarong/VirSorter2](https://github.com/jiarong/VirSorter2) (also mirrored at [bitbucket.org/MAVERICLab/virsorter2](https://bitbucket.org/MAVERICLab/virsorter2))

### CheckV

CheckV assesses the quality and completeness of viral genomes assembled from metagenomes — the direct viral counterpart to CheckM/CheckM2 (Chapter 4) — by comparing query sequences against a large database of complete viral reference genomes to estimate completeness, identifying closed genomes via terminal repeat signatures, and specifically trimming flanking bacterial host DNA from integrated proviruses before that contamination can distort downstream analysis.

**Prerequisites.** CheckV is a Python package (requiring Python 3.10 or later) installable via Conda, pip, or Docker, depending on DIAMOND, HMMER, and a specialized viral-aware version of Prodigal (Prodigal-gv) for its underlying gene predictions.

**Popularity.** CheckV has become the standard quality-control step immediately following viral identification tools like VirSorter2, and its own FAQ candidly redirects users elsewhere for a task it does not perform: for actual virus *prediction* (as opposed to quality assessment of sequences already believed to be viral), CheckV's documentation explicitly recommends geNomad, a related tool from an overlapping set of developers — a direct, documented pointer to a complementary rather than competing tool.

**Efficiency.** The original publication reports that, applying CheckV to large, diverse viral sequence collections including IMG/VR and the Global Ocean Virome, 44,652 high-quality (>90% complete) viral genomes were identified, while also revealing that the vast majority of assembled viral sequences remain small fragments — a finding that speaks less to CheckV's own runtime performance than to the genuine difficulty of assembling complete viral genomes from short-read metagenomic data in the first place, a limitation CheckV exists to measure and report accurately rather than to fix.

- **Citation:** Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S, Kyrpides NC. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nat Biotechnol.* 2021;39(5):578–585. PMID: [33349699](https://pubmed.ncbi.nlm.nih.gov/33349699/)
- **Repository:** No canonical GitHub repository; source hosted at [bitbucket.org/berkeleylab/checkv](https://bitbucket.org/berkeleylab/checkv)

### geNomad

geNomad identifies and classifies both viral and plasmid sequences in assembled genomes or metagenomes — the specific virus-*prediction* task that CheckV's own documentation, cited above, explicitly points users toward this tool for — combining gene-content marker profiles with a deep neural network, and additionally detecting proviruses integrated into host chromosomes using a conditional random field model.

**Prerequisites.** geNomad is a Python package installable via Conda, Mamba, or Pixi, depending on MMseqs2 for its protein-profile search step and a downloadable marker database of more than 200,000 profiles; it is also available as a web application through the Galaxy and NMDC EDGE platforms for users who prefer not to install it locally.

**Popularity.** geNomad is developed by an overlapping author group to CheckV (both trace to the DOE Joint Genome Institute), and its reference database directly underlies IMG/VR and IMG/PR, two major public repositories of viral and plasmid genomes respectively — meaning geNomad's classifications propagate into reference databases other tools throughout this series ultimately draw upon, not merely into individual users' own analyses.

**Efficiency.** The original publication reports high classification performance for both plasmids and viruses (Matthews correlation coefficients of 77.8% and 95.3% respectively), substantially outperforming prior tools in these benchmarks, and demonstrates this at a scale directly relevant to modern metagenomics: geNomad was applied to more than 2.7 trillion base pairs of sequencing data, leading to the discovery of millions of previously uncharacterized viruses and plasmids.

- **Citation:** Camargo AP, Roux S, Schulz F, Babinski M, Xu Y, Hu B, Chain PSG, Nayfach S, Kyrpides NC. Identification of mobile genetic elements with geNomad. *Nat Biotechnol.* 2024;42(8):1303–1312. PMID: [37735266](https://pubmed.ncbi.nlm.nih.gov/37735266/)
- **GitHub:** [github.com/apcamargo/genomad](https://github.com/apcamargo/genomad)

### vConTACT2

vConTACT2 assigns taxonomy to uncultivated viral genomes using a "guilt-by-contig-association" strategy: it builds a gene-sharing network across viral genomes (query sequences plus reference genomes), then applies network clustering and distance-based hierarchical clustering to group genomes into taxonomically meaningful clusters, addressing the absence of a universal, scalable taxonomic framework for viruses comparable to GTDB-Tk's role for bacteria and archaea (Chapter 4).

**Prerequisites.** vConTACT2 is a Python package with a substantial dependency chain (DIAMOND or BLAST for protein comparisons, MCL for clustering, ClusterONE for overlapping cluster detection) and hardware requirements its own documentation describes as potentially considerable, exceeding 48 GB of memory depending on dataset size and complexity.

**Popularity.** vConTACT2 has become a standard tool for viral taxonomic assignment at genus-level resolution specifically, distributed both as a downloadable package and as an app through iVirus, the same CyVerse-based viral ecology infrastructure that also hosts VirSorter2.

**Efficiency.** The original publication reports that vConTACT2 achieved near-identical (96%) replication of existing genus-level viral taxonomy assignments from the International Committee on Taxonomy of Viruses when benchmarked against NCBI RefSeq viral genomes, and successfully scaled to a real metagenomic dataset of 15,280 Global Ocean Virome genome fragments, providing taxonomic assignments for 31% of that data — a concrete demonstration that a method validated on curated reference genomes can carry over to noisier, real environmental sequences, though with a substantially lower assignment rate than on references, a gap worth noting rather than glossing over.

- **Citation:** Bin Jang H, Bolduc B, Zablocki O, Kuhn JH, Roux S, Adriaenssens EM, Brister JR, Kropinski AM, Krupovic M, Lavigne R, Turner D, Sullivan MB. Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. *Nat Biotechnol.* 2019;37(6):632–639. PMID: [31061483](https://pubmed.ncbi.nlm.nih.gov/31061483/)
- **Repository:** No canonical GitHub repository; source hosted at [bitbucket.org/MAVERICLab/vcontact2](https://bitbucket.org/MAVERICLab/vcontact2) (community-maintained GitHub mirrors, e.g. [github.com/Hocnonsense/vcontact2](https://github.com/Hocnonsense/vcontact2), also exist but are not the canonical source)

## Strain-Level Tracking

### inStrain

inStrain profiles intra-population genetic diversity ("microdiversity") directly from metagenomic paired-read alignments to a reference genome, quantifying nucleotide diversity and linkage disequilibrium and identifying single-nucleotide variants (SNVs, distinguishing synonymous from non-synonymous changes) — going one level of resolution finer than the species-representative MAGs that Chapter 4's dRep step deliberately collapses multiple similar genomes into.

**Prerequisites.** inStrain is a Python program distributed via GitHub, Conda, and Docker, developed in the Banfield Lab (the same lab and author, Matt Olm, behind dRep in Chapter 4), and depends on a standard genome-mapping and variant-calling toolchain.

**Popularity.** inStrain has become a standard tool for microdiversity-aware genomic comparison specifically because it improves accuracy over methods that compare genomes as if each were a single, homogeneous consensus sequence, and its application to a large infant-gut cohort in the original publication has become a frequently cited demonstration of what strain-resolved analysis can reveal that species-level analysis cannot.

**Efficiency.** Rather than a runtime benchmark, the original publication's central efficiency claim is methodological: inStrain performs strain-level genomic comparisons with higher accuracy and sensitivity than leading existing tools specifically because it explicitly accounts for microdiversity rather than ignoring it, applied in the paper to more than 1,000 fecal metagenomes from premature infants, where it found that siblings shared significantly more microbial strains than unrelated infants — but, notably, that identical twins shared no more strains than fraternal siblings, a genuinely surprising finding that a coarser, species-only analysis could not have distinguished at all.

- **Citation:** Olm MR, Crits-Christoph A, Bouma-Gregson K, Firek BA, Morowitz MJ, Banfield JF. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. *Nat Biotechnol.* 2021;39(6):727–736. PMID: [33462508](https://pubmed.ncbi.nlm.nih.gov/33462508/)
- **GitHub:** [github.com/MrOlm/inStrain](https://github.com/MrOlm/inStrain)

### StrainPhlAn

StrainPhlAn tracks individual microbial strains across large sets of metagenomic samples by reconstructing, for each species of interest, a strain-specific sequence directly from the reads mapped to that species' marker genes in the MetaPhlAn database (Chapter 1), then building a phylogenetic tree from a multiple-sequence alignment across all samples' reconstructed strains — a marker-gene-based route to strain resolution, contrasting with inStrain's whole-genome-alignment-based approach above.

**Prerequisites.** StrainPhlAn is bundled within, and requires, the MetaPhlAn Python package, and additionally depends on a phylogenetic tree-building tool (RAxML by default, though others can be substituted) to construct the final strain phylogeny from the reconstructed marker-gene alignments.

**Popularity.** StrainPhlAn is one of two named strain-tracking tools that grew directly out of the same bioBakery ecosystem responsible for HUMAnN and MetaPhlAn (Chapter 1), and has been applied at very large scale — the tool's own documentation describes a strain-sharing-inference tutorial built on 7,646 stool metagenomic samples, illustrating person-to-person strain transmission inference as a distinct, specific downstream use case.

**Efficiency.** The original publication demonstrates StrainPhlAn's core value is resolution rather than speed: because it works directly from a species' existing marker genes rather than requiring whole-genome assembly of each strain, it can recover strain-level population structure from datasets where a full genome-resolved (MAG-based) approach would fail to assemble a usable genome at all, trading some strain-discrimination completeness (limited to variation within marker-gene regions) for applicability to a far broader range of metagenomic datasets, including many too shallow or too complex for reliable MAG recovery.

- **Citation:** Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N. Microbial strain-level population structure and genetic diversity from metagenomes. *Genome Res.* 2017;27(4):626–638. PMID: [28167665](https://pubmed.ncbi.nlm.nih.gov/28167665/)
- **GitHub:** [github.com/biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn) (StrainPhlAn is bundled within the MetaPhlAn repository, not distributed separately)

## Comparative Assessment

**A striking, previously-unremarked pattern culminates in this chapter: two of the four viral tools have no canonical GitHub repository at all.** CheckV and vConTACT2 both live on Bitbucket, joining MetaBAT2 (Chapter 4) as tools whose authoritative source code sits outside the platform most readers instinctively check first; geNomad and VirSorter2, by contrast, are both properly hosted on GitHub, so this is a real split within the viral toolchain rather than a blanket rule. What is new here is the concentration: both Bitbucket-hosted viral tools trace back to the same MAVERICLab/Sullivan Lab ecosystem that also produced VirSorter2 (itself dual-hosted on both platforms), while MetaBAT2's Bitbucket hosting traces to a separate DOE/Berkeley Lab tradition — the same DOE Joint Genome Institute lineage, notably, that also produced geNomad and CheckV, meaning a single institutional hosting habit doesn't even hold consistently across that lab's own tools. This is not one coincidental holdout but two distinct labs' institutional software-hosting habits, each affecting multiple tools discussed across this series — a genuinely practical finding for anyone assembling a viral metagenomics pipeline expecting every dependency to be one `git clone` away from a GitHub search.

**Viral and strain-level tools both formalize a documented handoff to a sibling tool, rather than positioning themselves as complete solutions.** VirSorter2 to CheckV to DRAM-v is a published, named standard operating procedure from a single lab, not an emergent community convention; CheckV's own FAQ names geNomad by URL for the one task (virus prediction) it does not perform — a pointer this chapter now follows all the way through, rather than leaving geNomad as a name-dropped reference the reader has to look up elsewhere. This is a tighter, more explicit form of the tool-composition pattern seen throughout this series (DAS Tool consuming any binner's output in Chapter 4, DRAM consuming GTDB-Tk/CheckM output in Chapter 5) — here, the authors themselves publish the recommended pipeline order rather than leaving it to downstream users or review papers to establish by convention.

**inStrain and StrainPhlAn solve the same underlying problem (sub-species resolution) via genuinely different data requirements, echoing the Deblur/DADA2 and CheckV/CheckM pairings elsewhere in this series.** inStrain requires read alignment to an assembled reference genome (typically a MAG), making it dependent on successful genome-resolved metagenomics (Chapter 4) having already succeeded for the species of interest; StrainPhlAn requires only reads mapping to MetaPhlAn's existing marker genes, letting it recover strain signal in datasets too shallow or complex for reliable MAG assembly at all. Neither tool is a strict upgrade over the other — the choice between them is a direct function of which upstream data (a good MAG, or merely adequate marker-gene coverage) is actually available for a given species in a given study.

## Conclusion

This chapter's two themes — viral sequences that species-level bacterial and archaeal pipelines cannot see at all, and strain-level variation that those same pipelines deliberately average away — both represent genuine resolution boundaries of the Chapter 4–5 genome-resolved workflow, not incremental refinements of it. VirSorter2, CheckV, geNomad, and vConTACT2 together form a documented, lab-endorsed pipeline for viral discovery, quality assessment, and taxonomy that parallels Chapter 4's bacterial/archaeal pipeline step for step, while inStrain and StrainPhlAn offer two structurally different routes to the sub-species resolution that a dereplicated MAG collection cannot provide by construction. Combined with the five chapters preceding it, this chapter completes a genuinely comprehensive tour of modern microbiome bioinformatics: from raw reads, through community composition and statistics, through genome assembly and functional annotation, to the viral dark matter and strain-level population structure that lie at the edges of what a standard species-level MAG-based analysis can see. The Bitbucket-hosting pattern uncovered here is a fitting closing detail for a six-chapter series that has consistently found more value in how these tools relate to, depend on, and sometimes quietly redirect toward one another than in treating any single tool as a complete, standalone answer.

## References

1. Guo J, Bolduc B, Zayed AA, Varsani A, Dominguez-Huerta G, Delmont TO, Pratama AA, Gazitúa MC, Vik D, Sullivan MB, Roux S. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome.* 2021;9(1):37. PMID: 33522966. doi:10.1186/s40168-020-00990-y
2. Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S, Kyrpides NC. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nat Biotechnol.* 2021;39(5):578–585. PMID: 33349699. doi:10.1038/s41587-020-00774-7
3. Camargo AP, Roux S, Schulz F, Babinski M, Xu Y, Hu B, Chain PSG, Nayfach S, Kyrpides NC. Identification of mobile genetic elements with geNomad. *Nat Biotechnol.* 2024;42(8):1303–1312. PMID: 37735266. doi:10.1038/s41587-023-01953-y
4. Bin Jang H, Bolduc B, Zablocki O, Kuhn JH, Roux S, Adriaenssens EM, Brister JR, Kropinski AM, Krupovic M, Lavigne R, Turner D, Sullivan MB. Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. *Nat Biotechnol.* 2019;37(6):632–639. PMID: 31061483. doi:10.1038/s41587-019-0100-8
5. Olm MR, Crits-Christoph A, Bouma-Gregson K, Firek BA, Morowitz MJ, Banfield JF. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. *Nat Biotechnol.* 2021;39(6):727–736. PMID: 33462508. doi:10.1038/s41587-020-00797-0
6. Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N. Microbial strain-level population structure and genetic diversity from metagenomes. *Genome Res.* 2017;27(4):626–638. PMID: 28167665. doi:10.1101/gr.216242.116

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| VirSorter2 | Python | https://github.com/jiarong/VirSorter2 |
| CheckV | Python | No canonical GitHub repo; source at [bitbucket.org/berkeleylab/checkv](https://bitbucket.org/berkeleylab/checkv) |
| geNomad | Python | https://github.com/apcamargo/genomad |
| vConTACT2 | Python | No canonical GitHub repo; source at [bitbucket.org/MAVERICLab/vcontact2](https://bitbucket.org/MAVERICLab/vcontact2) |
| inStrain | Python | https://github.com/MrOlm/inStrain |
| StrainPhlAn | Python | https://github.com/biobakery/MetaPhlAn |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above.*
