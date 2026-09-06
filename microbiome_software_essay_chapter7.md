# Python and R Packages for Microbiome Analysis — Chapter 7: Pure Isolate WGS (Assembly, AMR, Serotyping, and Genome Comparison)

## Introduction

Every chapter so far has addressed mixed-community data: amplicon surveys, shotgun metagenomes, MAGs assembled from a soup of co-occurring organisms. This chapter turns to a related but methodologically distinct problem: whole-genome sequencing (WGS) of a single bacterial isolate grown in pure culture. Isolate WGS removes the binning and dereplication challenges central to Chapters 4 and 6 — there is exactly one genome to recover, not an unknown number mixed together — but introduces its own specialized toolchain, built around questions a metagenomic pipeline never has to ask: is this genome fully assembled into one contig per replicon, does it carry acquired antimicrobial resistance genes, what serotype is it, and how does its gene content compare to other isolates of the same species. This chapter covers seven tools across four requested categories: **assembly** (SPAdes, Unicycler), **antimicrobial resistance (AMR) detection** (AMRFinderPlus, CARD/RGI), **serotyping** (SeqSero2), and **genome comparison via pangenome analysis** (Roary, Panaroo). **Annotation** is deliberately not given new entries here: Prokka and Bakta, covered in Chapter 5, were in fact originally built for exactly this isolate-genomics use case (their MAG applicability came later), so this chapter cross-references those entries directly rather than duplicating them. As in every prior chapter, real patterns of tool succession, collaboration, and documented interdependency emerged directly from the research rather than being imposed in advance.

## Genome Assembly

### SPAdes

SPAdes assembles genomic sequencing reads into contigs using a de Bruijn graph approach, originally developed to handle the highly non-uniform coverage and elevated error rates of single-cell sequencing data, and subsequently extended into the default general-purpose short-read assembler for standard (multicell) bacterial isolates as well.

**Prerequisites.** SPAdes is a C++/Python toolkit distributed via GitHub and Conda, primarily designed for Illumina (and IonTorrent) short reads, with hybrid-assembly support for supplementing short reads with PacBio, Oxford Nanopore, or Sanger long reads.

**Popularity.** SPAdes is among the most heavily cited tools in this entire six-chapter series: independent tracking lists over 22,700 citations for the original 2012 publication, reflecting its status as the default starting point for bacterial isolate genome assembly across essentially the entire field, not merely a niche or specialized option.

**Efficiency.** The original publication demonstrated that SPAdes improved on the previously best-available single-cell assembler (E+V-SC) as well as popular multicell assemblers (Velvet, SoapDeNovo) of the time; its continued relevance more than a decade later is itself a notable efficiency data point, since few bioinformatics tools of any kind remain the field's default choice across such a long span without being displaced by a categorically different approach.

- **Citation:** Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, Nikolenko SI, Pham S, Prjibelski AD, Pyshkin AV, Sirotkin AV, Vyahhi N, Tesler G, Alekseyev MA, Pevzner PA. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *J Comput Biol.* 2012;19(5):455–477. PMID: [22506599](https://pubmed.ncbi.nlm.nih.gov/22506599/)
- **GitHub:** [github.com/ablab/spades](https://github.com/ablab/spades)

### Unicycler

Unicycler produces complete bacterial genome assemblies by combining short and long sequencing reads: it builds an initial, highly accurate assembly graph from short reads using SPAdes, then uses long reads and a custom semi-global aligner to resolve the graph's repeat regions and scaffold it toward a single contig per replicon — a "short-read-first" hybrid strategy specifically designed to make effective use of long reads even when their depth and accuracy are low.

**Prerequisites.** Unicycler is a Python/C++ package distributed via GitHub and Conda, depending directly on SPAdes (version 3.6.2 or later) as its short-read assembly engine, plus a long-read aligner for the scaffolding step.

**Popularity.** Unicycler became a standard tool for bacterial hybrid assembly following its 2017 publication, and — in a direct parallel to Prokka's author-endorsed succession by Bakta in Chapter 5 — Unicycler's own author now states in its GitHub README that, while the tool "is not completely out-of-date," he has since developed two newer tools, Trycycler and Polypolish, specifically because higher-depth, higher-accuracy modern long-read sequencing has made a different strategy ("long-read-first" assembly, i.e., long-read assembly followed by short-read polishing) often preferable to Unicycler's original short-read-first design.

**Efficiency.** The original publication reports that, on both synthetic and real sequencing read tests, Unicycler assembled larger contigs with fewer misassemblies than competing hybrid assemblers, even under low long-read depth and accuracy conditions — the specific scenario (sparse, noisy 2016-era Nanopore data) the tool was originally designed around, which its own author's more recent guidance suggests is now less universally the scenario researchers face.

- **Citation:** Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. *PLoS Comput Biol.* 2017;13(6):e1005595. PMID: [28594827](https://pubmed.ncbi.nlm.nih.gov/28594827/)
- **GitHub:** [github.com/rrwick/Unicycler](https://github.com/rrwick/Unicycler)

## Annotation (Cross-Reference)

Genome annotation for isolate WGS uses the same tools introduced in Chapter 5 — **Prokka** and **Bakta** — and is not repeated in full here. It is worth stating plainly, though, that isolate genomics was these tools' original and primary use case: Prokka's 2014 publication frames its motivation entirely around annotating "a draft bacterial genome" from an isolate, and Bakta's own README explicitly advertises annotation of "bacterial genomes, MAGs & plasmids" — MAGs listed second, isolate genomes first. Readers focused specifically on isolate WGS should treat Chapter 5's Prokka/Bakta discussion, including the direct, author-endorsed succession from one to the other, as this chapter's annotation section by reference.

## Antimicrobial Resistance Detection

### AMRFinderPlus

AMRFinderPlus identifies acquired antimicrobial resistance genes, resistance-conferring point mutations, and related stress-response and virulence genes in assembled bacterial genomes, drawing on NCBI's own curated Reference Gene Catalog — a database built as part of a multi-agency collaboration specifically to provide a comprehensive, quality-controlled alternative to the fragmented AMR gene databases that existed previously.

**Prerequisites.** AMRFinderPlus is distributed via GitHub and Conda (as the `ncbi/amr` package), accepts either assembled nucleotide sequences or predicted protein sequences as input, and applies taxon-specific filtering to exclude AMR genes that are so widespread within certain taxa (intrinsic genes) that reporting them would not be informative.

**Popularity.** AMRFinderPlus is used directly within NCBI's own Pathogen Detection Project, run against every bacterial isolate genome submitted to that system, giving it an unusually direct line to large-scale public health surveillance infrastructure rather than remaining solely an academic research tool.

**Efficiency.** The original publication reports that AMRFinderPlus was vetted against two well-characterized datasets and found to have high concordance with previously published results for both AMR and stress-response genes; notably, NCBI's own documentation describes an active, ongoing collaboration specifically with CARD (below) "to resolve issues and communicate updates and new genes" — a genuinely collaborative rather than competitive relationship between what might otherwise appear to be two rival AMR detection tools.

- **Citation:** Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep.* 2021;11(1):12728. PMID: [34135355](https://pubmed.ncbi.nlm.nih.gov/34135355/)
- **GitHub:** [github.com/ncbi/amr](https://github.com/ncbi/amr)

### CARD / RGI

CARD (the Comprehensive Antibiotic Resistance Database) provides a manually curated, ontology-driven knowledgebase of antimicrobial resistance determinants, paired with RGI (the Resistance Gene Identifier), the software that applies CARD's curated detection models to a query genome and reports matches under a tiered "Perfect / Strict / Loose" confidence paradigm reflecting how closely a match resembles a fully characterized reference sequence.

**Prerequisites.** RGI is a Python package distributed via GitHub and Conda, and can be run against genomic assemblies, proteomes, or raw metagenomic sequencing reads; it depends on the separately downloaded CARD reference database, which is updated on its own release cycle.

**Popularity.** CARD and RGI have been used to generate a large public resistome surveillance resource ("WildCARD"), in which RGI was applied to complete chromosome, plasmid, genomic island, and whole-genome-shotgun sequences across hundreds of pathogens available in NCBI Genomes, producing prevalence statistics for AMR gene distribution across pathogens and plasmids.

**Efficiency.** CARD's own 2023 publication notes explicitly that curation of the underlying database routinely runs ahead of RGI's software implementation, such that not every model type or parameter curated in CARD is yet supported by the RGI software itself — an unusually candid, publication-level acknowledgment of a real gap between a database's curated content and its companion tool's current capabilities, rather than a claim of complete feature parity.

- **Citation:** Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* 2023;51(D1):D690–D699. PMID: [36263822](https://pubmed.ncbi.nlm.nih.gov/36263822/)
- **GitHub:** [github.com/arpcard/rgi](https://github.com/arpcard/rgi)

## Serotyping

### SeqSero2

SeqSero2 predicts *Salmonella* serotype directly from whole-genome sequencing data — either raw reads or assembled genomes — by identifying serotype-determining O- and H-antigen sequence markers, reconstructing them via one of three interchangeable workflows (raw-read k-mer matching, targeted allele micro-assembly, or full genome assembly), and mapping the result back onto the traditional Kauffmann-White serotyping scheme still used in public health reporting.

**Prerequisites.** SeqSero2 is a Python package distributed via GitHub and Conda, and its allele micro-assembly and genome-assembly workflows depend on standard read-mapping and assembly tools (including, optionally, SPAdes above) as intermediate steps.

**Popularity.** Serotyping remains, per SeqSero2's own publication, a first-line *Salmonella* subtyping method even as public health surveillance has moved to whole-genome sequencing, and SeqSero (the original tool SeqSero2 replaced) had already seen routine use in food safety and public health laboratories in the United States and other countries prior to this update.

**Efficiency.** SeqSero2's own publication reports 96.1% serotype concordance across all three of its workflows on a large clinical isolate collection from CDC's NARMS surveillance system; however, an independent head-to-head comparison of four in silico *Salmonella* typing tools on a separate collection of 1,624 isolates found a different rank order — SISTR performed best at 94% correctly typed isolates, ahead of SeqSero2 at 87% and the original SeqSero at 81% — a useful, independently sourced reminder that a tool's own validation figures and its performance in an unaffiliated comparative study are not always the same number, and that consulting both is worthwhile before choosing a serotyping tool for a specific pathogen and use case.

- **Citation:** Zhang S, den Bakker HC, Li S, Chen J, Dinsmore BA, Lane C, Lauer AC, Fields PI, Deng X. SeqSero2: rapid and improved *Salmonella* serotype determination using whole-genome sequencing data. *Appl Environ Microbiol.* 2019;85(23):e01746-19. PMID: [31540993](https://pubmed.ncbi.nlm.nih.gov/31540993/)
- **GitHub:** [github.com/denglab/SeqSero2](https://github.com/denglab/SeqSero2)

## Genome Comparison (Pangenome Analysis)

### Roary

Roary calculates the pangenome of a set of related bacterial isolates — partitioning genes into core (present in all or nearly all genomes) and accessory (present in only a subset) categories — from annotated genome assemblies in GFF3 format, and was specifically designed to remain computationally tractable as population sequencing studies scaled from tens to hundreds or thousands of isolates.

**Prerequisites.** Roary is implemented in Perl and distributed via GitHub, and its documentation states explicitly that it expects GFF3 input files containing embedded nucleotide sequence, of exactly the kind Prokka (Chapter 5) produces by default — a direct, format-level dependency rather than a loose recommendation, and one more entry in this series' recurring pattern of tools designed to chain directly into a specific predecessor's output.

**Popularity.** Roary has been a standard tool for bacterial pangenome analysis since 2015, and remains widely used and taught (including in dedicated pangenome tutorials) despite the emergence of documented successor tools, including Panaroo below.

**Efficiency.** The original publication reports that Roary can analyze datasets with thousands of samples on a standard desktop computer, a scale its authors describe as computationally infeasible for existing pangenome methods at the time of publication — the same "make a previously impractical scale tractable" efficiency story told about DADA2, MEGAHIT, and Kraken 2 elsewhere in this series, here applied to gene-content comparison rather than sequence assembly or classification.

- **Citation:** Page AJ, Cummins CA, Hunt M, Wong VK, Reuter S, Holden MTG, Fookes M, Falush D, Keane JA, Parkhill J. Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics.* 2015;31(22):3691–3693. PMID: [26198102](https://pubmed.ncbi.nlm.nih.gov/26198102/)
- **GitHub:** [github.com/sanger-pathogens/Roary](https://github.com/sanger-pathogens/Roary)

### Panaroo

Panaroo re-approaches pangenome inference using a graph-based algorithm that explicitly models and corrects for errors introduced during genome annotation — fragmented gene calls, contamination, and diverse or mis-clustered gene families — by sharing structural information across all genomes in a dataset simultaneously, rather than treating each genome's annotation as independently correct input, as Roary and comparable "gold standard" methods do.

**Prerequisites.** Panaroo is a Python package distributed via GitHub under an MIT license, and like Roary consumes GFF3-format annotated assemblies (again, typically from Prokka or Bakta) as its primary input, additionally depending on CD-HIT for its initial gene-clustering step.

**Popularity.** Panaroo interfaces directly with downstream association-testing tools such as pyseer for linking pangenome gene presence/absence to phenotypes, and its authors explicitly position it as an alternative to, not merely an incremental update of, the "previous gold standard methods" it benchmarks against, of which Roary is the most prominent example discussed in this chapter.

**Efficiency.** The original publication's central demonstration is a highly clonal *Mycobacterium tuberculosis* outbreak dataset used as a deliberate negative control: because the outbreak's maximum pairwise SNP distance was only 9, essentially no true accessory genome variation should exist, yet every tool compared except Panaroo found in excess of 2,500 spurious accessory genes attributable purely to annotation error, while Panaroo alone maintained reasonable error control on this dataset — a controlled demonstration that a pangenome tool's naive gene-count outputs can be dominated by annotation artifacts rather than genuine biological variation if the tool does not explicitly correct for them.

- **Citation:** Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S, Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J. Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biol.* 2020;21(1):180. PMID: [32698896](https://pubmed.ncbi.nlm.nih.gov/32698896/)
- **GitHub:** [github.com/gtonkinhill/panaroo](https://github.com/gtonkinhill/panaroo)

## Comparative Assessment

**Isolate WGS makes the tool-chaining pattern observed throughout this series into an explicit, format-level contract rather than an inferred convention.** Roary's documentation does not merely suggest Prokka as a good annotation choice — it specifies that input must be GFF3 with embedded nucleotide sequence, which is precisely what Prokka (and, by extension, Bakta) produces, and precisely what NCBI's own GFF3 exports lack. This is more binding than the DAS Tool/any-binner or DRAM/GTDB-Tk-CheckM relationships from Chapters 4–5, which accept any tool producing a compatible output format; here, the compatible format is close enough to one specific tool's default behavior that deviating from it is a documented source of user error, not a hypothetical one.

**Two more author-endorsed successions appeared in this chapter's research, reinforcing that this is a recurring norm in bacterial genomics tooling specifically, not an isolated incident from Chapter 5.** Unicycler's own author now points users toward Trycycler and Polypolish for long-read-first assembly, using almost the same phrasing pattern as Prokka's author pointing toward Bakta; Panaroo positions itself explicitly against Roary as the corrected alternative for annotation-error-prone datasets, backed by a designed negative-control experiment rather than a general performance claim. Combined with Prokka→Bakta (Chapter 5) and DEICODE→Gemelli, ANCOM→ANCOM-BC, and SourceTracker→SourceTracker2/FEAST (Chapters 1–3), bacterial and microbiome genomics tooling appears to treat public, author-sanctioned succession as a normal and expected part of a tool's lifecycle, rather than something that only happens when a tool is abandoned or fails.

**AMR detection is the one category in this entire series where two prominent tools are described as actively collaborating rather than competing or succeeding one another.** NCBI's own AMRFinderPlus documentation states a direct, ongoing collaboration with CARD to reconcile database content; CARD's own publication candidly notes RGI's software support lags its own curation. Rather than one tool being positioned as authoritative and the other as legacy (the pattern everywhere else in this series), AMR detection appears to be treated by its own maintainers as a domain where two independently maintained, cross-referencing databases serve the field better than a single consolidated one — worth noting as a genuine structural difference from every other tool relationship surveyed across seven chapters.

## Conclusion

Isolate WGS is not a smaller or simpler version of the metagenomic workflows covered in Chapters 1–6; it is a parallel toolchain built around a different unit of analysis (one genome, known to be a single organism) and a different set of questions (is it fully assembled, does it carry resistance genes, what serotype is it, how does its gene content compare to its relatives). SPAdes and Unicycler recur in this chapter as assembly tools this series has not previously needed, precisely because Chapters 1–6 never needed to assemble a single, complete bacterial genome from a mixed sample; AMRFinderPlus and CARD/RGI address a question (acquired resistance) with no real analog in the community-composition-focused chapters before this one; SeqSero2 addresses subspecies classification via surface antigens, a concept that simply does not apply to a MAG assembled from an environmental sample of unknown taxonomic composition; and Roary/Panaroo's pangenome comparison presupposes exactly the kind of confident, per-isolate annotation that Chapter 5's tools (now revealed as originally isolate-focused) were built to provide. Across all seven chapters of this series, the throughline has held: no tool stands alone, most tools exist in documented relation to a predecessor, a successor, or a collaborator, and understanding those relationships is at least as valuable as knowing any single tool's feature list.

## References

1. Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, Nikolenko SI, Pham S, Prjibelski AD, Pyshkin AV, Sirotkin AV, Vyahhi N, Tesler G, Alekseyev MA, Pevzner PA. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *J Comput Biol.* 2012;19(5):455–477. PMID: 22506599. doi:10.1089/cmb.2012.0021
2. Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. *PLoS Comput Biol.* 2017;13(6):e1005595. PMID: 28594827. doi:10.1371/journal.pcbi.1005595
3. Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep.* 2021;11(1):12728. PMID: 34135355. doi:10.1038/s41598-021-91456-0
4. Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* 2023;51(D1):D690–D699. PMID: 36263822. doi:10.1093/nar/gkac920
5. Zhang S, den Bakker HC, Li S, Chen J, Dinsmore BA, Lane C, Lauer AC, Fields PI, Deng X. SeqSero2: rapid and improved *Salmonella* serotype determination using whole-genome sequencing data. *Appl Environ Microbiol.* 2019;85(23):e01746-19. PMID: 31540993. doi:10.1128/AEM.01746-19
6. Page AJ, Cummins CA, Hunt M, Wong VK, Reuter S, Holden MTG, Fookes M, Falush D, Keane JA, Parkhill J. Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics.* 2015;31(22):3691–3693. PMID: 26198102. doi:10.1093/bioinformatics/btv421
7. Tonkin-Hill G, MacAlasdair N, Ruis C, Weimann A, Horesh G, Lees JA, Gladstone RA, Lo S, Beaudoin C, Floto RA, Frost SDW, Corander J, Bentley SD, Parkhill J. Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biol.* 2020;21(1):180. PMID: 32698896. doi:10.1186/s13059-020-02090-4

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| SPAdes | C++/Python | https://github.com/ablab/spades |
| Unicycler | Python/C++ | https://github.com/rrwick/Unicycler |
| AMRFinderPlus | C++ | https://github.com/ncbi/amr |
| CARD / RGI | Python | https://github.com/arpcard/rgi |
| SeqSero2 | Python | https://github.com/denglab/SeqSero2 |
| Roary | Perl | https://github.com/sanger-pathogens/Roary |
| Panaroo | Python | https://github.com/gtonkinhill/panaroo |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above.*
