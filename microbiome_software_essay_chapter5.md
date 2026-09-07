# Python and R Packages for Microbiome Analysis — Chapter 5: Functional Annotation of MAGs and Metagenomes

## Introduction

Chapter 4 ended with a collection of high-quality, dereplicated, taxonomically classified metagenome-assembled genomes (MAGs) — a genuine scientific accomplishment, but not yet a biological answer. Knowing that a genome belongs to a particular GTDB taxon does not say what that organism actually *does*: what metabolic pathways it can run, what carbohydrates it can degrade, what genes it shares by common descent with organisms in other environments. Functional annotation closes this gap, and this chapter covers eight tools that do so at different levels of resolution — from the foundational gene-calling step every other tool here depends on, through whole-genome general-purpose annotation and narrowly targeted database searches, to comprehensive, multi-database metabolic synthesis. **Prodigal** performs the gene prediction that Prokka, Bakta, CheckM, CheckV, and eggNOG-mapper all depend on directly or through a derivative, and is added here specifically because this chapter's own text names it as a dependency without ever giving it its own treatment; **Prokka** and its explicitly designated successor **Bakta** provide rapid, general-purpose genome annotation; **eggNOG-mapper** and **InterProScan** assign fine-grained orthology and protein-domain annotations at genome-to-metagenome scale; **KofamScan** and **dbCAN2** perform targeted annotation against two specific, biologically important databases (KEGG Orthology and carbohydrate-active enzymes, respectively); and **DRAM** integrates several of these databases into a single, MAG-oriented metabolic summary — one that, notably, accepts GTDB-Tk taxonomy and genome-completeness estimates as direct input, closing the loop back to Chapter 4's pipeline explicitly rather than merely by convention.

## Foundational Gene Prediction: Prodigal

Before any of the annotation tools below can assign a function to a gene, something has to first identify where the genes actually are — a task Prodigal performs by finding protein-coding sequences and their translation initiation sites in bacterial and archaeal genomes using an unsupervised dynamic-programming algorithm that requires no training data, learning a genome's own coding statistics, start-codon usage, and ribosome-binding-site motifs directly from the input sequence. Prodigal is not itself an annotation tool — it identifies genes without saying what they do — but it is the specific dependency named, directly or through a derivative, by nearly every other tool in this chapter and several tools in Chapters 4 and 6: Prokka calls it directly, Bakta uses Pyrodigal (a faster reimplementation), eggNOG-mapper v2 uses it for de novo gene prediction from raw contigs, and CheckM/CheckV (Chapters 4 and 6) both depend on it or its viral-aware variant, Prodigal-gv, for their own underlying gene calls.

**Prerequisites.** Prodigal is implemented in C and distributed via GitHub as source code or precompiled binaries, with no external dependencies of its own — a deliberately lightweight design consistent with its role as a component embedded inside larger pipelines rather than a standalone analysis destination.

**Popularity.** Prodigal's ubiquity is best measured architecturally rather than by its own citation count alone: it is the gene-calling engine inside Prokka, Bakta (via Pyrodigal), CheckM, CheckM2, CheckV, and eggNOG-mapper, all covered elsewhere in this series, meaning a very large share of this survey's readers have already run Prodigal, whether or not they invoked it by name.

**Efficiency.** The original publication reports that Prodigal achieved good results compared to existing gene-prediction methods available at the time, with a specific emphasis on improving translation initiation site identification and reducing false positives — two goals the authors frame as at least as important as raw gene-count sensitivity, since an incorrectly placed start codon propagates errors into every downstream annotation step that depends on the resulting protein sequence.

- **Citation:** Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics.* 2010;11:119. PMID: [20211023](https://pubmed.ncbi.nlm.nih.gov/20211023/)
- **GitHub:** [github.com/hyattpd/Prodigal](https://github.com/hyattpd/Prodigal)

## General-Purpose Genome Annotation: Prokka and Bakta

Prokka annotates a draft bacterial, archaeal, or viral genome by coordinating a suite of existing feature-prediction tools (Prodigal for coding sequences, plus dedicated tools for rRNA, tRNA, and other non-coding features) into one command-line pipeline, producing standards-compliant output files without requiring a web- or email-based submission system.

**Prerequisites.** Prokka is a Perl-based command-line tool, installed via Conda or from source, that depends on BioPerl and BLAST+ alongside its bundled feature-prediction tools; it requires only a FASTA file of assembled contigs as mandatory input.

**Popularity.** Prokka's publication has accumulated over 5,400 citations, reflecting more than a decade as one of the most widely used bacterial genome annotation tools available — but its own GitHub repository now carries an unusually direct statement from its author: "I can no longer maintain it... I recommend replacing Prokka with Bakta in your analysis pipelines going forward," an explicit, author-sanctioned succession rather than the more common pattern (seen throughout this series) of a newer tool simply outcompeting an older one without the original author's endorsement.

**Efficiency.** The original publication reports that Prokka fully annotates a draft bacterial genome in about 10 minutes on a typical desktop computer, exploiting multiple processing cores where available — a design explicitly aimed at integration into automated genomic pipelines rather than one-off manual analysis.

Bakta, Prokka's designated successor, replaces several of Prokka's homology-search steps with an alignment-free sequence identification (AFSI) approach that exactly matches known protein sequences from RefSeq and UniProt before falling back to alignment-based search, and is explicitly designed to remain accurate on both well-studied species and the taxonomically novel genomes typical of MAGs.

**Prerequisites.** Bakta is a Python package installed via Conda or pip, depending on Pyrodigal (a faster reimplementation of Prodigal), DIAMOND, and several of the same non-coding-feature tools Prokka uses, alongside a substantially larger reference database that must be downloaded separately.

**Popularity.** Bakta has rapidly become one of the most established tools for bacterial genome annotation since its 2021 publication, and a dedicated web-based version, Bakta Web, was published in 2025 specifically to lower the technical barrier of installing and running the command-line tool for users less comfortable with bioinformatics infrastructure.

**Efficiency.** The original publication reports that Bakta achieves favorable annotation sensitivity and specificity across the full continuum from well-studied species to unknown MAG-derived genomes, outperforming other tools in the assignment of functional categories and database cross-references while providing wall-clock runtimes comparable to existing alternatives — the alignment-free identification step exists specifically to preserve Prokka-like speed while improving on Prokka's accuracy for novel lineages.

- **Citations:** Seemann T. Prokka: rapid prokaryotic genome annotation. *Bioinformatics.* 2014;30(14):2068–2069. PMID: [24642063](https://pubmed.ncbi.nlm.nih.gov/24642063/). Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom.* 2021;7(11):000685. PMID: [34739369](https://pubmed.ncbi.nlm.nih.gov/34739369/)
- **GitHub:** [github.com/tseemann/prokka](https://github.com/tseemann/prokka) and [github.com/oschwengers/bakta](https://github.com/oschwengers/bakta)

## Orthology and Domain Annotation: eggNOG-mapper and InterProScan

eggNOG-mapper assigns functional annotations to protein sequences by fast orthology assignment against precomputed eggNOG clusters and phylogenies, rather than performing a fresh homology search and functional inference for every query — a design choice that lets it scale specifically to metagenomic gene catalogs containing millions of predicted proteins.

**Prerequisites.** eggNOG-mapper is a Python package distributed via Conda and GitHub, and depends on a large precomputed eggNOG database plus a chosen search backend (DIAMOND, MMseqs2, or the slower but more sensitive HMMER); version 2 added the ability to perform de novo gene prediction directly from raw assembled contigs via Prodigal, removing a previously required separate gene-calling step.

**Popularity.** eggNOG-mapper is one of the most widely used tools for metagenomic functional annotation, and version 2's own publication demonstrates its scale specifically on microbiome data: reannotating 1.75 million proteins subsampled from a human gut metagenomic gene catalog yielded a 3.23% increase in annotation coverage over the previous version, on a dataset large enough that a fractional-percentage improvement still represents tens of thousands of newly annotated proteins.

**Efficiency.** The original publication reports that eggNOG-mapper v2 improves annotation rate (queries annotated per second) by 16% on average compared to version 1, despite the underlying reference databases having doubled in size, and benchmarks its annotation speed directly against Prokka on matched sets of input genomes — a rare case in this series of two tools from entirely different lineages being compared head-to-head on the same efficiency metric.

InterProScan classifies protein sequences by scanning them against a diverse collection of predictive models spanning multiple member databases (Pfam, PANTHER, PRINTS, and others) that InterPro integrates into a single combined output, providing broader functional coverage than any single member database alone.

**Prerequisites.** InterProScan is a Java-based tool distributed via GitHub, designed from its version 5 rewrite onward to exploit both multiprocessor machines and conventional compute clusters for scalable distributed analysis; it depends on a substantial collection of member-database data files that must be installed alongside the software itself.

**Popularity.** InterProScan 5's 2014 publication remains a standard citation for protein function classification at genome scale, and — in a detail unusually current for a tool this well-established — a successor, InterProScan 6, was published in *Bioinformatics Advances* in 2026, explicitly to modernize the pipeline's implementation while reproducing version 5's annotation results with near-identical precision and sensitivity.

**Efficiency.** InterProScan 5's central contribution, per its original publication, was a complete reimplementation of the software's underlying Java architecture specifically to address scalability for the "many millions of sequences" biologists were beginning to routinely generate; InterProScan 6's 2026 publication reports that this modernized pipeline reproduces the same annotations more efficiently and reproducibly, without sacrificing the accuracy of the version 5 results it replaces.

- **Citations:** Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol.* 2021;38(12):5825–5829. PMID: [34597405](https://pubmed.ncbi.nlm.nih.gov/34597405/). Jones P, Binns D, Chang HY, et al. InterProScan 5: genome-scale protein function classification. *Bioinformatics.* 2014;30(9):1236–1240. PMID: [24451626](https://pubmed.ncbi.nlm.nih.gov/24451626/)
- **GitHub:** [github.com/eggnogdb/eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) and [github.com/ebi-pf-team/interproscan6](https://github.com/ebi-pf-team/interproscan6)

## Targeted Database Annotation: KofamScan and dbCAN2

KofamScan assigns KEGG Orthology (KO) identifiers to protein sequences using a database of profile hidden Markov models (KOfam) paired with precomputed, family-specific adaptive score thresholds — a design intended to make the score cutoff for a positive match appropriately strict or lenient depending on the specific KO family being searched, rather than applying one fixed threshold across every family.

**Prerequisites.** KofamScan is a Perl-based command-line tool (with a newer, faster Python-based community reimplementation, PyKofamSearch, built on PyHMMER for high-memory systems) distributed alongside its KOfam database from GenomeNet, requiring HMMER as its underlying search engine.

**Popularity.** KofamScan is commonly paired with dbCAN2 (below) in modern MAG-annotation workflows specifically to jointly cover KEGG-pathway-level metabolic potential and carbohydrate-degradation capacity in the same analysis, a pairing reflected directly in several published microbiome studies' methods sections.

**Efficiency.** The original publication reports that KofamKOALA (the web-server front end to KofamScan) is faster than existing KO assignment tools (BlastKOALA, GhostKOALA, KAAS) while achieving accuracy comparable to the best-performing among them, benchmarked on a test set of 40 diverse genomes spanning both eukaryotes and prokaryotes.

dbCAN2 automates the annotation of carbohydrate-active enzymes (CAZymes) — the enzyme families responsible for synthesizing, modifying, and degrading complex carbohydrates — by combining three independent search strategies (HMMER against the dbCAN HMM database, DIAMOND against curated CAZyme sequences, and a short-conserved-motif search) and retaining only CAZyme calls supported by more than one method, directly relevant to microbiome function given how central dietary fiber and complex carbohydrate metabolism are to gut microbial ecology.

**Prerequisites.** The standalone version (`run_dbcan`, now in its third major iteration as dbCAN3) is a Python package distributed via GitHub and Conda, bundling HMMER and DIAMOND as its underlying search engines alongside the dbCAN-specific reference databases.

**Popularity.** dbCAN2 is the standard reference tool for CAZyme annotation in microbiome studies concerned with carbohydrate metabolism, and its associated database has been directly incorporated into other bioinformatics pipelines' workflows, including MOCAT2 and proGenomes, rather than remaining a standalone destination tool only.

**Efficiency.** The original publication reports that combining the three independent search methods and requiring agreement between at least two significantly improves CAZome annotation accuracy relative to any single method used alone — an explicit design decision to trade some raw sensitivity (a CAZyme called by only one method is discarded) for higher-confidence calls, particularly relevant given how error-prone single-method carbohydrate-enzyme classification has historically been.

- **Citations:** Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K, Kanehisa M, Goto S, Ogata H. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics.* 2020;36(7):2251–2252. PMID: [31742321](https://pubmed.ncbi.nlm.nih.gov/31742321/). Zhang H, Yohe T, Huang L, Entwistle S, Wu P, Yang Z, Busk PK, Xu Y, Yin Y. dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Res.* 2018;46(W1):W95–W101. PMID: [29771380](https://pubmed.ncbi.nlm.nih.gov/29771380/)
- **GitHub:** [github.com/takaram/kofam_scan](https://github.com/takaram/kofam_scan) and [github.com/linnabrown/run_dbcan](https://github.com/linnabrown/run_dbcan)

## Comprehensive Metabolic Synthesis: DRAM

DRAM (Distilled and Refined Annotation of Metabolism) integrates annotations from multiple databases at once — KEGG, Pfam, dbCAN, MEROPS, and UniRef among them — into a single, MAG-oriented output with three levels of detail: a raw per-gene annotation table, a "distillate" summarizing metabolic pathway completeness, and a "product" visualization intended to make a genome's overall metabolic potential interpretable at a glance, addressing what its authors describe as the absence of a scalable, metabolically resolved annotation framework prior to its development.

**Prerequisites.** DRAM is a Python package distributed via Conda and GitHub, and explicitly accepts externally computed taxonomy (e.g., from GTDB-Tk) and genome-completeness estimates (e.g., from CheckM) as direct input — a deliberate design choice connecting DRAM's output directly to the Chapter 4 pipeline's endpoint rather than requiring these steps to be redone; its reference database collection requires substantial disk space (commonly cited at 50 GB or more) and is best run on high-performance computing infrastructure.

**Popularity.** DRAM has become a standard tool specifically for MAG-level metabolic characterization in environmental and host-associated microbiome studies, and its companion mode, DRAM-v, extends the same framework to identify virally encoded auxiliary metabolic genes (AMGs) in viral MAGs — genes a bacteriophage carries that can augment its host's metabolism during infection, a distinct and increasingly studied class of finding in viral metagenomics.

**Efficiency.** Rather than a runtime benchmark, DRAM's original publication demonstrates its value through breadth of applicability, evaluating the tool across metabolically diverse genomes and explicitly designing its three-tiered output (raw, distillate, product) to let a user choose between full annotation detail and a rapidly interpretable metabolic summary depending on the scale of the analysis — a structural approach to managing the tension between annotation completeness and interpretability, rather than a claim about processing speed.

- **Citation:** Shaffer M, Borton MA, McGivern BB, Zayed AA, La Rosa SL, Solden LM, Liu P, Narrowe AB, Rodríguez-Ramos J, Bolduc B, Gazitúa MC, Daly RA, Smith GJ, Vik DR, Pope PB, Sullivan MB, Roux S, Wrighton KC. DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Res.* 2020;48(16):8883–8900. PMID: [32766782](https://pubmed.ncbi.nlm.nih.gov/32766782/)
- **GitHub:** [github.com/WrightonLabCSU/DRAM](https://github.com/WrightonLabCSU/DRAM)

## Comparative Assessment

**This chapter contains the clearest, most explicit tool-succession statement anywhere in this five-chapter series.** Every prior succession discussed (ANCOM→ANCOM-BC, DEICODE→Gemelli, SourceTracker→SourceTracker2/FEAST, GTDB-Tk v1→v2) has been a case of a newer tool superseding an older one through adoption, benchmarking, or documented technical improvement — but Prokka's own maintainer has written, in the tool's own README, an explicit recommendation to switch to Bakta. This is a qualitatively different and rarer event: an original author publicly redirecting the community to a different team's successor tool, rather than simply stepping back from maintenance and letting adoption patterns run their course.

**Annotation tools increasingly compose with each other's specific outputs, not just their general file formats.** DRAM's direct acceptance of GTDB-Tk taxonomy and CheckM completeness scores as input is a tighter form of integration than the "same feature-table format" compatibility seen in earlier chapters (e.g., DAS Tool consuming any binner's contig-to-bin mapping) — DRAM is explicitly designed around the assumption that a GTDB-Tk/CheckM(2) pipeline has already run, making Chapters 4 and 5 of this series less like two independent surveys and more like two halves of one continuous, tool-aware workflow.

**Redundant, disagreement-tolerant search strategies recur as a distinct efficiency pattern in this chapter.** dbCAN2's requirement that at least two of three independent search methods agree before calling a CAZyme is conceptually the same strategy DAS Tool (Chapter 4) uses across multiple binners and the same principle underlying ensemble methods generally — trading some sensitivity for higher-confidence output by requiring independent agreement, rather than trusting any single method's calls at face value. This is a third distinct occurrence of the same design idea across two chapters, suggesting it is a genuinely recurring solution to a genuinely recurring problem (any single computational prediction method has failure modes; requiring independent corroboration is a reliable, if conservative, mitigation) rather than a coincidence.

## Conclusion

Functional annotation is where genome-resolved metagenomics stops being a genome-assembly exercise and starts answering biological questions: what can this organism do, and does that capability help explain the community or host phenotype under study? The tools in this chapter span a full range of resolution and integration, from Prokka/Bakta's rapid whole-genome sweep, through eggNOG-mapper and InterProScan's fine-grained orthology and domain assignment, to KofamScan and dbCAN2's targeted, database-specific searches, to DRAM's explicit synthesis of several of these signals into one MAG-level metabolic picture. Combined with the amplicon-processing and statistical tools of Chapters 1–3 and the assembly-through-classification pipeline of Chapter 4, this chapter completes a genuinely end-to-end map of modern microbiome bioinformatics — from raw sequencing reads, through community composition and diversity statistics, to reconstructed genomes, to the metabolic functions those genomes actually encode. Consistent with every prior chapter's closing observation, and perhaps most explicitly here given Prokka's own author's words, these tools are best understood not as a menu of independent choices but as a connected, evolving pipeline in which today's standard tool is openly expected to be tomorrow's superseded one.

## References

1. Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics.* 2010;11:119. PMID: 20211023. doi:10.1186/1471-2105-11-119
2. Seemann T. Prokka: rapid prokaryotic genome annotation. *Bioinformatics.* 2014;30(14):2068–2069. PMID: 24642063. doi:10.1093/bioinformatics/btu153
3. Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom.* 2021;7(11):000685. PMID: 34739369. doi:10.1099/mgen.0.000685
4. Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol.* 2021;38(12):5825–5829. PMID: 34597405. doi:10.1093/molbev/msab293
5. Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, McWilliam H, Maslen J, Mitchell A, Nuka G, Pesseat S, Quinn AF, Sangrador-Vegas A, Scheremetjew M, Yong SY, Lopez R, Hunter S. InterProScan 5: genome-scale protein function classification. *Bioinformatics.* 2014;30(9):1236–1240. PMID: 24451626. doi:10.1093/bioinformatics/btu031
6. Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K, Kanehisa M, Goto S, Ogata H. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics.* 2020;36(7):2251–2252. PMID: 31742321. doi:10.1093/bioinformatics/btz859
7. Zhang H, Yohe T, Huang L, Entwistle S, Wu P, Yang Z, Busk PK, Xu Y, Yin Y. dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Res.* 2018;46(W1):W95–W101. PMID: 29771380. doi:10.1093/nar/gky418
8. Shaffer M, Borton MA, McGivern BB, Zayed AA, La Rosa SL, Solden LM, Liu P, Narrowe AB, Rodríguez-Ramos J, Bolduc B, Gazitúa MC, Daly RA, Smith GJ, Vik DR, Pope PB, Sullivan MB, Roux S, Wrighton KC. DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Res.* 2020;48(16):8883–8900. PMID: 32766782. doi:10.1093/nar/gkaa621

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| Prodigal | C | https://github.com/hyattpd/Prodigal |
| Prokka | Perl | https://github.com/tseemann/prokka |
| Bakta | Python | https://github.com/oschwengers/bakta |
| eggNOG-mapper | Python | https://github.com/eggnogdb/eggnog-mapper |
| InterProScan | Java | https://github.com/ebi-pf-team/interproscan6 |
| KofamScan | Perl | https://github.com/takaram/kofam_scan |
| dbCAN2 / run_dbcan | Python | https://github.com/linnabrown/run_dbcan |
| DRAM | Python | https://github.com/WrightonLabCSU/DRAM |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above.*
