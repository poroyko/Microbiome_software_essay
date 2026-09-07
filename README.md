# Python and R Packages for Microbiome Analysis

An eight-chapter series surveying the software ecosystem for microbiome and microbial genomics research — from amplicon sequencing statistics through genome-resolved metagenomics, functional annotation, viral and strain-level resolution, isolate whole-genome sequencing, and outbreak phylogenetics. Every tool entry is grounded in its peer-reviewed description (PubMed/NCBI, cited by PMID) and links to its GitHub repository (or, where one doesn't exist, its actual canonical source).

**77 tools across 8 chapters.** Each chapter follows the same structure: per-tool assessment of **prerequisites**, **popularity**, and **efficiency**, a **Comparative Assessment** synthesizing patterns across that chapter's tools, and a full **References** + **GitHub Repositories** section.

This is a revised edition following a systematic gap-analysis pass: 12 tools (mothur, plus 11 more identified by checking every chapter for major omissions and for tools named as dependencies in the text but never given their own entry) were added after the original release. See **Revision Note** at the end of this file for the full account.

## Table of Contents

| Ch. | Title | File | Tools |
|---|---|---|---|
| 1 | Foundational packages. Prerequisites, Popularity, and Efficiency | [`microbiome_software_essay_chapter1.md`](microbiome_software_essay_chapter1.md) | 20 |
| 2 | Compositional Inference, Diversity Networks, and Source Tracking | [`microbiome_software_essay_chapter2.md`](microbiome_software_essay_chapter2.md) | 8 |
| 3 | Alternative Pipelines, Ordination Visualization, and Longitudinal Statistics | [`microbiome_software_essay_chapter3.md`](microbiome_software_essay_chapter3.md) | 6 |
| 4 | Genome-Resolved Metagenomics (MAG Assembly, Binning, and Quality) | [`microbiome_software_essay_chapter4.md`](microbiome_software_essay_chapter4.md) | 10 |
| 5 | Functional Annotation of MAGs and Metagenomes | [`microbiome_software_essay_chapter5.md`](microbiome_software_essay_chapter5.md) | 8 |
| 6 | Beyond Species-Level MAGs (Viral Metagenomics and Strain Tracking) | [`microbiome_software_essay_chapter6.md`](microbiome_software_essay_chapter6.md) | 6 |
| 7 | Pure Isolate WGS (Assembly, AMR, Serotyping, and Genome Comparison) | [`microbiome_software_essay_chapter7.md`](microbiome_software_essay_chapter7.md) | 8 |
| 8 | Outbreak Phylogenetics and Long-Read Assembly Refinement | [`microbiome_software_essay_chapter8.md`](microbiome_software_essay_chapter8.md) | 7 |

*Tool counts above are primary entries (each with full prerequisites/popularity/efficiency treatment). The master index below additionally lists successor/port tools (e.g. `mia`, `gemelli`, `lefser`) that appear in GitHub tables without a full separate write-up — hence 77 total rows against 73 primary entries.*

## What Each Chapter Covers

**Chapter 1 — Foundational packages.** The R (Bioconductor/CRAN) and Python (QIIME 2 / scikit-bio) ecosystems' core tools for amplicon processing (DADA2, mothur), community data structures and ordination (phyloseq, vegan), differential abundance (ANCOM-BC, MaAsLin2, ALDEx2, metagenomeSeq), multi-omics integration (mixOmics), the QIIME 2 platform, functional prediction (PICRUSt2), compositional ordination (DEICODE), biomarker discovery (LEfSe), shotgun functional/taxonomic profiling (HUMAnN, MetaPhlAn, Kraken 2), and the sequence-processing/classification tools (VSEARCH, RDP Classifier) that sit underneath several of the above.

**Chapter 2 — Compositional inference, diversity networks, and source tracking.** ANCOM-BC's predecessor (ANCOM) and philosophical alternatives (Songbird, BIRDMAn, corncob); diversity estimation under ecological network structure (breakaway, DivNet); microbial association network inference (SpiecEasi); and the original Bayesian microbial source-tracking method (SourceTracker2).

**Chapter 3 — Alternative pipelines and visualization.** A DADA2 alternative (Deblur); the compositional balance-tree method that foreshadowed Songbird (gneiss); a faster source-tracking alternative to SourceTracker2 (FEAST); paired ordination and phylogenetic tree visualization (Emperor, Empress); and dedicated longitudinal/paired-sample statistics (q2-longitudinal).

**Chapter 4 — Genome-resolved metagenomics.** The complete MAG pipeline: assembly (MEGAHIT, metaSPAdes), binning (MetaBAT2, MaxBin2, CONCOCT), multi-binner refinement (DAS Tool), quality assessment (CheckM, CheckM2), dereplication (dRep), and taxonomic classification (GTDB-Tk).

**Chapter 5 — Functional annotation.** Foundational gene prediction (Prodigal), rapid general-purpose annotation (Prokka, Bakta), fine-grained orthology and domain assignment (eggNOG-mapper, InterProScan), targeted database search (KofamScan for KEGG, dbCAN2 for CAZymes), and comprehensive multi-database metabolic synthesis (DRAM).

**Chapter 6 — Beyond species-level MAGs.** Viral sequence identification, quality assessment, prediction, and taxonomy (VirSorter2, CheckV, geNomad, vConTACT2) alongside two structurally different approaches to strain-level resolution (inStrain, StrainPhlAn).

**Chapter 7 — Isolate whole-genome sequencing.** A parallel toolchain for single-organism genomics: assembly (SPAdes, Unicycler), antimicrobial resistance detection (AMRFinderPlus, CARD/RGI), serotyping (SeqSero2), sequence typing (mlst), and pangenome-based genome comparison (Roary, Panaroo).

**Chapter 8 — Outbreak phylogenetics and long-read refinement.** SNP calling, recombination correction, and tree building for outbreak investigation (Snippy, Gubbins, RAxML, IQ-TREE), plus consensus long-read assembly and polishing (Trycycler, Medaka, Polypolish).

## The Complete Tool Index

All 77 tools indexed across all 8 chapters, in the order they appear. For prerequisites, popularity, and efficiency details on any tool, follow the chapter link above.

| Ch. | Package | Language | Repository |
|---|---|---|---|
| 1 | DADA2 | R | https://github.com/benjjneb/dada2 |
| 1 | phyloseq | R | https://github.com/joey711/phyloseq |
| 1 | vegan | R | https://github.com/vegandevs/vegan |
| 1 | microbiome | R | https://github.com/microbiome/microbiome |
| 1 | mia (successor to microbiome) | R | https://github.com/microbiome/mia |
| 1 | ANCOMBC | R | https://github.com/FrederickHuangLin/ANCOMBC |
| 1 | MaAsLin2 | R | https://github.com/biobakery/Maaslin2 |
| 1 | mixOmics | R | https://github.com/mixOmicsTeam/mixOmics |
| 1 | ALDEx2 | R | https://github.com/ggloor/ALDEx2_dev |
| 1 | metagenomeSeq | R | https://github.com/HCBravoLab/metagenomeSeq |
| 1 | QIIME 2 | Python | https://github.com/qiime2/qiime2 |
| 1 | scikit-bio | Python | https://github.com/scikit-bio/scikit-bio |
| 1 | PICRUSt2 | Python | https://github.com/picrust/picrust2 |
| 1 | DEICODE | Python | https://github.com/biocore/DEICODE |
| 1 | gemelli (successor to DEICODE) | Python | https://github.com/biocore/gemelli |
| 1 | LEfSe | Python (+ R core deps) | https://github.com/SegataLab/lefse |
| 1 | lefser (R port of LEfSe) | R | https://bioconductor.org/packages/lefser |
| 1 | HUMAnN | Python | https://github.com/biobakery/humann |
| 1 | MetaPhlAn | Python | https://github.com/biobakery/MetaPhlAn |
| 1 | Kraken 2 | C++ | https://github.com/DerrickWood/kraken2 |
| 1 | mothur | C++ | https://github.com/mothur/mothur |
| 1 | VSEARCH | C/C++ | https://github.com/torognes/vsearch |
| 1 | RDP Classifier | Java | https://github.com/rdpstaff/classifier |
| 2 | ANCOM | R (no canonical repo; via QIIME 2 or ports) | [mortonjt/ancomP](https://github.com/mortonjt/ancomP) (Python), [ZRChao/fastANCOM](https://github.com/ZRChao/fastANCOM) (R) |
| 2 | corncob | R | https://github.com/statdivlab/corncob |
| 2 | pycorncob (Python port of corncob) | Python | https://github.com/jgolob/pycorncob |
| 2 | breakaway | R | https://github.com/adw96/breakaway |
| 2 | DivNet | R | https://github.com/adw96/DivNet |
| 2 | SpiecEasi | R | https://github.com/zdk123/SpiecEasi |
| 2 | Songbird | Python | https://github.com/biocore/songbird |
| 2 | BIRDMAn | Python | https://github.com/biocore/BIRDMAn |
| 2 | SourceTracker2 | Python | https://github.com/caporaso-lab/sourcetracker2 |
| 3 | Deblur | Python | https://github.com/biocore/deblur |
| 3 | gneiss | Python | https://github.com/biocore/gneiss |
| 3 | FEAST | R | https://github.com/cozygene/FEAST |
| 3 | Emperor | Python/JavaScript | https://github.com/biocore/emperor |
| 3 | Empress | Python/JavaScript | https://github.com/biocore/empress |
| 3 | q2-longitudinal | Python (QIIME 2 plugin) | https://github.com/qiime2/q2-longitudinal |
| 4 | MEGAHIT | C++ | https://github.com/voutcn/megahit |
| 4 | metaSPAdes | C++/Python | https://github.com/ablab/spades |
| 4 | MetaBAT2 | C++ | No canonical GitHub repo; source at [bitbucket.org/berkeleylab/metabat](https://bitbucket.org/berkeleylab/metabat) |
| 4 | MaxBin2 | Perl | No canonical GitHub repo; source at [sourceforge.net/projects/maxbin2](https://sourceforge.net/projects/maxbin2) |
| 4 | CONCOCT | Python | https://github.com/BinPro/CONCOCT |
| 4 | DAS Tool | R (CLI wrapper) | https://github.com/cmks/DAS_Tool |
| 4 | CheckM | Python | https://github.com/Ecogenomics/CheckM |
| 4 | CheckM2 | Python | https://github.com/chklovski/CheckM2 |
| 4 | dRep | Python | https://github.com/MrOlm/drep |
| 4 | GTDB-Tk | Python | https://github.com/Ecogenomics/GTDBTk |
| 5 | Prodigal | C | https://github.com/hyattpd/Prodigal |
| 5 | Prokka | Perl | https://github.com/tseemann/prokka |
| 5 | Bakta | Python | https://github.com/oschwengers/bakta |
| 5 | eggNOG-mapper | Python | https://github.com/eggnogdb/eggnog-mapper |
| 5 | InterProScan | Java | https://github.com/ebi-pf-team/interproscan6 |
| 5 | KofamScan | Perl | https://github.com/takaram/kofam_scan |
| 5 | dbCAN2 / run_dbcan | Python | https://github.com/linnabrown/run_dbcan |
| 5 | DRAM | Python | https://github.com/WrightonLabCSU/DRAM |
| 6 | VirSorter2 | Python | https://github.com/jiarong/VirSorter2 |
| 6 | CheckV | Python | No canonical GitHub repo; source at [bitbucket.org/berkeleylab/checkv](https://bitbucket.org/berkeleylab/checkv) |
| 6 | geNomad | Python | https://github.com/apcamargo/genomad |
| 6 | vConTACT2 | Python | No canonical GitHub repo; source at [bitbucket.org/MAVERICLab/vcontact2](https://bitbucket.org/MAVERICLab/vcontact2) |
| 6 | inStrain | Python | https://github.com/MrOlm/inStrain |
| 6 | StrainPhlAn | Python | https://github.com/biobakery/MetaPhlAn |
| 7 | SPAdes | C++/Python | https://github.com/ablab/spades |
| 7 | Unicycler | Python/C++ | https://github.com/rrwick/Unicycler |
| 7 | AMRFinderPlus | C++ | https://github.com/ncbi/amr |
| 7 | CARD / RGI | Python | https://github.com/arpcard/rgi |
| 7 | SeqSero2 | Python | https://github.com/denglab/SeqSero2 |
| 7 | mlst | Perl | https://github.com/tseemann/mlst |
| 7 | Roary | Perl | https://github.com/sanger-pathogens/Roary |
| 7 | Panaroo | Python | https://github.com/gtonkinhill/panaroo |
| 8 | Snippy | Perl | https://github.com/tseemann/snippy |
| 8 | Gubbins | Python/C | https://github.com/nickjcroucher/gubbins |
| 8 | RAxML | C | https://github.com/stamatak/standard-RAxML |
| 8 | IQ-TREE | C++ | https://github.com/iqtree/iqtree2 |
| 8 | Trycycler | Python | https://github.com/rrwick/Trycycler |
| 8 | Medaka | Python | https://github.com/nanoporetech/medaka |
| 8 | Polypolish | Rust | https://github.com/rrwick/Polypolish |

## Cross-Cutting Themes

Several patterns recurred across independent chapters, discovered through the research rather than planned in advance:

- **Explicit, author-endorsed tool succession is the norm, not the exception**, in this field. Original authors have publicly redirected users to newer tools in at least six documented cases across this series: Prokka→Bakta (Ch. 5), ANCOM→ANCOM-BC (Ch. 2), DEICODE→Gemelli (Ch. 1), SourceTracker→SourceTracker2/FEAST (Ch. 2–3), Unicycler→Trycycler→Autocycler (Ch. 7–8), and RAxML→RAxML-NG (Ch. 8). Two individual research groups — Torsten Seemann's and Ryan Wick/Kathryn Holt's — account for an unusually large share of these.
- **Not every canonical repository is on GitHub, and there are at least three different hosting platforms involved.** ANCOM (Ch. 2), MetaBAT2 and CheckV/vConTACT2 (Ch. 4, 6) live on Bitbucket; MaxBin2 (Ch. 4) lives on SourceForge. Two distinct lab hosting traditions are responsible (DOE/Berkeley Lab; Sullivan Lab/MAVERICLab), and notably, geNomad and CheckV — both from the same DOE Joint Genome Institute lineage — don't even share a consistent hosting platform with each other.
- **Not every widely used package has a peer-reviewed paper.** vegan, the `microbiome` R package, scikit-bio (Ch. 1), Snippy and `mlst` (Ch. 7–8), and Medaka (Ch. 8) are all cited as software directly — a real distinction for anyone compiling a methods section, since "cite the paper" isn't an option for these tools. A consistent pattern emerges specifically within Torsten Seemann's own body of work: his genome-scale platform tools (Prokka) get papers; his narrower single-purpose utilities (Snippy, `mlst`) do not.
- **Tool chaining is often a documented, format-level contract, not an inference.** Roary requires Prokka-style GFF3 with embedded sequence (Ch. 7); DRAM explicitly accepts GTDB-Tk and CheckM(2) output as input (Ch. 5); VirSorter2→CheckV→DRAM-v is a published standard operating procedure from one lab (Ch. 6); CheckM and CheckV both depend directly on Prodigal or its derivatives (Ch. 4–6) for their own underlying gene calls.
- **AMR detection is the one domain with documented cross-tool collaboration rather than succession.** NCBI's AMRFinderPlus and CARD actively coordinate database content (Ch. 7) — a structurally different relationship from every tool-succession case elsewhere in this series.
- **"Combine, don't choose" recurs as an explicit design philosophy**, most concentrated in Chapter 8: Snippy → Gubbins → RAxML/IQ-TREE and Trycycler → Medaka → Polypolish → further polishing are both pipelines their own authors describe as multi-stage by design, not single tools competing for the same job.
- **Institutional origin is usually academic, but not always** — Medaka (Ch. 8) is the one tool in this entire series maintained directly by a commercial sequencing-platform vendor (Oxford Nanopore Technologies) rather than an academic or public-health lab, a genuinely distinct funding and incentive structure worth flagging on its own.

## Citation Practice

Every tool description across all 8 chapters is a paraphrased summary of the cited source's findings and documentation — no verbatim quotation beyond short, appropriately attributed phrases. Each chapter's closing note states this explicitly. PMIDs link directly to PubMed; GitHub links are the actual current canonical repository (or an explicit note where none exists).

## License

MIT, @ Valeriy Poroyko.

## Revision Note

This edition reflects two rounds of correction following a systematic gap-analysis review:

**Round 1 (Chapter 1):** mothur (Schloss lab) was identified as a significant omission — one of the two dominant amplicon-analysis platforms alongside QIIME/QIIME 2, and by citation count one of the most influential tools in this entire survey. It was a gap in initial tool selection, not a principled exclusion.

**Round 2 (all chapters):** A full pass was made across all 8 chapters checking for two categories of gap: (a) major, well-known tools absent from the field entirely, and (b) tools named directly in each chapter's own text as dependencies, comparators, or recommendations, but never given their own full entry. This surfaced 11 further additions:

- **Chapter 1:** VSEARCH, RDP Classifier
- **Chapter 4:** metaSPAdes, MaxBin2, CONCOCT, CheckM (original)
- **Chapter 5:** Prodigal
- **Chapter 6:** geNomad
- **Chapter 7:** mlst (a missing sequence-typing category)
- **Chapter 8:** RAxML, Medaka

All 12 additions across both rounds were researched to the same standard as the original tools — verified PMID and GitHub repository (or an explicit note where neither exists) — and integrated into their chapters' existing prose, Comparative Assessment sections, References, and GitHub tables, not merely appended as afterthoughts.
