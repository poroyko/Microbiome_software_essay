# Python and R Packages for Microbiome Analysis — Chapter 2: Compositional Inference, Diversity Networks, and Source Tracking

## Introduction

Chapter 1 of this series surveyed the general-purpose backbone of microbiome bioinformatics: denoising and classification tools (DADA2, QIIME 2, Kraken 2), data containers and ordination frameworks (phyloseq, scikit-bio), and the most widely adopted differential-abundance methods (ANCOM-BC, MaAsLin2, ALDEx2, metagenomeSeq). This chapter goes one level deeper, examining eight additional tools that either preceded, extended, or directly challenge the assumptions of that first generation of methods: **ANCOM**, the original compositional-analysis framework that ANCOM-BC itself was built to improve upon; **Songbird** and its Bayesian successor **BIRDMAn**, which reframe differential abundance as a regression and inference problem respectively rather than a classical hypothesis test; **corncob**, which models both differential abundance and differential *variability* simultaneously; **SourceTracker2**, which asks a different question entirely — not "what changed" but "where did this community's members come from"; and a pair of tools, **breakaway** and **DivNet**, alongside the network-inference method **SpiecEasi**, which address diversity estimation and microbial association inference under the reality that taxa do not vary independently of one another. As in Chapter 1, every tool is linked to its GitHub repository and, where one exists, its peer-reviewed description indexed in PubMed/NCBI — and, notably, two of the eight tools discussed here remain formally unpublished preprints as of the most recent evidence available, a detail worth flagging rather than glossing over.

## The R Ecosystem

### ANCOM

ANCOM (Analysis of Composition of Microbiomes) was the first widely adopted method to explicitly model microbiome data using compositional log-ratios rather than treating taxon counts as independent quantities, predating and directly motivating its bias-corrected successor, ANCOM-BC, covered in Chapter 1.

**Prerequisites.** Unlike most tools discussed across both chapters, ANCOM was not released as a single, continuously maintained software package; its original implementation was distributed as R code in the paper's supplementary materials. In current practice, researchers most commonly access ANCOM either through QIIME 2's `q2-composition` plugin (the same plugin that now also hosts ANCOM-BC) or through independent third-party reimplementations, such as a Python port (`ancomP`) or a performance-oriented reimplementation (`fastANCOM`) built around the Mann-Whitney U statistic.

**Popularity.** ANCOM's influence is best measured by citation of the underlying methodology rather than downloads of a single canonical package: it remains a standard baseline comparator in nearly every subsequent differential-abundance methods paper, including the ANCOM-BC, corncob, and BIRDMAn papers discussed in this chapter.

**Efficiency.** The original publication reports that ANCOM scales well to comparisons involving thousands of taxa and controls the false discovery rate at the nominal level in settings where a standard t-test or the Zero-Inflated Gaussian (ZIG) method showed substantially inflated false discovery rates, in some instances as high as 68% and 60% respectively — a result that established compositional false-discovery control, rather than raw computational speed, as ANCOM's primary contribution.

- **Citation:** Mandal S, Van Treuren W, White RA, Eggesbø M, Knight R, Peddada SD. Analysis of composition of microbiomes: a novel method for studying microbial composition. *Microb Ecol Health Dis.* 2015;26:27663. PMID: [26028277](https://pubmed.ncbi.nlm.nih.gov/26028277/)
- **GitHub:** No single maintained canonical repository; commonly accessed via QIIME 2's `q2-composition` plugin, or third-party ports such as [github.com/mortonjt/ancomP](https://github.com/mortonjt/ancomP) (Python) and [github.com/ZRChao/fastANCOM](https://github.com/ZRChao/fastANCOM) (R)

### corncob

corncob (COunt Regression for Correlated Observations with the Beta-binomial) models microbial relative abundance using a beta-binomial distribution, which lets it test not only whether a taxon's mean abundance differs between conditions (differential abundance) but also whether its abundance becomes more or less *variable* between conditions (differential variability) — a question most other methods in this survey do not address at all.

**Prerequisites.** corncob is distributed via CRAN and GitHub, and can operate either directly on phyloseq objects or on plain count matrices, requiring no additional non-R dependencies beyond standard regression infrastructure.

**Popularity.** corncob is maintained by the StatDivLab at the University of Washington, the same group behind breakaway and DivNet (below), and is frequently included as a comparator in independent benchmarking studies of differential-abundance methods; a community-contributed Python reimplementation, `pycorncob`, mirrors the LEfSe/`lefser` cross-language porting pattern noted in Chapter 1, this time moving from R into Python rather than the reverse.

**Efficiency.** The original publication demonstrates via simulation that corncob's hypothesis-testing framework maintains valid Type I error control even at the small sample sizes typical of microbiome studies, a property the authors specifically contrast with methods that only perform well asymptotically at large sample sizes — making it a practical choice for the modest cohort sizes common in pilot and exploratory studies.

- **Citation:** Martin BD, Witten D, Willis AD. Modeling microbial abundances and dysbiosis with beta-binomial regression. *Ann Appl Stat.* 2020;14(1):94–115. PMID: [32983313](https://pubmed.ncbi.nlm.nih.gov/32983313/)
- **GitHub:** [github.com/statdivlab/corncob](https://github.com/statdivlab/corncob) (Python port: [github.com/jgolob/pycorncob](https://github.com/jgolob/pycorncob))

### breakaway

breakaway estimates species richness — the total number of taxa present in a community, including those never observed in a given sample — by fitting a nonlinear regression model to the ratios of consecutive frequency counts, a departure from the classical mixed-Poisson models that dominated richness estimation before it.

**Prerequisites.** breakaway is distributed via GitHub (its authors note the CRAN version is no longer actively maintained, and recommend installing directly from GitHub instead) and depends only on standard R infrastructure.

**Popularity.** breakaway is described by its own documentation as intended to be the primary tool for statistical analysis of microbial diversity in R, and its `betta` function for hypothesis testing on richness estimates is, in the authors' own words, comparatively underused relative to its usefulness — a candid, unusual admission worth noting for readers who may not be aware the function exists.

**Efficiency.** breakaway's regression-based approach is specifically designed to remain numerically stable in high-diversity settings where classical estimators struggle; a companion function, `breakaway_nof1`, further addresses the common practical problem of spurious singleton counts introduced by sequencing error, without requiring the singleton count as an input at all.

- **Citation:** Willis A, Bunge J. Estimating diversity via frequency ratios. *Biometrics.* 2015;71(4):1042–1049. PMID: [26038228](https://pubmed.ncbi.nlm.nih.gov/26038228/)
- **GitHub:** [github.com/adw96/breakaway](https://github.com/adw96/breakaway)

### DivNet

DivNet extends diversity estimation to explicitly account for ecological networks — the reality that microbial taxa positively and negatively co-occur, an assumption violated by the standard multinomial model underlying most "plug-in" diversity estimators, including simple Shannon and Simpson index calculations.

**Prerequisites.** DivNet is an R package, developed by the same StatDivLab group as breakaway and corncob, and is designed to integrate directly with breakaway's `betta` function for downstream hypothesis testing on the resulting diversity estimates.

**Popularity.** DivNet is commonly paired with breakaway in applied microbiome studies specifically because the two tools' outputs are designed to interoperate (DivNet's diversity estimates feed directly into breakaway's `betta` testing framework), making the pairing a de facto standard workflow within the R diversity-estimation ecosystem rather than two independently competing tools.

**Efficiency.** The original publication reports that DivNet is fast, accurate, and precise, and performs particularly well with large numbers of taxa and strongly networked communities — conditions under which classical multinomial-based diversity estimators are shown to be most biased — though the authors also note a specific limitation: with a very large number of latent (unobserved) taxa, DivNet may still miss the effects of those taxa, a weakness it shares with several competing estimators.

- **Citation:** Willis AD, Martin BD. Estimating diversity in networked ecological communities. *Biostatistics.* 2022;23(1):207–222. PMID: [32432696](https://pubmed.ncbi.nlm.nih.gov/32432696/)
- **GitHub:** [github.com/adw96/DivNet](https://github.com/adw96/DivNet)

### SpiecEasi

SpiecEasi (Sparse InversE Covariance estimation for Ecological Association and Statistical Inference, implementing the SPIEC-EASI method) infers microbial association networks — which taxa tend to co-occur or mutually exclude one another — by combining a compositional data transformation with sparse graphical model estimation, addressing both the compositionality and the severe underpowering (many more taxa than samples) that undermine naive correlation-based network inference.

**Prerequisites.** SpiecEasi is distributed via CRAN and GitHub and, since version 1.0, requires compiling source code (a Fortran toolchain in particular, which has historically caused installation friction on macOS); the package now depends on a companion package, `pulsar`, for its stability-based model selection procedure.

**Popularity.** SpiecEasi is one of the most widely used tools specifically for microbial network inference (as distinct from simple pairwise correlation, the naive approach discussed and critiqued in Chapter 1's co-occurrence network section), and its underlying compositional data transformation has itself been reused in later network-inference tools built on the same statistical foundation.

**Efficiency.** The original publication reports that SPIEC-EASI outperforms competing network-inference methods on synthetic benchmark data across a variety of underlying network topologies, and specifically identifies the degree distribution of the true underlying network — rather than any property of the inference algorithm itself — as having the largest effect on any method's practical performance, a nuanced finding about the limits of network inference in general, not just this tool.

- **Citation:** Kurtz ZD, Müller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA. Sparse and compositionally robust inference of microbial ecological networks. *PLoS Comput Biol.* 2015;11(5):e1004226. PMID: [25950956](https://pubmed.ncbi.nlm.nih.gov/25950956/)
- **GitHub:** [github.com/zdk123/SpiecEasi](https://github.com/zdk123/SpiecEasi)

## The Python Ecosystem

### Songbird

Songbird performs differential abundance analysis via multinomial regression, framing the problem explicitly around the concept of "reference frames" — a way of choosing a stable subset of taxa against which to measure the relative change of others, addressing common pitfalls the authors identify in how relative abundance comparisons are typically interpreted across samples.

**Prerequisites.** Songbird is a Python package built on TensorFlow, usable either as a standalone command-line tool or as a QIIME 2 plugin; the standalone version integrates with TensorBoard, letting users visually compare the effects of different model parameters across multiple runs.

**Popularity.** Songbird is maintained within the biocore GitHub organization alongside several other tools in this survey (DEICODE/Gemelli from Chapter 1) and integrates with Qurro, a companion visualization tool for interactively exploring the differential rankings it produces — reflecting the same ecosystem-of-interoperating-tools pattern seen with breakaway/DivNet in the R ecosystem above.

**Efficiency.** Songbird's underlying reference-frames approach was validated in the original publication using an oral-microbiome time-series experiment, where it was shown to reduce false positives relative to naive relative-abundance comparison and to produce results consistent across both raw sequencing data and independently cell-count-normalized data — a direct empirical check that the method's compositional corrections track real biological signal rather than a normalization artifact.

- **Citation:** Morton JT, Marotz C, Washburne A, Silverman J, Zaramela LS, Edlund A, Zengler K, Knight R. Establishing microbial composition measurement standards with reference frames. *Nat Commun.* 2019;10(1):2719. PMID: [31222023](https://pubmed.ncbi.nlm.nih.gov/31222023/)
- **GitHub:** [github.com/biocore/songbird](https://github.com/biocore/songbird)

### BIRDMAn

BIRDMAn (Bayesian Inferential Regression for Differential Microbiome Analysis) generalizes Songbird's regression-based approach into a fully Bayesian framework implemented via the Stan probabilistic programming language, allowing users to incorporate prior information and obtain full posterior uncertainty estimates rather than point estimates and p-values alone.

**Prerequisites.** BIRDMAn is a Python package (requiring Python 3.8 or later) that serves as an interface to Stan via `cmdstanpy`, meaning a working Stan installation is a required dependency in addition to the standard Python scientific stack (`numpy`, `scipy`, `arviz`, `xarray`, `biom-format`).

**Popularity.** BIRDMAn is a notably recent addition to this space, and — unlike every other tool discussed in both chapters of this series — remains, per the most recent evidence available, a peer-reviewed-pending bioRxiv preprint rather than a published journal article, despite already being indexed in PubMed and PubMed Central via the NIH's preprint-inclusion pilot program; a related university dissertation lists the manuscript as "submitted for publication" to *Nature Microbiology*, but independent confirmation of that publication was not found. This distinction is worth stating plainly rather than passed over, since it affects how the tool's claims should currently be weighed relative to fully peer-reviewed alternatives.

**Efficiency.** The preprint reports that BIRDMAn models are robust to uneven sequencing depth and provide more than a 20-fold improvement in statistical power over existing differential-abundance methods in the authors' simulations, a striking figure that — precisely because the underlying manuscript has not yet completed peer review — is best treated as a promising but not yet fully independently vetted result.

- **Citation:** Rahman G, Morton JT, Martino C, et al. BIRDMAn: A Bayesian differential abundance framework that enables robust inference of host-microbe associations. bioRxiv [Preprint]. 2023. PMID: [36778470](https://pubmed.ncbi.nlm.nih.gov/36778470/) (preprint; not yet confirmed as peer-reviewed journal publication)
- **GitHub:** [github.com/biocore/BIRDMAn](https://github.com/biocore/BIRDMAn)

### SourceTracker2

SourceTracker2 estimates what proportion of a given "sink" microbial community's composition can be explained by contributions from a set of candidate "source" environments, using Bayesian inference (a Gibbs sampler) — a fundamentally different question from every other tool in this survey, which addresses contamination tracing and environmental mixing rather than differential abundance.

**Prerequisites.** SourceTracker2 is Python 3 software, and is explicitly described by its own documentation as replicating and extending the functionality of the original SourceTracker, which was distributed as an R package; the Python rewrite's principal advantage is parallel execution, with runtime scaling approximately linearly with the number of parallel jobs specified, up to the number of sink samples being analyzed.

**Popularity.** The original SourceTracker method has, per an independent citing publication, been cited over 1,400 times and remains widely recognized as an effective tool for predicting microbial sources across applications ranging from environmental contamination studies to forensic and clinical research; notably, SourceTracker2 itself has never had its own dedicated peer-reviewed publication — its own README instructs users to cite the original 2011 SourceTracker paper "pending publication of SourceTracker 2," a state of affairs that, per available evidence, still holds.

**Efficiency.** SourceTracker2's Gibbs sampling procedure is computationally the most demanding step, and the tool's parallel `--jobs` option is specifically provided to address this, with the important caveat noted in its own documentation that a single sink sample cannot be split across multiple jobs — so the practical speedup ceiling is set by the number of sink samples in a given analysis, not by available compute alone.

- **Citation:** Knights D, Kuczynski J, Charlson ES, Zaneveld J, Mozer MC, Collman RG, Bushman FD, Knight R, Kelley ST. Bayesian community-wide culture-independent microbial source tracking. *Nat Methods.* 2011;8(9):761–763. PMID: [21765408](https://pubmed.ncbi.nlm.nih.gov/21765408/) (describes original SourceTracker; SourceTracker2 itself remains without a dedicated peer-reviewed publication)
- **GitHub:** [github.com/caporaso-lab/sourcetracker2](https://github.com/caporaso-lab/sourcetracker2)

## Comparative Assessment

**Two of eight tools in this chapter are not fully peer-reviewed, and this is a meaningfully different situation from Chapter 1's software-citation cases.** vegan, the `microbiome` R package, and scikit-bio (Chapter 1) have no single paper because they grew as general-purpose toolkits documented instead by their software citation — a citation-practice distinction, not a validation gap. BIRDMAn and SourceTracker2 are different: both are single-method tools whose core claims (a "20-fold improvement in statistical power," in BIRDMAn's case) rest on a manuscript or citation practice that has not completed, or has permanently skipped, formal peer review. This does not mean the tools are unreliable — both come from established, reputable research groups with strong publication track records elsewhere in this survey — but it does mean their specific performance claims warrant the same "promising but independently unverified" framing a careful reader would apply to any preprint finding, a distinction this essay tries to make explicit rather than silently smoothing over by citing a PMID as though it guaranteed peer review.

**A second, smaller software ecosystem has formed around statistical rather than institutional lines.** Chapter 1 showed bioBakery (bioBakery, HUMAnN, MetaPhlAn, LEfSe) and biocore (QIIME 2, DEICODE/Gemelli, scikit-bio) as the two dominant institutional hubs. This chapter reveals a third pattern: the StatDivLab ecosystem (breakaway, DivNet, corncob), where three separately citable tools are explicitly designed to interoperate — DivNet's diversity estimates are built to feed into breakaway's `betta` hypothesis-testing function, and corncob shares both its authorship and its beta-binomial statistical foundation with the group's broader diversity-estimation work. This is a materially different form of "ecosystem" than a shared software platform like QIIME 2: it is a shared statistical philosophy (explicit modeling of sampling variability and, where relevant, ecological network structure) applied consistently across otherwise-independent tools.

**ANCOM's lack of a single canonical repository is itself an instructive data point.** Every other tool across both chapters of this survey has a clear, actively maintained GitHub home. ANCOM does not — its practical availability today is mediated entirely through QIIME 2's plugin system or through third-party reimplementations that postdate the original paper by years. This is a useful reminder that a method's citation count and a method's software maintenance status are not the same thing, and that "how do I actually run this in 2026" is sometimes a materially harder question than "what should I cite."

## Conclusion

The eight tools in this chapter extend, rather than replace, the differential-abundance and network-analysis landscape mapped in Chapter 1. ANCOM's compositional log-ratio framework directly motivated ANCOM-BC's bias correction; corncob's beta-binomial regression and BIRDMAn's Bayesian regression both push differential abundance testing toward richer, more flexible statistical models than the classical hypothesis tests surveyed previously; breakaway and DivNet address diversity estimation with the same seriousness about sampling uncertainty and ecological network structure that SpiecEasi brings to association inference; and SourceTracker2 stands apart entirely, answering a source-attribution question no other tool in either chapter addresses. Taken together with Chapter 1, the fifteen R tools and ten Python tools surveyed across this two-part series make one point consistently: no tool, however well-cited, should be treated as a complete or final answer to a compositional data analysis question — each embeds specific, statable assumptions (independence, network structure, sampling depth, prior information) that a careful analyst should be able to name before choosing which method's output to trust.

## References

1. Mandal S, Van Treuren W, White RA, Eggesbø M, Knight R, Peddada SD. Analysis of composition of microbiomes: a novel method for studying microbial composition. *Microb Ecol Health Dis.* 2015;26:27663. PMID: 26028277. doi:10.3402/mehd.v26.27663
2. Martin BD, Witten D, Willis AD. Modeling microbial abundances and dysbiosis with beta-binomial regression. *Ann Appl Stat.* 2020;14(1):94–115. PMID: 32983313. doi:10.1214/19-AOAS1283
3. Willis A, Bunge J. Estimating diversity via frequency ratios. *Biometrics.* 2015;71(4):1042–1049. PMID: 26038228. doi:10.1111/biom.12332
4. Willis AD, Martin BD. Estimating diversity in networked ecological communities. *Biostatistics.* 2022;23(1):207–222. PMID: 32432696. doi:10.1093/biostatistics/kxaa015
5. Kurtz ZD, Müller CL, Miraldi ER, Littman DR, Blaser MJ, Bonneau RA. Sparse and compositionally robust inference of microbial ecological networks. *PLoS Comput Biol.* 2015;11(5):e1004226. PMID: 25950956. doi:10.1371/journal.pcbi.1004226
6. Morton JT, Marotz C, Washburne A, Silverman J, Zaramela LS, Edlund A, Zengler K, Knight R. Establishing microbial composition measurement standards with reference frames. *Nat Commun.* 2019;10(1):2719. PMID: 31222023. doi:10.1038/s41467-019-10656-5
7. Rahman G, Morton JT, Martino C, et al. BIRDMAn: A Bayesian differential abundance framework that enables robust inference of host-microbe associations. bioRxiv [Preprint]. 2023. PMID: 36778470. doi:10.1101/2023.01.30.526328 (preprint; peer-reviewed publication not confirmed as of latest available evidence)
8. Knights D, Kuczynski J, Charlson ES, Zaneveld J, Mozer MC, Collman RG, Bushman FD, Knight R, Kelley ST. Bayesian community-wide culture-independent microbial source tracking. *Nat Methods.* 2011;8(9):761–763. PMID: 21765408. doi:10.1038/nmeth.1650 (describes original SourceTracker; SourceTracker2 has no dedicated peer-reviewed publication)

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| ANCOM | R (no canonical repo; via QIIME 2 or ports) | [mortonjt/ancomP](https://github.com/mortonjt/ancomP) (Python), [ZRChao/fastANCOM](https://github.com/ZRChao/fastANCOM) (R) |
| corncob | R | https://github.com/statdivlab/corncob |
| pycorncob (Python port of corncob) | Python | https://github.com/jgolob/pycorncob |
| breakaway | R | https://github.com/adw96/breakaway |
| DivNet | R | https://github.com/adw96/DivNet |
| SpiecEasi | R | https://github.com/zdk123/SpiecEasi |
| Songbird | Python | https://github.com/biocore/songbird |
| BIRDMAn | Python | https://github.com/biocore/BIRDMAn |
| SourceTracker2 | Python | https://github.com/caporaso-lab/sourcetracker2 |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above. Where a tool's underlying manuscript remains an unpublished preprint (BIRDMAn) or lacks any dedicated peer-reviewed publication (SourceTracker2, ANCOM's canonical software), this is stated explicitly rather than implied otherwise by the presence of a PMID.*
