# Python and R Packages for Microbiome Analysis — Chapter 3: Alternative Pipelines, Ordination Visualization, and Longitudinal Statistics

## Introduction

Chapters 1 and 2 of this series established a working map of the microbiome bioinformatics landscape: sequence processing, differential abundance, diversity estimation, network inference, and source tracking. This chapter examines six tools that relate to that map not as new categories but as **alternatives to, visual companions of, or temporal extensions of** methods already covered — a deliberately different organizing principle from the previous two chapters' R-ecosystem/Python-ecosystem split, since five of the six tools here are Python-based and the sixth (FEAST) is better understood by its relationship to SourceTracker2 than by its implementation language alone. Specifically: **Deblur** offers a different algorithmic route than DADA2 to the same sub-OTU resolution goal; **gneiss** introduced the compositional balance-tree concept that directly foreshadowed Songbird's reference-frames approach from Chapter 2; **FEAST** answers the same source-tracking question as SourceTracker2, dramatically faster; and **Emperor**, **Empress**, and **q2-longitudinal** extend the ordination and statistical-modeling work of Chapters 1 and 2 into interactive 3D visualization, phylogenetically informed tree-and-ordination display, and dedicated longitudinal-study statistics, respectively — the last of these a particularly natural fit given this series' recurring attention to temporal microbiome designs. As before, every tool is linked to its GitHub repository and its peer-reviewed PubMed/NCBI description.

## Alternative Sequence Processing: Deblur

Deblur resolves sub-OTU (single-nucleotide) sequence variants from amplicon data — the same goal DADA2 pursues — but through a different algorithmic strategy: a greedy deconvolution that uses a known, pre-characterized read-error profile to separate true biological sequences from sequencing error, rather than DADA2's approach of inferring a sample-specific error model directly from the data.

**Prerequisites.** Deblur is a Python package distributed via GitHub and Conda, and is most commonly used through its QIIME 2 plugin (`q2-deblur`), which serves as a drop-in replacement for the `q2-dada2` step in a standard QIIME 2 workflow.

**Popularity.** Deblur was used to process the 16S data underlying the Earth Microbiome Project, one of the largest standardized microbiome meta-analyses ever assembled (spanning roughly 100 studies), giving it an outsized influence on the reference datasets many other tools in this series are subsequently benchmarked against.

**Efficiency.** The original publication reports that Deblur substantially reduces computational demands relative to comparable sub-OTU methods while achieving similar or better sensitivity and specificity in simulations, mock community mixtures, and real datasets, crediting this partly to Deblur's per-sample operation, which lets it scale to large meta-analyses without requiring joint processing across all samples simultaneously — the same practical benefit DADA2 achieves through a different route (Chapter 1).

- **Citation:** Amir A, McDonald D, Navas-Molina JA, Kopylova E, Morton JT, Zech Xu Z, Kightley EP, Thompson LR, Hyde ER, Gonzalez A, Knight R. Deblur rapidly resolves single-nucleotide community sequence patterns. *mSystems.* 2017;2(2):e00191-16. PMID: [28289731](https://pubmed.ncbi.nlm.nih.gov/28289731/)
- **GitHub:** [github.com/biocore/deblur](https://github.com/biocore/deblur)

## Compositional Balance Analysis: gneiss

gneiss introduced **balance trees** to microbiome analysis: rather than asking how any single taxon's abundance changes, it uses isometric log-ratio (ILR) transforms to ask how the *balance* between two subsets of a hierarchically organized community shifts — sidestepping the compositional-dependence problem (Chapter 1's Introduction) by design rather than by post hoc correction. This concept directly foreshadows the "reference frames" idea underlying Songbird (Chapter 2), published by an overlapping set of authors two years later.

**Prerequisites.** gneiss is a Python package, installable via Conda, whose core compositional statistics and tree data structures were folded into scikit-bio (version 0.4.1 and later) — meaning any environment with a sufficiently current scikit-bio already carries part of gneiss's underlying infrastructure.

**Popularity.** gneiss's balance-tree concept has been directly cited as motivation in a range of subsequent compositional methods beyond Songbird, and its central claim — that conventional univariate tests such as the t-test and Mann-Whitney test can flag nearly 100% of taxa as significantly different across environments purely as a compositional artifact — is one of the more widely quoted cautionary findings in the compositional-data-analysis literature.

**Efficiency.** The original publication demonstrates gneiss's approach using real soil and cystic fibrosis lung sputum datasets, showing that balances coupled with linear mixed-effects models produce more statistically robust results than naive per-taxon testing, without requiring the correction step that a method like ANCOM-BC must apply after the fact — the compositional robustness is structural to the balance-tree representation itself.

- **Citation:** Morton JT, Sanders J, Quinn RA, McDonald D, Gonzalez A, Vázquez-Baeza Y, Navas-Molina JA, Song SJ, Metcalf JL, Hyde ER, Lladser M, Dorrestein PC, Knight R. Balance trees reveal microbial niche differentiation. *mSystems.* 2017;2(1):e00162-16. PMID: [28144630](https://pubmed.ncbi.nlm.nih.gov/28144630/)
- **GitHub:** [github.com/biocore/gneiss](https://github.com/biocore/gneiss)

## Fast Source Tracking: FEAST

FEAST (Fast Expectation-mAximization microbial Source Tracking) answers the same question as SourceTracker2 (Chapter 2) — what proportion of a sink community's composition originates from each of a set of candidate source environments — using an expectation-maximization algorithm instead of SourceTracker2's Gibbs sampling, trading some of the Bayesian framework's flexibility for substantially faster computation.

**Prerequisites.** Unlike SourceTracker2, FEAST is an R package (requiring R 3.4.4 or later), with dependencies on `Rcpp` and `RcppArmadillo` for its performance-critical routines alongside `vegan` — a direct dependency link back to the R ecosystem's ordination and diversity infrastructure covered in Chapter 1.

**Popularity.** FEAST's own publication figures include a direct running-time comparison against "current state-of-the-art" methods (implicitly including SourceTracker2), and subsequent literature has come to describe SourceTracker and FEAST together as the two leading tools specifically developed for microbial source tracking, with newer methods such as STENSL (developed by an overlapping author group) building directly on FEAST's codebase.

**Efficiency.** The original publication's central efficiency claim is that FEAST can simultaneously estimate the contribution of thousands of potential source environments in a timely manner — a scale that Gibbs-sampling-based approaches like SourceTracker2 handle far more slowly — making FEAST the more practical choice specifically for studies with a very large number of candidate sources, such as large-scale environmental contamination surveys.

- **Citation:** Shenhav L, Thompson M, Joseph TA, Briscoe L, Furman O, Bogumil D, Mizrahi I, Pe'er I, Halperin E. FEAST: fast expectation-maximization for microbial source tracking. *Nat Methods.* 2019;16(7):627–632. PMID: [31182859](https://pubmed.ncbi.nlm.nih.gov/31182859/)
- **GitHub:** [github.com/cozygene/FEAST](https://github.com/cozygene/FEAST)

## Ordination and Phylogenetic Visualization: Emperor and Empress

Emperor and Empress are companion visualization tools from the same research group: Emperor renders interactive 3D principal coordinates plots (the direct visual output of the PCoA methods covered throughout Chapters 1 and 2), while Empress renders large phylogenetic trees and can display them side by side with, and interactively linked to, an Emperor ordination — selecting a tip in Empress's tree highlights the corresponding samples in Emperor's ordination, and vice versa.

**Prerequisites.** Both are Python packages usable either as standalone command-line tools producing self-contained HTML/JavaScript output or as QIIME 2 plugins producing `.qzv` visualization files; Empress additionally requires a phylogenetic tree as input, which is not required for a basic Emperor ordination plot.

**Popularity.** Emperor has been in continuous use since 2013 as the default 3D ordination viewer within the QIIME/QIIME 2 ecosystem, while Empress, released in 2021, has already been applied in large-scale contexts including a paired phylogenetic tree and ordination of the Earth Microbiome Project data (756,377 tree nodes alongside an unweighted UniFrac ordination of 26,035 samples) and in COVID-19 epidemiological surveillance work at UC San Diego.

**Efficiency.** Both tools report specific, scale-oriented performance claims rather than general speed comparisons: Emperor's original publication emphasizes a reduced loading time and memory footprint relative to prior visualization approaches, achieved in part through its small output file size (enabling sharing via email or embedded web pages without additional plugins), while Empress's publication demonstrates scalability to phylogenetic trees with well over 500,000 nodes, a scale the authors position explicitly against the limitations of earlier tree-viewing approaches.

- **Citations:** Vázquez-Baeza Y, Pirrung M, Gonzalez A, Knight R. EMPeror: a tool for visualizing high-throughput microbial community data. *GigaScience.* 2013;2(1):16. PMID: [24280061](https://pubmed.ncbi.nlm.nih.gov/24280061/). Cantrell K, Fedarko MW, Rahman G, et al. EMPress enables tree-guided, interactive, and exploratory analyses of multi-omic data sets. *mSystems.* 2021;6(2):e01216-20. PMID: [33727399](https://pubmed.ncbi.nlm.nih.gov/33727399/)
- **GitHub:** [github.com/biocore/emperor](https://github.com/biocore/emperor) and [github.com/biocore/empress](https://github.com/biocore/empress)

## Longitudinal Statistics: q2-longitudinal

q2-longitudinal is a QIIME 2 plugin purpose-built for the temporal and paired-sample designs this series has returned to repeatedly — it provides paired differences and distances testing, linear mixed-effects models, a nonparametric microbial interdependence test, first-differencing, and volatility analysis, wrapping several of these in a way specifically intended to make appropriate longitudinal statistics accessible to microbiome researchers without requiring them to implement mixed-effects modeling from scratch.

**Prerequisites.** As a QIIME 2 plugin, q2-longitudinal inherits QIIME 2's full installation and dependency footprint (Chapter 1); its feature-volatility analysis additionally wraps `q2-sample-classifier`, giving users access to multiple scikit-learn supervised-learning regressors (a random forest regressor by default) for predicting a continuous variable such as age or time from community composition.

**Popularity.** The plugin's authors frame its motivation directly around a documented gap: appropriate longitudinal statistical methods for microbiome data existed before q2-longitudinal, but were difficult for non-specialists to implement, requiring programming skills or statistical training that many microbiome researchers do not have, a barrier the plugin's QIIME 2 integration (including a graphical interface) is explicitly designed to remove.

**Efficiency.** Rather than reporting a runtime benchmark, the original publication demonstrates q2-longitudinal's statistical value directly on a real infant-gut-microbiome dataset (the ECAM study), where a linear mixed-effects model available through the plugin detected a statistically significant effect of diet on *Bifidobacterium* relative abundance at six months of age, alongside a significant diet-by-delivery-mode interaction — the kind of finding that requires exactly the repeated-measures modeling this series' own mock-community-generator work (and the plugin's `betta`-adjacent StatDivLab counterparts in Chapter 2) was built to support.

- **Citation:** Bokulich NA, Dillon MR, Zhang Y, Rideout JR, Bolyen E, Li H, Albert PS, Caporaso JG. q2-longitudinal: longitudinal and paired-sample analyses of microbiome data. *mSystems.* 2018;3(6):e00219-18. PMID: [30505944](https://pubmed.ncbi.nlm.nih.gov/30505944/)
- **GitHub:** [github.com/qiime2/q2-longitudinal](https://github.com/qiime2/q2-longitudinal)

## Comparative Assessment

**Alternatives are not always drop-in replacements, and this chapter's tools illustrate why in three different ways.** Deblur and DADA2 pursue the same sub-OTU resolution goal via genuinely different algorithms (pre-characterized error profiles versus data-inferred error models), so a study's choice between them is a methodological decision, not merely a preference; FEAST and SourceTracker2 answer the identical scientific question with a real, measured speed-versus-flexibility trade-off (expectation-maximization point estimates versus full Bayesian posteriors); and gneiss's balance trees versus ANCOM/ANCOM-BC's bias-corrected hypothesis testing represent two entirely different philosophies for handling compositionality — one restructures the question (asking about balances between taxon subsets) while the other corrects the answer (adjusting a per-taxon test for compositional bias). None of these three pairs is a simple upgrade path; each requires the analyst to understand what specifically changed.

**Visualization tools accumulate real, cross-referenced integration rather than existing as one-off outputs.** Emperor (2013) and Empress (2021) are separated by eight years, yet Empress was explicitly designed to interoperate with Emperor's existing ordination displays rather than to replace them — a single "Empire plot" now links a phylogenetic tree to a PCoA ordination interactively. This mirrors the StatDivLab pattern from Chapter 2 (breakaway feeding into DivNet) and the biocore pattern from Chapter 1 (DEICODE's output format designed for QIIME 2 plugin compatibility): the most durable tools in this space are increasingly the ones built to compose with siblings, not the ones that stand alone.

**q2-longitudinal's stated motivation is a useful check on this whole series.** Its authors justify the plugin's existence by pointing out that the *statistical methods* (mixed-effects models, paired testing) already existed — the gap was accessibility, not novel statistics. This is a useful reminder when evaluating any tool in this three-chapter survey: a new package is not always solving a new mathematical problem; sometimes, as here, its genuine contribution is packaging an existing solution so a working microbiologist can use it correctly without a statistics degree, which is a real and undervalued form of scientific contribution in its own right.

## Conclusion

This chapter closes the loop on several threads left open in Chapters 1 and 2: Deblur shows that DADA2's sub-OTU resolution goal was independently reachable by a different algorithmic path; gneiss reveals the compositional-balance idea that Songbird would later generalize into reference frames; FEAST demonstrates that SourceTracker2's Bayesian source-tracking framework has a faster, if less flexible, expectation-maximization counterpart; and Emperor, Empress, and q2-longitudinal show that ordination, phylogenetics, and longitudinal statistics — three separate technical concerns spanning both earlier chapters — increasingly interoperate as a connected toolchain rather than as isolated outputs. Across all three chapters of this series, the most consistent finding is not that any one language or single tool wins, but that the field's genuine progress lies in tools built with explicit awareness of, and connection to, the other tools around them — a pattern this closing chapter makes especially visible, since half of its six tools exist specifically because of, and in continuous dialogue with, a tool covered earlier in this survey.

## References

1. Amir A, McDonald D, Navas-Molina JA, Kopylova E, Morton JT, Zech Xu Z, Kightley EP, Thompson LR, Hyde ER, Gonzalez A, Knight R. Deblur rapidly resolves single-nucleotide community sequence patterns. *mSystems.* 2017;2(2):e00191-16. PMID: 28289731. doi:10.1128/mSystems.00191-16
2. Morton JT, Sanders J, Quinn RA, McDonald D, Gonzalez A, Vázquez-Baeza Y, Navas-Molina JA, Song SJ, Metcalf JL, Hyde ER, Lladser M, Dorrestein PC, Knight R. Balance trees reveal microbial niche differentiation. *mSystems.* 2017;2(1):e00162-16. PMID: 28144630. doi:10.1128/mSystems.00162-16
3. Shenhav L, Thompson M, Joseph TA, Briscoe L, Furman O, Bogumil D, Mizrahi I, Pe'er I, Halperin E. FEAST: fast expectation-maximization for microbial source tracking. *Nat Methods.* 2019;16(7):627–632. PMID: 31182859. doi:10.1038/s41592-019-0431-x
4. Vázquez-Baeza Y, Pirrung M, Gonzalez A, Knight R. EMPeror: a tool for visualizing high-throughput microbial community data. *GigaScience.* 2013;2(1):16. PMID: 24280061. doi:10.1186/2047-217X-2-16
5. Cantrell K, Fedarko MW, Rahman G, McDonald D, Yang Y, Zaw T, Gonzalez A, Janssen S, Estaki M, Haiminen N, Beck KL, Zhu Q, Sayyari E, Morton JT, Armstrong G, Tripathi A, Gauglitz JM, Marotz C, Matteson NL, Martino C, Sanders JG, Carrieri AP, Song SJ, Swafford AD, Dorrestein PC, Andersen KG, Parida L, Kim HC, Vázquez-Baeza Y, Knight R. EMPress enables tree-guided, interactive, and exploratory analyses of multi-omic data sets. *mSystems.* 2021;6(2):e01216-20. PMID: 33727399. doi:10.1128/mSystems.01216-20
6. Bokulich NA, Dillon MR, Zhang Y, Rideout JR, Bolyen E, Li H, Albert PS, Caporaso JG. q2-longitudinal: longitudinal and paired-sample analyses of microbiome data. *mSystems.* 2018;3(6):e00219-18. PMID: 30505944. doi:10.1128/mSystems.00219-18

## GitHub Repositories

| Package | Language | Repository |
|---|---|---|
| Deblur | Python | https://github.com/biocore/deblur |
| gneiss | Python | https://github.com/biocore/gneiss |
| FEAST | R | https://github.com/cozygene/FEAST |
| Emperor | Python/JavaScript | https://github.com/biocore/emperor |
| Empress | Python/JavaScript | https://github.com/biocore/empress |
| q2-longitudinal | Python (QIIME 2 plugin) | https://github.com/qiime2/q2-longitudinal |

*Note on citation practice: per copyright and academic integrity standards, all descriptions above are paraphrased summaries of the cited sources' findings and documentation. Readers wishing to quote these works directly should consult the original publications and software documentation linked above.*
