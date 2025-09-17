# Machine Learning Enhancement Blueprint for ATLAS

## Purpose and Scope
This document outlines how machine learning (ML) techniques can augment the Analysis Tools for Ancient and Low-depth Samples (ATLAS) toolkit across its major workflows. The guide summarizes existing capabilities, identifies ML-ready data assets, recommends high-impact enhancements, and proposes an implementation roadmap aimed at improving accuracy, robustness, and automation for ancient DNA analysis.

## Current Workflow Highlights
ATLAS organizes its command-line tasks into themed groups covering read-level processing, per-site inference, population analyses, variant call format (VCF) utilities, and simulation support. The main entry point registers these capabilities, offering a clear view of where ML modules could connect to existing pipelines.【F:atlas.cpp†L5-L138】

### Read-Level Tasks
Read-focused commands include filtering, merging, diagnostics, soft-clipping assessments, quality transformations, downsampling, and Post Mortem Damage (PMD) scoring. These routines operate on BAM inputs and emphasize curated heuristics for cleaning and characterizing ancient DNA reads.【F:atlas.cpp†L62-L71】

### Site-Level and Window Analyses
Site-oriented commands estimate sequencing and PMD error models, build genomic masks, compute genotype likelihood files (GLF), generate PSMC inputs, perform genotype calling, evaluate heterozygosity, quantify mutation load, and derive pileups and BED masks. They provide rich per-position statistics that can feed supervised or unsupervised models.【F:atlas.cpp†L75-L88】

### Population Genetics Toolkit
Population modules support GLF visualization, major/minor allele inference, genetic distance estimation, allele counting, population frequency estimation, inbreeding inference, F2 calculations, ancestral allele reconstruction, and site allele frequency (SAF) estimation, creating summary features suited to downstream ML models for demographic inference.【F:atlas.cpp†L92-L102】

### VCF Utilities and Simulations
VCF tools cover diagnostics, format conversion, pairwise comparisons, and Hardy-Weinberg equilibrium testing, all of which expose curated metadata for model-based quality control.【F:atlas.cpp†L104-L111】 Simulation commands generate BAM or VCF datasets via configurable haplotype and read simulators, which are valuable for synthetic data generation and model benchmarking.【F:atlas.cpp†L113-L117】【F:Simulations/TSimulator.h†L15-L92】

## Data Structures and Algorithmic Baseline
Understanding existing algorithms clarifies which components are candidates for ML augmentation.

### Read Filtering Architecture
The `TFilterBam` task inherits from a waiting-list BAM traverser to coordinate mate-pair handling, orphan detection, and per-read decisions before writing filtered outputs.【F:GenomeTasks/TFilterBam.h†L15-L27】 These deterministic handlers create labels (pass/fail, orphan) that can bootstrap ML classifiers for more nuanced authenticity scoring.

### Sequencing Error and Damage Estimation
`TErrorEstimator` orchestrates an expectation-maximization (EM) routine across genomic windows, maintaining per-region genotype distributions, per-site likelihood tables, recalibration maps, and PMD models while iteratively updating epsilon, rho, and psi parameters until convergence.【F:GenotypeLikelihoods/TErrorEstimator.h†L30-L95】 The explicit storage of reference bases, read data, and learned model parameters offers a structured dataset for training predictive or amortized inference models.

### Genotype Calling Stack
The genotype caller framework exposes a polymorphic hierarchy that writes richly annotated VCF records (GT, DP, GQ, AD, AP, GL, PL, GP, AB, AI) and supports multiple strategies such as random base selection, majority voting, downsampled consensus, allele presence, maximum likelihood, and Bayesian posteriors.【F:GenomeTasks/TCaller.h†L30-L318】 The modular design and prior interfaces provide clear insertion points for ML-driven scoring or consensus modules.

### Population Distance Estimation
Population distance inference relies on EM over genotype likelihood vectors to infer base frequencies, mismatch profiles, and distances, configurable with different weighting schemes and windowed analyses.【F:PopulationTools/TDistanceEstimator.h†L32-L177】 The intermediate statistics (phi, pi, per-window likelihoods) can serve as features for ML models estimating relatedness or population structure.

### Simulation Infrastructure
The simulation subsystem manages haplotype generation, read synthesis, BAM/VCF writing, and includes a roadmap item to model cross-contamination—an ideal use case for generative ML.【F:Simulations/TSimulator.h†L24-L92】 Synthetic data streams can underpin supervised training, reinforcement learning for pipeline optimization, or generative adversarial approaches to emulate authentic damage patterns.

## ML Improvement Opportunities
The following opportunities align ML techniques with ATLAS components while respecting existing APIs and data products.

### 1. Intelligent Read Filtering and Authenticity Scoring
* **Motivation:** Current filtering orchestrates mate handling and PMD-based heuristics but relies on hand-tuned thresholds.【F:atlas.cpp†L62-L71】【F:GenomeTasks/TFilterBam.h†L15-L27】
* **Approach:** Train supervised classifiers (gradient boosting, convolutional or transformer models on base-quality sequences) using labels from `filterBAM`, `PMDS`, and downstream authenticity assessments. Incorporate probabilistic outputs to complement existing pass/fail logic, e.g., by weighting reads instead of discarding them.
* **Integration:** Wrap the classifier inside `TFilterBam` or as a pre-filter that annotates reads via auxiliary BAM tags. Provide fallback paths to legacy heuristics for compatibility.
* **Data:** Use simulation outputs (`simulate`) with controlled contamination and damage to generate labeled corpora, supplemented by empirical truth sets when available.【F:atlas.cpp†L113-L117】【F:Simulations/TSimulator.h†L24-L92】

### 2. Sequencing Error and Damage Modeling with Learned Surrogates
* **Motivation:** The EM loop in `TErrorEstimator` iteratively adjusts epsilon/rho/psi parameters but can be slow on deep data or complex models.【F:GenotypeLikelihoods/TErrorEstimator.h†L30-L95】
* **Approach:** Develop neural amortized inference models that map summary statistics (per-site likelihood tables, quality histograms) to parameter estimates, accelerating convergence. Alternatively, use variational autoencoders to model joint distributions of sequencing errors and PMD patterns, providing uncertainty-aware predictions.
* **Integration:** Embed ML predictions as initializations for the EM routine, or expose a new task flag that substitutes EM updates with learned outputs when confidence is high. Persist learned parameters to ensure traceability.
* **Evaluation:** Compare log-likelihood trajectories and downstream heterozygosity estimates with and without ML initialization to quantify gains.

### 3. Adaptive Genotype Calling and Prior Calibration
* **Motivation:** Genotype calling currently supports deterministic, likelihood-based, and Bayesian variants with configurable priors and output annotations.【F:GenomeTasks/TCaller.h†L30-L318】
* **Approach:**
  * Train discriminative models that consume genotype likelihood vectors, read depth, and allele balance metrics to predict genotypes or posterior recalibrations, similar to modern deep callers.
  * Use meta-learning to adapt priors (`TGenotypePrior`) based on sample-specific characteristics inferred from GLF or population summaries.【F:atlas.cpp†L75-L102】
  * Introduce calibration networks that adjust quality scores (GQ, GP) to improve concordance with truth sets.
* **Integration:** Expose ML-enhanced callers as new subclasses of `TCaller`, leveraging existing VCF field writers to keep output consistent. Provide CLI flags for selecting the ML caller while maintaining compatibility with existing pipelines.
* **Data:** Generate high-confidence genotype truth via simulations or orthogonal datasets (e.g., high-coverage modern genomes) processed through the ATLAS pipeline to create aligned training inputs.

### 4. Population-Level Inference and Demographic Modeling
* **Motivation:** EM-based distance estimation and related statistics provide interpretable metrics but may not capture complex population structure or admixture patterns.【F:PopulationTools/TDistanceEstimator.h†L32-L177】
* **Approach:**
  * Train graph-based or embedding models on GLF-derived features to infer relatedness, population assignment, or admixture proportions.
  * Employ unsupervised clustering or variational inference to summarize SAF, allele counts, and F2 metrics into latent representations for downstream tasks.
  * Use ML regression to predict demographic parameters directly from per-window summaries, calibrated against coalescent simulations.
* **Integration:** Add optional ML-backed reporting modules that augment or replace existing outputs in `geneticDist`, `inbreeding`, or `calculateF2`, ensuring that classic metrics remain available for validation.

### 5. Automated Quality Control for VCF Diagnostics
* **Motivation:** VCF utilities currently report deterministic diagnostics, conversions, and comparisons.【F:atlas.cpp†L104-L111】
* **Approach:** Apply anomaly detection (autoencoders, isolation forests) to VCF-derived metrics (e.g., depth distributions, Hardy-Weinberg p-values) to flag outliers or suspect batches. Train models to predict mismatch patterns between VCF files for rapid regression testing.
* **Integration:** Extend `VCFDiagnostics` to output ML-derived quality flags or risk scores alongside existing summaries, enabling automated triage in pipelines.

### 6. Generative Simulation Enhancements
* **Motivation:** Simulations drive method testing and include a roadmap item for modeling cross-contamination, a phenomenon well-suited to ML-based generative models.【F:Simulations/TSimulator.h†L24-L92】
* **Approach:** Develop generative models (e.g., GANs, diffusion models) conditioned on reference haplotypes and sequencing settings to produce realistic read data, including contamination signatures. Use reinforcement learning to tune simulator parameters for matching empirical summary statistics.
* **Integration:** Offer ML-driven simulation modes within `simulate`, optionally writing parameter traces for reproducibility. The synthetic data can feed back into training for filtering, error modeling, and genotype calling modules.

## Implementation Strategy
1. **Pilot Data Collection:** Instrument existing tasks to emit structured training data (e.g., JSON summaries, parquet tables) alongside standard outputs. Focus on read-level features from `filterBAM`/`PMDS`, per-site likelihood arrays from `estimateErrors`, and genotype likelihood vectors from `call`/`GLF` tasks.【F:atlas.cpp†L62-L88】
2. **Model Prototyping:** Use notebooks or scripted experiments to benchmark candidate models on simulated and empirical datasets. Leverage the simulation framework to generate diverse training corpora with known truth.【F:Simulations/TSimulator.h†L32-L88】
3. **API Integration:** Implement ML wrappers as new task variants or optional flags to existing tasks (e.g., `--ml-filter`, `--ml-caller`). Reuse the polymorphic design of callers and estimators to keep code modular.【F:GenomeTasks/TCaller.h†L30-L318】【F:GenotypeLikelihoods/TErrorEstimator.h†L30-L95】
4. **Evaluation Harness:** Add regression tests comparing ML-enhanced outputs against baselines using known truth sets or established benchmarks. Track metrics such as concordance, calibration error, runtime, and memory usage.
5. **Deployment and Monitoring:** Ship pre-trained models with versioned metadata, and log predictions, confidence intervals, and fallbacks to detect drift. Use debug tasks (e.g., `thetaLLSurface`, `averageDepth`) as sanity checks when ML modules are active.【F:atlas.cpp†L126-L136】

## Data Management and Model Training Considerations
- **Label Generation:** Use deterministic outputs (pass/fail decisions, EM parameter estimates, genotype calls) as initial labels, then refine via expert review or cross-validation against orthogonal datasets.
- **Feature Engineering:** Combine raw read signals (qualities, patterns) with derived summaries (GLF likelihoods, allele balances, window-based distances) to provide ML models with both low-level and high-level context.【F:GenomeTasks/TCaller.h†L30-L318】【F:PopulationTools/TDistanceEstimator.h†L32-L177】
- **Model Versioning:** Track training datasets, hyperparameters, and evaluation scores alongside the models to ensure reproducibility and regulatory compliance for ancient DNA studies.
- **Compute Requirements:** Anticipate GPU resources for deep sequence models and ensure fallbacks for CPU-only environments.

## Evaluation and Monitoring Framework
- **Offline Benchmarks:** Compare ML outputs against simulated truth and curated empirical truth sets to quantify precision, recall, calibration, and runtime speedups.
- **Online Validation:** When ML modules run in production, periodically sample outputs for manual review, cross-checking with legacy pipelines to detect drift.
- **Ablation Studies:** Evaluate the contribution of ML initializations versus full replacements for existing EM or heuristic steps to justify integration complexity.

## Roadmap and Prioritization
1. **Short-Term (0–3 months):**
   - Collect training data during routine pipeline runs and finalize problem formulations for read filtering and error estimation.
   - Prototype lightweight classifiers/regressors that can be embedded without major dependency changes.
2. **Mid-Term (3–9 months):**
   - Integrate ML-augmented error estimation and genotype calling modules with optional CLI flags.
   - Develop automated QC scoring for VCF diagnostics and population summaries.
3. **Long-Term (9+ months):**
   - Launch ML-driven simulation modes with contamination modeling.
   - Explore end-to-end models that jointly optimize read authenticity, error calibration, and genotype calling, potentially leveraging multi-task learning.

## Risks and Mitigations
- **Data Scarcity:** Ancient DNA truth sets are limited; mitigate by combining simulations with semi-supervised learning and uncertainty-aware models.【F:Simulations/TSimulator.h†L24-L92】
- **Model Drift:** Monitor prediction quality and maintain easy fallbacks to deterministic routines to avoid regressions.【F:atlas.cpp†L62-L88】
- **Reproducibility:** Store model artifacts alongside ATLAS releases and document training procedures to meet transparency requirements in paleogenomics.

By aligning ML enhancements with the existing modular architecture, ATLAS can progressively incorporate data-driven intelligence while preserving the reproducibility and interpretability expected by the ancient DNA community.
