# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing pipeline for untargeted LC-MS metabolomics. Converts a
# MSDIAL alignment export or an XCMS feature table + sample sheet to notame
# format, applies pre-correction feature filters, and runs one or more
# drift/batch correction methods in parallel. QC metrics are computed per
# method for comparison.
#
# Daniel Brunnsåker, 2026-03-20
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(notame)
  library(notameViz)
  library(notameStats)
  library(openxlsx)
  library(doParallel)
  library(yaml)
  library(jsonlite)
})

source("R/config.R")
source("R/method_spec.R")
source("R/preflight.R")
source("R/notame_format.R")
source("R/msdial_to_notame.R")
source("R/xcms_to_notame.R")
source("R/qc_metrics.R")
source("R/drift_correction.R")
source("R/correction_methods.R")
source("R/serrf.R")
source("R/run_log.R")
source("R/report.R")

# Optional YAML config file — supplies default parameter values. Env vars
# still take precedence over the config file (see get_env() below), so
# Docker one-off overrides keep working unchanged.
config_file <- Sys.getenv("CONFIG_FILE", unset = "")
config <- if (config_file != "") load_config(config_file) else list()

get_env <- function(var, default) {
  val <- Sys.getenv(var, unset = NA)
  if (!is.na(val) && val != "") return(val)
  if (!is.null(config[[var]]) && nchar(config[[var]]) > 0) return(config[[var]])
  default
}

n_cores_env <- get_env("N_CORES", "")
registerDoParallel(cores = if (n_cores_env == "") parallel::detectCores() - 1 else as.integer(n_cores_env))

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# All parameters can be overridden via environment variables, or supplied in
# a YAML config file (CONFIG_FILE) — env vars take precedence over the
# config file when both are set.
# ─────────────────────────────────────────────────────────────────────────────

if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("
Usage: Rscript notame-workflow.r [--help]

Environment variables (required variables are marked; all others are optional with defaults shown).
Any of these may instead be supplied via a YAML file pointed to by CONFIG_FILE — env vars win
if both are set for the same parameter.

  CONFIG_FILE            Path to an optional YAML config file supplying default values for any
                        parameter below, plus an optional 'sample_type_rules' list overriding
                        the QC/ltQC/Blank/Wash/etc. filename-keyword classification. See
                        config.example.yaml.
                        Default: (none)

  IN_XLSX               Path to the MSDIAL alignment export (.xlsx). Required, unless
                        IN_FEATURE_TABLE + IN_SAMPLE_SHEET are set instead.

  IN_FEATURE_TABLE      Path to an XCMS-based feature table (.csv) — alternative to IN_XLSX.
                        Must be set together with IN_SAMPLE_SHEET.

  IN_SAMPLE_SHEET       Path to the XCMS pipeline's sample sheet (.xlsx) — alternative to
                        IN_XLSX. Must be set together with IN_FEATURE_TABLE.

  PROJECT_FOLDER        Root output directory (intermediates/ and output/ written here). Required.

  FORCE_RECONVERT       Force re-running the conversion step even if a cached
                        notame-formatted file from a previous run is found.
                        Default: FALSE
                        Values:  TRUE | FALSE

  COLUMN                Chromatographic column type. Required. Used together with POLARITY
                        to namespace output folders (e.g. RP_POS, HILIC_NEG).
                        Examples: RP | HILIC

  POLARITY              Ionisation polarity. Required.
                        Values:  POS | NEG

  CORRECTION_METHODS    Comma-separated list of correction 'recipes' to run. Each recipe is a
                        'drift:basis:batch' triple (e.g. 'loess:hybrid:combat'), not a single
                        opaque method name -- see the three axes below. Each recipe is saved to
                        its own output subfolder (colons replaced with dashes in the folder name).
                        Default: none:none:none,notame_spline:qc:ruv_s

                        DRIFT_METHOD (1st field) -- how within-batch drift is removed, before any
                        batch step:
                          none           no drift correction
                          loess          LOESS drift correction (LOESS_QC_SPAN / LOESS_SAMPLE_SPAN)
                          huber          Huber robust regression (HUBER_QC_K / HUBER_SAMPLE_K)
                          auto           auto-selected via cross-validation from LOESS spans
                                         (AUTO_LOESS_SPANS/AUTO_SAMPLE_LOESS_SPANS), Huber k's
                                         (AUTO_HUBER_KS/AUTO_SAMPLE_HUBER_KS), and a flat/no-op
                                         baseline -- ONE method is chosen for the whole run (per
                                         basis), not a different one per batch
                          notame_spline  notame's own per-feature cubic smoothing spline
                                         (notame::correct_drift()); only supports basis=qc

                        BASIS (2nd field) -- which data the drift method is fit against. Must be
                        'none' if and only if DRIFT_METHOD is 'none':
                          qc       fit only on QC samples, only for batches with
                                   >= DRIFT_MIN_QC_PER_BATCH QC; batches below that are left
                                   uncorrected, no fallback
                          samples  fit only on biological samples, for every batch, regardless of
                                   QC availability
                          hybrid   QC-based if the batch has enough QC (DRIFT_MIN_QC_PER_BATCH);
                                   otherwise a samples-based trial validated against ltQC/Sample
                                   D-ratio (DRIFT_MIN_LTQC_VALIDATE, DRIFT_HYBRID_VALIDATE), kept
                                   only if it helps; otherwise left uncorrected

                        BATCH_METHOD (3rd field) -- how between-batch differences are removed:
                          none            no between-batch step
                          combat          ComBat (COMBAT_MEAN_ONLY, COMBAT_PAR_PRIOR)
                          sva             Surrogate Variable Analysis (SVA_N_SV)
                          limma           limma::removeBatchEffect()
                          feature_median  per-feature median ratio normalisation
                          global_median   single global median ratio normalisation
                          ruv_s           notame's RUV-S (RUV_K), QC-anchored
                          cordbat         CordBat GGM-based alignment (CORDBAT_REF_BATCH)
                          batchcorr       batchCorr cluster spline + normalizeBatches() -- couples
                                          drift+batch internally, requires drift_method=none
                          waveica         WaveICA2.0 (WAVEICA_*) -- couples drift+batch, requires
                                          drift_method=none
                          waveica_v1      original WaveICA (WAVEICA_V1_*) -- couples drift+batch,
                                          requires drift_method=none
                          pmp_qcrsc       pmp QC-RSC spline -- couples drift+batch, requires
                                          drift_method=none
                          serrf           SERRF random-forest correction (SERRF_NUM) -- couples
                                          drift+batch, requires drift_method=none

                        Examples: 'loess:hybrid:combat', 'huber:qc:feature_median',
                        'none:none:cordbat', 'auto:hybrid:sva'

  QC_DETECTION_LIMIT    Min fraction of QC samples a feature must be detected in
                        Default: 0.60

  SAMPLE_DETECTION_LIMIT  Min fraction of biological samples a feature must be detected in
                        Default: 0.20

  BLANK_RATIO           Remove features where mean(Sample) <= BLANK_RATIO * mean(SolvBlank).
                        Set to 'none' to disable.
                        Default: none

  LOW_INT_FILTER_FRAC   Data-driven low-intensity filter. Removes features whose p80
                        intensity is below this fraction of the mean p80 across all features.
                        Ignored if LOW_INT_FILTER is set.
                        Default: 0.10

  LOW_INT_FILTER        Absolute low-intensity threshold (overrides LOW_INT_FILTER_FRAC).
                        Default: (disabled)

  LOW_INT_PERCENTILE    Percentile used for the low-intensity filter (0-1).
                        Default: 0.8

  MIN_QC_SAMPLE_DETECTION  Minimum fraction of features that must be detected in a QC or ltQC
                        sample for it to be used as a reference. Samples below this threshold
                        are removed before processing (empty injections, failed runs).
                        Default: 0.50

  QC_OUTLIER_MAD_K      Multivariate outlier check for QC/ltQC samples (PCA-distance-based),
                        run separately for QC and ltQC, per batch, after
                        MIN_QC_SAMPLE_DETECTION's removal (so badly-detected samples don't
                        distort the reference centroid this compares everyone else against).
                        Complements MIN_QC_SAMPLE_DETECTION: that catches a sample that failed
                        to detect most features; this catches one that detects fine but has an
                        anomalous intensity profile (contamination, carryover, a
                        degrading/recovering column) -- a failure mode detection-rate is
                        structurally blind to. A sample is flagged if its Euclidean distance
                        (in PCA score space, top 5 PCs, unit-variance scaled) from its batch
                        group's median score exceeds median(distances) + K * mad(distances).
                        Lower K = more aggressive (flags more samples); higher K = more
                        conservative. Batches with fewer than QC_OUTLIER_MIN_N samples of a
                        group are skipped (not enough points for a meaningful check). Set to 0
                        to disable.
                        Default: 5

  QC_OUTLIER_MIN_N      Minimum samples of a group (QC or ltQC) a batch must have for
                        QC_OUTLIER_MAD_K's check to run there at all; batches below this are
                        skipped for that group. 3 is the practical floor, not lower -- PCA
                        correctly caps at 2 dimensions with 3 points and MAD of 3 distances is
                        still well-defined (weak power, not meaningless), but with only 2
                        points there's no third reference to judge 'typical spread' against, so
                        neither could ever be meaningfully called the outlier. Default is 3
                        rather than a more conventional floor like 4-5 specifically because
                        ltQC groups are often exactly 3 per batch by design (a long-term QC
                        injected a few times per batch, not throughout it) -- a higher default
                        would silently skip ltQC's check in every batch for a setup like that.
                        Default: 3

  MIN_BATCH_DETECTION   Minimum number of detections a feature must have in every batch.
                        Features absent from any entire batch are removed — they have no
                        real measurements there and would be all-imputed placeholders.
                        Set to 0 to disable.
                        Default: 1

  MIN_BATCH_DETECTION_FRAC  Minimum fraction (0-1) of each batch's samples a feature must be
                        detected in. Unlike MIN_BATCH_DETECTION (an absolute count), this
                        scales with batch size, so batches of very different sizes get a
                        consistent relative bar instead of a fixed count that's stringent for
                        a small batch and lax for a large one. Applied in addition to
                        MIN_BATCH_DETECTION, not instead of it -- a feature must pass both.
                        Relevant for batch_method=waveica/waveica_v1 (and other methods whose
                        LoD/2-imputed placeholders participate directly in fitting the
                        correction): a batch with disproportionately more missingness than
                        others risks its placeholder pattern getting mistaken for real signal.
                        Set to 0 to disable.
                        Default: 0 (disabled)

  QC_RSD_FILTER         Pre-correction QC-RSD threshold (robust: MAD/median). Features are
                        removed if they fail the threshold in >= 50% of evaluable batches.
                        Set to 'none' to disable.
                        Default: none (disabled)

  RSD_THRESHOLD         RSD threshold used for post-correction output filtering.
                        Controls both the global (_rsdXX) and per-batch (_batchrsdXX) outputs.
                        Default: 0.30

  RUV_K                 Number of unwanted variation factors for batch_method=ruv_s (notame's RUV-S).
                        Default: 3

  SERRF_NUM             Number of correlated features used as RF predictors per feature model
                        in batch_method=serrf (SERRF). Lower = less overfitting risk; recommended
                        ≤ floor(QC_n/2). Automatically capped at floor(min_QC_per_batch / 2) at runtime.
                        Default: 5

  LOESS_QC_SPAN         LOESS smoothing span for drift_method=loess's QC-based fit (basis=qc, or
                        basis=hybrid's QC branch). Higher = smoother, more conservative correction.
                        Default: 0.75

  LOESS_SAMPLE_SPAN     LOESS smoothing span for drift_method=loess's samples-based fit
                        (basis=samples, or basis=hybrid's samples branch), fit on biological
                        samples instead of QC. Wider than LOESS_QC_SPAN by default since each point
                        is a unique biological measurement, not a technical replicate — a tighter
                        span risks fitting individual-sample noise as drift.
                        Default: 0.9

  DRIFT_SAMPLE_MIN_OBS  Minimum finite sample observations per feature required to attempt a
                        samples-based drift fit (basis=samples, or basis=hybrid's samples branch;
                        applies to both drift_method=loess and drift_method=huber). Deliberately
                        higher than the QC-based fit's threshold of 4, since sample points are far
                        noisier. Features below this are left uncorrected for that batch.
                        Default: 10

  DRIFT_MIN_QC_PER_BATCH  Minimum QC samples a batch must have to use the QC-based fit under
                        basis=hybrid (applies to drift_method=loess and drift_method=huber alike).
                        Batches meeting this threshold get the QC-based fit (LOESS_QC_SPAN /
                        HUBER_QC_K). Batches below it try the samples-based fit
                        (LOESS_SAMPLE_SPAN/HUBER_SAMPLE_K, DRIFT_SAMPLE_MIN_OBS) as a trial, kept
                        only if it improves the ltQC/Sample D-ratio — see DRIFT_MIN_LTQC_VALIDATE.
                        QC-based fitting is preferred whenever there's enough QC to support it,
                        since it doesn't risk removing real biological signal along with drift.
                        Also used by basis=qc as a hard cutoff (no samples-based fallback there —
                        batches below this threshold are simply left uncorrected).
                        Default: 4

  DRIFT_MIN_LTQC_VALIDATE  Minimum ltQC samples a batch must have to validate the samples-based
                        trial correction under basis=hybrid (used only for batches below
                        DRIFT_MIN_QC_PER_BATCH). The trial correction is compared against the
                        batch's uncorrected values by ltQC/Sample D-ratio (MAD(ltQC)/MAD(Sample) —
                        lower is better); it is kept if the D-ratio improves, discarded otherwise.
                        D-ratio (not raw ltQC RSD) is used because any real drift correction
                        shrinks measured sample variance somewhat, so a test that just required
                        sample variance not to shrink would reject working corrections too —
                        D-ratio instead only credits a disproportionate improvement in ltQC
                        relative to Sample. Batches with fewer ltQC than this have no way to
                        validate the trial and are left uncorrected.
                        Default: 3

  DRIFT_HYBRID_VALIDATE  Set to FALSE to skip the ltQC validation above entirely under
                        basis=hybrid: any batch below DRIFT_MIN_QC_PER_BATCH then always gets the
                        samples-based correction, regardless of ltQC availability or what it
                        shows. This reintroduces the risk the validation step exists to catch (the
                        trial can look fine on ltQC while still compressing real biological
                        signal) — use deliberately, not as a default.
                        Default: TRUE

  AUTO_LOESS_SPANS      Comma-separated LOESS spans drift_method=auto evaluates as candidates for
                        the QC-based selection (leave-one-out CV on QC; used when basis=qc, or
                        basis=hybrid's QC-tier — set together with AUTO_HUBER_KS as one QC-based
                        candidate pool; see AUTO_SAMPLE_LOESS_SPANS for the separate samples-only
                        pool).
                        Default: 0.5,0.75,0.9

  AUTO_HUBER_KS         Comma-separated Huber regression k values (MASS::rlm, psi.huber)
                        drift_method=auto evaluates as candidates for the QC-based selection.
                        Lower k = more robust to outlier QC points but less statistically
                        efficient; 1.345 is MASS::rlm's own default (~95% efficiency under
                        Gaussian errors).
                        Default: 1.0,1.345,2.0

  AUTO_SAMPLE_LOESS_SPANS  Comma-separated LOESS spans drift_method=auto evaluates for the
                        samples-based selection (basis=samples, or basis=hybrid's samples-tier;
                        fit on biological samples, validated against ltQC). Separate range from
                        AUTO_LOESS_SPANS — see the LOESS_SAMPLE_SPAN discussion above for why a
                        QC-free fit needs a different span range than a QC-anchored one.
                        Default: 0.3,0.6,0.9

  AUTO_SAMPLE_HUBER_KS  Comma-separated Huber k values for drift_method=auto's samples-based
                        candidate pool.
                        Default: 1.0,1.345,2.0

  AUTO_MIN_QC_PER_BATCH  Same role as DRIFT_MIN_QC_PER_BATCH, for drift_method=auto: minimum QC
                        samples a batch needs to contribute to (and, under basis=hybrid, receive)
                        the QC-based (leave-one-out CV) candidate selection instead of the
                        samples-only (ltQC-validated) one.
                        Default: 4

  AUTO_MIN_LTQC_VALIDATE  Same role as DRIFT_MIN_LTQC_VALIDATE, for drift_method=auto: minimum
                        ltQC samples a batch needs to contribute to the samples-based candidate
                        selection at all. Under basis=samples, once a winning candidate is chosen
                        it is applied to every batch regardless of this threshold — it only
                        affects which batches help pick the winner. Under basis=hybrid (and
                        basis=qc, where this pool isn't used), batches with fewer are left
                        uncorrected.
                        Default: 3

  AUTO_MIN_CV_OBS       Minimum finite training observations (QC, or samples for the ltQC-validated
                        pool) a feature needs before drift_method=auto attempts to fit any
                        candidate for it. Features below this are left uncorrected for that batch.
                        Default: 4

  HUBER_QC_K            Huber regression tuning constant (MASS::rlm, psi.huber) for
                        drift_method=huber's QC-based fit (basis=qc, or basis=hybrid's QC branch).
                        Fixed (shared across every feature), not auto-searched — see HUBER_QC_CV_KS
                        below for per-feature CV selection within basis=qc/hybrid, or
                        AUTO_HUBER_KS/drift_method=auto for a dataset-wide CV-chosen k instead.
                        Lower = more robust to outlier QC points but less statistically efficient;
                        1.345 is MASS::rlm's own default (~95% efficiency under Gaussian errors).
                        Default: 1.345

  HUBER_SAMPLE_K        Huber tuning constant for drift_method=huber's samples-based fit
                        (basis=samples, or basis=hybrid's samples branch). Separate setting from
                        HUBER_QC_K, same reasoning as LOESS_SAMPLE_SPAN vs LOESS_QC_SPAN — a
                        samples-only fit may warrant a different robustness/efficiency trade-off
                        than a QC-anchored one.
                        Default: 1.345

  LOESS_QC_CV_SPANS     Comma-separated LOESS span candidates for drift_method=loess's QC-based
                        step (basis=qc, or basis=hybrid's QC branch). When set, the span is chosen
                        per FEATURE via leave-one-out CV on QC (mirroring how notame::correct_drift()
                        's smooth.spline() step auto-selects its own smoothing parameter per feature,
                        rather than sharing one span dataset-wide) instead of using the fixed
                        LOESS_QC_SPAN for every feature. Safe to do per-feature specifically because
                        this only ever evaluates against QC (pure technical replicates, nothing
                        biological to overfit) — deliberately not offered for basis=samples, where
                        a flexible per-feature fit could overfit real biological variation.
                        Leave empty (default) to keep using the fixed LOESS_QC_SPAN for every feature.
                        Default: (disabled — uses LOESS_QC_SPAN)

  HUBER_QC_CV_KS        Same idea as LOESS_QC_CV_SPANS, for drift_method=huber's QC-based step:
                        per-feature CV selection from this comma-separated k grid instead of the
                        fixed HUBER_QC_K, when set.
                        Default: (disabled — uses HUBER_QC_K)

  SVA_N_SV              Number of surrogate variables SVA (sva package) estimates for
                        batch_method=sva's correction step. These are latent factors representing
                        systematic structure in the data not already explained by known Batch;
                        regressed out together with Batch via limma::removeBatchEffect(). Leave
                        unset to auto-estimate via sva::num.sv(..., method = 'be')
                        (Buja-Eyuboglu permutation test); set to 0 to disable SVA and correct for
                        known Batch only.
                        Default: (auto)

  CORDBAT_REF_BATCH     Reference batch ID for batch_method=cordbat. All other batches are
                        corrected onto this batch. Leave unset to auto-select the batch with the
                        lowest median feature RSD.
                        Default: (auto)

  WAVEICA_ALPHA         Comma-separated trade-off value(s) (0-1) for WaveICA2.0's internal ICA step
                        (unbiased_stICA()): 0 = spatial ICA, 1 = temporal ICA, balancing
                        independence across samples vs. across features in the decomposition. NOT
                        a significance/flagging threshold and does not itself control how many
                        components get removed -- see WAVEICA_CUTOFF for that. A single value keeps
                        that value fixed (one WaveICA_2.0() call); more than one value overall
                        across WAVEICA_ALPHA/WAVEICA_CUTOFF/WAVEICA_K triggers a search over the
                        full cross-product grid -- see WAVEICA_EVAL_GROUP below.
                        Default: 0.05

  WAVEICA_CUTOFF        Comma-separated threshold(s) (0-1): the minimum R² (against injection
                        order, via a GAM fit) an individual ICA component must reach to be treated
                        as injection-order-technical and subtracted out. This is the actual
                        'how aggressive' dial -- every wavelet decomposition level is always
                        ICA-decomposed regardless of this setting; Cutoff decides which of the
                        resulting *components* (not which levels) get removed. Lower = more
                        components qualify as technical (more removed, more aggressive); higher =
                        fewer qualify (more conservative). Same single-value-fixed /
                        multi-value-searched convention as WAVEICA_ALPHA.
                        Default: 0.10

  WAVEICA_K             Comma-separated number(s) of independent components batch_method=waveica
                        (WaveICA2.0) decomposes the data into. Each entry is either a number or the
                        literal 'auto' (2 x number of batches, resolved per run -- quite small for
                        typical feature counts, since a coarse decomposition gives ICA less room to
                        isolate narrow technical components from broad biological ones). Try
                        including a larger value (e.g. 10-20) if correction looks too aggressive.
                        Same single-value-fixed / multi-value-searched convention as WAVEICA_ALPHA.
                        Default: auto (= 2 x n_batches)

  WAVEICA_WF            Wavelet family used for the decomposition step (passed to WaveICA2.0's
                        wf argument, e.g. haar or a Daubechies family recognised by the
                        underlying wavelet package). Not part of the search grid -- kept fixed,
                        since it's a categorical choice rather than a more/less-aggressive dial and
                        including it would multiply the grid size for a dimension with no clear
                        prior on which value helps.
                        Default: haar

  WAVEICA_EVAL_GROUP    Which group -- 'ltQC' or 'QC' -- the WAVEICA_ALPHA/WAVEICA_CUTOFF/WAVEICA_K
                        search (when triggered) evaluates candidates against, via D-ratio vs.
                        Sample (MAD(group)/MAD(Sample), lower is better; the same candidate that
                        wins by D-ratio also gets its PCA-space distance ratio and PERMANOVA
                        R²(Batch) printed alongside every other candidate's, as an independent
                        cross-check -- WaveICA2.0 is an ICA-based method operating jointly across
                        features, so a purely per-feature metric like D-ratio alone could miss
                        damage to that joint structure). WaveICA2.0 never fits on QC or ltQC -- it
                        corrects using only injection order -- so either is a genuine held-out
                        reference regardless of which is chosen; 'QC' is worth trying if it has
                        more samples than ltQC in your data. Candidates are evaluated in parallel
                        via this pipeline's existing N_CORES/foreach setup, with each worker pinning
                        mc.cores to 1 for its own WaveICA_2.0() call to avoid nested-parallelism
                        oversubscription against that function's own internal use of
                        parallel::mclapply().
                        Values: ltQC | QC
                        Default: ltQC

                        Note: WAVEICA_ALPHA/WAVEICA_CUTOFF/WAVEICA_K/WAVEICA_WF are exposed as-is
                        from the WaveICA2.0 package (github.com/dengkuistat/WaveICA_2.0); the
                        direction-of-effect guidance above is based on the published method rather
                        than inspection of this specific package version's source, so confirm
                        empirically on your data.

  WAVEICA_V1_WF         Wavelet family for batch_method=waveica_v1 (the original WaveICA, not
                        WaveICA2.0). Same meaning as WAVEICA_WF, separate setting since the two
                        methods are independent packages. Not part of the search grid, same as
                        WAVEICA_WF.
                        Default: haar

  WAVEICA_V1_K          Comma-separated number(s) of components batch_method=waveica_v1's ICA step
                        decomposes into. Each entry is either a number or the literal 'auto' (2 x
                        number of batches -- not something the original package defines itself,
                        kept purely for parity with WAVEICA_K's convention). A single value keeps
                        it fixed (one WaveICA() call); more than one value overall across
                        WAVEICA_V1_ALPHA/WAVEICA_V1_T/WAVEICA_V1_K triggers a search over the full
                        cross-product grid, same mechanism as WAVEICA_ALPHA/WAVEICA_CUTOFF/WAVEICA_K
                        -- see WAVEICA_V1_EVAL_GROUP below.
                        Default: 20

  WAVEICA_V1_T          Comma-separated threshold(s) (0-1) for considering an ICA component
                        associated with batch in batch_method=waveica_v1. Unlike WAVEICA_CUTOFF
                        (WaveICA2.0's GAM-R²-against-injection-order test), this tests each
                        component's p-value directly against the real batch labels (a component is
                        removed if its p-value is below this threshold), since waveica_v1 uses
                        actual batch identity rather than injection order as a proxy for it. Higher
                        = a looser significance bar, so more components qualify as batch-associated
                        (more removed, more aggressive); lower = stricter, fewer qualify (more
                        conservative) -- opposite direction from WAVEICA_CUTOFF, which thresholds an
                        R² a component must exceed rather than a p-value it must fall under. Same
                        single-value-fixed / multi-value-searched convention as WAVEICA_V1_K.
                        Default: 0.05

  WAVEICA_V1_T2         Threshold (0-1) for considering an ICA component associated with a
                        biological comparison group in batch_method=waveica_v1. Not currently used
                        in practice -- this pipeline has no biological-group column to supply
                        waveica_v1's optional `group` argument, so that protection is inactive
                        regardless of this setting. Kept for parity with the package's own
                        parameters; not part of the search grid since it has no effect here.
                        Default: 0.05

  WAVEICA_V1_ALPHA      Comma-separated trade-off value(s) (0-1) between sample-wise and
                        variable-wise independence in batch_method=waveica_v1's ICA step. The same
                        KIND of parameter as WAVEICA_ALPHA (WaveICA2.0) -- both are ICA
                        spatial/temporal independence trade-offs, not a significance/flagging
                        threshold -- kept as a separate setting since they're independent packages
                        with separately-tuned defaults, not because the concept differs. Same
                        single-value-fixed / multi-value-searched convention as WAVEICA_V1_K.
                        Default: 0

  WAVEICA_V1_EVAL_GROUP Which group -- 'ltQC' or 'QC' -- the WAVEICA_V1_ALPHA/WAVEICA_V1_T/
                        WAVEICA_V1_K search evaluates candidates against (D-ratio vs. Sample). Same
                        mechanism as WAVEICA_EVAL_GROUP.
                        Default: ltQC

                        Note: waveica_v1 uses real batch labels directly rather than injection
                        order as a proxy for batch structure, which may make it a better fit
                        when batch labels are known and reliable. Defaults above are the
                        package's own defaults, not tuned for this pipeline specifically.

  COMBAT_MEAN_ONLY      Whether ComBat (batch_method=combat, with any drift_method) adjusts only
                        each feature's per-batch mean (TRUE) or also forces every batch's variance
                        to match a common value (FALSE, ComBat's own default). Forcing variance
                        equal across batches is the usual cause of PCA looking artificially
                        'flattened' after correction, if batches genuinely differ in spread (e.g.
                        different biological composition). 'auto' tries both and keeps whichever
                        gives the better ltQC/Sample D-ratio (see combat_correct() in
                        R/correction_methods.R) -- set TRUE or FALSE directly to skip the search
                        and force a specific behaviour.
                        Values: auto, TRUE, FALSE
                        Default: auto

  COMBAT_PAR_PRIOR      Whether ComBat estimates its empirical Bayes prior parametrically
                        (TRUE, assumes a Normal/Inverse-Gamma shape for batch effects -- faster)
                        or non-parametrically (FALSE, a more flexible density estimate -- slower,
                        more robust if batch effects are non-Gaussian). 'auto' tries both and keeps
                        whichever gives the better ltQC/Sample D-ratio, same mechanism as
                        COMBAT_MEAN_ONLY.
                        Values: auto, TRUE, FALSE
                        Default: auto

  NORMALIZATION         Post-correction normalisation method. Uses pooled QC samples as
                        reference when available, otherwise median of biological samples.
                        Default: none
                        Values:  none | pqn

  SAVE_PRE_CORRECTION_PLOTS  Save QC plots before correction (slow on large datasets).
                        Default: TRUE
                        Values:  TRUE | FALSE

  N_CORES               Number of CPU cores for parallelisation.
                        Default: all available cores minus one

")
  quit(status = 0)
}

project_folder <- get_env("PROJECT_FOLDER", "")
if (project_folder == "") stop("PROJECT_FOLDER is required. Set it via env var or CONFIG_FILE.")

# Input: an MSDIAL alignment export (IN_XLSX), or an XCMS-based feature
# table + sample sheet (IN_FEATURE_TABLE + IN_SAMPLE_SHEET) — mutually
# exclusive, auto-detected from which are set.
in_xlsx          <- get_env("IN_XLSX", "")
in_feature_table <- get_env("IN_FEATURE_TABLE", "")
in_sample_sheet  <- get_env("IN_SAMPLE_SHEET", "")

has_msdial <- in_xlsx != ""
has_xcms   <- in_feature_table != "" || in_sample_sheet != ""

if (has_msdial && has_xcms) {
  stop("Set either IN_XLSX (MSDIAL) or IN_FEATURE_TABLE + IN_SAMPLE_SHEET (XCMS), not both.")
} else if (has_xcms) {
  if (in_feature_table == "" || in_sample_sheet == "")
    stop("XCMS input requires both IN_FEATURE_TABLE and IN_SAMPLE_SHEET to be set.")
  input_mode <- "xcms"
} else if (has_msdial) {
  input_mode <- "msdial"
} else {
  stop("An input is required: set IN_XLSX (MSDIAL) or IN_FEATURE_TABLE + IN_SAMPLE_SHEET (XCMS), via env var or CONFIG_FILE.")
}

polarity <- get_env("POLARITY", "")
if (polarity == "") stop("POLARITY is required. Set it to 'POS' or 'NEG' (env var or CONFIG_FILE).")

column <- get_env("COLUMN", "")
if (column == "") stop("COLUMN is required. Set it to the chromatographic column type, e.g. 'RP' or 'HILIC' (env var or CONFIG_FILE).")

mode_label <- paste0(column, "_", polarity)  # e.g. RP_POS — used to namespace all output folders

out_xlsx   <- file.path(project_folder, "intermediates", mode_label, "notame_rev.xlsx")
interdir   <- file.path(project_folder, "intermediates", mode_label)
output_dir <- file.path(project_folder, "output",        mode_label)

QC_DETECTION_LIMIT     <- as.numeric(get_env("QC_DETECTION_LIMIT",     "0.60"))
SAMPLE_DETECTION_LIMIT <- as.numeric(get_env("SAMPLE_DETECTION_LIMIT", "0.20"))

blank_ratio_env <- get_env("BLANK_RATIO", "none")
BLANK_RATIO <- if (blank_ratio_env %in% c("none", "skip", "")) NA_real_ else
               as.numeric(blank_ratio_env)

LOW_INT_FILTER      <- suppressWarnings(as.numeric(get_env("LOW_INT_FILTER", "")))
LOW_INT_FILTER_FRAC <- as.numeric(get_env("LOW_INT_FILTER_FRAC", "0.10"))
LOW_INT_PERCENTILE  <- as.numeric(get_env("LOW_INT_PERCENTILE",  "0.8"))
qc_rsd_env    <- get_env("QC_RSD_FILTER", "none")
QC_RSD_FILTER           <- if (qc_rsd_env %in% c("none", "")) NA_real_ else as.numeric(qc_rsd_env)
MIN_QC_SAMPLE_DETECTION <- as.numeric(get_env("MIN_QC_SAMPLE_DETECTION", "0.50"))
QC_OUTLIER_MAD_K        <- as.numeric(get_env("QC_OUTLIER_MAD_K", "5"))
QC_OUTLIER_MIN_N        <- as.integer(get_env("QC_OUTLIER_MIN_N", "3"))
MIN_BATCH_DETECTION     <- as.integer(get_env("MIN_BATCH_DETECTION", "1"))
MIN_BATCH_DETECTION_FRAC <- as.numeric(get_env("MIN_BATCH_DETECTION_FRAC", "0"))
RSD_THRESHOLD <- as.numeric(get_env("RSD_THRESHOLD", "0.30"))

# Correction methods: each CORRECTION_METHODS entry is a "drift:basis:batch"
# spec, not one opaque name per combination -- see the --help text above
# (searchable: "DRIFT_METHOD (1st field)") for the full vocabulary and
# legality rules, and R/method_spec.R for where that vocabulary is defined.
CORRECTION_METHODS <- strsplit(get_env("CORRECTION_METHODS", "none:none:none,notame_spline:qc:ruv_s"), ",")[[1]]

RUV_K      <- as.integer(get_env("RUV_K",       "3"))
SERRF_NUM  <- as.integer(get_env("SERRF_NUM",   "5"))
LOESS_QC_SPAN             <- as.numeric(get_env("LOESS_QC_SPAN", "0.75"))
LOESS_SAMPLE_SPAN         <- as.numeric(get_env("LOESS_SAMPLE_SPAN", "0.9"))
DRIFT_SAMPLE_MIN_OBS      <- as.integer(get_env("DRIFT_SAMPLE_MIN_OBS", "10"))
DRIFT_MIN_QC_PER_BATCH    <- as.integer(get_env("DRIFT_MIN_QC_PER_BATCH", "4"))
DRIFT_MIN_LTQC_VALIDATE   <- as.integer(get_env("DRIFT_MIN_LTQC_VALIDATE", "3"))
DRIFT_HYBRID_VALIDATE     <- as.logical(get_env("DRIFT_HYBRID_VALIDATE", "TRUE"))
cordbat_ref_env   <- get_env("CORDBAT_REF_BATCH", "")
CORDBAT_REF_BATCH <- if (cordbat_ref_env == "") NULL else cordbat_ref_env
parse_num_list <- function(s) as.numeric(strsplit(s, ",")[[1]])
# WAVEICA_ALPHA/WAVEICA_CUTOFF: comma-separated lists -- a single value keeps
# today's fixed behaviour (one WaveICA_2.0() call); more than one candidate
# overall (across alpha/cutoff/K together) triggers a search, evaluated
# against WAVEICA_EVAL_GROUP/Sample D-ratio (see select_waveica_params() in
# R/correction_methods.R).
WAVEICA_ALPHA   <- parse_num_list(get_env("WAVEICA_ALPHA",  "0.05"))
WAVEICA_CUTOFF  <- parse_num_list(get_env("WAVEICA_CUTOFF", "0.10"))
# WAVEICA_K: same list convention, but "auto" (or an empty token) means
# "2 x n_batches", represented internally as NA -- resolved at call time
# since it depends on the data, not at parse time.
parse_waveica_k_list <- function(s) {
  vapply(strsplit(s, ",")[[1]], function(tok) {
    tok <- trimws(tok)
    if (tok == "" || tolower(tok) == "auto") NA_real_ else as.numeric(tok)
  }, numeric(1), USE.NAMES = FALSE)
}
WAVEICA_K          <- parse_waveica_k_list(get_env("WAVEICA_K", "auto"))
WAVEICA_WF         <- get_env("WAVEICA_WF", "haar")
WAVEICA_EVAL_GROUP <- get_env("WAVEICA_EVAL_GROUP", "ltQC")
WAVEICA_V1_WF         <- get_env("WAVEICA_V1_WF", "haar")
WAVEICA_V1_K          <- parse_waveica_k_list(get_env("WAVEICA_V1_K", "20"))
WAVEICA_V1_T          <- parse_num_list(get_env("WAVEICA_V1_T", "0.05"))
WAVEICA_V1_T2         <- as.numeric(get_env("WAVEICA_V1_T2", "0.05"))
WAVEICA_V1_ALPHA      <- parse_num_list(get_env("WAVEICA_V1_ALPHA", "0"))
WAVEICA_V1_EVAL_GROUP <- get_env("WAVEICA_V1_EVAL_GROUP", "ltQC")
AUTO_LOESS_SPANS        <- parse_num_list(get_env("AUTO_LOESS_SPANS",        "0.5,0.75,0.9"))
AUTO_HUBER_KS           <- parse_num_list(get_env("AUTO_HUBER_KS",           "1.0,1.345,2.0"))
AUTO_SAMPLE_LOESS_SPANS <- parse_num_list(get_env("AUTO_SAMPLE_LOESS_SPANS", "0.3,0.6,0.9"))
AUTO_SAMPLE_HUBER_KS    <- parse_num_list(get_env("AUTO_SAMPLE_HUBER_KS",    "1.0,1.345,2.0"))
AUTO_MIN_QC_PER_BATCH   <- as.integer(get_env("AUTO_MIN_QC_PER_BATCH",  "4"))
AUTO_MIN_LTQC_VALIDATE  <- as.integer(get_env("AUTO_MIN_LTQC_VALIDATE", "3"))
AUTO_MIN_CV_OBS         <- as.integer(get_env("AUTO_MIN_CV_OBS",       "4"))
HUBER_QC_K     <- as.numeric(get_env("HUBER_QC_K",     "1.345"))
HUBER_SAMPLE_K <- as.numeric(get_env("HUBER_SAMPLE_K", "1.345"))
sva_n_sv_env <- get_env("SVA_N_SV", "")
SVA_N_SV     <- if (sva_n_sv_env == "") NULL else as.integer(sva_n_sv_env)
LOESS_QC_CV_SPANS <- parse_num_list(get_env("LOESS_QC_CV_SPANS", ""))
HUBER_QC_CV_KS    <- parse_num_list(get_env("HUBER_QC_CV_KS",    ""))
COMBAT_MEAN_ONLY <- get_env("COMBAT_MEAN_ONLY", "auto")
COMBAT_PAR_PRIOR <- get_env("COMBAT_PAR_PRIOR", "auto")
NORMALIZATION              <- get_env("NORMALIZATION",              "none")
SAVE_PRE_CORRECTION_PLOTS  <- as.logical(get_env("SAVE_PRE_CORRECTION_PLOTS", "TRUE"))
FORCE_RECONVERT            <- as.logical(get_env("FORCE_RECONVERT", "FALSE"))

# ─────────────────────────────────────────────────────────────────────────────
# PREFLIGHT CHECKS
# ─────────────────────────────────────────────────────────────────────────────

run_preflight_checks(
  input_mode = input_mode, in_xlsx = in_xlsx,
  in_feature_table = in_feature_table, in_sample_sheet = in_sample_sheet,
  project_folder = project_folder, column = column, polarity = polarity,
  correction_methods = CORRECTION_METHODS, normalization = NORMALIZATION,
  qc_detection_limit = QC_DETECTION_LIMIT, sample_detection_limit = SAMPLE_DETECTION_LIMIT,
  low_int_filter_frac = LOW_INT_FILTER_FRAC, low_int_percentile = LOW_INT_PERCENTILE,
  min_qc_sample_detection = MIN_QC_SAMPLE_DETECTION, min_batch_detection = MIN_BATCH_DETECTION,
  min_batch_detection_frac = MIN_BATCH_DETECTION_FRAC, qc_outlier_mad_k = QC_OUTLIER_MAD_K,
  qc_outlier_min_n = QC_OUTLIER_MIN_N,
  rsd_threshold = RSD_THRESHOLD, ruv_k = RUV_K, serrf_num = SERRF_NUM, loess_qc_span = LOESS_QC_SPAN,
  loess_sample_span = LOESS_SAMPLE_SPAN, drift_sample_min_obs = DRIFT_SAMPLE_MIN_OBS,
  drift_min_qc_per_batch = DRIFT_MIN_QC_PER_BATCH,
  drift_min_ltqc_validate = DRIFT_MIN_LTQC_VALIDATE,
  auto_loess_spans = AUTO_LOESS_SPANS, auto_huber_ks = AUTO_HUBER_KS,
  auto_sample_loess_spans = AUTO_SAMPLE_LOESS_SPANS, auto_sample_huber_ks = AUTO_SAMPLE_HUBER_KS,
  auto_min_qc_per_batch = AUTO_MIN_QC_PER_BATCH, auto_min_ltqc_validate = AUTO_MIN_LTQC_VALIDATE,
  auto_min_cv_obs = AUTO_MIN_CV_OBS,
  huber_qc_k = HUBER_QC_K, huber_sample_k = HUBER_SAMPLE_K,
  sva_n_sv = if (is.null(SVA_N_SV)) NA_integer_ else SVA_N_SV,
  loess_qc_cv_spans = LOESS_QC_CV_SPANS, huber_qc_cv_ks = HUBER_QC_CV_KS,
  waveica_alpha = WAVEICA_ALPHA, waveica_cutoff = WAVEICA_CUTOFF, waveica_k = WAVEICA_K,
  waveica_eval_group = WAVEICA_EVAL_GROUP,
  waveica_v1_k = WAVEICA_V1_K, waveica_v1_t = WAVEICA_V1_T,
  waveica_v1_t2 = WAVEICA_V1_T2, waveica_v1_alpha = WAVEICA_V1_ALPHA,
  waveica_v1_eval_group = WAVEICA_V1_EVAL_GROUP,
  combat_mean_only = COMBAT_MEAN_ONLY, combat_par_prior = COMBAT_PAR_PRIOR,
  blank_ratio = BLANK_RATIO, low_int_filter = LOW_INT_FILTER, qc_rsd_filter = QC_RSD_FILTER,
  save_pre_correction_plots = SAVE_PRE_CORRECTION_PLOTS,
  config_file = config_file, raw_sample_type_rules = config$sample_type_rules
)

# Resolved after preflight has validated it (pattern/type present, pattern a
# valid regex) so all config problems are still reported together up front.
sample_type_rules <- resolve_sample_type_rules(config)

# ─────────────────────────────────────────────────────────────────────────────
# 1) CONVERT & IMPORT
# ─────────────────────────────────────────────────────────────────────────────

dir.create(interdir, showWarnings = FALSE, recursive = TRUE)

# Cache the conversion step: re-parsing the source export(s) is one of the
# slower steps and only needs to be redone when the input(s) change. The
# annotations table (returned in memory, not part of out_xlsx) is cached
# alongside it so a cache hit can skip the converter call entirely.
input_files     <- if (input_mode == "xcms") c(in_feature_table, in_sample_sheet) else in_xlsx
annotations_rds <- file.path(interdir, "conversion_annotations.rds")
newest_input    <- max(file.mtime(input_files))
cache_valid <- !FORCE_RECONVERT &&
  file.exists(out_xlsx) && file.exists(annotations_rds) &&
  file.mtime(out_xlsx)        >= newest_input &&
  file.mtime(annotations_rds) >= newest_input &&
  (config_file == "" ||
     (file.mtime(out_xlsx) >= file.mtime(config_file) &&
      file.mtime(annotations_rds) >= file.mtime(config_file)))  # invalidate on sample_type_rules changes too

if (cache_valid) {
  message("==> Using cached ", toupper(input_mode), " conversion: ", out_xlsx,
          " (set FORCE_RECONVERT=TRUE to re-run)")
  convert_result <- readRDS(annotations_rds)
} else if (input_mode == "xcms") {
  convert_result <- xcms_to_notame(in_feature_table, in_sample_sheet, out_xlsx, column, polarity)
  saveRDS(convert_result, annotations_rds)
} else {
  convert_result <- msdial_to_notame(in_xlsx, out_xlsx, sample_type_rules = sample_type_rules)
  saveRDS(convert_result, annotations_rds)
}

mode_name           <- convert_result$mode
feature_annotations <- convert_result$annotations

message("==> Importing")
data <- import_from_excel(file = out_xlsx, sheet = 1, name = mode_name)
names(assays(data)) <- "abundances"
data <- fix_object(data, assay.type = "abundances")

cat("Imported:", nrow(data), "features,", ncol(data), "samples\n")
print(table(colData(data)$QC))

# ─────────────────────────────────────────────────────────────────────────────
# 2) PRE-CORRECTION FEATURE FILTERING
# ─────────────────────────────────────────────────────────────────────────────

message("==> Pre-correction feature filtering")

data    <- mark_nas(data, value = 0)
n_before <- nrow(data)

# Blank filter
if (!is.na(BLANK_RATIO)) {
  blank_idx <- which(colData(data)$QC == "Blank")
  if (length(blank_idx) > 0) {
    sample_idx   <- which(colData(data)$QC == "Sample")
    blank_means  <- rowMeans(assay(data)[, blank_idx,  drop = FALSE], na.rm = TRUE)
    sample_means <- rowMeans(assay(data)[, sample_idx, drop = FALSE], na.rm = TRUE)
    keep_blank   <- is.na(blank_means) | blank_means == 0 | sample_means > BLANK_RATIO * blank_means
    keep_blank[is.na(keep_blank)] <- FALSE
    data <- data[keep_blank, ]
  } else {
    cat("No blank samples found — skipping blank filter\n")
  }
}
n_after_blank <- nrow(data)

# Remove non-analytical sample types; retain ltQC for downstream evaluation
data <- data[, !colData(data)$QC %in% c("Blank", "Wash", "Cond", "MSe", "MS2", "SST", "MatrixBlank")]

# Low-intensity filter: remove features below a fraction of the mean pN intensity
{
  mat_int      <- assay(data); mat_int[is.na(mat_int)] <- 0
  int_quantile <- apply(mat_int, 1, quantile, probs = LOW_INT_PERCENTILE)
  mean_p       <- mean(int_quantile[int_quantile > 0])

  low_int_cutoff <- if (!is.na(LOW_INT_FILTER)) {
    LOW_INT_FILTER
  } else if (!is.na(LOW_INT_FILTER_FRAC)) {
    LOW_INT_FILTER_FRAC * mean_p
  } else {
    NA_real_
  }

  if (!is.na(low_int_cutoff)) {
    data <- data[int_quantile >= low_int_cutoff, ]
    cat(sprintf("Low-intensity filter (p%.0f < %.4g): removed %d features\n",
                LOW_INT_PERCENTILE * 100, low_int_cutoff, nrow(data) - n_after_blank))
  }
}
n_after_lowint <- nrow(data)

# Remove QC/ltQC samples with insufficient feature detection (empty injections,
# failed runs). Runs after feature filters so detection rate is assessed on
# meaningful features only. ltQC is checked here too, using the same threshold
# as QC — it's the held-out group the loess_samples_combat/limma validation
# trial and the ltqc_permanova_* metrics rely on, and with typically only a
# few ltQC per batch, one bad injection would badly distort both.
for (qc_group in c("QC", "ltQC")) {
  qc_cols <- which(colData(data)$QC == qc_group)
  if (length(qc_cols) == 0) next
  mat_qc    <- assay(data)[, qc_cols, drop = FALSE]
  detect_qc <- colMeans(mat_qc > 0 & !is.na(mat_qc))
  bad_qc    <- qc_cols[detect_qc < MIN_QC_SAMPLE_DETECTION]
  if (length(bad_qc) > 0) {
    bad_names <- colData(data)$Sample_ID[bad_qc]
    bad_batch <- colData(data)$Batch[bad_qc]
    message("==> Removing ", length(bad_qc), " ", qc_group, " sample(s) with detection rate < ",
            round(MIN_QC_SAMPLE_DETECTION * 100), "%:")
    for (k in seq_along(bad_names))
      message("    ", bad_names[k], " (batch: ", bad_batch[k], ", detection: ",
              round(detect_qc[detect_qc < MIN_QC_SAMPLE_DETECTION][k] * 100, 1), "%)")
    data <- data[, -bad_qc]
  }
}

# Multivariate QC/ltQC outlier removal (PCA-distance-based, per batch).
# Runs after the detection-rate removal above so badly-detected samples
# (which would distort the reference centroid) are already gone; catches a
# different failure mode -- a sample that detects fine but has an anomalous
# intensity profile (contamination, carryover, a degrading/recovering
# column). QC and ltQC checked separately (see detect_qc_outliers()).
if (QC_OUTLIER_MAD_K > 0) {
  for (qc_group in c("QC", "ltQC")) {
    outlier_idx <- detect_qc_outliers(data, group = qc_group, mad_k = QC_OUTLIER_MAD_K,
                                       min_n = QC_OUTLIER_MIN_N)
    if (length(outlier_idx) > 0) data <- data[, -outlier_idx]
  }
}

# QC detection filter
data <- flag_detection(data, qc_limit = QC_DETECTION_LIMIT)
data <- drop_flagged(data)
n_after_qc <- nrow(data)

# Biological sample detection filter
sample_idx  <- which(colData(data)$QC == "Sample")
detect_rate <- rowMeans(!is.na(assay(data)[, sample_idx, drop = FALSE]))
data        <- data[detect_rate >= SAMPLE_DETECTION_LIMIT, ]
n_after_sample <- nrow(data)

# Zero-variance filter
zero_var <- apply(assay(data), 1, function(x) {
  v <- var(x, na.rm = TRUE)
  !is.na(v) && v < .Machine$double.eps
})
data <- data[!zero_var, ]
n_after_zerovar <- nrow(data)

# Per-batch detection filter: feature must have >= MIN_BATCH_DETECTION observations
# in every batch. Features absent from a whole batch have no real measurements
# there — all values would be LoD/2 placeholders.
if (MIN_BATCH_DETECTION > 0) {
  cd_batch  <- as.data.frame(colData(data))
  batches   <- unique(cd_batch$Batch)
  batch_det <- do.call(cbind, lapply(batches, function(b) {
    idx <- which(cd_batch$Batch == b)
    rowSums(!is.na(assay(data)[, idx, drop = FALSE]))
  }))
  pass_all_batches <- apply(batch_det, 1, function(x) all(x >= MIN_BATCH_DETECTION))
  data <- data[pass_all_batches, ]
}
n_after_batchdet <- nrow(data)

# Per-batch detection filter (fraction-based): feature must have detection in
# >= MIN_BATCH_DETECTION_FRAC of each batch's samples. Unlike
# MIN_BATCH_DETECTION (an absolute count), this scales with batch size, so
# batches of very different sizes get a consistent relative bar instead of a
# fixed count that's stringent for a small batch and lax for a large one --
# relevant for methods (e.g. batch_method=waveica/waveica_v1) whose
# LoD/2-imputed placeholders participate directly in fitting the correction,
# where a batch with disproportionately more missingness risks its
# placeholder pattern getting mistaken for real signal.
if (MIN_BATCH_DETECTION_FRAC > 0) {
  cd_batch_frac  <- as.data.frame(colData(data))
  batches_frac   <- unique(cd_batch_frac$Batch)
  batch_det_frac <- do.call(cbind, lapply(batches_frac, function(b) {
    idx <- which(cd_batch_frac$Batch == b)
    rowMeans(!is.na(assay(data)[, idx, drop = FALSE]))
  }))
  pass_all_batches_frac <- apply(batch_det_frac, 1, function(x) all(x >= MIN_BATCH_DETECTION_FRAC))
  data <- data[pass_all_batches_frac, ]
}
n_after_batchdetfrac <- nrow(data)

# Pre-imputation QC-RSD filter (per batch; pass = acceptable RSD in >= 1 batch)
if (!is.na(QC_RSD_FILTER)) {
  cd_pre  <- as.data.frame(colData(data))
  batches <- unique(cd_pre$Batch)
  rsd_mat <- do.call(cbind, lapply(batches, function(b) {
    idx <- which(cd_pre$QC == "QC" & cd_pre$Batch == b)
    if (length(idx) < 2) return(rep(NA_real_, nrow(data)))
    apply(assay(data)[, idx, drop = FALSE], 1, function(x) {
      x <- x[!is.na(x)]
      med <- median(x)
      if (length(x) < 2 || med <= 0) return(NA_real_)
      mad(x) / med
    })
  }))
  colnames(rsd_mat) <- as.character(batches)
  passes_any <- apply(rsd_mat, 1, function(r) any(is.na(r) | r <= QC_RSD_FILTER))
  data <- data[passes_any, ]
}
n_after_qcrsd <- nrow(data)

# Print and save filter summary
filter_log <- data.frame(
  step = c(
    "Blank filter",
    if (!is.na(low_int_cutoff))
      sprintf("Low-intensity filter (p%.0f >= %.4g)", LOW_INT_PERCENTILE * 100, low_int_cutoff)
    else
      "Low-intensity filter (disabled)",
    sprintf("QC detection (>= %.0f%%)", QC_DETECTION_LIMIT * 100),
    sprintf("Sample detection (>= %.0f%%)", SAMPLE_DETECTION_LIMIT * 100),
    "Zero variance",
    if (MIN_BATCH_DETECTION > 0) sprintf("Per-batch detection (>= %d per batch)", MIN_BATCH_DETECTION) else "Per-batch detection (disabled)",
    if (MIN_BATCH_DETECTION_FRAC > 0) sprintf("Per-batch detection (>= %.0f%% per batch)", MIN_BATCH_DETECTION_FRAC * 100) else "Per-batch detection, fraction (disabled)",
    if (!is.na(QC_RSD_FILTER)) sprintf("QC-RSD filter (<= %.0f%% in >= 1 batch)", QC_RSD_FILTER * 100) else "QC-RSD filter (disabled)"
  ),
  features_removed = c(
    n_before             - n_after_blank,
    n_after_blank        - n_after_lowint,
    n_after_lowint       - n_after_qc,
    n_after_qc           - n_after_sample,
    n_after_sample       - n_after_zerovar,
    n_after_zerovar      - n_after_batchdet,
    n_after_batchdet     - n_after_batchdetfrac,
    n_after_batchdetfrac - n_after_qcrsd
  ),
  features_remaining = c(
    n_after_blank, n_after_lowint,
    n_after_qc, n_after_sample, n_after_zerovar, n_after_batchdet, n_after_batchdetfrac, n_after_qcrsd
  )
)
cat("\n--- Pre-filtering summary ---\n")
print(filter_log, row.names = FALSE)
write.csv(filter_log, file.path(interdir, "prefilter_log.csv"), row.names = FALSE)

report_batch_summary(data, file = file.path(interdir, "batch_summary.csv"))

# Run parameters, kept as a data.frame so the same values can be written to
# run_parameters.txt and embedded in each method's Settings sheet.
run_params <- list(
  "Run timestamp"           = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "CONFIG_FILE"             = if (config_file == "") "(none)" else config_file,
  "INPUT_MODE"              = input_mode,
  "IN_XLSX"                 = if (input_mode == "msdial") in_xlsx else "(n/a)",
  "IN_FEATURE_TABLE"        = if (input_mode == "xcms") in_feature_table else "(n/a)",
  "IN_SAMPLE_SHEET"         = if (input_mode == "xcms") in_sample_sheet else "(n/a)",
  "COLUMN"                  = column,
  "POLARITY"                = polarity,
  "CORRECTION_METHODS"      = paste(CORRECTION_METHODS, collapse = ", "),
  "QC_DETECTION_LIMIT"      = QC_DETECTION_LIMIT,
  "SAMPLE_DETECTION_LIMIT"  = SAMPLE_DETECTION_LIMIT,
  "BLANK_RATIO"             = BLANK_RATIO,
  "LOW_INT_FILTER"          = if (!is.na(LOW_INT_FILTER)) LOW_INT_FILTER else "(disabled)",
  "LOW_INT_FILTER_FRAC"     = LOW_INT_FILTER_FRAC,
  "LOW_INT_PERCENTILE"      = LOW_INT_PERCENTILE,
  "LOW_INT_CUTOFF"          = if (!is.na(low_int_cutoff)) low_int_cutoff else "(disabled)",
  "MIN_QC_SAMPLE_DETECTION" = MIN_QC_SAMPLE_DETECTION,
  "QC_OUTLIER_MAD_K"        = QC_OUTLIER_MAD_K,
  "QC_OUTLIER_MIN_N"        = QC_OUTLIER_MIN_N,
  "MIN_BATCH_DETECTION"     = MIN_BATCH_DETECTION,
  "MIN_BATCH_DETECTION_FRAC" = MIN_BATCH_DETECTION_FRAC,
  "QC_RSD_FILTER"           = if (!is.na(QC_RSD_FILTER)) QC_RSD_FILTER else "(disabled)",
  "RSD_THRESHOLD"           = RSD_THRESHOLD,
  "RUV_K"                   = RUV_K,
  "SERRF_NUM"               = SERRF_NUM,
  "LOESS_QC_SPAN"           = LOESS_QC_SPAN,
  "LOESS_SAMPLE_SPAN"       = LOESS_SAMPLE_SPAN,
  "DRIFT_SAMPLE_MIN_OBS"    = DRIFT_SAMPLE_MIN_OBS,
  "CORDBAT_REF_BATCH"       = if (is.null(CORDBAT_REF_BATCH)) "(auto)" else CORDBAT_REF_BATCH,
  "NORMALIZATION"           = NORMALIZATION,
  "N_CORES"                 = if (n_cores_env == "") paste(parallel::detectCores() - 1, "(auto)") else n_cores_env
)
run_params_df <- data.frame(Parameter = names(run_params),
                            Value     = unlist(run_params, use.names = FALSE),
                            stringsAsFactors = FALSE)

writeLines(sprintf("%-24s%s", paste0(run_params_df$Parameter, ":"), run_params_df$Value),
           file.path(interdir, "run_parameters.txt"))

# Sanity check: injection order must be finite for all samples
bad_inj <- !is.finite(colData(data)$Injection_order)
if (any(bad_inj)) {
  cat("WARNING: samples with non-finite Injection_order:\n")
  print(as.data.frame(colData(data))[bad_inj,
    intersect(c("Sample_ID", "Original_name", "QC", "Batch", "Injection_order"),
              colnames(colData(data)))])
  stop("Non-finite injection orders found — fix the ", input_mode, "_to_notame() conversion before proceeding.")
} else {
  cat("Injection_order: OK\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) CORRECTION & OUTPUT (per method)
# ─────────────────────────────────────────────────────────────────────────────

# Pre-correction QC plots
if (SAVE_PRE_CORRECTION_PLOTS) {
  dir.create(file.path(output_dir, "pre_correction"), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    save_QC_plots(data, prefix = file.path(output_dir, "pre_correction/"), id = "Sample_ID",
                  perplexity = safe_perplexity(ncol(data))),
    error = function(e) message("WARNING: save_QC_plots failed (pre-correction): ", conditionMessage(e))
  )
}

old_summaries <- list.files(interdir, pattern = "^qc_summary_.+\\.csv$", full.names = TRUE)
if (length(old_summaries) > 0) file.remove(old_summaries)

# Build raw reference for signal-preservation metric (biological samples only).
# Intentionally NOT imputed — NAs are kept so that signal_preservation_r is
# computed only over originally-observed positions. Positions missing in the
# raw data are excluded from the correlation naturally via is.finite() checks.
message("==> Building raw reference for signal-preservation metric")
raw_ref <- tryCatch({
  samp_idx <- which(colData(data)$QC == "Sample")
  mat      <- assay(data, 1)[, samp_idx, drop = FALSE]
  colnames(mat) <- colData(data)$Sample_ID[samp_idx]
  write.csv(as.data.frame(mat), file.path(output_dir, "raw_reference.csv"))
  mat
}, error = function(e) {
  message("WARNING: could not build raw reference (signal_preservation_r will be NA): ", conditionMessage(e))
  NULL
})

# Uncorrected baseline QC metrics
save_correction_summary(assess_quality(data), method = "uncorrected", interdir = interdir, raw_ref = raw_ref)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(colData(data)), file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)

run_log <- new_run_log()

# Single shared params list read by run_correction()'s drift/batch dispatch
# (R/correction_methods.R) for every CORRECTION_METHODS entry. serrf_num_eff
# (SERRF_NUM capped at half the smallest per-batch QC count, to avoid
# overfitting) is computed once here from `data`, same value every
# batch_method=serrf run would have computed inline before.
serrf_min_qc_per_batch <- if (any(colData(data)$QC == "QC"))
  min(table(colData(data)$Batch[colData(data)$QC == "QC"])) else 0L
serrf_num_eff <- max(1L, min(SERRF_NUM, floor(serrf_min_qc_per_batch / 2L)))
if (serrf_num_eff < SERRF_NUM)
  message("==> SERRF: reducing num from ", SERRF_NUM, " to ", serrf_num_eff,
          " (min QC per batch = ", serrf_min_qc_per_batch, ")")

params <- list(
  loess_qc_span = LOESS_QC_SPAN, loess_sample_span = LOESS_SAMPLE_SPAN,
  huber_qc_k = HUBER_QC_K, huber_sample_k = HUBER_SAMPLE_K,
  drift_sample_min_obs = DRIFT_SAMPLE_MIN_OBS,
  min_qc_per_batch = DRIFT_MIN_QC_PER_BATCH, min_ltqc_validate = DRIFT_MIN_LTQC_VALIDATE,
  drift_hybrid_validate = DRIFT_HYBRID_VALIDATE,
  loess_qc_cv_spans = LOESS_QC_CV_SPANS, huber_qc_cv_ks = HUBER_QC_CV_KS,
  auto_loess_spans = AUTO_LOESS_SPANS, auto_huber_ks = AUTO_HUBER_KS,
  auto_sample_loess_spans = AUTO_SAMPLE_LOESS_SPANS, auto_sample_huber_ks = AUTO_SAMPLE_HUBER_KS,
  auto_min_qc_per_batch = AUTO_MIN_QC_PER_BATCH, auto_min_ltqc_validate = AUTO_MIN_LTQC_VALIDATE,
  auto_min_cv_obs = AUTO_MIN_CV_OBS,
  combat_mean_only = COMBAT_MEAN_ONLY, combat_par_prior = COMBAT_PAR_PRIOR, sva_n_sv = SVA_N_SV,
  ruv_k = RUV_K, cordbat_ref_batch = CORDBAT_REF_BATCH,
  waveica_alpha = WAVEICA_ALPHA, waveica_cutoff = WAVEICA_CUTOFF, waveica_k = WAVEICA_K, waveica_wf = WAVEICA_WF,
  waveica_eval_group = WAVEICA_EVAL_GROUP,
  waveica_v1_wf = WAVEICA_V1_WF, waveica_v1_k = WAVEICA_V1_K, waveica_v1_t = WAVEICA_V1_T,
  waveica_v1_t2 = WAVEICA_V1_T2, waveica_v1_alpha = WAVEICA_V1_ALPHA,
  waveica_v1_eval_group = WAVEICA_V1_EVAL_GROUP,
  serrf_num_eff = serrf_num_eff
)

for (method in CORRECTION_METHODS) {
  message("\n############################################################")
  message("# METHOD: ", toupper(method))
  message("############################################################")

  method_out   <- file.path(output_dir, sanitize_method_id(method))
  method_start <- Sys.time()
  dir.create(file.path(method_out, "QC_plots"), showWarnings = FALSE, recursive = TRUE)

  result <- tryCatch({
    spec <- parse_correction_method_spec(method)
    if (!is.list(spec)) stop(spec)  # already preflight-checked; re-validated here defensively
    run_correction(data, drift_method = spec$drift, basis = spec$basis,
                   batch_method = spec$batch, params = params)
  }, error = function(e) {
    message("ERROR in method '", method, "': ", conditionMessage(e))
    message("Skipping.")
    run_log <<- log_method_result(run_log, method, "failed",
                                   as.numeric(difftime(Sys.time(), method_start, units = "secs")),
                                   error = conditionMessage(e))
    NULL
  })

  if (is.null(result)) next

  combined <- result$post
  obs_mask <- result$obs_mask

  # Optional post-correction normalisation for dilution effects
  if (NORMALIZATION == "pqn") {
    tryCatch({
      suppressPackageStartupMessages(library(pmp))
      classes    <- colData(combined)$QC
      qc_present <- any(classes == "QC")
      combined   <- pqn_normalisation(df = combined, classes = classes,
                                      qc_label = if (qc_present) "QC" else NULL)
      message("==> PQN normalisation applied (reference: ",
              if (qc_present) "pooled QC" else "median of samples", ")")
    }, error = function(e) {
      message("WARNING: PQN normalisation failed: ", conditionMessage(e))
    })
  }

  tryCatch(
    save_QC_plots(combined, prefix = file.path(method_out, "QC_plots/post_correction_"), id = "Sample_ID",
                  perplexity = safe_perplexity(ncol(combined))),
    error = function(e) message("WARNING: save_QC_plots failed (", method, "): ", conditionMessage(e))
  )

  combined <- assess_quality(combined)
  save_correction_summary(combined, method = method, interdir = interdir, obs_mask = obs_mask, raw_ref = raw_ref)
  report_batch_summary(combined, file = file.path(method_out, "batch_summary_post_correction.csv"))

  # Drop RUV W-factor columns from final output
  w_cols <- grep("^W_", colnames(colData(combined)), value = TRUE)
  if (length(w_cols) > 0)
    colData(combined) <- colData(combined)[, !colnames(colData(combined)) %in% w_cols]

  message("==> Writing output: ", method)
  combined <- add_batch_qc_metrics(combined)

  # Settings sheet content for this method: global run parameters plus the
  # method and feature-set context for this particular workbook.
  method_settings_df <- function(feature_set_label) {
    rbind(run_params_df,
          data.frame(Parameter = c("CORRECTION_METHOD", "FEATURE_SET"),
                     Value     = c(method, feature_set_label),
                     stringsAsFactors = FALSE))
  }

  # Helper: write the full (uncompressed) and clustered (compressed) results
  # workbooks for a given SE.
  write_outputs <- function(se, suffix, feature_set_label) {
    settings_df <- method_settings_df(feature_set_label)
    write_results_workbook(se, feature_annotations, settings_df,
                            file = file.path(method_out, paste0("results_full", suffix, ".xlsx")))
    se_c <- tryCatch(compress_clusters(cluster_features(se, all_features = TRUE)),
                     error = function(e) { message("WARNING: clustering failed: ", conditionMessage(e)); NULL })
    if (!is.null(se_c)) {
      write_results_workbook(se_c, feature_annotations, settings_df,
                              file = file.path(method_out, paste0("results_clustered", suffix, ".xlsx")))
    }
  }

  # Full feature set
  write_outputs(combined, "", "All features")

  # Global QC-RSD filter
  rsd_suffix <- paste0("_rsd", round(RSD_THRESHOLD * 100))
  tryCatch({
    global_keep <- !is.na(rowData(combined)$RSD_r) & rowData(combined)$RSD_r < RSD_THRESHOLD
    write_outputs(combined[global_keep, ], rsd_suffix,
                  sprintf("QC RSD < %.0f%%", RSD_THRESHOLD * 100))
  }, error = function(e) message("WARNING: global RSD filter export failed (", method, "): ", conditionMessage(e)))

  # Batchwise QC-RSD filter (pass in >= 50% of batches)
  batch_rsd_cols <- grep("^RSD_r_", colnames(rowData(combined)), value = TRUE)
  if (length(batch_rsd_cols) > 0) {
    tryCatch({
      batch_rsd_mat <- as.matrix(as.data.frame(rowData(combined))[, batch_rsd_cols, drop = FALSE])
      frac_passing  <- rowMeans(batch_rsd_mat < RSD_THRESHOLD, na.rm = TRUE)
      write_outputs(combined[!is.na(frac_passing) & frac_passing >= 0.5, ],
                    paste0("_batchrsd", round(RSD_THRESHOLD * 100)),
                    sprintf("Per-batch QC RSD < %.0f%% in >= 50%% of batches", RSD_THRESHOLD * 100))
    }, error = function(e) message("WARNING: batchwise RSD filter export failed (", method, "): ", conditionMessage(e)))
  }

  run_log <- log_method_result(run_log, method, "success",
                                as.numeric(difftime(Sys.time(), method_start, units = "secs")),
                                n_features = nrow(combined))
}

run_log_df <- write_run_log(run_log, interdir)

compare_corrections(interdir, output_dir)

primary_input <- if (input_mode == "xcms") in_feature_table else in_xlsx
tryCatch(
  write_run_report(interdir, output_dir, mode_label, primary_input, run_params_df, run_log_df),
  error = function(e) message("WARNING: report generation failed: ", conditionMessage(e))
)

message("==> FINISHED. Output at: ", output_dir)
