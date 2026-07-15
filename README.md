# Protein-Language-Model Features vs. Gene Expression — TCGA-BRCA

A methods project on a well-characterised breast-cancer cohort (TCGA-BRCA), in two acts:

1. **Unsupervised stratification as a positive control** — starting only from molecular
   data, can we recover biologically distinct patient subgroups, and do they differ in
   survival? *(notebooks 01–03)*
2. **Do protein-language-model (ESM2) features add predictive signal over expression?**
   Tested two independent ways — reference-sequence embeddings and mutation
   variant-effect scores. *(notebooks 04–05)*

The cohort is used as a **controlled testbed**, not a discovery set: if a method works,
it should recover the known biology; if an added feature helps, it should beat a plain
expression baseline under cross-validation. The headline finding of Act 2 is an honest,
well-powered **negative result**, reported and explained rather than tuned away.

---

## Act 1 — Unsupervised stratification (01–03)

### Question
From molecular data alone (no labels), can we recover biologically distinct patient
subgroups, and do those subgroups differ in survival?

### Data
Source: [UCSC Xena](https://xenabrowser.net/), cohort *TCGA Breast Cancer (BRCA)*.

| Modality | Dataset | Shape (genes × samples) |
|---|---|---|
| Gene expression | RNA-seq, IlluminaHiSeq, `log2(norm+1)` | ~20k × 1218 |
| Copy number | gene-level GISTIC2 (continuous) | ~24k × 1080 |
| Clinical | phenotype matrix (survival + PAM50) | 1080 × 193 |

After keeping primary tumours only and intersecting the three tables: **~1,070 patients**
measured in all modalities, of which **801** have complete follow-up and ~820 carry a
PAM50 label. Raw data is **not** committed — `01_download.ipynb` fetches it from Xena.

### Method
```
01_download      →  02_preprocess           →  03_stratify_survival
fetch from Xena     align + integrate           cluster + validate
```

**Preprocessing (per modality, independently):** transpose to patients × genes and keep
primary tumours (`-01` barcodes); intersect samples across all three tables (the honest
cost of multimodal integration — fewer patients); select the 2,000 most-variable genes
**before** z-scoring (variance is scale-dependent, so selecting after scaling would be
meaningless); z-score; reduce each modality with PCA (50 components).

**Integration:** concatenate the two 50-component PCA blocks. Reducing each modality
*first* and joining *after* ("intermediate" integration) keeps neither omic from
dominating just because it has more genes, and stays interpretable.

**Stratification:** cluster with **HDBSCAN** on a 10-dimensional UMAP of the integrated
space — not on the 2D map (HDBSCAN degrades in ~100 dimensions; the 2D map is only a
picture, and clustering on a picture invents structure). HDBSCAN is density-based, so the
number of clusters is not chosen in advance and patients in no dense region are reported
as noise rather than forced into a group.

**Validation:** *clinical* — Kaplan–Meier per cluster + log-rank test; *biological* —
overlap with PAM50 (used only afterwards to interpret clusters, never fed into clustering).

### Results
Two clusters, no noise:

| Cluster | n (with survival) | PAM50 composition |
|---|---|---|
| 0 | 129 | **98% Basal** |
| 1 | 672 | 59% LumA · 28% LumB · 10% Her2 · 3% Normal |

**Biological validation — strong.** The pipeline cleanly isolates the basal-like
(triple-negative) tumours into Cluster 0 without ever seeing PAM50. Unsupervised method
and known biology agree on this population — the positive control passes.

**Survival — no significant difference (log-rank p = 0.46).** The two clusters do not
separate overall survival robustly. This is expected and interpretable: the split the
method found is *Basal vs. everything else*, but the survival contrast in this cohort
lives among the **non-basal** subtypes (LumA good vs. LumB/Her2 worse) — all of which sit
inside Cluster 1, averaged together. A two-way split cannot express a signal that needs
distinguishing ≥3 groups. The Kaplan–Meier curves also **cross**, violating the
proportional-hazards assumption the log-rank test relies on, which further reduces its
power. So the high p-value means "no clean, monotone two-group separation", not "nothing
is happening".

Reported as-is: a real biological signal (basal vs. luminal) that does **not** translate
into a significant survival split is a legitimate, complete result. Robustness: the split
is stable across `min_cluster_size` ∈ {20, …, 40}, indicating real structure rather than
a single-parameter artefact.

![UMAP](outputs/figures/umap_clusters.png)
![Kaplan–Meier](outputs/figures/kaplan_meier.png)
![PAM50 overlap](outputs/figures/cluster_pam50_overlap.png)

---

## Act 2 — Does ESM2 add predictive signal over expression? (04–05)

[ESM2](https://www.science.org/doi/10.1126/science.ade2574) is a protein language model:
pass it a protein sequence and it returns an embedding capturing structure and function
learned from evolution. The question: does adding ESM2-derived features improve prediction
over gene expression alone? Two independent constructions were tested — both with the same
honest protocol (5-fold cross-validation; PCA + scaler + model fit **inside** each training
fold; report mean ± std).

### 04 — Reference-sequence embeddings
For each of the top-500 variable genes, fetch the canonical human protein sequence from
UniProt, embed it with ESM2, and build a **per-patient** embedding as an
expression-weighted average of those gene embeddings.

| Task | Expression | ESM2 | Combined |
|---|---|---|---|
| PAM50 subtype (5-fold CV accuracy) | 0.854 | 0.827 | 0.853 |
| Overall survival (5-fold CV C-index) | 0.584 ± 0.047 | 0.555 ± 0.131 | 0.614 ± 0.054 |

**Verdict: ESM2 adds nothing.** The construction is the reason — every patient receives the
*same reference sequence* per gene, so the only thing that varies patient-to-patient is the
expression weight. The patient embedding is therefore a **deterministic linear transform of
expression** (`weights.T @ fixed_embedding_matrix`); it cannot contain information expression
doesn't already have. ESM2-alone is slightly *worse* (a lossy re-encoding), and Combined ≈
Expression. The small survival "bump" for Combined was traced to **PCA truncation** (each
feature set gets its own top-20 PCs): raising `n_components` erases the gap entirely,
confirming it was an artifact, not signal.

### 05 — Mutation variant-effect scores
The fix for 04's flaw: give ESM2 sequences that **differ between patients** — their somatic
mutations. For each missense mutation, score it with ESM2's masked-marginal
**log-likelihood ratio** `logP(mutant) − logP(wildtype)` at that position (more negative =
more likely damaging; Meier et al., 2021). Build a patient × driver-gene "damage" matrix and
re-run the comparison.

Cohort after aligning expression + survival (all patients, not just mutated ones):
**801 patients, 90 deaths, 190 with a panel mutation.**

| Overall survival (5-fold CV C-index) | Expression | Mutation-ESM2 | Combined |
|---|---|---|---|
| mean ± std | 0.582 ± 0.048 | 0.474 ± 0.027 | 0.587 ± 0.042 |

**Verdict: mutation features add nothing either.** Because the folds are paired (shared
seed), Combined vs. Expression can be compared directly: Combined wins only 2 of 5 folds,
mean difference **+0.005** — a wash. Mutation-ESM2 alone is reliably ~random (tight std).
Three reasons, all biological rather than technical:

- **Sparsity** — only 190/801 patients carry a panel mutation, so 76% have an all-zero
  mutation vector; a feature block blank for three-quarters of the cohort can't carry much.
- **Redundancy** — by the time a driver mutation matters to a tumour, its consequences are
  downstream in the transcriptome (a TP53 mutation yields a p53-loss expression signature,
  etc.). Expression already encodes most of that, so mutation status is largely redundant
  *for outcome prediction*.
- **Hard endpoint** — 90 events with opposing-direction drivers (PIK3CA mutation trends
  good-prognosis, TP53 bad) leave little net overall-survival signal for anything to capture.

> **A note on rigor:** the first run of 05 gave uninterpretable results (351 patients, 39
> events, C-index std ≈ 0.17). The cause was a **cohort-selection bug** — the mutation matrix
> was built only over mutated patients, silently dropping everyone with a quiet genome and
> biasing the cohort. Fixing it (align on expression + survival, fill absent mutation rows
> with zeros) restored 801 patients / 90 events and produced the clean result above.

### What Act 2 shows
Across two independent constructions — reference-sequence embeddings and mutation
variant-effect scores — **ESM2 features do not add prognostic signal over gene expression in
TCGA-BRCA overall survival.** Expression is close to a sufficient statistic for this endpoint.
This is the kind of negative result worth reporting well: a plausible, popular hypothesis,
tested fairly, with a mechanistic explanation for *why* it fails — which is more informative
than a fragile positive.

Where ESM2 genuinely *does* carry signal is its intended use, **variant-effect prediction**
(pathogenic vs. benign), which does not require it to beat an expression baseline. That is the
natural next project rather than a rescue of this one.

---

## Repository structure
```
.
├── 01_download.ipynb            # fetch RNA-seq + CNV + clinical from UCSC Xena
├── 02_preprocess.ipynb          # transpose, align, select, scale, PCA, integrate
├── 03_stratify_survival.ipynb   # UMAP + HDBSCAN + Kaplan–Meier + PAM50 overlap
├── 04_esm2_embeddings.ipynb     # reference-sequence ESM2 features vs. expression
├── 05_mutation_esm2.ipynb       # mutation variant-effect ESM2 features vs. expression
├── outputs/                     # generated figures (committed)
├── data/                        # raw + processed data (gitignored)
├── requirements.txt
└── README.md
```

## Quickstart
```bash
pip install -r requirements.txt
# run the notebooks in order: 01 → 02 → 03  (Act 1),  then 04 → 05  (Act 2)
```

`requirements.txt`:
```
pandas
numpy
scikit-learn
umap-learn
hdbscan
lifelines
matplotlib
requests
transformers
torch
tqdm
```

## Limitations & honest notes
- **TCGA-BRCA is heavily characterised.** The point is not discovery but a controlled test
  of whether methods recover known biology and whether added features beat a baseline.
- **Overall survival is a hard, low-event endpoint** in BRCA; a modest C-index (~0.58 for
  expression) is expected, and the comparison — not the absolute number — is what matters.
- **ESM2 was the smallest model** (`esm2_t6_8M`) and mutation scoring was **missense-only**,
  so truncating drivers (many TP53/PTEN hits) are excluded and UniProt-vs-transcript isoform
  mismatches drop some variants. These limit sensitivity but do not explain the negative,
  which is driven by redundancy with expression and feature sparsity.
- **UMAP axes are not interpretable**; clusters are read by returning to the data (PAM50, PCA
  loadings), not by map geometry.
- **The cohort is imbalanced** (LumA dominates; Normal n≈19–22), so small-subtype proportions
  are noisier.
