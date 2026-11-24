# Analysis plan

The rules below were fixed before any code was run or any data downloaded. Any later change will be listed at the end of this file, with the date and the reason, as in the [AIDev pilot](https://github.com/laroshan/aidev-agent-pr-pilot).

## 1. Question

The original benchmark ([Abdelnaby & Moussa, bioRxiv 2025](https://doi.org/10.1101/2025.02.04.636499); [J. Comput. Biol. 2026](https://doi.org/10.1177/15578666261436821)) compares PCA with random projection (RP) for single-cell RNA-seq, on datasets of about 2,900 to 10,500 cells. Its scalability test covers 1,000 to 10,000 cells.

This project asks how the same methods behave at up to **250,000 cells**, with the data kept **sparse**:
1. How do run time and peak memory grow with the number of cells?
2. Does agreement with the authors' cell-type labels hold up?
3. Are rare cell types (under 1% of cells, including rare T-cell states) still found?

It is a short benchmark on one dataset and one laptop, not a paper. It claims no new method.

## 2. Bridge: reproduce one of the lab's results first

- **Code:** [moussa-lab/BenchmarkingPCA-RP](https://github.com/moussa-lab/BenchmarkingPCA-RP), pinned to commit `e6d0f92` (MIT licence).
- **Data:** the repo's sorted-PBMC set, `Datasets/PBMC-Zheng2017/PBMC_SC1.csv` (2,882 cells × 7,174 genes, 7 labelled populations).
- **Target:** the committed table `output/test/Labeled_PBMC/Labeled_PBMC_evaluation.csv`. It was made by reducing to 500 components, then a 3-D UMAP (`n_neighbors=5`, `min_dist=0.3`), then spherical k-means with k = 7.
- **Also run:**
  - the journal version's Matching Sparsity RP through the same route;
  - every method with five more seeds, to measure the run-to-run spread;
  - this project's own clustering route (§6) on the same data.
- **Success:** the same ranking of methods, with values inside the measured spread. Differences caused by newer library versions are reported, not hidden.

## 3. Data

- **Source:** CELLxGENE Census, release `2025-11-08` (long-term-support build), accessed with the `cellxgene_census` Python API.
- **Dataset:** one single-study blood dataset with at least 250,000 primary cells and author cell-type labels. Candidates, in order of preference:
  - **A.** OneK1K (Yazar et al., *Science* 2022; 1,248,980 cells, 29 labels, healthy donors, 10x 3′ v2).
  - **B.** AIDA Data Freeze v2 (Kock et al., *Cell* 2025; 1,265,624 cells, 32 labels).
  - **C.** COVID-19 PBMC (Stephenson et al., *Nat. Med.* 2021; 647,366 cells, 46 labels). Used only if A and B fail.
- **Choice rule:** read only the cell metadata (no expression data) for A and B. Pick the one with more rare (<1%) T-cell labels that still have at least 200 cells in a 250,000-cell sample. If they tie, pick A.
- **Labels:** `cell_type`, the authors' annotation mapped to the Cell Ontology. Only cells whose label is not a population (for example "peripheral blood mononuclear cell" or "unknown") are dropped. The exact list is recorded in `results/data_manifest.json`.
- **Subsets:** the cell IDs are shuffled once (seed 2026). The first 10,000 cells form the 10k set, and so on, giving nested sets of 10k, 25k, 50k, 100k and 250k cells.
- **Genes:** 2,000 highly variable genes are chosen once, with Scanpy's default method, on the 50k set with all genes. The same genes are used at every size.
- **Preprocessing:** each cell is divided by its total counts over all genes (`raw_sum` from the Census) and multiplied by 10,000, then transformed with log(1 + x). There is no scaling to unit variance, and the matrix stays sparse. No quality control beyond the authors' own filtering.

## 4. Methods

| # | Method | Implementation | Input |
|---|---|---|---|
| 1 | PCA, full SVD | scikit-learn `PCA(svd_solver="full")` | dense copy |
| 2 | PCA, randomized SVD | `PCA(svd_solver="randomized")` | dense copy |
| 3 | PCA, covariance method | `PCA(svd_solver="covariance_eigh")` | sparse |
| 4 | Gaussian RP | `GaussianRandomProjection` | sparse |
| 5 | Sparse RP | `SparseRandomProjection` (density 1/√d) | sparse |
| 6 | Matching Sparsity RP | the lab's `MatchingSparseRandomProjection` | sparse |

- Methods 1, 2, 4, 5 and 6 are the lab's settings.
- **Method 3 is an addition, made for fairness.** scikit-learn's full and randomized solvers don't accept sparse input. The covariance solver computes the same PCA without a dense copy.
- **Dense copy:** for methods 1 and 2, the time and memory to make the dense copy count as part of the method. They are also recorded separately.
- **Components:** 50 at every size, and 500 (the lab's evaluation setting) wherever it fits.
- **The lab's code:** `MatchingSparseRandomProjection`, `SphericalKMeans`, `calculate_accuracy`, `dunn_index` and `gap_statistic` are copied unchanged from `code_updated_journal_version.py` (commit `e6d0f92`), with the lab's MIT licence notice.

## 5. Sizes, seeds and runs

| Cells | 10,000 | 25,000 | 50,000 | 100,000 | 250,000 |
|---|---:|---:|---:|---:|---:|
| Seeds | 5 | 5 | 5 | 5 | 3 |

- **Grid:** 6 methods × 2 component settings × 23 (size, seed) pairs = 276 reduction runs.
- **Seeds:** the seed sets every random choice, i.e. the projection matrix, the randomized-SVD start and the clustering starts. Runs are executed in a shuffled order (seed 2026).

## 6. Measurements

**Time**
- Wall time and CPU time of the reduction step only.

**Memory**
- Each run is a fresh process, wrapped in `/usr/bin/time -l` (macOS).
- Memory is also sampled every 5 ms with `psutil` during the reduction step.
- **Headline:** the peak during the step minus the memory just before it. The whole-process peak is also reported.
- **Limits:** a run that passes 20 GiB is stopped and recorded as "did not fit"; each run also has a 60-minute limit.

**Agreement with labels**
- **Clustering:** two routes on every embedding:
  - the lab's spherical k-means, with k = number of labels;
  - Leiden on a 15-nearest-neighbour graph (Scanpy, `flavor="igraph"`), with the resolution tuned by bisection to give k clusters.
- **Metrics:** ARI, NMI, the lab's Hungarian accuracy, and raw mutual information.
- **Headline metric (fixed now):** ARI with Leiden.

**Rare populations**
- **Rare** means under 1% of the 250k cells. The list is fixed in the manifest before any method runs.
- **Best-match F1:** for each cell type, the best F1 over all clusters. It is reported only for types with at least 20 cells at that size.
- **Neighbour purity:** the share of each cell's 15 nearest neighbours with the same label, averaged per type. It uses no clustering.

**Distance preservation**
- The Spearman correlation between cell-to-cell distances before and after reduction, on 1,000 fixed cells.

**Summaries**
- Median and min–max over seeds. No p-values.
- Differences smaller than the seed range are reported as "no clear difference".

## 7. Steps that don't scale, and what replaces them

| Step | Memory at 10k | 25k | 50k | 250k |
|---|---:|---:|---:|---:|
| Pairwise distances (`pdist`, used by hierarchical clustering) | 0.4 GB | 2.5 GB | 10 GB | 250 GB |
| Full square matrix (`squareform`, used by the Dunn index), on top | 0.8 GB | 5 GB | 20 GB | 500 GB |

- **At every size:** spherical k-means and Leiden.
- **At 10k and 25k:** the lab's hierarchical clustering (Ward on cosine distances) and the lab's UMAP route, to check that both routes rank the methods the same way.
- **At 10k only (seed 0):** the lab's Dunn index and gap statistic, unchanged.

## 8. Reporting rules

- All results are published, including failed runs and results that favour PCA.
- Every number in the README comes from a file in `results/` written by a script.
- If a method doesn't fit at 250,000 cells, its results are reported up to the largest size that ran, and this is stated plainly.

## 9. Related work

- Tsuyuzaki et al. 2020, *Genome Biology* 21:9 ([doi](https://doi.org/10.1186/s13059-019-1900-3)): benchmarks PCA algorithms for large-scale scRNA-seq.
- Jindal et al. 2018 (FiRE), *Nature Communications* 9:4719 ([doi](https://doi.org/10.1038/s41467-018-07234-6)): sketching to find rare cells.
- Wan et al. 2020 (SHARP), *Genome Research* 30:205–213 ([doi](https://doi.org/10.1101/gr.254557.119)): ensemble random projection for large scRNA-seq.
- CZ CELLxGENE Discover, *Nucleic Acids Research* 53:D886–D900 ([doi](https://doi.org/10.1093/nar/gkae1142)). Census data is reusable under CC BY 4.0.

## Changes after this plan

None yet.
