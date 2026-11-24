# Random projection vs. PCA for single-cell RNA-seq at up to 250,000 cells

**Status: plan only.** No code has been written and no data has been downloaded yet. The rules are fixed in advance in [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md).

**What this is:** I extend the Moussa Lab's public benchmark of random projection and PCA for single-cell RNA-seq ([Abdelnaby & Moussa, J. Comput. Biol. 2026](https://doi.org/10.1177/15578666261436821); code: [moussa-lab/BenchmarkingPCA-RP](https://github.com/moussa-lab/BenchmarkingPCA-RP)). The original benchmark goes up to 10,000 cells. This one uses one public blood dataset from the CELLxGENE Census, at 10,000 to 250,000 cells, with the data kept sparse.

For each method and size it measures:
- run time;
- peak memory;
- agreement with the authors' cell-type labels;
- whether rare cell types are still found.

**Methods:**
- PCA (full SVD);
- PCA (randomized SVD);
- PCA (covariance method, which works on sparse data);
- Gaussian random projection;
- sparse random projection;
- the lab's Matching Sparsity random projection.

**Credit:** the Matching Sparsity random projection class and the lab's clustering and metric functions will be copied unchanged from the lab's repository (MIT licence, © 2024 moussa-lab).

Results, figures and limitations will be added here once the runs are done.
