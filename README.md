# Case-control effective N for LDSC meta-analysis

This code-only module illustrates how to replace the participant-count `N` column in an already munged LDSC summary-statistics file with a SNP-specific sum of contributing cohorts' effective sample sizes. It contains **no HA GWAS data, actual cohort counts, results, credentials, or cluster paths**. The included tests use invented examples only.

For SNP *j* and its verified contributing cohorts *Cj*:

`Neff(j) = sum[k in Cj] 4 * cases(k,j) * controls(k,j) / (cases(k,j) + controls(k,j))`.

The provided implementation assumes each cohort's post-QC case/control split is constant across its contributed SNPs. If per-SNP counts vary, this implementation must be extended; do not silently substitute cohort-wide counts. It also requires an independently validated, SNP-aligned 0/1 membership file. Historical meta-analysis participant N or an `effects` string alone is not accepted as proof of cohort membership.

## What LDSC estimates

LD Score regression relates each SNP's GWAS chi-squared statistic (the square of its Z score) to its LD score, a summary of how much variation that SNP tags in an ancestry-appropriate reference panel. Under the model, polygenic signal tends to increase with LD score; the regression intercept captures inflation not explained by that LD-score relationship. The slope informs SNP-based heritability, with uncertainty estimated by block jackknife. For a binary phenotype, converting the observed-scale estimate to a liability scale requires an assumed population prevalence and a sample-prevalence parameter. These quantities are distinct. See [Bulik-Sullivan et al.](https://www.nature.com/articles/ng.3211) and the [LDSC documentation](https://github.com/bulik/ldsc).

## Inputs and safety gates

`counts.tsv` has `cohort`, `cases`, `controls`; `membership.tsv.gz` has `SNP` followed by one 0/1 column per cohort, in the same order as `counts.tsv`; and an existing LDSC file has `SNP`, `A1`, `A2`, `N`, `Z`. The code validates matching SNP order, nonmissing fields, positive counts, and a contributor at every SNP. It preserves SNP, alleles and Z exactly, changes only N, and refuses to overwrite an output. These inputs are **not included** here and must not be pushed if restricted.

Run the synthetic checks with `python3 -m unittest -v`. A permitted local analysis can then use:

`python3 compute_effective_n.py --sumstats original.sumstats.gz --membership verified_membership.tsv.gz --counts verified_counts.tsv --out effective_n.sumstats.gz`

For the effective-N LDSC analysis, use `--samp-prev 0.5` and separate runs at justified population prevalences such as `--pop-prev 0.005` and `0.01`, holding SNPs, Z scores, LD reference and other settings fixed. Compare with the original participant-N parameterization on the same SNP set. The 0.5 sample prevalence is a balanced-sample convention paired with effective N, **not** the population prevalence. See [Grotzinger et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10066905/) for the statistical rationale and the [LDSC repository](https://github.com/bulik/ldsc) for software usage.

Example command after installing LDSC separately and validating the input (replace paths; never commit those inputs or outputs):

```bash
python /path/to/ldsc/ldsc.py --h2 effective_n.sumstats.gz \
  --ref-ld-chr /path/to/ancestry_matched_ld/ \
  --w-ld-chr /path/to/ancestry_matched_ld/ \
  --samp-prev 0.5 --pop-prev 0.005 \
  --out /private/output/h2_effective_n_pop005
```

For multiple assumed population prevalences, `bash run_ldsc.sh /path/to/ldsc.py effective_n.sumstats.gz /path/to/ancestry_matched_ld/ /private/output/ 0.005 0.01` writes separate logs without overwriting existing logs. Run the original participant-N input with its documented original sample prevalence separately; this code does not manufacture a pooled prevalence.

## What is required to reproduce a real result

An investigator must provide authorized access to the exact meta-analysis summary statistics, post-QC cohort case/control counts, a validated SNP-by-cohort membership map, an ancestry-appropriate LD reference, software versions, and any justified population-prevalence assumptions. Run logs and checksums should be retained privately. Because these research inputs are not public in this repository, the tests reproduce the **algorithm**, not any unpublished cohort-specific h² estimate. Before publishing a result, compare retained SNP count, alleles and Z statistics with the original input and independently verify every count and cohort assignment.

The helper script invokes a separately installed LDSC; this repository does **not** bundle LDSC, infer cohort membership, recover case/control counts, or certify a heritability result. Do not use an SE-implied N estimate as a substitute for the documented cohort case/control formula. No license is asserted here; ownership and redistribution permissions must be confirmed before publication.
