# cacoa example notebooks (dev_lm)

Example analyses with the model interface of the `dev_lm` branch of [cacoa](https://github.com/kharchenkolab/cacoa):
model formulas with covariates and grouping terms `(1 | g)`, expression shifts (gsd significance, distance-model
magnitude), composition (lmCoda), density and expression-shift maps, differential expression, and checks of the model
(design check, residual sample maps, fit diagnostics, sensitivity to the covariate set, jackknife intervals).

Each folder holds the notebook (`<name>_cacoa.Rmd`), its rendered output (`<name>_cacoa.md` with the figures in
`<name>_cacoa_files/`) and the scripts that prepare its input (`prepare_*.R`; `prepare_harmony_common.R` is shared). The
input paths are set in each notebook's `params` and point to the lab server.

| folder | data | what the model comparisons show |
|---|---|---|
| `alz` | Alzheimer's disease nuclei, 6 AD and 6 control samples | the 10x batch is nested in the condition; exchanging whole batches `(1 | batch)` leaves 20 arrangements and no significant result |
| `ms` | multiple sclerosis lesions | sequencing batch as a grouping term; lesion stage as a three-level factor (whole-factor test and two contrasts) |
| `pf` | pulmonary fibrosis (IPF) lungs | age is nearly aliased with the diagnosis; what an age adjustment would rest on |
| `scc` | squamous cell carcinoma, tumor and normal skin per patient | pairing by patient `(1 | patient)`: relabeling within patients, random intercept, patient as a fixed factor |
| `natri` | pulmonary fibrosis lungs of two centers (Natri et al. 2024) | `(1 | cohort)`; age adjustment with and without restriction to the common age range; IPF against other ILD |
| `sle` | lupus PBMC (Perez et al. 2022) | `(1 | cohort)`; flare samples confined to one processing cohort; the SLE effect within each ancestry |
| `onek1k` | OneK1K PBMC, 200 donors | age as a continuous tested variable per decade, `(1 | pool)`; continuous age against a median split; pool as a random intercept against a fixed factor in differential expression |

Package state of the rendered outputs: `scc`, `natri`, `sle`, `onek1k` on commit 7ae483b of `dev_lm`; `alz`, `ms`, `pf`
on 22ca2bf. The later commit changes two things: the stated smallest p-value of designs with as many samples of each
compared group (in `alz` the model report prints 0.0011 where the attainable value is 0.0022, as its text says), and
differential expression with a grouping term, which none of these three uses. The Harmony companion notebooks
(`*_cacoa_harmony.Rmd`: the same analyses on a Harmony kNN graph of the same cells) were rendered earlier, on a4a2f9f.
