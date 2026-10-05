# INCOMMON tumour-type merges: a detailed check

Oct 5, 2026 · Ran Tene

## Summary

The INCOMMON study on MSK-MET 2021 joins 32 OncoTree codes into 9 tumour groups before its survival analysis. In four groups the merge changes published results, because the joined tumour types differ in both mutation rate and survival.

| Group | Model | Published | Main code only | What the merge did |
| --- | --- | --- | --- | --- |
| PAAD + PANET | KRAS, high dosage | HR 3.19, p 1.4e-17 | HR 1.73, p 0.0006 | WT group is 62% neuroendocrine tumours |
| OV (4 histologies) | TP53, all classes | HR 1.20–1.26, not significant | HR 2.60–2.79, p 0.03–0.04 | WT group is 81% non-high-grade serous |
| MEL (9 codes) | NRAS, balanced | HR 1.30, p 0.073 | HR 1.85, p 0.00053 | 42% of samples are not skin melanoma |
| BRCA (IDC + ILC) | CDH1 | Balanced HR 1.37, p 0.042 | Not testable (28 mutants) | Mutants are 91% lobular: compares histologies |

For the first three, a control that removes the same number of random samples does not reproduce the change. So they come from which samples were merged, not from sample size. The other merges (CRC, BLCA, LUAD, UCEC) change results little.

Two more points:

- 220 of 227 published models in merged groups reproduce exactly from the authors' released files.
- The preprint does not mention the merges. Its breast PIK3CA result (HR 1.36) matches a ductal-only analysis, while the released code merges ductal and lobular cases (HR 1.21).

## Study and data

The study is INCOMMON (Calonaci et al., Nature Genetics 2026; preprint medRxiv 2024.05.13.24307238). It classifies each mutation in a tumour sample by gene dosage, then tests whether dosage predicts survival.

**Data.** MSK-MET 2021 (cBioPortal `msk_met_2021`): 25,775 samples, one per patient, sequenced with MSK-IMPACT. Each sample has an OncoTree code (`ONCOTREE_CODE`), survival time (`OS_MONTHS`) and status (`OS_STATUS`).

**Dosage classes.** INCOMMON uses the variant allele frequency and tumour purity to give each mutant sample one of three classes:

| Class | Short name | Meaning |
| --- | --- | --- |
| Low dosage | PWK | Few mutant copies relative to wild-type copies |
| Balanced | BK | Equal mutant and wild-type copies |
| High dosage | PMK | More mutant copies, often by loss of the wild-type copy |

The reference group (WT) is every sample of the same tumour group that has no mutation in the gene.

**Survival model.** For each tumour group and gene with at least 100 mutant samples, the authors fit a Cox model (Efron ties):

```
Surv(OS_MONTHS, OS_STATUS) ~ class + age + TMB (<10 / >=10) + sample type + FGA (+ sex)
```

Sex is left out for sex-specific cancers. The published result file `survival_analysis_msk_met.rds` holds 250 such models (fit type `general`). 227 of them are in tumour groups made by a merge (next section).

## The merge code

One step in `1.1.prepare_msk_met.R` (lines 84–98) replaces the OncoTree code with a broader group. The code comment calls this "Collapse related OncoTree codes into broader groups". All later steps see only the new group name.

| New group | OncoTree codes joined (MSK-MET 2021 samples) | What is joined |
| --- | --- | --- |
| PAAD | PAAD 1,779 · PANET 211 | Pancreatic adenocarcinoma with pancreatic neuroendocrine tumour |
| MEL | SKCM 699 · MUP 165 · CSCC 105 · UM 103 · ARMM 71 · HNMUCM 48 · VMM 48 · URMM 5 · ESMM 3 | Skin melanoma with melanoma of unknown primary, cutaneous squamous cell carcinoma, uveal and mucosal melanomas |
| HCC | HCC 205 · IHCH 407 · GBC 121 · EHCH 119 · GBAD 50 | Liver cancer with intra- and extrahepatic cholangiocarcinoma and gallbladder cancers |
| OV | HGSOC 997 · CCOV 96 · LGSOC 89 · CCBOV 1 | High-grade serous with clear cell and low-grade serous ovarian cancer |
| BRCA | IDC 2,236 · ILC 373 | Ductal with lobular breast cancer |
| UCEC | UEC 835 · USC 288 · UCS 192 | Endometrioid with serous endometrial cancer and carcinosarcoma |
| CRC | COAD 2,312 · READ 781 · COADREAD 455 | Colon with rectal cancer |
| LUAD | LUAD 4,064 · LUNE 84 | Lung adenocarcinoma with lung neuroendocrine tumour |
| BLCA | BLCA 961 · UTUC 197 | Bladder with upper tract urothelial cancer |

The survival script `4.survival_analysis.R` later splits three groups again by `SUBTYPE_ABBREVIATION`: BRCA by receptor status, UCEC and CRC by molecular subtype. The UCEC split puts the serous and carcinosarcoma cases back into their own strata, so it undoes most of that merge. The BRCA split does not: it removes the IDC/ILC label first, so ductal and lobular cases stay together in each receptor subtype. LUAD is split by KRAS and EGFR status, which does not separate LUNE.

The same script also joins mutations: all KRAS G12 variants become one group `p.G12`, and all G13 variants become `p.G13`.

## Main case: KRAS in PAAD+PANET

In MSK-MET 2021, the published KRAS model for pancreatic cancer compares mostly PANET patients with PAAD patients. Its high-dosage hazard ratio falls from 3.19 to 1.73 when only PAAD samples are used.

**Reproduction.** I took the published model object and its own input data from `survival_analysis_msk_met.rds`. A refit gives the same coefficients, covariance, log-likelihood, N and deaths (tolerance 1e-10). I then joined each sample to its original OncoTree code in the MSK-MET 2021 sample file.

**Who is in each class** (samples in the fit):

| Class | PAAD | PANET | PANET share |
| --- | --- | --- | --- |
| WT (reference) | 82 | 132 | 62% |
| Low dosage | 533 | 5 | 0.9% |
| Balanced | 663 | 3 | 0.5% |
| High dosage | 346 | 4 | 1.1% |

About 90% of pancreatic adenocarcinomas carry a KRAS mutation. Most pancreatic neuroendocrine tumours do not. So the KRAS wild-type group in the merged data is mostly PANET, a different disease with longer survival. The model then reads "not PANET" as "KRAS mutant".

**Result** (MSK-MET 2021, KRAS, PAAD group; HR vs WT, 95% CI):

| Term | Published (PAAD+PANET) | PAAD only |
| --- | --- | --- |
| N (deaths) | 1,768 (1,088) | 1,624 (1,035) |
| Low dosage | 2.05 (1.57–2.67), p 1.5e-7 | 1.19 (0.88–1.62), p 0.26 |
| Balanced | 2.29 (1.77–2.97), p 4.5e-10 | 1.32 (0.98–1.79), p 0.069 |
| High dosage | 3.19 (2.45–4.17), p 1.4e-17 | 1.73 (1.26–2.37), p 0.0006 |
| High vs balanced | 1.39, p 4.9e-5 | 1.31, p 0.0013 |
| FGA (covariate) | 0.91 (0.64–1.29), p 0.59 | 2.46 (1.57–3.85), p 8e-5 |

With PAAD only, low and balanced dosage are no longer significant. High dosage stays significant, with about half the effect. The dosage trend (high vs balanced) remains. The FGA covariate changes most. A likely reason: PANET has fewer copy-number changes than PAAD, so in the merged data FGA partly tracked the tumour type.

**Is it only the smaller sample?** No. I removed 144 random samples (the same number as PANET) from the merged data, 300 times. The high-dosage HR stayed at 3.20 (90% range 3.02–3.43), and every term stayed significant in all 300 runs. The drop to 1.73 comes from who was removed, not from how many.

**Same pattern in other PAAD models.** KRAS `p.G12` (N 1,631) goes from 3.27 to 1.73 for high dosage, with the same WT group (132 of 214 PANET). TP53 in PAAD shrinks for all three classes (high dosage 1.78 → 1.27, p 1.6e-9 → 0.019); there the WT group is 26% PANET against 2% of mutants. The random-drop control kept TP53 at 1.78 (1.69–1.88).

## Other merges

I refit all 227 published MSK-MET 2021 models in merged groups. Each refit keeps only the main OncoTree code of its group. 220 models reproduce exactly; the 7 that do not are all in the UCEC hypermutated stratum (see Limits).

A refit on fewer samples loses power, so a lost p < 0.05 alone means little. The test that matters is the random-drop control: remove the same number of random samples and see if the result moves the same way.

**Effect on published significance** (class terms, p < 0.05, reproduced models):

| Group | Models | Median N | Median removed | Significant terms (published) | Lost | Gained |
| --- | --- | --- | --- | --- | --- | --- |
| PAAD | 7 | 1,768 | 144 | 13 | 6 | 0 |
| MEL | 50 | 1,137 | 478 | 16 | 10 | 7 |
| BRCA | 16 | 1,567 | 339 | 9 | 5 | 1 |
| HCC | 4 | 752 | 577 | 3 | 2 | 0 |
| OV | 1 | 978 | 124 | 0 | 0 | 3 |
| CRC | 72 | 346 | 73 | 31 | 9 | 4 |
| BLCA | 24 | 1,083 | 179 | 7 | 1 | 4 |
| LUAD | 33 | 1,294 | 61 | 31 | 2 | 0 |
| UCEC | 6 | 511 | 0 | 1 | 0 | 0 |

UCEC counts only the endometrioid and hypermutated strata. Its serous and carcinosarcoma strata hold no UEC samples, so they have no refit.

### OV: TP53

In MSK-MET 2021, the merged OV group hides a TP53 effect. Almost all high-grade serous cancers carry TP53 mutations. So the WT group is mostly clear cell and low-grade serous cancer: 110 of 136 (81%), against 14 of 842 mutants (2%).

| Term | Published (OV merged) | HGSOC only |
| --- | --- | --- |
| Low dosage | 1.26 (0.84–1.87), p 0.26 | 2.79 (1.12–6.98), p 0.028 |
| Balanced | 1.20 (0.83–1.73), p 0.32 | 2.60 (1.06–6.42), p 0.038 |
| High dosage | 1.22 (0.84–1.78), p 0.29 | 2.66 (1.09–6.54), p 0.032 |

Random drop of 124 samples: the HR stays near 1.22, and at most 0.7% of 300 runs reach p < 0.05. So the published null result comes from the mix of histologies. The HGSOC-only result is weak: only 26 HGSOC samples are TP53 wild-type.

### BRCA: CDH1

In MSK-MET 2021, CDH1 mutants in the HR+/HER2− BRCA stratum are 91% lobular (296 of 324). WT samples are 3.5% lobular. CDH1 loss is the defining event of lobular cancer. So this model compares lobular with ductal cancer more than it compares CDH1 dosage. The published balanced term is HR 1.37 (1.01–1.87), p 0.042. A ductal-only refit keeps only 28 mutants (HR 1.57, 0.70–3.53). The question cannot be answered inside one histology with these data.

### MEL: NRAS and TP53

In MSK-MET 2021, 478 of 1,137 samples in the MEL models are not skin melanoma. They are 98 cutaneous squamous cell carcinomas, 85 uveal and 146 mucosal melanomas, and 149 melanomas of unknown primary.

- **NRAS.** WT is 47% non-SKCM, mutants 25%. Balanced: 1.30 (0.98–1.74), p 0.073 → 1.85 (1.31–2.62), p 0.00053 with SKCM only. High dosage: 1.65 → 2.16. Random drop: balanced HR 1.29 (90% range 1.06–1.61), significant in 20% of runs. The SKCM value is outside that range.
- **TP53.** 86 of 326 TP53 mutants are squamous cell carcinomas, not melanomas. Balanced: 0.93, p 0.65 → 0.58 (0.35–0.95), p 0.03 with SKCM only. Random drop: HR 0.93, significant in 0.7% of runs.
- **TERT, BRAF V600E.** WT is 71% and 52% non-SKCM, against 27% and 20% of mutants. No significance changed in these models.

Six MEL refit terms cannot be estimated (no deaths in one class after the restriction).

### HCC: TERT and TP53

In MSK-MET 2021, 577 of 752 samples in the HCC models are bile duct or gallbladder cancers, not liver cancer. TERT promoter mutants are 78% liver cancer, while the WT group is 89% biliary. TP53 mutants are 78% biliary. TP53 loses significance with HCC only (high dosage 1.66, p 0.0018 → 1.41, p 0.37). But only 175 samples remain, and the random-drop control loses significance in 66–74% of runs. So the HCC data are too small to tell a composition effect from a power loss.

### CRC, BLCA, LUAD and UCEC

These merges join close diseases, and the results change little. CRC joins colon and rectal cancer; most changed p-values sit near 0.05. For BRAF V600E in CRC MSS, high dosage goes 2.23, p 0.0042 → 1.88, p 0.061, and the random drop is significant in only 70% of runs. LUAD removes 61 LUNE samples; the largest HR change is 0.90 → 1.11 (ARID1A, high dosage, not significant either way). BLCA removes 179 UTUC samples; 4 terms gain significance and 1 loses it, all near 0.05. The UCEC subtype split undoes that merge.

## What the paper says

The preprint (medRxiv 2024.05.13.24307238) does not mention the merges. It has no word for PANET, OncoTree or a grouping of tumour types. It names merged groups by their main member:

- "KRAS in pancreatic adenocarcinoma (HR = 3.19, CI = 2.45–4.17)". This is the PAAD+PANET model; PAAD alone gives 1.73.
- "PIK3CA (HR = 1.36, CI = 1.07–1.74) in HR+HER2− invasive ductal breast cancer". The released result file gives 1.21 (0.97–1.51) for this model, which includes 175 lobular PIK3CA mutants. My ductal-only refit gives 1.36 (1.07–1.74), the same as the preprint. So the preprint text matches a ductal-only analysis, while the released code and results merge ductal and lobular cases.
- "TP53 … in hepatocellular carcinomas" and "in melanoma" (metastasis analysis). In the MSK-MET 2021 survival models, the HCC group is 77% biliary cancers and the MEL group is 42% not skin melanoma.

Other preprint numbers match the released merged file: BRAF 2.55 and NRAS 2.02 in MSS colorectal cancer, PIK3CA 1.58 and TERT 1.55 in bladder cancer.

I have not checked the Nature Genetics 2026 version or its supplementary tables. They may describe the grouping.

## Limits

- **Classes are not recomputed.** I kept each sample's published dosage class. INCOMMON sets its priors per tumour type, and the merge changes those priors too. A full fix would rerun the classification on the original codes. For KRAS in PAAD, PANET adds only 12 of about 1,560 mutant samples (<1%), so the prior effect there is likely small. I did not check the other groups.
- **One code per group.** Each refit keeps only the main code (MEL → SKCM, OV → HGSOC, BRCA → IDC, UCEC → UEC, CRC → COAD, HCC → HCC). Other choices are valid, for example a separate model per code, or tumour type as a covariate. I did not test them.
- **Power.** Every refit has fewer samples. The random-drop control was run for 7 models only (PAAD KRAS and TP53, OV TP53, MEL NRAS and TP53, HCC TP53, CRC MSS BRAF V600E). For the other models, a lost or gained p < 0.05 may be power alone.
- **Seven models not reproduced.** All are in the UCEC hypermutated stratum (CTCF, FAT1, KMT2D, PIK3CA, POLE, PTEN, ZFHX3). That stratum has only UEC samples, so the merge does not touch it. I did not find the cause.
- **Uncorrected p-values.** All p-values here are per term, without the paper's FDR correction.
- **Only survival models.** I did not test the metastasis and organotropism analyses, which use the same merged groups.

## Methods and reproduction files

**Inputs** (all public):

- `survival_analysis_msk_met.rds`: the authors' released survival results (250 general models). Each model's own input data sit inside the model object and are read with `eval(fit$call$data, envir = environment(fit$terms))`.
- `data_clinical_sample.txt` from cBioPortal `msk_met_2021`: the original `ONCOTREE_CODE`, joined to each model's samples by `SAMPLE_ID`.
- The authors' scripts `1.1.prepare_msk_met.R` and `4.survival_analysis.R`.

**Steps** (R, `survival` package):

1. Refit each published model on its own data with the same formula and ties method. Count it as reproduced when coefficients match to 1e-10 and N and deaths are equal.
2. Add the original OncoTree code to each sample. Count codes by dosage class among the samples used in the fit.
3. Refit the same formula on samples with the group's main code only.
4. For 7 models, remove the same number of random samples 300 times (seed 1) and refit each time.

**Files** (in this repository):

| File | Content |
| --- | --- |
| [`incommon-exact-sensitivity.R`](incommon-exact-sensitivity.R) | KRAS in PAAD: exact reproduction and PAAD-only refit |
| [`allmerges2.R`](allmerges2.R) | All 227 merged-group models: reproduction, composition, main-code refit |
| [`randdrop.R`](randdrop.R) | Random-drop control |
| [`incommon-allmerges.csv`](incommon-allmerges.csv) | One row per model: N, composition, HR, CI and p, published and refit |
| [`incommon-allmerges-composition.csv`](incommon-allmerges-composition.csv) | Samples per original code and class, per model |
| [`incommon-randdrop.csv`](incommon-randdrop.csv) | Random-drop results |
| [`incommon-allmerges-scored.csv`](incommon-allmerges-scored.csv) | Same as `incommon-allmerges.csv`, one column per class term (published and refit) |
| [`incommon-exact-sensitivity.csv`](incommon-exact-sensitivity.csv) | KRAS in PAAD: published and PAAD-only terms |
| [`incommon-exact-sensitivity-input.csv`](incommon-exact-sensitivity-input.csv) | KRAS in PAAD: model input with the original OncoTree code |

The scripts read their inputs from `/private/tmp`. `incommon-survival.rds` is the released `survival_analysis_msk_met.rds`. `incommon-paper-sample.txt` is the cBioPortal `data_clinical_sample.txt`. `incommon-exact-kras-model.rds` and `incommon-exact-kras-input.rds` are the KRAS PAAD model and its data, taken from the released file. Change the paths to run them elsewhere.

**Sources**

1. Calonaci N. et al. INCOMMON. Nature Genetics 2026. Preprint: [medRxiv 2024.05.13.24307238](https://www.medrxiv.org/content/10.1101/2024.05.13.24307238).
2. Nguyen B. et al. Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients. Cell 2022. Data: [cBioPortal msk_met_2021](https://www.cbioportal.org/study/summary?id=msk_met_2021).
