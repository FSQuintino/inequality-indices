## Overview

This page provides the codes for estimation and data fitting used in 
Vila et al. (2026). Unbiased estimation of normalized scale-invariant indices under the Gamma distribution. Preprint.

## Author Information

- **Roberto Vila**: [Scopus](https://www.scopus.com/authid/detail.uri?authorId=56924856000&origin=AuthorProfile), [Orcid](https://orcid.org/0000-0003-1073-0114),
  [Google Scholar](https://scholar.google.com/citations?hl=en&user=-u38Si8AAAAJ&view_op=list_works&sortby=pubdate)
- **Helton Saulo**: [Scopus](https://www.scopus.com/authid/detail.uri?authorId=51061286700&origin=AuthorProfile), [Orcid](https://orcid.org/0000-0002-4467-8652), [Google Scholar](https://scholar.google.com/citations?hl=en&user=okYv8sQAAAAJ&view_op=list_works&sortby=pubdate)
- **Felipe Quintino** [Scopus](https://www.scopus.com/authid/detail.uri?authorId=57604163400), [Orcid](https://orcid.org/0000-0003-0286-0541), [Google Scholar](https://scholar.google.com/citations?user=BY-KytEAAAAJ&hl=en&oi=ao)



## Repository Structure

📄 [01_plots](./01-Dataset) # rfam08 and risfam08 datasets
 

📄 [02_auxiliar_functions](./02-Functions) # Basic codes for PDF, CDF, random sample generation and parameter estimation methods for EB model

📄 [03_modeling](./03-Modelling) # Data modeling scripts (requires 01 and 02)

## Requirements

- R (version ≥ 4.5 recommended)
- Required R packages: gsl, MASS, latex2exp, AdequacyModel, goftest

## Citation

If you use these codes in your work, please cite the original paper: 
Vila, R., Saulo, H. and Quintino, F. Unbiased estimation of normalized scale-invariant indices under the Gamma distribution. arXiv (2026). 
https://doi.org/10.48550/arXiv.2606.22712

