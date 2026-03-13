# Fama-French 5 Factors: PCA & MDS Geometry (Monthly)

This project analyzes the internal geometry of the Fama-French 5-factor model using dimensionality reduction, then tests whether latent factors are useful in portfolio-level asset pricing. The goal is not only statistical fit, but a practical answer to a real question: when does factor compression remain economically meaningful?

## Why this matters

In empirical finance, factors are often treated as independent building blocks, while in practice they are correlated and partially redundant. This project shows how to map that hidden structure and convert it into a compact, interpretable representation for risk analysis.

## Methodology

- Monthly Fama-French 5-factor data (Mkt-RF, SMB, HML, RMW, CMA)
- PCA and rotated PCA (Varimax) for latent linear structure
- Classical MDS and SMACOF MDS for distance-preserving geometry
- Procrustes analysis for PCA vs MDS alignment
- Portfolio validation on 25 size-value portfolios (5x5):
  - Model A: latent PCA factors (PC1-PC3)
  - Model B: original FF5 factors
  - Comparison via R2 and adjusted R2

## Key findings

- The first 3 principal components explain about 82% of total variance.
- 2D SMACOF stress is 0.165, indicating useful but not perfect 2D representation.
- The PCA latent-factor model retains about 91.9% of FF5 R2 on average.
- In 60% of portfolios, PCA reaches at least 95% of FF5 adjusted R2.
- Practical conclusion: PCA is a strong compact layer for structure discovery and exposure mapping, while full FF5 remains preferable when maximizing in-sample fit is the priority.

## Report and data

- RPubs report: https://rpubs.com/timothypawelczyk/1409900
- Data source (Kenneth R. French Data Library): https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
