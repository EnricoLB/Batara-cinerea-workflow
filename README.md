# Batara-cinerea-workflow
Workflow of the codes used to perform distinction between Batara cinerea species using BioEncoder and BirdNET

# *Batara cinerea* — Morphological and Vocal Analysis

Code supplement for the manuscript:
**"Morphological and vocal differentiation among *Batara cinerea* subspecies"**
Enrico Lopes Brevi — UNESP, Instituto de Biociências de Botucatu

---

## Description

This repository contains the R code used to perform morphological and vocal
analyses comparing three subspecies of *Batara cinerea*:
- *B. c. cinerea*
- *B. c. argentina*
- *B. c. excubitor*

Analyses include sex comparisons, multivariate tests, hierarchical linear
models, and linear discriminant analysis (LDA) with biplots.

---

## Repository structure
Batara-cinerea-workflow/
├── Batara_cinerea_analysis.R # Main analysis script
├── data/
│ ├── morfo_batara.csv # Museum morphological measurements
│ └── batara.csv # Vocal measurements (Raven Pro output)
└── output/
 └── plots/ # Generated figures saved here


---

## Dependencies

R >= 4.2.0. Install all required packages with:

```r
install.packages(c("MASS", "ggplot2", "ggpubr", "dplyr", "emmeans",
                   "nlme", "tidyr", "multcompView", "multcomp",
                   "ggExtra", "ggrepel", "cowplot"))

Analyses performed
Sex comparisons — Shapiro-Wilk normality tests, Welch t-tests,
Mann-Whitney U tests
MANOVA — Multivariate tests by subspecies, split by sex
Hierarchical LME models — Linear mixed-effects models with museum
as random effect (males only)
Tukey post-hoc tests — Pairwise comparisons for all significant
variables
LDA — Morphology — Linear discriminant analysis on wing, tail,
beak, and tarsus
LDA — Vocal — Linear discriminant analysis on seven acoustic
parameters from Raven Pro

Data sources
Morphological data were collected from museum specimens held at MZUSP
and other institutions. Vocal data were measured in Raven Pro from
recordings sourced from xeno-canto and museum collections.

Contact
Enrico Lopes Breviglieri
enricolopesbrevi@ufl.edu
enrico.breviglieri@unesp.br
