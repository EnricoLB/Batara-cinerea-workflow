
```markdown
# *Batara cinerea* — Integrative Taxonomic Analysis

Code supplement for:
**"Avian visual model and vocal machine learning help clarify species limits
in the Giant Antshrike (Aves, Thamnophilidae)"**
Enrico L. Breviglieri & Vagner Aparecido Cavarzere Júnior
UNESP — Instituto de Biociências de Botucatu, São Paulo, Brazil

---

## Description

This repository contains all code and workflows used in the morphological,
vocal, and machine-learning analyses comparing three subspecies of
*Batara cinerea*: *B. c. cinerea*, *B. c. argentina*, and *B. c. excubitor*.

Analyses integrate traditional bioacoustic and morphometric approaches with
deep metric learning (BioEncoder) and neural network audio embeddings
(BirdNET), providing a reproducible framework for integrative taxonomic
assessment in Neotropical birds.

---

## Repository structure

```
Batara-cinerea-workflow/
│
├── R/
│   └── Batara_cinerea_analysis.R       # Morphometric and vocal statistics
│
├── Python/
│   └── BirdNET_UMAP_SVM_analysis.ipynb # BirdNET embeddings, UMAP, SVM
│
├── HiPerGator/
│   └── SAM3_BioEncoder_workflow.md     # SAM3 segmentation + BioEncoder
│                                         training (platform-agnostic)
│
├── data/
│   ├── morfo_batara.csv                # Museum morphometric measurements
│   └── batara.csv                      # Vocal measurements (Raven Pro)
│
└── output/
    └── plots/                          # Generated figures
```

---

## Analyses

### R script — `R/Batara_cinerea_analysis.R`
Statistical analyses of morphometric and vocal data:
- Shapiro-Wilk normality tests, Welch t-tests, Mann-Whitney U tests
  (sex comparisons)
- MANOVA by subspecies, split by sex
- Hierarchical linear mixed-effects models (LME) with museum/observer
  as random effect
- Tukey post-hoc pairwise comparisons
- Linear Discriminant Analysis (LDA) with biplots — morphology and vocal

### Python notebook — `Python/BirdNET_UMAP_SVM_analysis.ipynb`
Machine-learning vocal analysis using BirdNET embeddings:
- Spectrogram generation and visualization
- BirdNET feature embedding extraction (1024-dimensional)
- UMAP dimensionality reduction (2D and 3D projections)
- SVM subspecies classifier (median accuracy: 0.99)
- Dummy classifier baseline (accuracy: 0.51)
- Confusion matrix visualization

### SAM3 + BioEncoder workflow — `HiPerGator/SAM3_BioEncoder_workflow.md`
Deep metric learning of plumage phenotype from images:
- Automated background segmentation using SAM3 (Segment Anything Model)
  via the `autodistill` framework
- BioEncoder two-stage training pipeline (EfficientNet-B5 backbone,
  Supervised Contrastive Loss)
- t-SNE projection of plumage embeddings
- GradCAM and Contrastive GradCAM visualization of discriminative
  plumage regions
- Compatible with HiPerGator HPC, local GPU machines, or Google Colab

---

## Requirements

### R (>= 4.2.0)
```r
install.packages(c("MASS", "ggplot2", "ggpubr", "dplyr", "emmeans",
                   "nlme", "tidyr", "multcompView", "multcomp",
                   "ggExtra", "ggrepel", "cowplot"))
```

### Python — BirdNET analysis (Google Colab recommended)
```python
pip install umap-learn librosa scikit-learn seaborn
```

### Python — SAM3 + BioEncoder (GPU required)
```python
pip install autodistill autodistill-sam3 roboflow inference \
            opencv-python numpy scikit-learn sam3 bioencoder
```

---

## Data

- `data/morfo_batara.csv` — morphometric measurements from museum
  specimens (MZUSP and other institutions)
- `data/batara.csv` — vocal measurements extracted in Raven Pro from
  recordings sourced from the Macaulay Library and Xeno-Canto
- Audio `.wav` files and BirdNET embedding `.txt` files are available
  upon request (file sizes exceed GitHub limits)
- Plumage images were sourced from the Macaulay Library and are not
  redistributed here

---

## Citation

If you use this code, please cite:

> Breviglieri, E. L. & Cavarzere Júnior, V. A. (2025). Avian visual model
> and vocal machine learning help clarify species limits in the Giant
> Antshrike (Aves, Thamnophilidae). *[Journal name]*.

---

## Authors

**Enrico L. Breviglieri**
enricolopesbrevi@ufl.edu
enrico.breviglieri@unesp.br
University of Florida - Florida Museum of Natural History

**Vagner Aparecido Cavarzere Júnior**
vagner.cavarzere@unesp.br
UNESP — Universidade Estadual Paulista, Instituto de Biociências de Botucatu

---

## Funding

This study was funded by the São Paulo Research Foundation (FAPESP),
Grant Nos. 2023/09512-6 and 2022/04384-7, and the Edward W. Rose
Postdoctoral Fellowship through the Cornell Lab of Ornithology.
```

