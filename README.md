# Intra-individual variability in functional brain organization

### This is the repository for the publication:
Bianca Serio, Deniz Yilmaz, Laura Pritschet, Hannah Grotzinger, Emily G. Jacobs, Simon B. Eickhoff, & Sofie L. Valk (2025). **Exploring neuroendocrine influences on the sensorimotor-association axis in a female and male individual**. _Imaging Neuroscience_. https://doi.org/10.1162/imag_a_00474.

Preprint version available [here](https://www.biorxiv.org/content/10.1101/2024.05.04.592501v1).

## Scripts

**1. Main analyses**
- `main.ipynb` computes and visualizes main analyses
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_He_AM.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Male testosterone cortisol)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_He_AM_estr_test.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Male estrogen testosterone)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_He_AM_stress.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Male stress) 
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female estrogen progesterone)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me_estr_test.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female estrogen testosterone)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me_stress.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female stress)

**2. Supplementary analyses**
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_He_suppl_test_cort_pss.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Male testosterone cortisol stress)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me_estr_test_suppl_15.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female estrogen testosterone n=15)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me_stress_suppl_20.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female stress n=20)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me_suppl_20.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female estrogen progresterone n=20)
- `p28Me_spin_permutation_lmer_within_between_network_dispersion_Me_suppl_estr_prog_pss.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses (Female estrogen progresterone stress)
- `p28_suppl_VAR.R` Estimate first- and second-order vector autoregressive (VAR) models of S-A axis loadings vs estradiol using parallel processing
- ´permTestEdgeVAR.R´ Temporal permutation testing for edgewise VAR models
- ´p28_suppl_power.R´ Power analysis

**3. Functions**
- `p28Me_myfunctions.ipynb` contains functions used for main analyses


## Data
- Open-access 28&Me sample (Pritschet et al., 2020), available on [OpenNeuro](https://openneuro.org/datasets/ds002674/versions/1.0.5).
- Open-access the 28&He sample (Grotzinger et al., 2024), available on [OpenNeuro](https://openneuro.org/datasets/ds005115/versions/1.0.0). 


## Support
Please address any questions about the analyses or code to [Bianca Serio](mailto:serio@cbs.mpg.de)

---

### Research poster presented at:
- Annual Meeting of the Organization for Human Brain Mapping (OHBM), Seoul 2024

![alt text](https://github.com/biancaserio/MC_gradients/tree/master/Poster.png?raw=true)
