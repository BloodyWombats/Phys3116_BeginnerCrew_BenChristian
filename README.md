# Faber-Jackson relation with SAMI DR3

This repository explores the **Faber-Jackson relation** using data from the **SAMI Galaxy Survey (DR3)** and the GAMA input catalogues. The goal is to reproduce and compare the scaling relation between luminosity and stellar velocity dispersion, and identify the difference between this and the classic literature value

## Data
## Data
You will need the following data tables from SAMI DR3 (download links available from the [SAMI DR3 data release page](https://datacentral.org.au/services/download/)):

- `samiDR3VisualMorphology.csv` – visual morphological classifications  
- `samiDR3StelKin.csv` – stellar kinematics, including velocity dispersions  
- `InputCatGAMADR3.csv` – host photometry and stellar masses from GAMA  
- `samiDR3InputCatClusters.csv` – cluster membership and environment info  

Put these files into the `data/` directory.

---

## How to run

1. Clone this repository:
   ```bash
   git clone https://github.com/<your-username>/sami-faber-jackson.git
   cd sami-faber-jackson
