Sex-Specific Atrial Modeling Framework

This repository contains the code and data supporting the article:

Herrera N, Ni H, Smith CER, Wu Y, Dobrev D, Morotti S, Grandi E.
Mechanistic insights into sex differences in atrial electrophysiology and arrhythmia vulnerability through sex‐specific computational models.
Journal of Physiology. 2025. doi: 10.1113/JP289425

This work builds on and extends:
H. Ni, S. Morotti, X. Zhang, D. Dobrev, E. Grandi. Integrative human atrial modelling unravels interactive protein kinase A and Ca²⁺/calmodulin-dependent protein kinase II signalling as key determinants of atrial arrhythmogenesis. Cardiovasc Res. 2023;119(13):2294–2311. doi: 10.1093/cvr/cvad118

Repository Contents
🔹 Single Cell Framework
    Code to simulate single atrial myocyte behavior.
    Main entry point: Main_Platform.m
        Runs the NH_Single_Cell model.

🔹 Population of Models
    Framework for generating and analyzing populations of atrial models.
    Includes both:
       Pre-generated data for loading/analysis.
       Scripts for creating new populations.

Key scripts:
A0_generate_perturbations.m → generates parameter perturbation matrices.
A1_Main_platform_population.m → runs the population simulations.
Linear_Regression_Analysis_APD_CaT_Features.m → plots regression coefficients for action potential duration (APD) and CaT amplitude features.
Linear_Regression_Analysis_DAD_Alternans.m → plots regression coefficients for DADs and alternans BCL thresholds.

Getting Started:
  Clone or download this repository.
  For single cell simulations, start with Main_Platform.m.
  For population studies:
      Generate perturbations with A0_generate_perturbations.m.
      Run simulations using A1_Main_platform_population.m.
      Perform regression analyses with the provided analysis scripts.

Citation

If you use this code in your work, please cite:
Herrera N, Ni H, Smith CER, Wu Y, Dobrev D, Morotti S, Grandi E. Mechanistic insights into sex differences in atrial electrophysiology and arrhythmia vulnerability through sex‐specific computational models. Journal of Physiology. 2025. doi: 10.1113/JP289425

Ni H, Morotti S, Zhang X, Dobrev D, Grandi E. Integrative human atrial modelling unravels interactive protein kinase A and Ca²⁺/calmodulin-dependent protein kinase II signalling as key determinants of atrial arrhythmogenesis. Cardiovasc Res. 2023;119(13):2294–2311. doi: 10.1093/cvr/cvad118
