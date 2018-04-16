


# R package *cntmf*
# Prioritisation of single nucleotide variants using Corrected non-negative matrix tri-factorisation


A data fusion framework for prioritising reliable associations between single nucleotide variants (SNVs) and traits.

This algorithm allows for studying the effect of SNVs on categorical traits, thanks to its main features :

1. It captures the interrelatedness between variants data, the SNVs deleteriousness effect and the protein-protein interactions (PPIs) that might be disrupted.  

2. It simultaneously accounts for the patient's outcome and ancestry by means of kernels functions, minimizing the confounding for population structures.


The cnmtf package provides four categories of functions for preprocessing data, clustering, scoring SNVs and comparing results.

* **Preprocessing functions:**

These functions will help you to create the inputs for the algorithm.

* **Factorisation functions:**

Main functions to integrate the input data, generate the low-dimmensional matrices and find consensus clusters.

* **Scoring functions:**

A set of functions to score SNVs and prioritise significant SNV-trait associations from the low-dimmensional matrices.

* **Comparing functions:**

Auxiliary functions to compare your results across different settings of the algorithm.


Author: *Luis G. Leal, lgl15@imperial.ac.uk*


