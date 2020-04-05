# SARS-CoV2-Bloom-Filter
Scripts for simulating and plotting an approach for population-scale SARS-CoV2 testing, accompanying preprint. 
Most model details are in Supplementary Note 1 (PDF available in repo). 

An easy-to-use interactive graphing calculator for plotting/calculating error rates only accounting for barcode loss is available here: https://www.desmos.com/calculator/okym5pyunh. You can vary any parameter in scenario 2.

Files:
- SARS_CoV_2_Testing_as_Bloom_Filter_s_.pdf is a PDF of Supplementary Note 1, describing the models. 
- For numerical simulation:
    - SARSCoV2barcoding.ipynb contains original numerical simulations for skewing errors across samples due to viral titer load, including plots that produced the numbers used in the preprint.
    - SARSCoV2barcoding_v2.ipynb is the most up-to-date script for numerical simulations, featuring parallelized computation, and the ability to model barcode loss, sample skewing, or both.
- For plotting figures: 
    - SARS_CoV2_Scalable_Testing_Simulations_v3.m contains code for making all of the figures in the text and in Supplementary Note 3.
    - SARS_CoV2_Scalable_Testing_Simulations_v3.mlx is the live notebook version of the previous file, with a table of contents and integrated inline plots.
