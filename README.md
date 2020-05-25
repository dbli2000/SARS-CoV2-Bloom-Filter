# SARS-CoV2-Bloom-Filter
Scripts for simulating and plotting an approach for population-scale SARS-CoV2 testing using LAMP-Seq, accompanying [preprint.](https://www.biorxiv.org/content/10.1101/2020.04.06.025635v1) 
Most model details are in Supplementary Note 1 and Supplementary Note 2 (PDFs available in repo). 

An easy-to-use interactive graphing calculator for plotting/calculating error rates only accounting for barcode loss in scenario 2 is available here: https://www.desmos.com/calculator/okym5pyunh. You can vary any parameter.

Files:
- SARS_CoV_2_Testing_as_Bloom_Filter_s_.pdf is a PDF of Supplementary Note 1, describing the single-indexed (FIP or BIP barcode) models. 
- Supplementary_Note_2.pdf is a PDF of Supplementary Note 2, describing the dual-indexed (FIP and BIP barcode) models. 
- For numerical simulation:
    - SARSCoV2barcoding.ipynb contains original numerical simulations for skewing errors across samples due to viral titer load, including plots that produced the numbers used in the preprint in supplementary note 1.
    - SARSCoV2barcoding_v2.ipynb is the most up-to-date script for numerical simulations of single-indexed models, featuring parallelized computation, and the ability to model barcode loss, sample skewing, or both.
    - SARSCoV2barcoding_v3.ipynb is the most up-to-date script for numerical simulations of dual-indexed models, featuring parallelized computation, and the ability to model barcode loss, sample skewing, and template switching.
- For plotting figures: 
    - SARS_CoV2_Scalable_Testing_Simulations_v3.m contains code for making all of the figures in the text and in Supplementary Note 1.
    - SARS_CoV2_Scalable_Testing_Simulations_v3.mlx is the live notebook version of the previous file, with a table of contents and integrated inline plots.
