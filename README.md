# Nucleotide spectrum analysis in T2T-primates XY paper
The directory includes all the scripts, inputs and outputs to conduct nucleotide spectrum analysis in T2T-primates XY paper.

# Step 1: 
count substitutions in ten branches separately for chr X and Y. 
## Script:
python MAF_to_counts.py
## Inputs: 
7species_new_wholeblock_PARfiltered_alignment.chrX.maf and 7species_new_wholeblock_PARfiltered_alignment.chrY.maf
## Outputs 
include PARfiltered.chrX_changed.variety.tc_0815.csv and PARfiltered.chrY_changed.variety.tc_0815.csv

# Step 2: 
compare the number of substitutions between X and Y. 
## Script:
Rscript compare_x_y_visualize.R
## Inputs: 
PARfiltered.chrX_changed.variety.tc_0815.csv and PARfiltered.chrY_changed.variety.tc_0815.csv
## Outputs: 
single_mutation_dot_box_0322_with_p_nopar.pdf


