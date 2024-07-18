from Bio import AlignIO
from Bio.AlignIO import MafIO
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyranges as pr
print(sys.version)
from pybedtools import BedTool

## summarize all the bases which remain identical in all of the 7 species and 6 reconstructed ancestral sequences, summarize all the bases that are different among any of the 13 sequences.
unchanged=[]
changed=[]
for alignment in list(AlignIO.parse("7species_new_wholeblock_PARfiltered_alignment.chrX.maf", "maf")):

    for r in range(1,len(alignment[11].seq)-1):
        if str.lower(alignment[:,r-1]) in('ggggggggggggg','aaaaaaaaaaaaa','ccccccccccccc','ttttttttttttt') and str.lower(alignment[:,r+1]) in('ggggggggggggg','aaaaaaaaaaaaa','ccccccccccccc','ttttttttttttt') and "-" not in alignment[:,r] and "-" not in alignment[:,r-1] and "-" not in alignment[:,r+1]:
            if str.lower(alignment[:,r]) in('ggggggggggggg','aaaaaaaaaaaaa','ccccccccccccc','ttttttttttttt'):
                unchanged.append([alignment[11].annotations["start"],r, alignment[0,r-1:r+2].seq])               
            else:
                changed.append([alignment[11].annotations["start"],r, alignment[0,r-1:r+2].seq,alignment[1,r-1:r+2].seq,alignment[2,r-1:r+2].seq,alignment[3,r-1:r+2].seq,alignment[4,r-1:r+2].seq,alignment[5,r-1:r+2].seq,alignment[6,r-1:r+2].seq,alignment[7,r-1:r+2].seq,alignment[8,r-1:r+2].seq,alignment[9,r-1:r+2].seq,alignment[10,r-1:r+2].seq,alignment[11,r-1:r+2].seq,alignment[12,r-1:r+2].seq])


with open('7species.chrX.wo.snp.unchanged.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(unchanged)


with open('/7species.chrX.wo.snp.changed.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(changed)

## translate sequence differences into substitutions in any of 12 branches
changed_df = pd.read_csv('7species.chrX.wo.snp.changed.csv', header=None)
print(changed_df)
for column in list(range(2,15)):
    changed_df[column] = changed_df[column].str.lower()
indexnames = changed_df[(changed_df.iloc[:,8]==changed_df.iloc[:,9]) & (changed_df.iloc[:,8]==changed_df.iloc[:,10]) & (changed_df.iloc[:,8]==changed_df.iloc[:,11]) & (changed_df.iloc[:,8]==changed_df.iloc[:,12]) & (changed_df.iloc[:,8]==changed_df.iloc[:,13]) & (changed_df.iloc[:,13]==changed_df.iloc[:,14]) & (changed_df.iloc[:,8]==changed_df.iloc[:,9]) & (changed_df.iloc[:,8]==changed_df.iloc[:,14]) & (changed_df.iloc[:,10]==changed_df.iloc[:,9]) & (changed_df.iloc[:,9]==changed_df.iloc[:,11]) & (changed_df.iloc[:,12]==changed_df.iloc[:,9]) & (changed_df.iloc[:,13]==changed_df.iloc[:,9]) & (changed_df.iloc[:,14]==changed_df.iloc[:,9]) & (changed_df.iloc[:,10]==changed_df.iloc[:,11]) &(changed_df.iloc[:,10]==changed_df.iloc[:,12]) & (changed_df.iloc[:,13]==changed_df.iloc[:,10]) &(changed_df.iloc[:,10]==changed_df.iloc[:,14]) & (changed_df.iloc[:,12]==changed_df.iloc[:,11]) &(changed_df.iloc[:,11]==changed_df.iloc[:,13]) & (changed_df.iloc[:,11]==changed_df.iloc[:,14]) &(changed_df.iloc[:,12]==changed_df.iloc[:,13]) & (changed_df.iloc[:,13]==changed_df.iloc[:,14]) & (changed_df.iloc[:,12]==changed_df.iloc[:,14])].index
changed_df=changed_df.drop(indexnames)
changed_df.columns = ["human_block_start", "snp_loci_in_block", "Ancestor.chrX","AncBC.chrX","AncGHBC.chrX","AncHBC.chrX","AncO.chrX","AncOGHBC.chrX","bonobo.chrX","borang.chrX","chimp.chrX","gibbon.chrX","gorilla.chrX","human.chrX","sorang.chrX"]
changed_df['ALL_OGHBC'] = changed_df['Ancestor.chrX'] + '_' + changed_df['AncOGHBC.chrX']
changed_df['ALL_gibbon'] = changed_df['Ancestor.chrX'] + '_' + changed_df['gibbon.chrX']
changed_df['OGHBC_GHBC'] = changed_df['AncOGHBC.chrX'] + '_' + changed_df['AncGHBC.chrX']
changed_df['OGHBC_O'] = changed_df['AncOGHBC.chrX'] + '_' + changed_df['AncO.chrX']
changed_df['GHBC_HBC'] = changed_df['AncGHBC.chrX'] + '_' + changed_df['AncHBC.chrX']
changed_df['GHBC_gorilla'] = changed_df['AncGHBC.chrX'] + '_' + changed_df['gorilla.chrX']
changed_df['HBC_BC'] = changed_df['AncHBC.chrX'] + '_' + changed_df['AncBC.chrX']
changed_df['HBC_human'] = changed_df['AncHBC.chrX'] + '_' + changed_df['human.chrX']
changed_df['BC_bonobo'] = changed_df['AncBC.chrX'] + '_' + changed_df['bonobo.chrX']
changed_df['BC_chimp'] = changed_df['AncBC.chrX'] + '_' + changed_df['chimp.chrX']
changed_df['O_sorang'] = changed_df['AncO.chrX'] + '_' + changed_df['sorang.chrX']
changed_df['O_borang'] = changed_df['AncO.chrX'] + '_' + changed_df['borang.chrX']
changed_transition = changed_df.iloc[:,15:27]
changed_transition.head
changed_counts = changed_transition.apply(pd.value_counts)
all_types = np.array(changed_counts.index)
all_types 
types=pd.DataFrame(all_types)

types.columns=['entries']
types[['anc','new']] = types['entries'].str.split("_",expand=True)


types[['anc0','anc1','anc2','anc3','anc4']] = types['anc'].str.split("",expand=True)
types[['new0','new1','new2','new3','new4']] = types['new'].str.split("",expand=True)
types
unchanged_from_changed_types_indexnames = types[(types['anc1']==types['new1']) & (types['anc2']==types['new2']) & (types['anc3']==types['new3'])].index
unchanged_from_changed_counts = changed_counts.take(unchanged_from_changed_types_indexnames)
real_changed_types_indexnames = types[(types['anc1']==types['new1']) & (types['anc2']!=types['new2']) & (types['anc3']==types['new3'])].index
real_changed_counts = changed_counts.take(real_changed_types_indexnames)
real_changed_counts
real_changed_counts_variety = real_changed_counts
real_changed_counts_variety
unchanged_df = pd.read_csv('7species.chrX.wo.snp.unchanged.csv', header=None)
print(unchanged_df)
unchanged_df.columns = ['blockstart','idinblock','genotype']
type(unchanged_df['genotype'])
unchanged_list = unchanged_df['genotype'].str.lower()
unchanged_counts = unchanged_list.value_counts().to_frame('counts')
new_list_unchanged = [item[0:3] for item in unchanged_from_changed_counts.index]
unchanged_from_changed_counts.index = new_list_unchanged


C = pd.DataFrame(index=unchanged_from_changed_counts.index)


for col in unchanged_from_changed_counts.columns:

    C[col] = unchanged_counts.values.flatten() + unchanged_from_changed_counts[col].values

print(C)
C.to_csv('10branches.chrX.wo.unchanged.csv',header=True, index=True)
real_changed_counts_variety.to_csv('10branches.chrX.wo.changed.variety.csv',header=True, index=True)

real_changed_counts_variety = pd.read_csv('10branches.chrX.changed.variety.csv',index_col=0)
C = pd.read_csv('10branches.chrX.unchanged.csv',index_col=0)

real_changed_counts_variety = real_changed_counts_variety.iloc[:,:12]

mapping = {idx: idx[1]+idx[5] + '_' + idx[:3] for idx in real_changed_counts_variety.index}

# Rename the index using the mapping dictionary
real_changed_counts_variety = real_changed_counts_variety.rename(index=mapping)

## consolidate 192 types of substitutions into 96 types based on the middle base

new_real_changed_counts_variety = real_changed_counts_variety.sort_index()
new_real_changed_counts_variety.iloc[:,2:12]






old_table = new_real_changed_counts_variety

# Initialize the new table with 96 rows and 10 columns
new_table = pd.DataFrame(columns=old_table.columns)

for i in range(96):
    old_row1 = i
    old_row2 = 191 - i
    
    new_row = old_table.iloc[old_row1] + old_table.iloc[old_row2]

    new_table = new_table.append(new_row, ignore_index=True)

reversed_list = old_table.index[144:192][::-1]
reversed_list
new_table.index = reversed_list.append(old_table.index[48:96])
real_changed_counts_variety_new = new_table


table_a = real_changed_counts_variety_new
table_b = C

groups = table_a.groupby(table_a.index.str[-3:])

# Calculate the sum for each group
group_sums = groups.sum()
## calculate the probability of substitutions in ten branches
# Create table C with the same dimensions as table A
table_c = pd.DataFrame(0, index=real_changed_counts_variety_new.index, columns=real_changed_counts_variety_new.columns)


# Iterate over the groups and compute the ratios
for group_name, group in groups:
    # Get the corresponding row sums from group_sums
    row_sum = group_sums.loc[group_name]
    
    # Calculate the ratios and assign them to table C
    table_c.loc[group.index] = group.div(row_sum + table_b.loc[group_name])


## save the counts of 96 types of substitutions in ten branches
real_changed_counts_variety_new.to_csv('PARfiltered.chrX_changed.variety.tc_0815.csv',header=True, index=True)

