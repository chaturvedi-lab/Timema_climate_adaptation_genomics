#!/usr/bin/env python
# coding: utf-8

# 

# In[1]:


# ! sudo pip3 install pandas numpy


# In[2]:


import numpy as np
import os
import math, time
import pandas as pd
import sys
#from IPython.core.display import HTML
import matplotlib.pyplot as plt
from statistics import median

data_dir = "./"


# ### Read Master scaffold and determine bin count for each
# - compute the number of bins per scaffold

# In[3]:


bin_steps = 100000
bin_error = 999999
scaffold_master_file = os.path.join(data_dir, "scaffold_length.txt")
scaffold_master_dict = {}
with open(scaffold_master_file) as f:
    for i, line in enumerate(f):
        if i == 0: continue
        (key, val) = line.split()
        bin_count = math.ceil(int(val)/bin_steps)
        scaffold_master_dict[str(key)] = [bin_count, int(val)]

print("Sample Scaff: ",scaffold_master_dict[str(3)])
print("Sample Scaff: ",scaffold_master_dict[str(12)])
print("Sample Scaff: ",scaffold_master_dict[str(112)])


# ### Read Species Bayes-Factor Means (BFMEANS)

# In[4]:


prefix = sys.argv[1]
species_list = ["bart","cali", "chum", "cris", "knul", "land", "podu", "popp"]
output_csv_dict = {}
debug_scaffold_id = 31.0


# In[5]:


def read_species_bfmeans(sp_fname):
    _df = None
    if(os.path.exists(sp_fname)):
        _df = pd.read_csv(sp_fname, comment='#', sep=',', encoding="utf-8", skipinitialspace=True)
    return _df

def compute_bin(x):
    #print(type(x), (str(x) == "nan"))
    if str(x) != "nan" and (type(x) == int or type(x) == float):
        return int(math.ceil(x/bin_steps))
    else: 
        return bin_error
        
# add bins to the dataframe
def bin_species_dataframe(sp_bfmeans, sp):
    # lambda x: True if x % 2 == 0 else False
    sp_bfmeans['species'] = sp
    sp_bfmeans['bin'] = sp_bfmeans['position'].apply(compute_bin)
    print(sp_bfmeans[sp_bfmeans.scaffold == debug_scaffold_id])
    return sp_bfmeans    

def fill_row_out_data(row, sp, prev_entry=None):
    entry_idx = species_list.index(sp)
    if prev_entry is None:
        entry = [float(row.lg),float(row.scaffold), int(row.length), int(row.bin), \
                                     "miss", "miss", "miss", "miss", "miss", "miss", "miss", "miss",  \
                                     "miss", "miss", "miss", "miss", "miss", "miss", "miss", "miss",  \
                                     "miss", "miss", "miss", "miss", "miss", "miss", "miss", "miss",  \
                                     "miss", "miss", "miss", "miss", "miss", "miss", "miss", "miss",  \
                                     "miss", "miss", "miss", "miss", "miss", "miss", "miss", "miss",  ]
    else:
        entry = prev_entry

    entry[4+ 0+entry_idx] = float(row.bfmeans['mean'])
    entry[4+ 8+entry_idx] = float(row.bfmeans['min'])
    entry[4+16+entry_idx] = float(row.bfmeans['max'])
    entry[4+24+entry_idx] = int(row.bfmeans['count'])
    entry[4+32+entry_idx] = float(row.bfmeans['median'])
    return entry


# In[6]:


fig, a = plt.subplots(8,1)

for i, sp in enumerate(species_list):
    sp_fname = "{0}{1}".format(prefix, sp)
    print("-------------------------------\n",i+1, ") Read filename: {0}".format(sp_fname))
    sp_bfmeans = None
    sp_bfmeans = read_species_bfmeans(os.path.join(data_dir, sp_fname))
    sp_df = bin_species_dataframe(sp_bfmeans, sp)
    out_df = sp_df.groupby(['lg', 'scaffold', 'length', 'bin']).agg({'bfmeans': ['mean','count','min', 'max','median']}).reset_index()
    #ax = out_df.xs('bfmean').plot(ax=a[0], kind='line', x='bin')
    #display(HTML(out_df.to_html()))
 
    # Per Species and Per Scaffold
    for index, row in out_df.iterrows():
        key = "{}-{:.1f}-{}".format(int(row.lg), float(row.scaffold), int(row.bin))
        if key in output_csv_dict:
            prev_entry = output_csv_dict[key]
            output_csv_dict[key] = fill_row_out_data(row, sp, prev_entry)
        else:
            output_csv_dict[key] = fill_row_out_data(row, sp)
        
        if(float(row.scaffold) == debug_scaffold_id and sp == 'bart'):
            print(key, "----> ", sp, "\n", row)
    #break # Debug Only
    
#plt.legend(loc='best')
#plt.show()


# 
# - For every species assign a bin id based on the position 

# ### Output 
# 
# Result Table 
# 
# |scaffold| length| bin | bart(mBF)| "cali"| "chum"| "cris"| "knul"| "land"| "podu"| "popp"|
# | ----- |-----:| -----:| --------:|  ----:| -----:| -----:| -----:| -----:| -----:| -----:|
# |12     | 791103|     8|       0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|
# 
# |scaffold| length| bin | bart(min)| "cali"| "chum"| "cris"| "knul"| "land"| "podu"| "popp"|
# | ----- |-----:| -----:| --------:|  ----:| -----:| -----:| -----:| -----:| -----:| -----:|
# |       |       |     "|       0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|
# 
# |scaffold| length| bin | bart(max)| "cali"| "chum"| "cris"| "knul"| "land"| "podu"| "popp"|
# | ----- |-----:| -----:| --------:|  ----:| -----:| -----:| -----:| -----:| -----:| -----:|
# |       |       |     "|       0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|
# 
# |scaffold| length| bin | bart(frq)| "cali"| "chum"| "cris"| "knul"| "land"| "podu"| "popp"|
# | ----- |-----:| -----:| --------:|  ----:| -----:| -----:| -----:| -----:| -----:| -----:|
# |       |       |     "|     0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|    0.0|
# 

# In[7]:


print(output_csv_dict.keys())             


# In[8]:


import csv
with open(os.path.join(data_dir, sys.argv[2]), 'w', newline='') as file:
    writer = csv.writer(file, delimiter=',')
    writer.writerow(["lg", "scaffold", "length", "bin", \
                    "mbf_bart", "mbf_cali", "mbf_chum", "mbf_cris", "mbf_knul", "mbf_land", "mbf_podu", "mbf_popp", \
                    "min_bart", "min_cali", "min_chum", "min_cris", "min_knul", "min_land", "min_podu", "min_popp", \
                    "max_bart", "max_cali", "max_chum", "max_cris", "max_knul", "max_land", "max_podu", "max_popp", \
                    "frq_bart", "frq_cali", "frq_chum", "frq_cris", "frq_knul", "frq_land", "frq_podu", "frq_popp", \
                    "med_bart", "med_cali", "med_chum", "med_cris", "med_knul", "med_land", "med_podu", "med_popp"])
    writer.writerows(output_csv_dict.values())


# ## Plots

# In[ ]:






# In[ ]:




