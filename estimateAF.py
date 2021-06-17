# import required packages for data preprocessing and statistical analysis

import os
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels import stats as stats2
import random


# import required packages for plotting

import matplotlib.pyplot as plt
import seaborn as sns


# Listing all VCF files within the current directory (from where the script has been executed)

CURRENT_DIR = os.getcwd()

SET_OF_VCF = [] 

for file in os.listdir(CURRENT_DIR):
    
    if file.endswith(".vcf"):
        
        all_path2_file=os.path.join(CURRENT_DIR, file)
        SET_OF_VCF.append(all_path2_file)

        

# Create a folder where all output files and results will be located
if not os.path.exists('Results'):
    os.makedirs('Results')
   

# Definition of two empty dictionaries that will index all bi-allelic sites common in all replicates ( 1) REF_ALT read counts, and 2) Total coverage)

hash_BA={} # This records REF and ALT reads for each allele
hash_COV={} # This records total coverage for each site


# Parse every replicate VCF file, and if condition applies, record bi-allelic sites in hash_BA

for file in SET_OF_VCF:
    
    print(file)
    
    with open(file,"r") as VCF:
        
        for line in VCF:
            
            if line.startswith("#"):
                continue
            else:
                
                ParsedLine=line.split("\t")
                GT_alleles=ParsedLine[9].split(":")[2].split(",")
                GT_reads=ParsedLine[9].split(":")[2]
                
                QUAL = float(ParsedLine[5]) #Phred-scaled quality score value
                
                if len(GT_alleles) is 2 and QUAL > 30: # Only if the variant has two alleles (ref and alt1) and QUAL>30
                  
                    REF=ParsedLine[3]
                    ALT=ParsedLine[4]
                    
                    BASES = ["A","T","C","G"]
    
                    if REF in BASES and ALT in BASES:
                            
                        # That means it is a single-base substitution
                        Site_ID = ParsedLine[0]+"_"+ParsedLine[1] #Scaffold_Position
                        Site_Alleles = REF+","+ALT
                        Site_reads = GT_reads.replace(",","_")
                            
                        if Site_ID in hash_BA:
                                
                            # Check out all conditions and if applies, record second replicate
                                
                            Site_info_0 = hash_BA[Site_ID][0] # Alleles found in the first replicate
                                
                            if Site_info_0 == Site_Alleles:
                                    
                                # Alleles are consistent among replicates, so all information is recalled in hash_BA and hash_COV
                                    
                                hash_BA[Site_ID].append(Site_reads)
                                
                                REF_int = int(Site_reads.split("_")[0])
                                ALT_int = int(Site_reads.split("_")[1])
                                COV_int = REF_int+ALT_int
                                
                                hash_COV[Site_ID].append(COV_int)
                                
                                
                        else:
                                
                            #  Bi-allelic site hasn't been yet entered in the hashes, so it is now recalled
                                
                            hash_BA[Site_ID] = []
                            hash_BA[Site_ID].append(Site_Alleles)
                            hash_BA[Site_ID].append(Site_reads)
                            
                            
                            # Total coverage at this site is recorded at the same time in hash_COV
                            
                            REF_int = int(Site_reads.split("_")[0])
                            ALT_int = int(Site_reads.split("_")[1])
                            COV_int = REF_int+ALT_int
                            
                            hash_COV[Site_ID] = []
                            hash_COV[Site_ID].append(COV_int)

                            

# Remove from the two dictionaries all the bi-allelic sites that hadn't been detected in all replicates

for key in list(hash_BA):
    
    if len(hash_BA[key]) < len(SET_OF_VCF)+1:
        del hash_BA[key]
        del hash_COV[key]
        
                                    
## hash_BA contains all the substitutions that were found to be common in all replicates and for each, its read counts REF/ALT. 
## hash_COV contains the total coverage for the same sites that have been recorded in hash_BA


# Load hash_BA and hash_COV as a Pandas dataframe: data and cov_data

# ColNames of the dataframes are defined accordingly to the names of the VCF input files
repnames = []

for file in SET_OF_VCF:
    parsedName = file.split("/")
    rep = parsedName[len(parsedName)-1].split(".")[0]
    repnames.append(rep)

colnames = ["Alleles"]+repnames

# The name of the AMF isolate being analized
sample_name = repnames[0].split("-")[0]

# Loading hash_BA information into a pd.DataFrame (data) and saving RAW_dataset.csv file (check-up point)
data = pd.DataFrame.from_dict(hash_BA, orient="index",columns=colnames)
data.to_csv('Results/RAW_dataset_'+sample_name+".csv",sep=";")


# Loading hash_COV information into a pd.DataFrame (cov_data)
cov_data = pd.DataFrame.from_dict(hash_COV,orient="index",columns=repnames)


#### Preprocessing of data ####

# Alleles column is not needed anymore (it was only used when building the hashes), so it is dropped from data
data = data.drop("Alleles",axis=1)

# Convert str from data into lists of int values [REF,ALT]
for column in data:
    data[column] = data[column].str.split("_")



## Filtering of sites according to read coverage criteria ##

# A site for a certain replicate is considered not covered sufficiently if total coverage is less than 20x #
# A site for a certain replicate is considered 'unusually' covered if they are outside the Q1-Q3 bounderies of that particular replicate #

## When a site can't pass the filter for at least 80% of the samples (replicates), it is disregarded ##



# --- Data is copied into a new pd.DataFrame where transformation of variables occurr without losing REF and ALT info --- #


# Definition of a function that determines whether the coverage of a site passes or not the filtering conditions

def score_function(value,q1,q3,minDP):
    
    score=0
    
    if value > q1 and value < q3 and value > minDP:
        score = 1
    else:
        score = 0
        
    return score



# A new np.ndarray (score_vector) of same lenght as sites in cov_data will recall the number of replicates that have passed the QC
score_vector = np.zeros((len(cov_data),),dtype=int)


# For each replicate we estimate the Q1-Q3 values and then we apply the score_function to each of their sites
for column in cov_data:
    
    Q1 = pd.Series(cov_data[column]).quantile(0.25)
    Q3 = pd.Series(cov_data[column]).quantile(0.75)
    min_reads = 20
    
    i=0
    
    for VALUE in cov_data[column]:
        
        score = score_function(VALUE,Q1,Q3,min_reads)
        
        score_vector[i]=score_vector[i]+score
        
        i=i+1
        
#print(score_vector)

# A new column 'score' is added to both DataFrames (data and cov_data) indicating number of replicates passing QC
cov_data['score']=score_vector
data['score']=score_vector

# Generation of a boxplot (cov_plot_SAMPLE.pdf) to display coverage distribution of all replicates
cov_plot = cov_data.boxplot(column=repnames, grid=False, figsize=[20,10], vert = False)
plt.savefig("Results/cov_plot_"+sample_name+".pdf")

## Only sites that have 80% of the samples with coverage > 20 reads and within Q1-Q3, are not droped from data

data = data[data.score > 0.8*(len(repnames))]
data = data.drop('score',axis=1)

# Perform a chi-squared test for each site to determine whether allele frequencies are consistent among replicates 

pvalues_prop = []
chi2_prop = []

for index,row in data.iterrows():
    
    REF_counts = []
    ALT_counts = []
    
    for rep in row:
        
        REF_counts.append(int(rep[0]))
        ALT_counts.append(int(rep[1]))
    
    if 0 in REF_counts:
        data = data.drop(index) # if a site has no reads at all for the REF allele is droped here
    else:
        OBSERVED_table = np.array([REF_counts,ALT_counts])
        chi2, p, dof, ex = stats.chi2_contingency(OBSERVED_table, correction=True, lambda_=None)
        pvalues_prop.append(p)
        chi2_prop.append(chi2)

data["pvalue"] = pvalues_prop
data["chi2"] = chi2_prop


## Generate a histogram (plot_pval_SAMPLE.pdf) of the distribution of p-values

plotA = data.hist(figsize=[20,10],grid=False, bins=20)
plt.savefig("Results/chi2_plot_"+sample_name+".pdf")


# Save data with pvalues and chi2 values into a csv file (chi2_dataset_SAMPLE.csv)
data.to_csv('Results/chi2_dataset_'+sample_name+".csv",sep=";")

# Filter out sites that have allele freq that aren't consistent among replicates
data = data[data.pvalue > 0.05]

# Drop from data the variables that are no longer needed (chi2 and pvalue)
data = data.drop('pvalue',axis=1)
data = data.drop('chi2',axis=1)

# Estimation of allele frequencies by pooling reads from all replicates and their intervals of confidence at alpha=0.05
FREQ = []
UPPER_CI = []
LOWER_CI = []

for index,row in data.iterrows():
    
    REF_counts = []
    ALT_counts = []
    
    for rep in row:
        
        REF_counts.append(int(rep[0]))
        ALT_counts.append(int(rep[1]))
    
    Total_REF = np.sum(REF_counts)
    Total_ALT = np.sum(ALT_counts)
    
    freq = Total_REF/(Total_REF+Total_ALT)
    lw_ci = stats2.proportion.proportion_confint(Total_REF,(Total_REF+Total_ALT))[0]
    up_ci = stats2.proportion.proportion_confint(Total_REF,(Total_REF+Total_ALT))[1]
    
    FREQ.append(freq)
    LOWER_CI.append(lw_ci)
    UPPER_CI.append(up_ci)

    
# A final pd.DataFrame (freq_data) is created to display Freq(REF);LowerCI;UpperCI;
freq_dic = {'Freq_REF':FREQ ,'LowerCI':LOWER_CI, 'UpperCI':UPPER_CI} 
freq_data = pd.DataFrame(freq_dic, index=data.index.tolist())

# Check if the difference between upper and lower CI is larger than 0.2 and then discard position

freq_data = freq_data[freq_data.UpperCI - freq_data.LowerCI < 0.2]
        
# Save freq_data into a csv file (Final_dataset_SAMPLE.csv)
freq_data.to_csv('Results/Final_dataset_'+sample_name+".csv",sep=";")

# Generation of a histogram + density function curve graph

freq_data = freq_data[freq_data.Freq_REF > 0.2] 
freq_data = freq_data[freq_data.Freq_REF < 0.8]
freq_data['Freq_ALT'] = 1-freq_data["Freq_REF"]
randomized_frequencies = []

for index,row in freq_data.iterrows():
    locus_freq = [row[0],row[3]]
    randomized_frequencies.append(random.sample(locus_freq,1))   

temp_df = pd.DataFrame(randomized_frequencies, columns=['Random_freq'])

fig, ax = plt.subplots(figsize=[20,10])
sns.distplot(temp_df['Random_freq'], hist=True, kde=True, 
             bins=30, color = 'black', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4,  'bw':0.02})

plt.savefig("Results/distplot_"+sample_name+".pdf")
