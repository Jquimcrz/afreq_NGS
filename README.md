# Allele frequencies analysis from NGS data ##

Joaquim Cruz Corella <br>
Université de Lausanne <br>
November 2019 <br>

mail: joaquim.cruzcorella@unil.ch

You can read our full research article in New Phytologist:

https://doi.org/10.1111/nph.17530

### Context and background

The analysis presented here was originally conceived to investigate a type of intra-specific genetic variation in *Rhizophagus irregularis*, a model species of arbuscular mycorrhizal fungi (AMF). *R. irregularis* isolates have coenocytic hyphae and harbour hundreds of haploid nuclei in a common cytoplasm, which is thought to increase its phenotypic plasticity to cope with rapidly changing environments. Some *R. irregularis* isolates have been found to harbour two genetically distinct populations of nuclei (dikaryon isolates), based on both the presence of two different alleles of a putative MAT-locus and the distribution of allele frequencies in their bi-allelic sites from WGS data. *See the related literature below for a more detailed explanation*.

Because it was suspected that the proportions of the two different nuclear genotypes in dikaryon isolates varied among clonally produced offspring, we expected these allele frequencies distributions to shift from the expected unimodal shape when nuclear ratios are 1:1, to produce bi-modal distributions in some cases. In order to investigate this hypothesis from a reduced set of bi-allelic sites obtained by ddRAD-sequencing we performed the following analysis, improving the accuracy of our estimates by computing the pooled allele frequencies of each bi-allelic site using replicates sequencing data.

This pipeline was designed to input a collection of variant files (VCFs 4.1) from different replicates of a biological sample; process the read count information recorded at each bi-allelic site and compute, when possible, a pooled estimation of the frequency of the reference allele (REF) and its interval of confidence at 95% (IC). From this pooled estimates, we produced the histograms displayed in the publication.

* Pawlowska, T., Taylor, J. *Organization of genetic variation in individuals of arbuscular mycorrhizal fungi*. **Nature 427, 733–737 (2004)** https://doi:10.1038/nature02290
* Hijri, M., Sanders, I. *Low gene copy number shows that arbuscular mycorrhizal fungi inherit genetically different nuclei*. **Nature 433, 160–163 (2005)** https://doi:10.1038/nature03069
* Wyss, T., Masclaux, F., Rosikiewicz, P. et al. *Population genomics reveals that within-fungus polymorphism is common and maintained in populations of the mycorrhizal fungus Rhizophagus irregularis*. **ISME J 10, 2514–2526 (2016)** https://doi:10.1038/ismej.2016.29
* Ropars, J., et al. *Evidence for the sexual origin of heterokaryosis in arbuscular mycorrhizal fungi*. **Nature Microbiology volume1, Article number: 16033 (2016)** https://doi.org/10.1038/nmicrobiol.2016.33
* Masclaux, F.G., Wyss, T., Mateus-Gonzalez, I.D. et al. *Variation in allele frequencies at the bg112 locus reveals unequal inheritance of nuclei in a dikaryotic isolate of the fungus Rhizophagus irregularis*. **Mycorrhiza (2018) 28: 369**. https://doi.org/10.1007/s00572-018-0834-z
* Kokkoris V, et al. *Nuclear Dynamics in the Arbuscular Mycorrhizal Fungi*. **Trends in Plant Science Vol 25, issue 8 (2020)*. https://doi.org/10.1016/j.tplants.2020.05.002


### Pre-processing of the data

Given that *R. irregularis* nuclei are haploid, most bi-allelic sites likely represent true differences between nuclear genotypes in a dikaryon isolate. To ensure that these bi-allelic sites belong to different nuclei though, it is necessary to filter out the ones located in repeated sequences of the genome. This filtering of bi-allelic sites is done before running this pipeline. We used the BEDtools suite to remove bi-allelic sites that intersected with annotated repeated sequences of the reference genome (repeats annotation was generated with RepeatModeler/RepeatMasker).

1. **Only consider bi-allelic sites that are located in non-repeated sequences of the genome.**

### Input data

A collection of VCF files v4.1 from replicates of a given sample (rep1.vcf, rep2.vcf, ...) are given as **input**. The script is executed from the directory where the vcf files are located (Sample1), and it will process sequentially all the VCF files present in that folder. A suggested directory organization is displayed below:

>Sample1<br>
|-- rep1.vcf<br>
|-- rep2.vcf<br>
|-- rep3.vcf<br>
|-- rep4.vcf<br>

### Processing of bi-allelic sites information from VCFs

The number of reads in a bi-allelic site that support each allele are used to estimate their frequencies. Thus, VCFs are scanned to retrieve this piece of information for each replicate and build up a table (pd.dataframe) containing only bi-allelic sites that are detected in all replicates.

2. **Only include in the dataset bi-allelic sites that are detected in all replicates.**

### QC of the bi-allelic sites and filtering according to coverage

Sufficient coverage is crutial to accurately detect variants in NGS data, and even more important to estimate allele frequencies. Thus, a minimum of 20 reads (20x) is set up as an absolute threshold below which, a site would be normally disregarded. However, in order to minimize the loss of sites (and information) in our dataset, only a site is droped if it doesn't pass the threshold in at least 80% of the total number of replicates (eg. if the site has good coverage in 5/6 replicates, it will not be disregarded).

In addition, in genomic data, sudden changes or unsual values of coverage may indicate the presence of repeats or complex genomic regions. For this reason, when analysing genomic data (ddRADseq), also the sites that have coverage outside the Q1-Q3 range in more than 20% of the replicates will be also disregarded.

3. **Only sites with a coverage > 20 reads in at least 80% of the replicates will be considered.**
4. **Only sites with a coverage within the IQR in at least 80% of the replicates will be considered.**

### Pooling of reads from replicates to estimate allele frequencies at each bi-allelic site

Variation in allele frequencies is not expected among replicates and is most likely explained by PCR biases during library preparation. (This is an important assumption we make, and it should be always tested with technical replication in the design). Thus, we are only considering sites in the analysis that do have consistent frequencies among the replicates. Each bi-allelic site is tested with a chi-squared test for equal proportions among the replicates and if it yields a *p-value* below 0.05 is, therefore, disregarded.

5. **Only sites that DO NOT show statistically significant differences among replicates in their allele frequencies (*p-value* > 0.05) are pooled together in order to estimate the frequency of the reference allale more accurately.**

### Plotting allele frequencies distribution

Pooled allele frequencies within the range [0.2 - 0.8] are plotted, only if their intervals of confidence at 95%(CI) are smaller than 0.2.

6. **Only sites with a frequency above 0.2 and below 0.8 will be plotted. Extreme frequencies are likely to be due to sequencing artifacts and/or rare mutations that occurred at a very low rate. Neither of these cases would reflect varying nuclear ratios and, therefore, they are not displayed.**

7. **Only sites with an IC smaller that 0.2 will be plotted. If the frequency cannot be estimated precisely, it is not displayed to reduce the noise in the final distribution.**

### Output files

All output files will be written into a new folder (Results) in the directory where the script run. The output files will include the following:

* RAW_dataset_*samplename*.csv - A dataset of all bi-allelic sites recorded in all VCFs, with the read counts per each.
* cov_plot_*samplename*.pdf - A boxplot of distributions of coverage in the recorded bi-alleic sites per each isolate.
* chi2_plot_*samplename*.pdf - A couple of histograms of (1) the distribution of the chi-squared statistic and (2) the p-value for the comparison of allele frequencies among replicates. 
* chi2_dataset_*samplename*.csv - The dataset with the chi-squared statistic and p-value for each bi-allelis site resulting from the comparison among replicates. 
* Final_dataset_*samplename*.csv - Curated dataset after all the filters have been applied.
* distplot_*samplename*.pdf - Histogram of the pooled allele frequencies. 

***
