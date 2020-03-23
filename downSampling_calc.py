#---------------------------------------------------------------------------
# CALIBRATED ChIP-seq ANALYSIS: DownSampling Calculations
# Emily Georgiades, March 2020.
#---------------------------------------------------------------------------

#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Load the textfile containing read counts:
read_counts = pd.read_csv('./readCounts.txt', sep=' ')
read_counts['REPLICATE'] = read_counts['SAMPLE'].str.extract('rep(\d+)', expand = False)
cols = read_counts.columns.tolist()
cols = cols[-1:] + cols[:-1]
read_counts = read_counts[cols]
# Filter for the rows containing ChIP samples (rather than inputs controls)
# Input samples must have "input" in sample name. 
sample_counts = read_counts[~read_counts['SAMPLE'].str.contains("input")]
sample_counts['RATIO_MM10_UNIQ'] = sample_counts['GENOME_READS']/sample_counts['TOTAL_READS']
sample_counts['RATIO_DM6_UNIQ'] = sample_counts['SPIKEIN_READS']/sample_counts['TOTAL_READS']
sample_counts['RATIO_DM6vMM10'] = sample_counts['SPIKEIN_READS']/sample_counts['GENOME_READS']
sample_counts['DM6_NORM'] = min(sample_counts['SPIKEIN_READS'])/sample_counts['SPIKEIN_READS'] # will need to change when have more data

# Filter for the rows containing input controls
input_counts = read_counts[read_counts['SAMPLE'].str.contains("input")]
input_counts['RATIO_MM10_UNIQ'] = input_counts['GENOME_READS']/input_counts['TOTAL_READS']
input_counts['RATIO_DM6_UNIQ'] = input_counts['SPIKEIN_READS']/input_counts['TOTAL_READS']
input_counts['RATIO_DM6vMM10'] = input_counts['SPIKEIN_READS']/input_counts['GENOME_READS']
input_counts['DM6_NORM'] = min(input_counts['SPIKEIN_READS'])/input_counts['SPIKEIN_READS'] # will need to change when have more data

# Re-merge the sample and input rows
merged = pd.merge(sample_counts, input_counts, on="REPLICATE")
new_merged = merged.iloc[:,0:9]
new_merged['INPUT_RATIO'] = merged['RATIO_DM6vMM10_y']
new_merged.columns = ['REPLICATE','SAMPLE','TOTAL_READS','GENOME_READS','SPIKEIN_READS','RATIO_MM10_UNIQ','RATIO_DM6_UNIQ','RATIO_DM6vMM10','DM6_NORM','INPUT_RATIO']

# Calculate the downsampling ratios
new_merged['DOWNSAMPLE_FRACTION'] = new_merged['INPUT_RATIO']*new_merged['DM6_NORM']
new_merged['DOWNSAMPLE_FACTOR'] = new_merged['DOWNSAMPLE_FRACTION']/max(new_merged['DOWNSAMPLE_FRACTION'])*0.9999
new_merged['READS_POST_DOWNSAMPLE'] = (new_merged['GENOME_READS']*new_merged['DOWNSAMPLE_FACTOR']).astype(int)

# Export results to textfile
new_merged.to_csv('./downsamplingCalculations.txt', sep='\t', index=False, header=None)
