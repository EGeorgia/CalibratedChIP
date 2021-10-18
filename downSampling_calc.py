#!/usr/bin/env python

#---------------------------------------------------------------------------
# CALIBRATED ChIP-seq ANALYSIS: DownSampling Calculations
# Emily Georgiades, March 2020.
# Updated October 2021
#---------------------------------------------------------------------------

# Import required modules:
import pandas as pd
import warnings
import re
import argparse

# Option to specify whether input samples were included
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputs", help="Project contains input samples? (yes/no)")
args = parser.parse_args()

warnings.filterwarnings('ignore')

# Load the textfile containing read counts:
read_counts = pd.read_csv('./readCounts.txt', sep=' ')

# Load the textfile containing read counts:
rep = []
cond = []
for i in read_counts['SAMPLE']:
    rep.append(re.sub('\.input$', '', i))
    if "input" in i:
        cond.append("input")
    else:
         cond.append("ChIP")
    
read_counts['REPLICATE'] = rep
read_counts['SAMPLE'] = cond
cols = read_counts.columns.tolist()
cols = cols[-1:] + cols[:-1]
read_counts = read_counts[cols]

if args.inputs == 'yes':
    print("Normalising by spike-in and input controls")
    # Filter for the rows containing ChIP samples (rather than inputs controls)
    # Input samples must have "input" in sample name.
    sample_counts = read_counts[read_counts['SAMPLE'] != "input"]

    # Calculate ratios:
    sample_counts['RATIO_GENOME_UNIQ'] = (sample_counts['GENOME_READS'])/(sample_counts['TOTAL_READS'])
    sample_counts['RATIO_SPIKEIN_UNIQ'] = sample_counts['SPIKEIN_READS']/sample_counts['TOTAL_READS']
    sample_counts['RATIO_SPIKEINvGENOME'] = sample_counts['SPIKEIN_READS']/sample_counts['GENOME_READS']
    sample_counts['SPIKEIN_NORM'] = min(sample_counts['SPIKEIN_READS'])/sample_counts['SPIKEIN_READS'] 

    # Filter for the rows containing input controls
    # Input samples must have "input" in sample name.
    input_counts = read_counts[read_counts['SAMPLE'] == "input"]

    # Calculate ratios:
    input_counts['RATIO_GENOME_UNIQ'] = input_counts['GENOME_READS']/input_counts['TOTAL_READS']
    input_counts['RATIO_SPIKEIN_UNIQ'] = input_counts['SPIKEIN_READS']/input_counts['TOTAL_READS']
    input_counts['RATIO_SPIKEINvGENOME'] = input_counts['SPIKEIN_READS']/input_counts['GENOME_READS']
    input_counts['SPIKEIN_NORM'] = min(input_counts['SPIKEIN_READS'])/input_counts['SPIKEIN_READS'] 

    # Re-merge the sample and input rows
    merged = pd.merge(sample_counts, input_counts, on="REPLICATE")
    new_merged = merged.iloc[:,0:9]
    new_merged['INPUT_RATIO'] = merged['RATIO_SPIKEINvGENOME_y']
    new_merged.columns = ['REPLICATE','SAMPLE','TOTAL_READS','GENOME_READS','SPIKEIN_READS','RATIO_GENOME_UNIQ','RATIO_SPIKEIN_UNIQ','RATIO_SPIKEINvGENOME','SPIKEIN_NORM','INPUT_RATIO']
    
    # Calculate the downsampling ratios
    new_merged['DOWNSAMPLE_FRACTION'] = new_merged['INPUT_RATIO']*new_merged['SPIKEIN_NORM']
    new_merged['DOWNSAMPLE_FACTOR'] = new_merged['DOWNSAMPLE_FRACTION']/max(new_merged['DOWNSAMPLE_FRACTION'])*0.9999
    new_merged['READS_POST_DOWNSAMPLE'] = (new_merged['GENOME_READS']*new_merged['DOWNSAMPLE_FACTOR']).astype(int)
    
    # Export results to textfile
    new_merged.to_csv('./downsamplingCalculations.txt', sep='\t', index=False, header=None)

else:
    print("Proceeding with no input normalisation")
    # Calculate ratios:
    sample_counts = read_counts[read_counts['SAMPLE'] != "input"]
    sample_counts['RATIO_GENOME_UNIQ'] = (sample_counts['GENOME_READS'])/(sample_counts['TOTAL_READS'])
    sample_counts['RATIO_SPIKEIN_UNIQ'] = sample_counts['SPIKEIN_READS']/sample_counts['TOTAL_READS']
    sample_counts['RATIO_SPIKEINvGENOME'] = sample_counts['SPIKEIN_READS']/sample_counts['GENOME_READS']
    sample_counts['SPIKEIN_NORM'] = min(sample_counts['SPIKEIN_READS'])/sample_counts['SPIKEIN_READS'] 
    sample_counts['DOWNSAMPLE_FACTOR'] = sample_counts['SPIKEIN_NORM']/max(sample_counts['SPIKEIN_NORM'])*0.9999
    sample_counts['READS_POST_DOWNSAMPLE'] = (sample_counts['GENOME_READS']*sample_counts['DOWNSAMPLE_FACTOR']).astype(int)
    
    # Export results to textfile
    sample_counts.to_csv('./downsamplingCalculations.txt', sep='\t', index=False, header=None)