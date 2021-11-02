#---------------------------------------------------------------------------
# CALIBRATED ChIP-seq ANALYSIS
# Re-written pipeline from Rob Klose's lab (courtesy of Nadya Fursova).

# Bowtie2 v2.1.0
# Sambamba v0.6.6
# Samtools v1.3 (using htslib 1.3)
# Python v3.7.4
# macs2 v2.0.10
# ucsctools v373

# Emily Georgiades, February 2020.
# Modifications: September 2021
#---------------------------------------------------------------------------
#!/bin/bash

# Specify parameters:
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage
while getopts ":g:s:b:p:i" arg; do
  case $arg in
    g) # Specify genome build (e.g. mm39).
      GENOME=${OPTARG};;
    s) # Specify spike-in genome build (e.g. hg38)
      SPIKEIN=${OPTARG};;
    b) # Give path to directory containing bt2 files include prefix (e.g. ./Bowtie2_mm39.hg38/mm39.hg38)
      BT2DIR=${OPTARG};;
    p) # Give path to public directory where bigwigs will be saved.
      public_dir=${OPTARG};;
    i) # Specify whether the project contains input samples (if so will be needed for each sample)
      input_sample=${OPTARG};;
    h) # Display help.
      usage
      exit 0
      ;;
  esac
done

INPUT_FASTQS="./paths_to_fastqs.txt"
start="$(date)"

module load ucsctools
module load bowtie2
module load samtools
module load sambamba
module load python-cbrg
module load deeptools

echo ""
echo "CALIBRATED CHIP-SEQ"
echo "Run started : " "$start"
echo ""
# ---------------------------------------------------------------------------
# Step 1:
# Align paired-end FASTQ data with bowtie2.
# Uniquely mapped reads are sorted, indexed and PCR duplicates are removed
# with sambamba. Number of reads aligning from each genome are counted.
# ---------------------------------------------------------------------------
# echo "---------------------------"
# echo "STEP 1"
# echo "---------------------------"
# echo "SAMPLE TOTAL_READS GENOME_READS SPIKEIN_READS" >> readCounts.txt

# mkdir outputBams

# for file in ${INPUT_FASTQS}
# do
#   columns=$(awk '{print NF}' $file | sort -nu | tail -n 1)
#   read_num=$(($columns-1))

#   if [[ $read_num = "2" ]]; then
#     READ_TYPE="paired"
#     while IFS=$'\t' read -r sampleName READ1 READ2
#     do
#       SAMPLE="${sampleName}"
#       FASTQ1="${READ1}"
#       FASTQ2="${READ2}"

#       echo "Processing sample: ${SAMPLE}."
#       echo ""
#       echo "FASTQ1 located here: ${FASTQ1}"
#       echo "FASTQ2 located here: ${FASTQ2}"

#       echo ""
#       echo "Aligning to concatenated genome file..."
#       bowtie2 -x ${BT2DIR} -1 ${FASTQ1} -2 ${FASTQ2} -p 56 --no-mixed --no-discordant| grep -v XS: - | samtools view -b -h -S -F4 - > ./outputBams/${SAMPLE}\_UniqMapped.bam
#       sambamba sort  -t 56 -m 60G -o ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped.bam
#       sambamba markdup -r -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
#       echo "Extracting reads aligning uniquely to ${GENOME}"
#       samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${SPIKEIN} | sed s/${GENOME}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
#       echo "Extracting reads aligning uniquely to ${SPIKEIN}."
#       samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${GENOME} | sed s/${SPIKEIN}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam

#       echo ""
#       echo "Indexing bam files..."
#       sambamba index -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
#       echo "${SAMPLE} is indexed"
#       sambamba index -t 56 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
#       echo "${SAMPLE} for ${GENOME} is indexed"
#       sambamba index -t 56 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam
#       echo "${SAMPLE} for ${SPIKEIN} is indexed"

#       echo ""
#       echo "Counting reads in bam files..."
#       COUNT=$(sambamba view -c -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam)
#       echo "Number of total uniquely aligning reads is ${COUNT}"
#       GENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam)
#       echo "Number of reads aligning to ${GENOME} is ${GENOMECOUNT}"
#       SPIKEGENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam)
#       echo "Number of reads aligning to ${SPIKEIN} is ${SPIKEGENOMECOUNT}"
#       echo "${SAMPLE} ${COUNT} ${GENOMECOUNT} ${SPIKEGENOMECOUNT}" >> readCounts.txt
#       echo "---------------------------"
#     done < "${INPUT_FASTQS}"

#   elif [[ $read_num = "1" ]]; then
#     READ_TYPE="single"
#     while IFS=$'\t' read -r sampleName READ
#     do
#       SAMPLE="${sampleName}"
#       FASTQ="${READ}"

#       echo "Processing sample: ${SAMPLE}."
#       echo ""
#       echo "FASTQ located here: ${FASTQ}"
#       echo ""
#       echo "Aligning to concatenated genome file..."
#       bowtie2 -x ${BT2DIR} -U ${FASTQ} -p 56 --no-mixed --no-discordant| grep -v XS: - | samtools view -b -h -S -F4 - > ./outputBams/${SAMPLE}\_UniqMapped.bam
#       sambamba sort  -t 56 -m 60G -o ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped.bam
#       sambamba markdup -r -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
#       echo "Extracting reads aligning uniquely to ${GENOME}"
#       samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${SPIKEIN} | sed s/${GENOME}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
#       echo "Extracting reads aligning uniquely to ${SPIKEIN}."
#       samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${GENOME} | sed s/${SPIKEIN}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam

#       echo ""
#       echo "Indexing bam files..."
#       sambamba index -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
#       echo "${SAMPLE} is indexed"
#       sambamba index -t 56 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
#       echo "${SAMPLE} for ${GENOME} is indexed"
#       sambamba index -t 56 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam
#       echo "${SAMPLE} for ${SPIKEIN} is indexed"

#       echo ""
#       echo "Counting reads in bam files..."
#       COUNT=$(sambamba view -c -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam)
#       echo "Number of total uniquely aligning reads is ${COUNT}"
#       GENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam)
#       echo "Number of reads aligning to ${GENOME} is ${GENOMECOUNT}"
#       SPIKEGENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam)
#       echo "Number of reads aligning to ${SPIKEIN} is ${SPIKEGENOMECOUNT}"
#       echo "${SAMPLE} ${COUNT} ${GENOMECOUNT} ${SPIKEGENOMECOUNT}" >> readCounts.txt
#       echo "Written ${SAMPLE} read counts to ./readCounts.txt"
#       echo "---------------------------"
#     done < "${INPUT_FASTQS}"
#   else
#     echo "Error! paths_to_fastqs.txt file incorrectly formatted"
#   fi
# done



# echo ""
# echo "Cleaning up..."
# echo ""
# rm ./outputBams/*\_UniqMapped.bam
# rm ./outputBams/*\_UniqMapped_sorted.bam
# rm ./outputBams/*\_UniqMapped_sorted.bam.bai

#---------------------------------------------------------------------------
# Step 2:
# Calculate downsampling factors and use these to downsample the bamfiles
# for each sample.
#---------------------------------------------------------------------------
echo "---------------------------"
echo "STEP 2"
echo "---------------------------"

echo ""
echo "Proceeding with downsampling calculations..."

# Need to specify whether input samples are included or not using -i.
python -i ${input_sample} ./downSampling_calc.py 

echo "Running ./downSampling_calc.py to create downsamplingCalculations.txt"
echo "Downsampling ratios calculated."

DOWNSAMPLING="./downsamplingCalculations.txt"

#TODO: Need to modify this for input sample example:
while IFS=$'\t' read -r REPLICATE SAMPLE TOTAL_READS GENOME_READS SPIKEIN_READS RATIO_GENOME_UNIQ RATIO_SPIKEIN_UNIQ RATIO_SPIKEINvGENOME SPIKEIN_NORM DOWNSAMPLE_FACTOR READS_POST_DOWNSAMPLE
do
  downsample_factor="${DOWNSAMPLE_FACTOR}"
  sample_name="${REPLICATE}"

  echo ""
  echo "Creating downsampled bam for ${sample_name} in ${GENOME} genome"
  echo "downsample factor is ${downsample_factor}"
  sambamba view -h -t 56 -f bam --subsampling-seed=123 -s ${downsample_factor} ./outputBams/${sample_name}\_${GENOME}.UniqMapped_sorted_rmdup.bam -o ./outputBams/${sample_name}\_${GENOME}\_downsampled.bam
  echo "Indexing downsampled bam..."
  sambamba index -t 56 ./outputBams/${sample_name}\_${GENOME}\_downsampled.bam

  DOWNSAMP_COUNT=$(sambamba view -c -t 56 ./outputBams/${sample_name}\_${GENOME}\_downsampled.bam)
  echo "Read count after downsampling is ${DOWNSAMP_COUNT}"

  echo ""
  echo "Creating downsampled bam for ${sample_name} in ${SPIKEIN} genome"
  sambamba view -h -t 56 -f bam --subsampling-seed=123 -s ${downsample_factor} ./outputBams/${sample_name}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam -o ./outputBams/${sample_name}\_${SPIKEIN}\_downsampled.bam
  echo "Indexing downsampled bam..."
  sambamba index -t 56 ./outputBams/${sample_name}\_${SPIKEIN}\_downsampled.bam

  DOWNSAMP_COUNT=$(sambamba view -c -t 56 ./outputBams/${sample_name}\_${SPIKEIN}\_downsampled.bam)
  echo "Read count after downsampling is ${DOWNSAMP_COUNT}"

done < "${DOWNSAMPLING}"
echo ""
echo "Downsampling complete!"

#---------------------------------------------------------------------------
# Step 3:
# Convert downsampled bams to bigwigs and view in UCSC.
#---------------------------------------------------------------------------
echo "---------------------------"
echo "STEP 3"
echo "---------------------------"
echo "Creating downsampled bigwigs for each sample..."
echo ""
cd outputBams

for bamfile in *_downsampled.bam
do
  echo "Creating bigwig for ${bamfile%%.*}"
  bamCoverage -b ${bamfile} -o ${bamfile%%.*}.bw
done
cd ..

# # echo "If happy with replicates: merge replicate bams and create new bigwigs"
# # echo ""
# # #echo "Cleaning up..."
# # # rm ./outputBams/*\_rmdup.bam
# # # rm ./outputBams/*\_rmdup.bam.bai
# # # rm ./outputBams/*\_downsampled.bam.bai
# # echo ""
# # echo "RUN SUCCESSFUL!"
# # end="$(date)"
# # echo "Run finished : " "$end"
