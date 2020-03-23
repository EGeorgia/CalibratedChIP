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
#---------------------------------------------------------------------------
#!/bin/bash

# Specify parameters:
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage
while getopts ":g:s:b:p:" arg; do
  case $arg in
    g) # Specify genome build (e.g. mm10).
      GENOME=${OPTARG};;
    s) # Specify spike-in genome build (e.g. dm6)
      SPIKEIN=${OPTARG};;
    b) # Give path to directory containing bt2 files include prefix (e.g. ./Bowtie2_mm10.dm6/mm10.dm6)
      BT2DIR=${OPTARG};;
    p) # Give path to public directory where bigwigs will be saved.
      public_dir=${OPTARG};;
    h) # Display help.
      usage
      exit 0
      ;;
  esac
done

INPUT_FASTQS="./paths_to_fastqs.txt"
start="$(date)"

module load ucsctools
module load bowtie2/2.3.5

echo ""
echo "CALIBRATED CHIP-SEQ"
echo "Run started : " "$start"
echo ""
#---------------------------------------------------------------------------
# Step 1:
# Align paired-end FASTQ data with bowtie2.
# Uniquely mapped reads are sorted, indexed and PCR duplicates are removed
# with sambamba. Number of reads aligning from each genome are counted.
#---------------------------------------------------------------------------
echo "---------------------------"
echo "STEP 1"
echo "---------------------------"
echo "SAMPLE TOTAL_READS GENOME_READS SPIKEIN_READS" >> readCounts.txt

mkdir outputBams

for file in ${INPUT_FASTQS}
do
  columns=$(awk '{print NF}' $file | sort -nu | tail -n 1)
  read_num=$(($columns-1))

  if [[ $read_num = "2" ]]; then
    READ_TYPE="paired"
    while IFS=$'\t' read -r sampleName READ1 READ2
    do
      SAMPLE="${sampleName}"
      FASTQ1="${READ1}"
      FASTQ2="${READ2}"

      echo "Processing sample: ${SAMPLE}."
      echo ""
      echo "FASTQ1 located here: ${FASTQ1}"
      echo "FASTQ2 located here: ${FASTQ2}"

      echo ""
      echo "Aligning to concatenated genome file..."
      bowtie2 -x ${BT2DIR} -1 ${FASTQ1} -2 ${FASTQ2} -p 56 --no-mixed --no-discordant| grep -v XS: - | samtools view -b -h -S -F4 - > ./outputBams/${SAMPLE}\_UniqMapped.bam
      sambamba sort  -t 56 -m 60G -o ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped.bam
      sambamba markdup -r -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
      echo "Extracting reads aligning uniquely to ${GENOME}"
      samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${SPIKEIN} | sed s/${GENOME}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
      echo "Extracting reads aligning uniquely to ${SPIKEIN}."
      samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${GENOME} | sed s/${SPIKEIN}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam

      echo ""
      echo "Indexing bam files..."
      sambamba index -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
      echo "${SAMPLE} is indexed"
      sambamba index -t 56 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
      echo "${SAMPLE} for ${GENOME} is indexed"
      sambamba index -t 56 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam
      echo "${SAMPLE} for ${SPIKEIN} is indexed"

      echo ""
      echo "Counting reads in bam files..."
      COUNT=$(sambamba view -c -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam)
      echo "Number of total uniquely aligning reads is ${COUNT}"
      GENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam)
      echo "Number of reads aligning to ${GENOME} is ${GENOMECOUNT}"
      SPIKEGENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam)
      echo "Number of reads aligning to ${SPIKEIN} is ${SPIKEGENOMECOUNT}"
      echo "${SAMPLE} ${COUNT} ${GENOMECOUNT} ${SPIKEGENOMECOUNT}" >> readCounts.txt
      echo "---------------------------"
    done < "${INPUT_FASTQS}"

  elif [[ $read_num = "1" ]]; then
    READ_TYPE="single"
    while IFS=$'\t' read -r sampleName READ
    do
      SAMPLE="${sampleName}"
      FASTQ="${READ}"

      echo "Processing sample: ${SAMPLE}."
      echo ""
      echo "FASTQ located here: ${FASTQ}"
      echo ""
      echo "Aligning to concatenated genome file..."
      bowtie2 -x ${BT2DIR} -U ${FASTQ} -p 56 --no-mixed --no-discordant| grep -v XS: - | samtools view -b -h -S -F4 - > ./outputBams/${SAMPLE}\_UniqMapped.bam
      sambamba sort  -t 56 -m 60G -o ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped.bam
      sambamba markdup -r -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted.bam ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
      echo "Extracting reads aligning uniquely to ${GENOME}"
      samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${SPIKEIN} | sed s/${GENOME}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
      echo "Extracting reads aligning uniquely to ${SPIKEIN}."
      samtools view -h ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam | grep -v ${GENOME} | sed s/${SPIKEIN}\_chr/chr/g | samtools view -bhS - > ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam

      echo ""
      echo "Indexing bam files..."
      sambamba index -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam
      echo "${SAMPLE} is indexed"
      sambamba index -t 56 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam
      echo "${SAMPLE} for ${GENOME} is indexed"
      sambamba index -t 56 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam
      echo "${SAMPLE} for ${SPIKEIN} is indexed"

      echo ""
      echo "Counting reads in bam files..."
      COUNT=$(sambamba view -c -t 56 ./outputBams/${SAMPLE}\_UniqMapped_sorted_rmdup.bam)
      echo "Number of total uniquely aligning reads is ${COUNT}"
      GENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${GENOME}.UniqMapped_sorted_rmdup.bam)
      echo "Number of reads aligning to ${GENOME} is ${GENOMECOUNT}"
      SPIKEGENOMECOUNT=$(sambamba view -c -t 50 ./outputBams/${SAMPLE}\_${SPIKEIN}.UniqMapped_sorted_rmdup.bam)
      echo "Number of reads aligning to ${SPIKEIN} is ${SPIKEGENOMECOUNT}"
      echo "${SAMPLE} ${COUNT} ${GENOMECOUNT} ${SPIKEGENOMECOUNT}" >> readCounts.txt
      echo "Written ${SAMPLE} read counts to ./readCounts.txt"
      echo "---------------------------"
    done < "${INPUT_FASTQS}"
  else
    echo "Error! paths_to_fastqs.txt file incorrectly formatted"
  fi
done



echo ""
echo "Cleaning up..."
echo ""
rm ./outputBams/*\_UniqMapped.bam
rm ./outputBams/*\_UniqMapped_sorted.bam
rm ./outputBams/*\_UniqMapped_sorted.bam.bai

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
python ./downSampling_calc.py
echo "Running ./downSampling_calc.py to create downsamplingCalculations.txt"
echo "Downsampling ratios calculated."

DOWNSAMPLING="./downsamplingCalculations.txt"

while IFS=$'\t' read -r REPLICATE SAMPLE_NAME TOTAL_READS GENOME_READS SPIKEIN_READS RATIO_MM10_UNIQ RATIO_DM6_UNIQ RATIO_DM6vMM10 DM6_NORM INPUT_RATIO DOWNSAMPLE_FRACTION DOWNSAMPLE_FACTOR READS_POST_DOWNSAMPLE
do
  downsample_factor="${DOWNSAMPLE_FACTOR}"
  sample_name="${SAMPLE_NAME}"

  echo ""
  echo "Creating downsampled bam for ${sample_name} in ${GENOME} genome"
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
  IFS="\_" read -r  chip celltype condition genome type  <<< "$bamfile"
  sample_name="${chip}_${celltype}"
  replicate="${condition}"
  genome="${genome}"
  type="${type}"
  echo "Sample: ${sample_name}_${replicate}_${genome}"
  if [[ $READ_TYPE = "single" ]]; then
    echo "Analysing single-end reads."
    macs2 pileup -i ${bamfile} -f BAM -o ${sample_name}\_${replicate}\_${genome}\_downsampled.bg > /dev/null 2>&1
    wigToBigWig -clip ${sample_name}\_${replicate}\_${genome}\_downsampled.bg ../sampleData/${genome}.chrom.sizes ${sample_name}\_${replicate}\_${genome}\_downsampled.bw > /dev/null 2>&1
    rm ${sample_name}\_${replicate}\_${genome}\_downsampled.bg

  elif [[ $READ_TYPE = "paired" ]];then
    echo "Analysing paired-end reads."
    macs2 callpeak -t ${bamfile} -f BAMPE -g mm --bdg -n ${sample_name}\_${replicate}\_${genome} --tempdir /data/tmp/ > /dev/null 2>&1
    wigToBigWig -clip ${sample_name}\_treat_pileup.bdg ../sampleData/${genome}.chrom.sizes > /dev/null 2>&1 ${sample_name}\_${replicate}\_${genome}\_downsampled.bw
  	rm ${sample_name}\_${replicate}\_${genome}\_control_lambda.bdg
  	rm ${sample_name}\_${replicate}\_${genome}\_peaks.narrowPeak
  	rm ${sample_name}\_${replicate}\_${genome}\_peaks.xls
  	rm ${sample_name}\_${replicate}\_${genome}\_summits.bed
  	rm ${sample_name}\_${replicate}\_${genome}\_treat_pileup.bdg

  else
    echo "Error! Unrecognised read type. Should be either 'single' or 'paired'"

  fi

  echo "Copying bigwigs to public folder"
  cp *.bw ${public_dir}
  echo "Generating UCSC custom track:"
  echo "track type=bigWig name=\"${sample_name}_${replicate}_${genome}\" description=\"Calib.ChIP ${sample_name}_${replicate}_${genome}\" bigDataUrl=http://sara.molbiol.ox.ac.uk/public/egeorgia/RAPID-release/CalibratedChIP//${sample_name}_${replicate}_${genome}_downsampled.bw"
  echo ""
  echo ""

done
cd ..

echo "Copy and paste UCSC urls to view bigwigs in UCSC"
echo "If happy with replicates: merge replicate bams and create new bigwigs"
echo ""
echo "Cleaning up..."
rm ./outputBams/*\_rmdup.bam
rm ./outputBams/*\_rmdup.bam.bai
rm ./outputBams/*\_downsampled.bam.bai
echo ""
echo "RUN SUCCESSFUL!"
end="$(date)"
echo "Run finished : " "$end"
