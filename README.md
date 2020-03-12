# CalibratedChIP
###Pipeline for analysing and calibrating ChIP-seq data with spike-in genome.

###Emily Georgiades, Hughes Group
###March 2020


This pipeline has been re-written using scripts from Nadya Fursova (Robert Klose Lab, Biochem)
and is designed to take fastq files as input and output downsampled bigwig files to view in UCSC.

Downsampling is achieved by using a 4% spike-in in ChIP.


### PRELIMINARY REQUIREMENTS:
1. Catenated genome

..* Take the two genomes of interest and rename chromosomes so that thet include species: 
...`sed 's/>chr/>mm10_chr/g' /databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa > ./mm10_genome.fa`
...`sed 's/>chr/>dm6_chr/g' /databank/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome.fa > ./dm6_genome.fa`

..* Catenate these two genomes:
...`cat /databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa /databank/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome.fa > catenated_mm10_dm6.fa &`

..* Then need to build bowtie2 index:
...`bowtie2-build /path/concatenated.fa output_prefix`

[These instructions are helpful:](http://homer.ucsd.edu/homer/basicTutorial/mapping.html)

2. paths_to_fastqs.txt

..* Needs to tab separated and without headers:
...sampleName  pathtoRead1 pathtoRead2




### SUMMARY
You should have a new directory containing the following:
* paths_to_fastqs.txt
* downSampling_calc.py
* calibratedChIP_pipeline.sh

Run: `$ bash calibratedChIP_pipeline.sh -g genome -s spike-in genome -b bt2_dir -p path/public_dir`
For help see: `$ bash -h calibratedChIP_pipeline.sh`
