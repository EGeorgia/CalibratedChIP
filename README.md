# CalibratedChIP
#### Pipeline for analysing and calibrating ChIP-seq data with spike-in genome.

##### Emily Georgiades, Hughes Group (March 2020)
***

This pipeline has been re-written using scripts from Nadya Fursova (Robert Klose Lab, Biochem)
and is designed to take fastq files as input and output downsampled bigwig files to view in UCSC.

Downsampling is achieved by using a 4% spike-in in ChIP.

***
### Version requirements :gear:
* Bowtie2 v2.1.0
* Sambamba v0.6.6
* Samtools v1.3 (using htslib 1.3)
* Python v3.7.4
* macs2 v2.0.10
* ucsctools v373

### Preliminary requirements 
#### 1. Catenated genome 

Take the two genomes of interest and rename chromosomes so that thet include species: 

```sed 's/>chr/>mm10_chr/g' /databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa > ./mm10_genome.fa```

```sed 's/>chr/>dm6_chr/g' /databank/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome.fa > ./dm6_genome.fa```

Catenate these two genomes:

```cat /databank/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa /databank/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome.fa > catenated_mm10_dm6.fa &```

Then need to build bowtie2 index:

```bowtie2-build /path/concatenated.fa output_prefix```

[See instructions on Homer webpage](http://homer.ucsd.edu/homer/basicTutorial/mapping.html)

#### 2. paths_to_fastqs.txt 

Needs to tab separated and without headers:
> sampleName  pathtoRead1 pathtoRead2

#### 3. chrom.size files

These need to be saved in the same directory as the bowtie2 indexing.

They can be downloaded from UCSC e.g. [hg19.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)

***

### SUMMARY
You should have a new directory containing the following (can sym link to scripts):
* paths_to_fastqs.txt
* downSampling_calc.py
* calibratedChIP_pipeline.sh

Run: ```$ bash calibratedChIP_pipeline.sh -g genome -s spike-in genome -b bt2_dir -i yes/no -p path/public_dir```

For help see: ```$ bash -h calibratedChIP_pipeline.sh```

