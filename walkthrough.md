# Run public sample guide

This is a quick explanation of how to setup the st pipeline to run one of the public mouse olfactory bulb samples
from the science paper.
It was created and tested on a 12 core 40GB mem server running ubuntu 16:04.

The following procedure will place all files in your current working directory.
Feel free to modify paths to keep a nice and tidy structure while working :)

## 1. Installing the pipeline
  
### 1.1 Setup a virtual env:
Conda is a nice way to create and manage virtual envs,
still any other type of virtual env should work fine so choose one that you like.
  
If you want to go for conda you can find it here: https://conda.io/miniconda.html

Copy pasting the following to your terminal window would start the installer for the latest linux version of miniconda:
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```
After the install completes you would need to close your terminal window and open a new one.

Next step is to create the virtual env.
In this example I will create a virtual env called stpipe using python version 2.7
```
conda create -n stpipe python=2.7 --yes;
```

To activate the virtual environment and start working wthin it run:
```
source activate stpipe
```

### 1.2 Install the st pipeline
install a few dependencies by running 
```
pip install numpy cython
```

now to install the pipeline simply run
```
pip install stpipeline
```

## 2. Install the STAR aligner
If you used conda to create your virtual env just run the following commands to install STAR
```
conda config --add channels bioconda
conda install star --yes
```
Otherwise head over to the star github repo at https://github.com/alexdobin/STAR
and follow the instruction to install the software.

## 3. Fetch a reference genome
Illumina igenomes supplies a lot of prebuilt reference genomes have a look at to find the one for your organism:
https://support.illumina.com/sequencing/sequencing_software/igenome.html

As we are currently analysing mouse olfactory bulb we will fetch the mouse genome by running the following command
```
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/Ensembl/GRCm38/Mus_musculus_Ensembl_GRCm38.tar.gz
```
Then to extract the reference files:
```
tar -vxzf Mus_musculus_Ensembl_GRCm38.tar.gz
```

## 4. Build a STAR index
```
mkdir Mus_musculus/Ensembl/GRCm38/Sequence/STAR_Index
STAR \
--runMode genomeGenerate \
--genomeDir Mus_musculus/Ensembl/GRCm38/Sequence/STAR_Index \
--genomeFastaFiles Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf \
--sjdbOverhang 100 \
--runThreadN 12
```
## 5. Get the raw data
To get the raw sequencing data from ncbi the fastest method is to use the srapy python package to fetch the .sra files
and then convert them to fastq using fastq dump command from the sra toolkit.

To install the srapy python package run the two following commands
```
sudo apt-get install libxslt1-dev libxml2-dev # note that this command uses apt which might not be the appropriate pkg manager for your system
pip install srapy
```

Then to fetch the .sra file for the smallest MOB replicate run
```
get-run.py -i SRR3382385 -F {acc}.sra
```
SRR3382385 is the accession id for the run on sra.
Head over to https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=316587
to find the accessions for other MOB samples

To get the sra toolkit go to https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
and download the version appropriate for your system.
for ubuntu this would be
```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```

and then extract the compressed arcvhive
```
tar -xvzf sratoolkit.current-ubuntu64.tar.gz
```

Convert the sra file to a pair of fastq files
```
sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --split-files --origfmt --gzip SRR3382385.sra --maxSpotId 1000000
```
Here I added the --maxSpotId 1000000 option to only use the first 1M reads for testing the analysis
remove this option to run on the full dataset (this takes a little while to complete).

## 6. Fetch the id files
```
git clone https://github.com/SpatialTranscriptomicsResearch/st_pipeline.git
```

## 7. Run the st pipeline
```
mkdir -p MOB5_testrun/temp

st_pipeline_run.py \
--ref-map Mus_musculus/Ensembl/GRCm38/Sequence/STAR_Index \
--ref-annotation Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf \
--expName MOB5_testrun \
SRR3382385_1.fastq.gz \
SRR3382385_2.fastq.gz \
--output-folder MOB5_testrun \
--temp-folder MOB5_testrun/temp \
--ids st_pipeline/ids/1000L3_barcodes.txt \
--log-file MOB5_testrun/log.txt \
--mapping-threads 12 \
--verbose \
--keep-discarded-files \
--no-clean-up;
```
