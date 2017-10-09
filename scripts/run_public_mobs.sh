#!/usr/bin/env bash

# to speed up the fetching of reads you can prefetch the sra cahces using the following command:
conda create -n SRApy python=2.7 --yes;
source activate SRApy
mkdir -p /data/ncbi/public/sra/
ln -s /data/ncbi ~/
cd /data/ncbi/public/sra/
pip install srapy
get-project-sras.py -p 316587 -F {acc}.sra &
source deactivate

# build the reference genome
genomeDir=~/references/Mus_musculus/Ensembl/GRCm38/Sequence/STAR_Index
gtfFile=~/references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
mkdir ~/references
cd ~/references
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/Ensembl/GRCm38/Mus_musculus_Ensembl_GRCm38.tar.gz
tar -xvzf Mus_musculus_Ensembl_GRCm38.tar.gz 
STAR --runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles ~/references/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile $gtfFile \
--sjdbOverhang 100 \
--runThreadN 4 &

testfolder=~/test_pipeline
cp `basename "$0"` $testfolder/the_script_ran_at_$(date +"%y%m%d_%H_%M").sh
# setup the analysis envs
mkdir -p $testfolder/repos
cd $testfolder/repos
git clone https://github.com/elhb/st_pipeline.git
cp -vr st_pipeline st_pipeline_parallel &
cp -vr st_pipeline st_pipeline_coordparse &
cp -vr st_pipeline st_pipeline_anno_par &
cd st_pipeline/repos

pipe_versions=( \
    pipe_master \
    pipe_coordparse \
    pipe_151 \
    pipe_parallel \
    pipe_anno_par )

# setup pipeline versions
for pipeline_version in "${pipe_versions[@]}"; do
    conda create -n $pipeline_version python=2.7 --yes;
    source activate $pipeline_version
    conda install STAR --yes
    pip install numpy
    pip install -r requirements.txt
    source deactivate
    done

source activate pipe_151
pip install stpipeline==1.5.1 #&& python tests/pipeline_run_test.py
source deactivate

source activate pipe_master
python setup.py build && python setup.py install #&& python tests/pipeline_run_test.py
source deactivate

cd ../st_pipeline_parallel
source activate pipe_parallel
git checkout parallelizing_filterInputReads
python setup.py build && python setup.py install #&& python tests/pipeline_run_test.py
source deactivate

cd ../st_pipeline_coordparse
source activate pipe_coordparse
git checkout proof_of_concept_coordinateParser
python setup.py build && python setup.py install #&& python tests/pipeline_run_test.py
source deactivate

cd ../st_pipeline_anno_par
source activate pipe_anno_par
git checkout test_anntotation_parallel
python setup.py build && python setup.py install #&& python tests/pipeline_run_test.py
source deactivate

wait

cd $testfolder

#
# Define the accessions, names and id file connections
#
accessions_by_size=( \
                   SRR3382385   # 0   45M \
                   SRR3382387   # 1  154M \
                   SRR3382389   # 2  156M \
                   SRR3382373   # 3  169M \
                   SRR3382386   # 4  172M \
                   SRR3382388   # 5  179M \
                   SRR3382374   # 6  181M \
                   SRR3382390   # 7  224M \
                   SRR3382384   # 8  386M \
                   SRR3382383   # 9  401M \
                   SRR3382372   # 10 424M \
                   SRR3382371 ) # 11 482M

declare -A names_by_accessions=( \
                               [SRR3382371]=MOB_REPLICATE_1 \
                               [SRR3382372]=MOB_REPLICATE_2 \
                               [SRR3382383]=MOB_REPLICATE_3 \
                               [SRR3382384]=MOB_REPLICATE_4 \
                               [SRR3382385]=MOB_REPLICATE_5 \
                               [SRR3382386]=MOB_REPLICATE_6 \
                               [SRR3382387]=MOB_REPLICATE_7 \
                               [SRR3382388]=MOB_REPLICATE_8 \
                               [SRR3382389]=MOB_REPLICATE_9 \
                               [SRR3382390]=MOB_REPLICATE_10 \
                               [SRR3382373]=MOB_REPLICATE_11 \
                               [SRR3382374]=MOB_REPLICATE_12 )

declare -A ids_by_accessions=( \
                             [SRR3382371]=1000L2 \
                             [SRR3382372]=1000L2 \
                             [SRR3382383]=1000L2 \
                             [SRR3382384]=1000L2 \
                             [SRR3382385]=1000L3 \
                             [SRR3382373]=1000L5 \
                             [SRR3382374]=1000L5 \
                             [SRR3382386]=1000L5 \
                             [SRR3382387]=1000L5 \
                             [SRR3382388]=1000L5 \
                             [SRR3382389]=1000L5 \
                             [SRR3382390]=1000L5 )

path_2_ids=$testfolder/st_pipeline/ids

#
# Loop throught the samples get files and start the pipe
#
for i in {0..11}; do

# get the info from the arrays
acc=${accessions_by_size[$i]}
name=${names_by_accessions[$acc]}
id=${ids_by_accessions[$acc]}

# fetch raw data
fastq-dump --accession $acc --split-files --origfmt #--minSpotId 1000000 --maxSpotId 11000000

for pipeline_version in "${pipe_versions[@]}"; do

source activate $pipeline_version

# make directories for results
echo $name $acc $id $pipeline_version
results_dir=results_$(date +"%y%m%d_%H_%M")/$name\_$pipeline_version
temp_dir=temp_$(date +"%y%m%d_%H_%M")/$name\_$pipeline_version
mkdir -p $results_dir
mkdir -p $temp_dir

dstat \
-Tmcd \
--top-cpu \
--noheader \
--output $results_dir/sys_stat.$name\_$pipeline_version.csv > $results_dir/sys_stat.$name\_$pipeline_version.txt & stat_logger_pid=$(echo $!)

# start pipeline
st_pipeline_run.py \
--ref-map $genomeDir \
--ref-annotation $gtfFile \
--expName $name\_$pipeline_version \
$acc\_1.fastq \
$acc\_2.fastq \
--output-folder $results_dir \
--temp-folder $temp_dir \
--ids $path_2_ids/$id\_barcodes.txt \
--log-file $results_dir/$name\_$pipeline_version.log.txt \
--mapping-threads 12 \
--verbose \
--htseq-no-ambiguous \
--disable-multimap;

kill $stat_logger_pid

source deactivate

done

# remove the fastqs
rm -v $acc\_*.fastq

done

for pipeline_version in "${pipe_versions[@]}"; do
    conda-env remove -n $pipeline_version --yes;
    done