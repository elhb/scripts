#!/usr/bin/env bash

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

declare -A pids

testfolder=/data/subsetTest
rm -rvf $testfolder
mkdir $testfolder
start_time=$(date +"%y%m%d_%H_%M")
cp -v "$0" $testfolder/the_script_ran_at_$start_time.sh
mkdir -p $testfolder/rawdata

## To speed up the fetching of reads you can prefetch the sra cahces using the following command:
conda create -n SRApy python=2.7 --yes;
source activate SRApy
mkdir -p /data/ncbi/public/sra/
ln -s /data/ncbi ~/
cd /data/ncbi/public/sra/
pip install srapy
get-project-sras.py -p 316587 -F {acc}.sra &
source deactivate

cd $testfolder/rawdata
for i in {0..11}; do
# get the info from the arrays
acc=${accessions_by_size[$i]}
name=${names_by_accessions[$acc]}
id=${ids_by_accessions[$acc]}
# fetch raw data
fastq-dump --accession $acc --split-files --origfmt --gzip & pids[$acc]=$(echo $!)
# could add the following to get a subset --minSpotId 1000000 --maxSpotId 1010000
done
cd $testfolder

## build the reference genome
genomeDir=~/references/Mus_musculus/Ensembl/GRCm38/Sequence/STAR_Index
gtfFile=~/references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf
mkdir ~/references
cd ~/references
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/Ensembl/GRCm38/Mus_musculus_Ensembl_GRCm38.tar.gz
tar -xvzf Mus_musculus_Ensembl_GRCm38.tar.gz 
mkdir -p $genomeDir
STAR --runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles ~/references/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile $gtfFile \
--sjdbOverhang 100 \
--runThreadN 4 &

# setup the analysis envs
mkdir -p $testfolder/repos
cd $testfolder/repos
git clone https://github.com/elhb/scripts.git
git clone https://github.com/elhb/st_pipeline.git
path_2_ids=$testfolder/repos/st_pipeline/ids
cp -vr st_pipeline st_pipeline_parallel &
cp -vr st_pipeline st_pipeline_151
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
git checkout ff732f10e7cd86fa2cd682cfbd6d6aee356a5e7e
make STAR
mkdir -p ~/bin
ln -s $testfolder/repos/STAR/source/STAR ~/bin/STAR
cd $testfolder/repos/st_pipeline/

pipe_versions=( \
    pipe_parallel \
    pipe_master \
    pipe_151 )

# setup pipeline versions
for pipeline_version in "${pipe_versions[@]}"; do
    conda create -n $pipeline_version python=2.7 --yes;
    source activate $pipeline_version
    pip install numpy
    pip install -r requirements.txt
    source deactivate
    done

source activate pipe_151
cd $testfolder/repos/st_pipeline_151
conda install STAR --yes
git checkout a7fe1222f4c835aecb9a922e0c5b5252ef340d83
python setup.py build && python setup.py install
#pip install stpipeline==1.5.1 #&& python tests/pipeline_run_test.py
source deactivate

source activate pipe_master
cd $testfolder/repos/st_pipeline
git checkout master
python setup.py build && python setup.py install #&& python tests/pipeline_run_test.py
source deactivate

cd $testfolder/repos/st_pipeline_parallel
source activate pipe_parallel
git checkout parallel_filterInput_4_v1.6.0
git pull
python setup.py build && python setup.py install #&& python tests/pipeline_run_test.py
source deactivate

cd $testfolder

#
# Loop throught the samples get files and start the pipe
#
for i in {0..11}; do

# get the info from the arrays
acc=${accessions_by_size[$i]}
name=${names_by_accessions[$acc]}
id=${ids_by_accessions[$acc]}

# uncomment next line to wait for all fastq-dump commands to finish before running pipeline
wait ${pids[@]}

for pipeline_version in "${pipe_versions[@]}"; do

source activate $pipeline_version

# make directories for results
echo $name $acc $id $pipeline_version
STAR --version
results_dir=results_$start_time/$name/$pipeline_version
temp_dir=$results_dir/temp
rm -vr $results_dir
mkdir -p $results_dir
mkdir -p $temp_dir

echo "waiting for fq.gz file to be created ..."
wait ${pids[$acc]}
echo "read file generation completed."

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
$testfolder/rawdata/$acc\_1.fastq.gz \
$testfolder/rawdata/$acc\_2.fastq.gz \
--output-folder $results_dir \
--temp-folder $temp_dir \
--ids $path_2_ids/$id\_barcodes.txt \
--log-file $results_dir/$name\_$pipeline_version.log.txt \
--mapping-threads 12 \
--verbose \
--keep-discarded-files \
--no-clean-up \
--overhang 0 \
--umi-cluster-algorithm hierarchical;

kill $stat_logger_pid

python $testfolder/repos/scripts/scripts/csv_converter_2.py $results_dir/sys_stat.$name\_$pipeline_version.csv > $results_dir/sys_stat.$name\_$pipeline_version.tsv
python $testfolder/repos/st_pipeline/scripts/convertEnsemblToNames.py --annotation $gtfFile $results_dir/$name\_$pipeline_version\_stdata.tsv --output $results_dir/$name\_$pipeline_version\_stdata.names_fixed.tsv

source deactivate

done

# remove the fastqs
rm -v rawdata/$acc\_*.fastq.gz

done

for pipeline_version in "${pipe_versions[@]}"; do
    conda-env remove -n $pipeline_version --yes;
    done