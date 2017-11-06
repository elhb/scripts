# create the test folder
test_folder=~/STARtest2
rm -vr $test_folder
mkdir $test_folder
cd $test_folder

# get the infiles
cp /tmp/st_pipeline_test_tempGmK4Pn/R2_quality_trimmed.fastq ./R2_quality_trimmed.v151.fastq
cp /tmp/st_pipeline_test_temprv5KEP/R2_quality_trimmed.bam ./R2_quality_trimmed.v160.bam

# convert bam file to fastq
samtools view R2_quality_trimmed.v160.bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > R2_quality_trimmed.v160.fastq

# convert fastq to sam
source activate bam_pipeline
picard FastqToSam F1=R2_quality_trimmed.v151.fastq O=R2_quality_trimmed.v151.sam RG=aa QUALITY_FORMAT=Standard SAMPLE_NAME=unknown_sample
source deactivate

# convert sam to bam and bam to sam
samtools view -Sb R2_quality_trimmed.v151.sam -o R2_quality_trimmed.v151.bam
samtools view -h R2_quality_trimmed.v160.bam -o R2_quality_trimmed.v160.sam

# The two STAR versions are installed in two conda envs loop through both
for pipeline_version in bam_pipeline pipe_151;
do
    source activate $pipeline_version

# for each combination of input-file and STAR version, make a folder for the output-files
    for version in v151 v160;
    do
        for file_type in bam sam fastq;
        do
            if [ $(STAR --version) == "STAR_2.5.3a" ] && [ file_type == "bam" ]
            then echo "skipping old version bam"
            elif [ $(STAR --version) == "STAR_2.5.3a" ] && [ file_type == "sam" ]
            then echo "skipping old version sam"
            else
                mkdir $(STAR --version)\_$version\_$file_type
            fi
        done

# run STAR for each combination of input-files and STAR version
        if [ $(STAR --version) == "STAR_2.5.3a_modified" ];
        then
            echo $(STAR --version)\_$version\_bam
            cd $(STAR --version)\_$version\_bam
            STAR --runThreadN 12 --alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outFilterMatchNmin 25 --genomeDir ~/code_repos/my_mods/st_pipeline/tests/config/contaminant_genomes/R45S5_R5S1/ --readFilesIn ~/STARtest2/R2_quality_trimmed.$version.*bam --readFilesType SAM SE --readFilesCommand samtools view -h
            cd $test_folder
        fi

        if [ $(STAR --version) == "STAR_2.5.3a_modified" ];
        then
            echo $(STAR --version)\_$version\_sam
            cd $(STAR --version)\_$version\_sam
            STAR --runThreadN 12 --alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outFilterMatchNmin 25 --genomeDir ~/code_repos/my_mods/st_pipeline/tests/config/contaminant_genomes/R45S5_R5S1/ --readFilesIn ~/STARtest2/R2_quality_trimmed.$version.*sam --readFilesType SAM SE --readFilesCommand samtools view -h
            cd $test_folder
        fi

        echo $(STAR --version)\_$version\_fastq
        cd $(STAR --version)\_$version\_fastq
        STAR --runThreadN 12 --alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outFilterMatchNmin 25 --genomeDir ~/code_repos/my_mods/st_pipeline/tests/config/contaminant_genomes/R45S5_R5S1/ --readFilesIn ~/STARtest2/R2_quality_trimmed.$version.*fastq
        cd $test_folder

    done
    source deactivate
done

# print stats to tsv
echo -e "\t"$(ls |grep STAR) | sed s/\ /\\t/g > stats.tsv
paste STAR_2.5.3*/Log.final.out | cut -f 1,2,4,6,8,10,12,14,16 >> stats.tsv