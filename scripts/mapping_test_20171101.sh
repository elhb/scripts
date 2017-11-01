rm -vr ~/STARtest
mkdir ~/STARtest
cd ~/STARtest

cp /tmp/st_pipeline_test_tempGmK4Pn/R2_quality_trimmed.fastq ./R2_quality_trimmed.v151.fastq
cp /tmp/st_pipeline_test_temprv5KEP/R2_quality_trimmed.bam ./R2_quality_trimmed.v160.bam

samtools view R2_quality_trimmed.v160.bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > R2_quality_trimmed.v160.fastq

source activate bam_pipeline
picard FastqToSam F1=R2_quality_trimmed.v151.fastq O=R2_quality_trimmed.v151.sam RG=aa QUALITY_FORMAT=Standard SAMPLE_NAME=unknown_sample
source deactivate

samtools view -Sb R2_quality_trimmed.v151.sam -o R2_quality_trimmed.v151.bam
samtools view -h R2_quality_trimmed.v160.bam -o R2_quality_trimmed.v160.sam

for pipeline_version in bam_pipeline pipe_151;
do
    source activate $pipeline_version
    for version in v151 v160;
    do
        for file_type in bam sam fastq;
        do
            mkdir $(STAR --version)\_$version\_$file_type
        done

        echo $(STAR --version)\_$version\_bam
        cd $(STAR --version)\_$version\_bam
        STAR --genomeDir ~/code_repos/my_mods/st_pipeline/tests/config/genomes/mouse_grcm38/ --readFilesIn ~/STARtest/R2_quality_trimmed.$version.*bam --readFilesType SAM SE --readFilesCommand samtools view -h
        cd ~/STARtest

        echo $(STAR --version)\_$version\_sam
        cd $(STAR --version)\_$version\_sam
        STAR --genomeDir ~/code_repos/my_mods/st_pipeline/tests/config/genomes/mouse_grcm38/ --readFilesIn ~/STARtest/R2_quality_trimmed.$version.*sam --readFilesType SAM SE --readFilesCommand samtools view -h
        cd ~/STARtest

        echo $(STAR --version)\_$version\_fastq
        cd $(STAR --version)\_$version\_fastq
        STAR --genomeDir ~/code_repos/my_mods/st_pipeline/tests/config/genomes/mouse_grcm38/ --readFilesIn ~/STARtest/R2_quality_trimmed.$version.*fastq
        cd ~/STARtest

    done
    source deactivate
done

echo $(ls |grep STAR)
paste STAR_2.5.3*/Log.final.out | cut -f 1,2,4,6,8,10,12,14,16

