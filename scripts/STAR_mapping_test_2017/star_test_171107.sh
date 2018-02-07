# Create the test folder
test_folder=~/STARtest
rm -vr $test_folder
mkdir $test_folder
cd $test_folder

# Get the infiles from the st pipeline runs of version 151 and 160
cp /tmp/st_pipeline_test_tempGmK4Pn/R2_quality_trimmed.fastq ./R2_quality_trimmed.v151.fastq
cp /tmp/st_pipeline_test_temprv5KEP/R2_quality_trimmed.bam ./R2_quality_trimmed.v160.bam

# Convert the version 160 bam file to fastq format
samtools view R2_quality_trimmed.v160.bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > R2_quality_trimmed.v160.fastq

# Convert the version 151 fastq to sam format
source activate bam_pipeline
picard FastqToSam F1=R2_quality_trimmed.v151.fastq O=R2_quality_trimmed.v151.sam RG=aa QUALITY_FORMAT=Standard SAMPLE_NAME=unknown_sample
source deactivate

# Convert newly converted sam and bam files to bam resp sam so all input versions will be present in all formats
samtools view -Sb R2_quality_trimmed.v151.sam -o R2_quality_trimmed.v151.bam
samtools view -h R2_quality_trimmed.v160.bam -o R2_quality_trimmed.v160.sam

# declare three different parameters combinations to test
declare -A parameter_strings=( \
            [1]="" \
            [2]="--alignIntronMin 1" \
            [3]="--alignIntronMin 1 --alignIntronMax 1" \
            [4]="--alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1" \
            [5]="--alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterMatchNmin 25" \
            [6]="--alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterMatchNmin 25 --alignEndsType Local" \
            [7]="--alignIntronMin 1 --alignIntronMax 1 --outFilterMultimapNmax 1 --outFilterMatchNmin 25 --alignEndsType EndToEnd" )

# Define which reference genome to map the reads to
reference=~/code_repos/my_mods/st_pipeline/tests/config/contaminant_genomes/R45S5_R5S1/
# Alternatively (larger reference longer runtime)
# reference=~/code_repos/my_mods/st_pipeline/tests/config/genomes/mouse_grcm38/

# The two STAR versions that will be used are installed in two conda envs loop through both and use $(STAR --version) to get the current version
for pipeline_version in bam_pipeline pipe_151;
do
    source activate $pipeline_version

    # loop through the different parameter combinations
    for parameters_index in {1..7};
    do

        # The st_pipeline_version variable determines the original input file before conversion v151=fastq and v160=bam
        for st_pipeline_version in v151 v160;
        do

            # The input_file_type varible determines the format of the input file used for mapping with STAR
            for input_file_type in bam sam fastq;
            do
                
                # For each combination of parameters, input-file, format and STAR version, make a folder for the output-files and run STAR
                # skip bam and sam input for STAR_2.5.3a as its not supported
                if [ $(STAR --version) == "STAR_2.5.3a" ] && [ $input_file_type == "bam" ]
                then echo "skipping old version bam"
                elif [ $(STAR --version) == "STAR_2.5.3a" ] && [ $input_file_type == "sam" ]
                then echo "skipping old version sam"
                else
                    
                    echo "Running: "$(STAR --version)", on input-files from st_pipeline-version="$st_pipeline_version", parameter-set="$parameters_index ", input-file-format="$input_file_type
                    mkdir $(STAR --version)\.$st_pipeline_version\.param_$parameters_index\.$input_file_type

                    parameters=${parameter_strings[$parameters_index]}
                    # Add parameters defining sam and bam input when appropriate
                    if [ $input_file_type == "bam" ] || [ $input_file_type == "sam" ]
                    then
                        parameters="$parameters --readFilesType SAM SE --readFilesCommand samtools view -h"
                    fi

                    # run STAR for each combination of input-files and STAR version
                    cd $(STAR --version)\.$st_pipeline_version\.param_$parameters_index\.$input_file_type
                    STAR \
                    --runThreadN 12 \
                    --genomeDir $reference \
                    --readFilesIn ~/STARtest/R2_quality_trimmed.$st_pipeline_version.*$input_file_type \
                    $parameters
                    grep "STAR --" Log*
                    echo -e "\t"$(STAR --version)\.$st_pipeline_version\.param_$parameters_index\.$input_file_type > Log.final.out.header_and_no_space
                    cat Log.final.out | sed s/\ \ /\ /g | sed s/\ \ /\ /g | sed s/\ \ /\ /g | sed s/\ \ /\ /g | sed s/\ \ /\ /g | sed s/\ \ /\ /g| sed s/\ \ /\ /g >> Log.final.out.header_and_no_space
                    cd $test_folder

                fi
            done
        done
    done
    source deactivate
done

# Print all the stats to tsv
number_of_testfolders=$(ls -lh STAR_2.5.3*/Log.final.out | wc -l | awk '{print $0*2}')
paste STAR_2.5.3*/Log.final.out.header_and_no_space | $( echo cut -f 1 $(for ((X = 2; X <= $number_of_testfolders; X += 2)); do echo ",$X" ;done) | sed s/\ \,/\,/g) > stats.tsv