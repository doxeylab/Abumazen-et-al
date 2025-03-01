#running Salmon

module load r python StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 salmon/1.7.0 #loading relevant modules

idx="gencode"
data="""" #replace with location of project folder
reference=gencode_v39_idx #location of gencode folder
output="""" #replace with location of data files
 
# data.txt is a metadata file containing sample ID, file name for r1 sequence file and filename for r2 sequence file for sample with that sampleID. Each sample on its own line
# loop to iterate over lines in file data.txt

while IFS=$'\t' read sample r1 r2
do
        echo "Processing sample $sample with reads $r1 and reads $r2"
        salmon quant -i $reference -l A\
                -1 $data/$r1 \
                -2 $data/$r2 \
                --validateMappings \
                --seqBias \
                --gcBias \
                -p 12 -o $output/$sample

done < $data/data.txt

# post processing step to merge salmon outputs for certain downstream analysis
salmon quantmerge --quants $output/* -o $output/merged_${idx}.quant
