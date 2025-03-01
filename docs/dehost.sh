module load bowtie2
module load samtools

dehost_out=### replace with dehost out folder
data=### replace with data folder
export BOWTIE2_INDEXES=/GRCh38_noalt_as/ #replace with GRCh38_noalt_as location

while IFS=$'\t' read sample r1 r2

do     

    echo "bowtie2 -p 32 -x GRCh38_noalt_as -q -1 $data/$r1 -2 $data/$r2 -S $dehost_out/${sample}_mapped_and_unmapped.sam"

    echo "samtools view -bS $dehost_out/${sample}_mapped_and_unmapped.sam > $dehost_out/${sample}_mapped_and_unmapped.bam"

    echo "samtools view -b -f 12 -F 256 $dehost_out/${sample}_mapped_and_unmapped.bam > $dehost_out/${sample}_bothReadsUnmapped.bam"

    echo "samtools sort -n -m 5G -@ 4 $dehost_out/${sample}_bothReadsUnmapped.bam  -o $dehost_out/${sample}_bothReadsUnmapped_sorted.bam"

    echo "samtools fastq -@ 16 $dehost_out/${sample}_bothReadsUnmapped_sorted.bam -1 $dehost_out/${sample}_host_removed_R1.fastq.gz -2 $dehost_out/${sample}_host_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n"

done < $data/samplelist.txt
