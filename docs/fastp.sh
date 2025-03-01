import os
import sys

directorylist = ### replace with list of all data folders
output_directory = ### replace with output directory

for inputdir in directorylist:
    for root, subdirs, files in os.walk(inputdir):
        for file in files:
            if file.endswith("R1_001.fastq.gz"):
                input1 = os.path.join(root, file)
                input2 = input1.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                
                output1 = output_directory + file.replace("fastq.gz", "fastp.gz")
                output2 = output_directory + file.replace("R1_001.fastq.gz", "R2_001.fastp.gz")
                
                if os.path.exists(output1) is not True:
                    cmd = "fastp --dedup -i %s -I %s -o %s -O %s"%(input1, input2, output1, output2)
                    print(cmd)
