#!/bin/bash

python="/home/avm27/anaconda3/envs/RNA_Seq/bin/python"
stringtiescript="/home/avm27/anaconda3/envs/RNA_Seq/bin/prepDE.py"
experimentname="CIVSA"
experimentnamefull="Cocaine Intravenous Self Administration"

echo "Final Code: CIVSA Experiment Pegasus Pre-Processing"

echo "Loading Environment"

cd /home/avm27/Documents/Raw_Sequencing_Data/CIVSA

echo "Combining Samples by Lanes"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/CIVSA -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
do 

    echo "Merging R1" ${samplename}

    cat ${samplename}_L00*_R1_001.fastq.gz > ${samplename}_R1.fastq.gz
   
    echo "Merging R2" ${samplename}

    cat ${samplename}_L00*_R2_001.fastq.gz > ${samplename}_R2.fastq.gz

done
wait

echo "Trimming Samples"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/CIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    trim_galore --cores 6 --dont_gzip --paired ${samplename}_R1.fastq.gz ${samplename}_R2.fastq.gz

wait
done
wait

echo "Trimming Completed"

echo "Begin Aligning to the mm10 Genome"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/CIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

 	STAR --runThreadN 8 \
    --runMode alignReads  \
    --genomeDir /home/avm27/genome/mm10/ncbi_STAR \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD XS \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNoverLmax 0.3 \
    --outFileNamePrefix ${samplename} \
    --quantMode TranscriptomeSAM GeneCounts \
    --readFilesIn ${samplename}_R1_val_1.fq ${samplename}_R2_val_2.fq \
    --outTmpDir ${samplename}
    
wait

gzip ${samplename}_R1_val_1.fq
gzip ${samplename}_R2_val_2.fq

wait
done
wait

echo "Obtain unique mapped reads"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/CIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do 

    samtools view ${samplename}Aligned.sortedByCoord.out.bam | grep -w 'NH:i:1' | samtools view -bt /home/avm27/genome/mm10/mm10.fa.fai > ${samplename}.bam 

wait
done
wait

echo "Index Samples"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/CIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do

	samtools index ${samplename}.bam 

wait
done
wait

echo "Create GTF Files"

for samplename in $(find /home/avm27/Documents/Raw_Sequencing_Data/CIVSA -type f -name "**.fastq.gz" | while read F; do basename $F | rev | cut -c 13- | rev; done | sort | uniq)
do

	stringtie -p 12 -G /home/avm27/genome/mm10/mm10.ncbiRefSeq.gtf -e -o ${samplename}.gtf -A ${samplename}_FPKM.txt -l ${samplename} ${samplename}.bam

wait
done


## Generate Sample Info List for Stringtie. Name the file "sample_lst.txt" and save to the working directory. Below is the code ran for this experiment.


## First Rename Sample Files

mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215755-01_S13.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT01DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215756-01_S14.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT02DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215757-01_S15.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT03DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215758-01_S16.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT04DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215759-01_S17.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_MN01DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215760-01_S19.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_MN02DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215761-01_S20.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_MN03DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215762-01_S18.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV01DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215763-01_S21.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV02DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215764-01_S19.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV03DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215765-01_S20.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV04DA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215768-01_S23.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT01nDA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215767-01_S24.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT02nDA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215769-01_S21.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT03nDA.gtf
mv /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/202215770-01_S25.bam /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT04nDA.gtf

echo "Generating sample info list for StringTie"

echo "FT01DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT01DA.gtf" > sample_lst_$experimentname.txt
echo "FT02DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT02DA.gtf" >> sample_lst_$experimentname.txt
echo "FT03DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT03DA.gtf" >> sample_lst_$experimentname.txt
echo "FT04DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT04DA.gtf" >> sample_lst_$experimentname.txt
echo "MN01DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_MN01DA.gtf" >> sample_lst_$experimentname.txt
echo "MN02DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_MN02DA.gtf" >> sample_lst_$experimentname.txt
echo "MN03DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_MN03DA.gtf" >> sample_lst_$experimentname.txt
echo "CRV01DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV01DA.gtf" >> sample_lst_$experimentname.txt
echo "CRV02DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV02DA.gtf" >> sample_lst_$experimentname.txt
echo "CRV03DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV03DA.gtf" >> sample_lst_$experimentname.txt
echo "CRV04DA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_CRV04DA.gtf" >> sample_lst_$experimentname.txt
echo "FT01nDA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT01nDA.gtf" >> sample_lst_$experimentname.txt
echo "FT02nDA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT02nDA.gtf" >> sample_lst_$experimentname.txt
echo "FT03nDA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT03nDA.gtf" >> sample_lst_$experimentname.txt
echo "FT04nDA /home/avm27/Documents/Raw_Sequencing_Data/CIVSA2/CIVSA_FT04nDA.gtf" >> sample_lst_$experimentname.txt

echo "Completed Sample Info File"

## Generate Gene Count matrix for DESEQ2. Name the file "gene_count_matrix_CIVSA.csv" and save to the working directory.

echo "Generating Gene Count Matrix File"

$python $stringtiescript -i sample_lst_$experimentname.txt -g gene_count_matrix_$experimentname.csv

echo $experimentnamefull "RNA-Sequencing Pre-Processing Analysis Complete"



