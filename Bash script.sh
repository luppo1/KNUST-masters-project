### Quality Check
fastqc -t 28 -o /home/abdul/albert/QC_before *.fastq.gz

#### Trimming with Trimmomatic
#!/bin/bash
output=/home/abdul/albert/Trimmed_Data
input=/home/abdul/albert/Data
for i in $input/*_1.fastq.gz;
do
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*.fastq.gz}"
sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'` 
trimmomatic PE -threads 28 -trimlog $output/"${base}".log.gz $input/"${base}"_1.fastq.gz $input/"${base}"_2.fastq.gz $output/"${base}"_1.trimmed_PE.fastq.gz $output/"${base}"_1.trimmed_SE.fastq.gz $output/"${base}"_2.trimmed_PE.fastq.gz $output/"${base}"_2.trimmed_SE.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


### Quality Check
fastqc -t 28 -o /home/abdul/albert/QC_after *PE.fastq.gz


### ALIGNMENT WITH KALLISTO
### BUILDING INDEX FILES
#!/bin/bash
kallisto index -i /home/abdul/albert/Reference/kallistoindex.idx -k 31 Homo_sapiens.GRCh38.cdna.all.fa.gz
done


## Actual Alignment
#!/bin/bash 
input=/home/abdul/albert/Trimmed_Data
for i in $input/*_1.trimmed_PE.fastq.gz; 
do 
withpath="${i}"
filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}" 
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
kallisto quant -i /home/abdul/albert/Reference/kallistoindex.idx -o /home/abdul/albert/Kallisto/"${base}" --genomebam --bias  --gtf /home/abdul/albert/Reference/Homo_sapiens.GRCh38.99.gtf.gz --chromosomes /home/abdul/albert/Reference/chrom.txt -b 50 -t 27 $input/"${base}"_1.trimmed_PE.fastq.gz $input/"${base}"_2.trimmed_PE.fastq.gz &> /home/abdul/albert/Kallisto/"${base}".log
done





## ALIGNMNET WITH SALMON
## BUILDING INDEX

grep "^>" <(gunzip -c /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat /home/abdul/albert/Reference/Homo_sapiens.GRCh38.cdna.all.fa.gz /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > gentrome.fa.gz
salmon index -p 27 -t /home/abdul/albert/Reference/gentrome.fa.gz -d /home/abdul/albert/Reference/decoys.txt -i salmon_index -k 31


## ALIGNMENT
#!/bin/bash
input=/home/abdul/albert/Trimmed_Data
for i in $input/*_1.trimmed_PE.fastq.gz;
do 
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*1.trimmed_PE.fastq.gz}"
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
salmon quant -i /home/abdul/albert/Salmon/salmon_index -l A -1 $input/"${base}"*_1.trimmed_PE.fastq.gz -2 $input/"${base}"*_2.trimmed_PE.fastq.gz --validateMappings -p 27 -o /home/abdul/albert/Salmon/"${base}" –-rangeFactorizationBins 4 --seqBias --numBiasSamples 2000000 --gcBias --numBootstraps 50
done




## ALIGNMENT WITH HISAT2
## BUILDING INDEX
hisat2-build -p 27 /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa hisatindex


# ALIGNMNET
#!/bin/bash
input=/home/abdul/albert/Trimmed_Data
output=/home/abdul/albert/Hisat2
for i in $input/*_1.trimmed_PE.fastq.gz;
do 
withpath="${i}" filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}"
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
hisat2 -p 27 -x /home/abdul/albert/Hisat2/Index/hisatindex -1 $input/"${base}"*_1.trimmed_PE.fastq.gz -2 $input/"${base}"*_2.trimmed_PE.fastq.gz -S $output/"${base}".hisat.sam --summary-file $output/"${base}".txt 
echo "$sample_name done!"
samtools view -@ 26 -m 26G -ub $output/"${base}".hisat.sam -o $output/"${base}".hisat.bam
echo "${sample_name} hisat.sam change to bam done!"
samtools sort -n -@ 16 -m 2G -T /tmp/ $output/"${base}".hisat.bam -o $output/"${base}".hisat.sorted.bam
rm $output/"${base}".hisat.sam $output/"${base}".hisat.bam
echo "${sample_name} hisat.sorted.bam sort done!"
echo "$base done!"
done


## Quantification step
featureCounts -p -T 27 -t exon -g gene_id --extraAttributes gene_id,gene_biotype -a ~/reference_genomes/ensembl/Homo_sapiens.GRCh38.92.gtf -o hisat2.txt *.sorted.bam
 


## ALIGNMENT WITH STAR

## INDEXING
#!/bin/bash 
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir ~/albert/STAR/Index --genomeFastaFiles /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ~/reference_genomes/ensembl/Homo_sapiens.GRCh38.92.gtf --genomeSAindexNbases 12 --genomeSAsparseD 2 --sjdbOverhang 100
done



## ACTUAL ALIGNMNENT
#!/bin/bash 
input=/home/abdul/albert/Trimmed_Data
for i in $input/*_1.trimmed_PE.fastq.gz; 
do 
withpath="${i}"
filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}" 
sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'`
mkdir /home/abdul/albert/STAR/"${base}"
STAR --runThreadN 25 --genomeDir /home/abdul/albert/STAR/Index --readFilesIn $input/"${base}"*1.trimmed_PE.fastq.gz $input/"${base}"*2.trimmed_PE.fastq.gz --outFileNamePrefix /home/abdul/albert/STAR/"${base}"/"${base}"_ --quantMode GeneCounts --sjdbOverhang 100 --genomeSAsparseD 2 --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted SortedByCoordinate --outSAMmapqUnique 60
done


## Quantification
featureCounts -p -T 27 -t exon -g gene_id --extraAttributes gene_id,gene_biotype -a ~/reference_genomes/ensembl/Homo_sapiens.GRCh38.92.gtf -o star.txt *sortedByCoord.out.bam





