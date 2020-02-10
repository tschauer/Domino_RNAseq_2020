#! /bin/bash
#
# STAR.sbatch
#
#SBATCH -J STAR_array   # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 8                        # 8 tasks
#SBATCH -N 1                         # one node
#SBATCH -t 0-12:00                    # Running time of 2 hours
#SBATCH --mem 20000                  # Memory request in megabytes
#SBATCH -o STAR_%A_%a.out          # Standard output
#SBATCH -e STAR_%A_%a.err          # Standard error

# grab filename
FILENAME=`ls *_1.txt.gz  | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
FILEBASE=`echo ${FILENAME} | sed -e "s/_1.txt.gz//g"`

module load ngs/STAR/2.5.3a
module load ngs/samtools






### change to your STAR index directory ###
STAR_INDEX="/work/project/becbec_002/genome/STAR_dmel-all-chromosome-r6.17"

### change to your GTF file ###
GTF_FILE="/work/project/becbec_002/genome/dmel-all-r6.17.gtf"


STAR \
--runThreadN 8 \
--readFilesCommand gunzip -c \
--quantMode GeneCounts \
--genomeDir ${STAR_INDEX} \
--sjdbGTFfile ${GTF_FILE} \
--readFilesIn ${FILEBASE}_1.txt.gz  \
--outFileNamePrefix ${FILEBASE}.dmel. \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 5000000000 \
--outFilterMultimapNmax 1







### change to your STAR index directory ###
STAR_INDEX="/work/project/becbec_002/genome/STAR_dvir-all-chromosome-r1.07"

### change to your GTF file ###
GTF_FILE="/work/project/becbec_002/genome/dvir-all-r1.07.gtf"


STAR \
--runThreadN 8 \
--readFilesCommand gunzip -c \
--quantMode GeneCounts \
--genomeDir ${STAR_INDEX} \
--sjdbGTFfile ${GTF_FILE} \
--readFilesIn ${FILEBASE}_1.txt.gz  \
--outFileNamePrefix ${FILEBASE}.dvir. \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 5000000000 \
--outFilterMultimapNmax 1
