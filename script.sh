#!/bin/bash

# === Logging setup ===
LOGFILE="/home/pls_trainee/adreeja/RNA/pipeline_log_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1
start_time=$(date +%s)

echo " RNA-Seq Pipeline Started at $(date) "

# === Directories ===
ref="/home/pls_trainee/adreeja/ref/hg38"
known_sites="/home/pls_trainee/adreeja/known_sites"
reads="/home/pls_trainee/adreeja/RNA/reads"
aligned_reads="/home/pls_trainee/adreeja/RNA/aligned_reads"
fastqc_data="/home/pls_trainee/adreeja/RNA/reads/fastqc_data"
trimmed_data="$reads/trimmed_data"
bam_files="/home/pls_trainee/adreeja/RNA/bam_files"
sam_files="/home/pls_trainee/adreeja/RNA/sam_files"
results="/home/pls_trainee/adreeja/RNA/results"
gtf="/home/pls_trainee/adreeja/gtf/gencode.v44.annotation.gtf"
multiqc_output="/home/pls_trainee/adreeja/RNA/multiqc_report"

# === Functions ===

run_fastp() {
    echo "Running fastp trimming..."
    for fq1 in "$reads"/*_1.fastq.gz; do
        fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
        sample=$(basename "$fq1" _1.fastq.gz)

        echo "Trimming $sample"
        fastp \
            -i "$fq1" -I "$fq2" \
            -o "$trimmed_data/${sample}_1.trimmed.fastq.gz" \
            -O "$trimmed_data/${sample}_2.trimmed.fastq.gz" \
            --detect_adapter_for_pe \
            --thread 10 \
            --html "$trimmed_data/${sample}_fastp.html" \
            --json "$trimmed_data/${sample}_fastp.json"
    done
}

run_fastqc() {
    echo "Running FastQC..."
    for fq in "$trimmed_data"/*.trimmed.fastq.gz; do
        echo "FastQC on $fq"
        fastqc "$fq" -t 10 -o "$fastqc_data"
    done
}

run_multiqc() {
    echo "Running MultiQC..."
    multiqc "$fastqc_data" "$trimmed_data" -o "$multiqc_output"
}

run_alignment_and_sorting() {
    echo "Running Alignment and Sorting..."
    for fq1 in "$trimmed_data"/*_1.trimmed.fastq.gz; do
        fq2="${fq1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
        sample=$(basename "$fq1" _1.trimmed.fastq.gz)

        echo "Aligning $sample"
        hisat2 -p 10 --dta -x "$ref" -1 "$fq1" -2 "$fq2" -S "$sam_files/${sample}.sam"

        echo "Converting SAM to BAM: $sample"
        samtools view -@ 10 -bS "$sam_files/${sample}.sam" -o "$bam_files/${sample}.bam"

        echo "Name-sorting BAM for fixmate"
        samtools sort -n -@ 10 "$bam_files/${sample}.bam" -o "$bam_files/${sample}.name_sorted.bam"

        echo "Fixing mate info"
        samtools fixmate -m "$bam_files/${sample}.name_sorted.bam" "$bam_files/${sample}.fixmate.bam"

        echo "Coordinate-sorting BAM"
        samtools sort -@ 10 "$bam_files/${sample}.fixmate.bam" -o "$bam_files/${sample}.sorted.bam"

        echo "Indexing final BAM"
        samtools index "$bam_files/${sample}.sorted.bam"

        # Optional cleanup
        rm "$bam_files/${sample}.bam" "$bam_files/${sample}.name_sorted.bam" "$bam_files/${sample}.fixmate.bam"
    done
}


run_featurecounts() {
    echo "Running featureCounts"
    featureCounts -T 10 \
        -a "$gtf" \
        -o "$bam_files/gene_counts.txt" \
        -t exon -g gene_id \
        "$bam_files"/*.sorted.bam
}

run_qualimap() {
    echo "Running Qualimap..."
    for bam in "$bam_files"/*.sorted.bam; do
        sample=$(basename "$bam" .sorted.bam)
        outdir="$bam_files/qualimap_${sample}"

        echo "Qualimap for $sample"
        [ -f "${bam}.bai" ] || samtools index "$bam"
        [ -d "$outdir" ] && rm -rf "$outdir"

        qualimap rnaseq \
            -bam "$bam" \
            -gtf "$gtf" \
            -outdir "$outdir" \
            -outformat PDF:HTML \
            -p non-strand-specific \
            --java-mem-size=16G

        echo "Qualimap completed for $sample"
    done
}

run_fastp
run_fastqc
run_multiqc
run_alignment_and_sorting
run_featurecounts
run_qualimap

# Calculate total duration
end_time=$(date +%s)
runtime=$((end_time - start_time))
printf " Pipeline completed at $(date). Total duration: %02dh:%02dm:%02ds \n " \
  $((runtime/3600)) $(((runtime%3600)/60)) $((runtime%60))
