rule MergeFastqFiles:
    input:
        lambda wildcards: expand("{RNASeqSample}", RNASeqSample=RNASeqSampleToFastq_dict[wildcards.RNASeqSample])
    output:
        config["temp_files_prefix"] + "RNASeqFastq/{RNASeqSample}.fastq.gz"
    shell:
        "cat {input} > {output}"

rule STAR_make_index:
    input:
        fasta=config["ref"]["genome"],
        gtf=config["ref"]["genomegtf"]
    output:
        config["ref"]["STAR_index_folder"] + "chrLength.txt"
    log:
        "logs/STAR/MakingIndex.log"
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN 4 --genomeDir MiscOutput/STAR_index/ --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """

rule STAR_alignment:
    input:
        index = config["ref"]["STAR_index_folder"] + "chrLength.txt",
        fastq = get_RNASeq_merged_fastq
    log:
        "logs/STAR/{RNASeqSample}.log"
    threads: 12
    output:
        config["RNASeqMapping"]["PathTo_STAR_alignment_Results"] + "{RNASeqSample}/SJ.out.tab",
        config["RNASeqMapping"]["PathTo_STAR_alignment_Results"] + "{RNASeqSample}/Aligned.sortedByCoord.out.bam",
        config["RNASeqMapping"]["PathTo_STAR_alignment_Results"] + "{RNASeqSample}/ReadsPerGene.out.tab"
    shell:
        "STAR --runThreadN {threads} --genomeDir MiscOutput/STAR_index/ --readFilesIn {input.fastq} --outSAMtype BAM SortedByCoordinate --outWigStrand Unstranded --outWigType wiggle --alignEndsType EndToEnd --quantMode GeneCounts --twopassMode Basic --readFilesCommand zcat --outFileNamePrefix RNASeq/STAR/{wildcards.RNASeqSample}/ &> {log}"

rule index_RNA_seq_bams:
    input:
        bam = lambda wildcards: "RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample)
    output:
        bai= "RNASeq/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule Gather_RNA_seq_FlowCellInfo:
    input: expand("{myfiles}", myfiles=RNASeqBasenameToFastq.values())
    output: "MiscOutput/RNA_seq_FlowCellInfo.txt"
    log: "logs/Gather_RNA_seq_FlowCellInfo"
    shell: "bash scripts/GetFastqIdentifierInfo.sh {input} > {output} 2> {log}"
