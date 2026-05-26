rule build_hisat2_index:
    input:  # 输入：配置文件中指定的基因组fasta文件
        genome=rules.unzip_genome_annotation.output.genome
    output:  # 输出：HISAT2索引文件
        idx=expand("../data/ref/genome.hisat2.{idx}.ht2", idx=range(1, 9))
    log: "logs/build_hisat2_index.log"
    conda: "rnaseq"
    threads: 10
    shell:  # hisat2-build构建索引
        """
        hisat2-build -p {threads} {input.genome} ../data/ref/genome.hisat2 1>{log} 2>&1
        """

rule align_hisat2:
    input:  # 输入：质控后fastp + HISAT2索引
        r1="01cleandata/{sample}_R1.clean.fq.gz",
        r2="01cleandata/{sample}_R2.clean.fq.gz",
        idx=expand("../data/ref/genome.hisat2.{idx}.ht2", idx=range(1, 9))
    output:  # 输出：SAM文件 + 比对日志
        sam=temp("03aligned/{sample}.sam"),  
        log="03aligned/{sample}.hisat2.log"
    threads: 10
    conda: "rnaseq"
    shell: 
        """
        hisat2 -x ../data/ref/genome.hisat2 -p {threads} -1 {input.r1} -2 {input.r2} \
            --new-summary --rna-strandness RF \
            -S {output.sam} \
            1> {output.log} 2>&1
        """

rule samtools_sort_idx:
    input:  # 输入：SAM文件
        sam="03aligned/{sample}.sam"
    output:  # 输出：排序后的BAM + BAM索引
        bam=protected("03aligned/{sample}.bam"),
        idx="03aligned/{sample}.bam.bai"
    threads: 10
    log: "logs/{sample}.sort.log"
    conda: "rnaseq"
    shell:  # samtools排序+建索引
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} 1>{log} 2>&1
        samtools index {output.bam}
        """
