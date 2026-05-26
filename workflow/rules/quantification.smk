rule quantification:
    input: 
        bam = "03aligned/{sample}.bam",
        gtf = rules.unzip_genome_annotation.output.annotation
    output: 
        expr="04quantification/{sample}.count",
        log="04quantification/{sample}.log"
    params: 
        outputprefix = "04quantification/{sample}"
    conda: "rnaseq"
    shell:
        "Rscript ../workflow/script/run-featurecounts.R -b {input.bam} -g {input.gtf} -f exon -a gene_id -i TRUE -s 2 -o {params.outputprefix}"

rule abundance_estimates_to_matrix:
    input:  # 输入：所有样本的count文件
        exprs=expand("04quantification/{sample}.count", sample=SAMPLES)
    output:  # 输出：合并后的计数矩阵 + TPM矩阵
        "05expression/genes.counts.matrix",
        "05expression/genes.FPKM.EXPR.matrix",
        "05expression/genes.TPM.EXPR.matrix",
        "05expression/genes.TMM.EXPR.matrix",
        "05expression/genes.TPM.TMM.EXPR.matrix"
    conda: "rnaseq"
    run:  
        unique_exprs = list(set(input.exprs))
        with open("05expression/quant_files.txt", "w") as f:
            for expr in unique_exprs:
                f.write(f"{expr}\n")
        shell("""
        Rscript ../workflow/script/abundance_to_matrices.R \
            --est_method featureCounts \
            --out_prefix 05expression/genes \
            --quant_files 05expression/quant_files.txt
        """)