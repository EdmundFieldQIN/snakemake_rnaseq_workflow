rule trim_fastq:
    input: 
        r1="../data/rawdata/{sample}_1.fastq.gz",
        r2="../data/rawdata/{sample}_2.fastq.gz"
    output: 
        r1="01cleandata/{sample}_R1.clean.fq.gz",
        r2="01cleandata/{sample}_R2.clean.fq.gz",
        report="01cleandata/{sample}.fastp.html",
        json="01cleandata/{sample}.fastp.json"
    threads: 4 # 4线程
    log: "logs/{sample}.fastp.log"# 日志文件
    conda: "rnaseq"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -w {threads} -h {output.report} -j {output.json} > {log} 2>&1"

rule fastqc_multiqc:
    input: 
        expand("01cleandata/{sample}.fastp.html", sample=SAMPLES),
        expand("01cleandata/{sample}.fastp.json", sample=SAMPLES)
    output: 
        "02fastqc/multiqc_report.html"
    threads: 10
    conda: "rnaseq"
    shell:
        "fastqc -t {threads} -q -o 02fastqc/ 01cleandata/*.gz && multiqc 02fastqc/* -o 02fastqc/"
