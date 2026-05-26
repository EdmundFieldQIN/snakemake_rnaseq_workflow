# 识别genome和注释文件是否为gz结尾，如果是gz结尾，则解压到data/ref目录下
rule unzip_genome_annotation:
    input: 
        genome=config["ref"]["genome"],
        annotation=config["ref"]["annotation"]
    output: 
        genome="../data/ref/" + os.path.basename(config["ref"]["genome"]).replace(".gz", ""),
        annotation="../data/ref/" + os.path.basename(config["ref"]["annotation"]).replace(".gz", "")
    log: "logs/unzip_genome_annotation.log"
    conda: "rnaseq"
    shell: 
        """
        if [[ {input.genome} == *.gz ]]; then
            gzip -d -c {input.genome} > {output.genome}
        else
            cp {input.genome} {output.genome}
        fi
        if [[ {input.annotation} == *.gz ]]; then
            gzip -d -c {input.annotation} > {output.annotation}
        else
            cp {input.annotation} {output.annotation}
        fi
        """

rule make_samples_file:
    input: 
        config["units"]
    output: 
        "../config/samples.tsv"
    shell: # 去除第一行，提取前两列，生成samples.tsv文件
        "awk 'NR>1 {{print $1\"\\t\"$2}}' {input} > {output}"
