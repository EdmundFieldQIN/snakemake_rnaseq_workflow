rule de_analysis:
    input: 
        martrix="05expression/genes.counts.matrix",
        contrast=config["contrast"],
        samples_file=rules.make_samples_file.output
    output: 
        directory("06de_analysis/DESeq2_Result/")
    conda: "rnaseq"
    shell: # 读取config文件的units.tsv忽略表头并提取前两列重新生成samples.tsv文件，作为DE分析--samples_files的输入
        """
        mkdir -p 06de_analysis
        run_DE_analysis.pl --matrix {input.martrix} --method DESeq2 --samples_file {input.samples_file} --contrast {input.contrast} --output {output} 
        """