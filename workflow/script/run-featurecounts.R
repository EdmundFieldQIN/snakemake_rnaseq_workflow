#!/usr/bin/env Rscript

# 导入必要的库
library(argparser, quietly=TRUE)
library(Rsubread)
library(limma)
library(edgeR)

# 创建参数解析器
p <- arg_parser("Run featureCounts and calculate FPKM/TPM")

# 添加参数
p <- add_argument(p, "--bam", help="input: bam file", type="character")
p <- add_argument(p, "--strandSpecific", help="Strand specificity (0, 1, 2)", type="numeric", default=0)
p <- add_argument(p, "--gtf", help="input: gtf file", type="character", default=NA)
p <- add_argument(p, "--saf", help="input: saf file", type="character", default=NA)
p <- add_argument(p, "--isPairedEnd", help="Paired-end reads, TRUE or FALSE", type="logical", default=TRUE)
p <- add_argument(p, "--featureType", help="Feature type in GTF annotation (default: exon)", type="character", default="exon")
p <- add_argument(p, "--attrType", help="Attribute type in GTF annotation (default: gene_id)", type="character", default="gene_id")
p <- add_argument(p, "--output", help="Output prefix", type="character")

# 解析命令行参数
argv <- parse_args(p)

# 检查 gtf 和 saf 参数是否互斥
if (!is.na(argv$gtf) && !is.na(argv$saf)) {
  stop("Cannot use both --gtf and --saf. Choose one.")
}

# 主要计算流程
bamFile <- argv$bam
annotFile <- if (!is.na(argv$gtf)) argv$gtf else argv$saf
isGTF <- !is.na(argv$gtf)
outFilePref <- argv$output

# 输出文件路径
outStatsFilePath  <- paste0(outFilePref, '.log')
outCountsFilePath <- paste0(outFilePref, '.count')

# 运行 featureCounts
fCountsArgs <- list(
  files = bamFile,
  annot.ext = annotFile,
  nthreads = 1,
  isPairedEnd = argv$isPairedEnd,
  strandSpecific = argv$strandSpecific
)

if (isGTF) {
  fCountsArgs$isGTFAnnotationFile <- TRUE
  fCountsArgs$GTF.featureType <- argv$featureType
  fCountsArgs$GTF.attrType <- argv$attrType
} else {
  fCountsArgs$isGTFAnnotationFile <- FALSE
}

fCountsList <- do.call(featureCounts, fCountsArgs)

# 创建 DGEList
dgeList <- DGEList(counts = fCountsList$counts, genes = fCountsList$annotation)

# 计算表达量
cpm <- cpm(dgeList)
fpkm <- rpkm(dgeList, dgeList$genes$Length)
tpm <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

# 写入输出文件
write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# 使用 --attrType 参数的值作为第一列的列名
firstColName <- argv$attrType
featureCounts <- cbind(fCountsList$annotation[,1], fCountsList$counts, fpkm, tpm, cpm)
colnames(featureCounts) <- c(firstColName, 'counts', 'fpkm', 'tpm', 'cpm')

write.table(featureCounts, outCountsFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


