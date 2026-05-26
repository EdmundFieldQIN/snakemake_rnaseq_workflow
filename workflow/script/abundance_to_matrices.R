suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(edgeR)
})

option_list <- list(
  make_option("--est_method", type = "character",
              help = "featureCounts|RSEM|eXpress|kallisto|salmon"),
  make_option("--quant_files", type = "character", default = "",
              help = "file containing one quant file path per line"),
  make_option("--name_sample_by_basedir", action = "store_true", default = FALSE,
              help = "use dirname component as sample name"),
  make_option("--basedir_index", type = "integer", default = -2,
              help = "0-based index; negative allowed, default: -2"),
  make_option("--out_prefix", type = "character", default = "matrix",
              help = "output prefix, default: matrix")
)

parser <- OptionParser(option_list = option_list)
parsed <- parse_args(parser, positional_arguments = TRUE)
opt <- parsed$options
argv_files <- parsed$args

valid_methods <- c("featurecounts", "rsem", "express", "kallisto", "salmon")
method <- tolower(opt$est_method %||% "")
if (!(method %in% valid_methods)) {
  stop("Error: --est_method must be one of featureCounts|RSEM|eXpress|kallisto|salmon")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

get_files <- function(quant_files, argv_files) {
  if (nzchar(quant_files)) {
    files <- readLines(quant_files, warn = FALSE)
    files <- files[nzchar(trimws(files))]
    return(files)
  }
  if (length(argv_files) > 0) {
    return(argv_files)
  }
  stop("Error: provide quant files via positional args or --quant_files")
}

get_colspec <- function(m) {
  if (m == "rsem") {
    list(acc = "transcript_id", count = "expected_count", fpkm = "FPKM", tpm = "TPM")
  } else if (m == "express") {
    list(acc = "target_id", count = "eff_counts", fpkm = "fpkm", tpm = "tpm")
  } else if (m == "kallisto") {
    list(acc = "target_id", count = "est_counts", fpkm = "tpm", tpm = "tpm")
  } else if (m == "salmon") {
    list(acc = "Name", count = "NumReads", fpkm = "TPM", tpm = "TPM")
  } else if (m == "featurecounts") {
    list(acc = "gene_id", count = "counts", fpkm = "fpkm", tpm = "tpm")
  } else {
    stop("Unsupported est_method")
  }
}

perl_index_pick <- function(parts, idx) {
  pos <- if (idx < 0) length(parts) + idx + 1 else idx + 1
  if (pos < 1 || pos > length(parts)) {
    stop("Error: basedir_index out of range for path: ", paste(parts, collapse = "/"))
  }
  parts[[pos]]
}

sample_name_from_path <- function(path, by_basedir, basedir_index) {
  if (by_basedir) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    nm <- perl_index_pick(parts, basedir_index)
  } else {
    nm <- basename(path)
  }
  nm <- sub("\\.(genes|isoforms)\\.results$", "", nm)
  nm <- sub("\\.count$", "", nm)
  nm
}

load_quant <- function(file, colspec) {
  dt <- fread(file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  need <- c(colspec$acc, colspec$count, colspec$fpkm, colspec$tpm)
  miss <- setdiff(need, colnames(dt))
  if (length(miss) > 0) {
    stop("Error in file ", file, ": missing columns: ", paste(miss, collapse = ", "))
  }

  out <- data.frame(
    gene = as.character(dt[[colspec$acc]]),
    count = suppressWarnings(as.numeric(dt[[colspec$count]])),
    fpkm = suppressWarnings(as.numeric(dt[[colspec$fpkm]])),
    tpm = suppressWarnings(as.numeric(dt[[colspec$tpm]])),
    stringsAsFactors = FALSE
  )
  out
}

merge_metric <- function(dfs, samples, metric_name) {
  pieces <- Map(
    function(df, s) {
      x <- df[, c("gene", metric_name), drop = FALSE]
      colnames(x)[2] <- s
      x
    },
    dfs, samples
  )
  merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE, sort = FALSE), pieces)
  merged
}

write_matrix <- function(df, file) {
  fwrite(df, file = file, sep = "\t", quote = FALSE, na = "NA")
}

files <- get_files(opt$quant_files, argv_files)
if (length(files) == 0) {
  stop("Error: no quant files found")
}

colspec <- get_colspec(method)

sample_names <- vapply(
  files,
  sample_name_from_path,
  character(1),
  by_basedir = isTRUE(opt$name_sample_by_basedir),
  basedir_index = opt$basedir_index
)

if (any(duplicated(sample_names))) {
  stop("Error: duplicated sample names: ", paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "))
}

quant_list <- lapply(files, load_quant, colspec = colspec)

counts_mat <- merge_metric(quant_list, sample_names, "count")
fpkm_mat <- merge_metric(quant_list, sample_names, "fpkm")
tpm_mat <- merge_metric(quant_list, sample_names, "tpm")

counts_file <- paste0(opt$out_prefix, ".counts.matrix")
fpkm_file <- paste0(opt$out_prefix, ".FPKM.EXPR.matrix")
tpm_file <- paste0(opt$out_prefix, ".TPM.EXPR.matrix")
tmm_file <- paste0(opt$out_prefix, ".TMM.EXPR.matrix")
tpm_tmm_file <- paste0(opt$out_prefix, ".TPM.TMM.EXPR.matrix")

write_matrix(counts_mat, counts_file)
write_matrix(fpkm_mat, fpkm_file)
write_matrix(tpm_mat, tpm_file)

# 对count矩阵做TMM标准化
count_numeric <- as.matrix(counts_mat[, -1, drop = FALSE])
mode(count_numeric) <- "numeric"
count_numeric[is.na(count_numeric)] <- 0
count_numeric <- round(count_numeric)

if (ncol(count_numeric) >= 2) {
  dge <- DGEList(counts = count_numeric, group = factor(colnames(count_numeric)))
  dge <- calcNormFactors(dge, method = "TMM")
  tmm_numeric <- t(t(dge$counts) / dge$samples$lib.size / dge$samples$norm.factors) * 1000000
  # tmm_numeric <- sweep(count_numeric, 2, dge$samples$norm.factors, "/")
} else {
  tmm_numeric <- count_numeric
}

tmm_mat <- data.frame(gene = counts_mat$gene, tmm_numeric, check.names = FALSE)
write_matrix(tmm_mat, tmm_file)

# 对TPM矩阵做TMM标准化
tpm_numeric <- as.matrix(tpm_mat[, -1, drop = FALSE])
mode(tpm_numeric) <- "numeric"
tpm_numeric[is.na(tpm_numeric)] <- 0

if (ncol(tpm_numeric) >= 2) {
  dge <- DGEList(counts = tpm_numeric, group = factor(colnames(tpm_numeric)))
  dge <- calcNormFactors(dge, method = "TMM")
  tmm_tpm_numeric <- t(t(dge$counts) / dge$samples$lib.size / dge$samples$norm.factors) * 1000000
} else {
  tmm_tpm_numeric <- tpm_numeric
}

tmm_tpm_mat <- data.frame(gene = tpm_mat$gene, tmm_tpm_numeric, check.names = FALSE)
write_matrix(tmm_tpm_mat, tpm_tmm_file)
