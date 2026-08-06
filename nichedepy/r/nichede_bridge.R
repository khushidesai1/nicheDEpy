suppressPackageStartupMessages({
  library(jsonlite)
  library(Matrix)
  library(NicheDE)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: nichede_bridge.R <command> <params.json>")
}

command <- args[[1]]
params <- fromJSON(args[[2]], simplifyVector = TRUE)

read_table_vector <- function(path, column) {
  data <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data[[column]]
}

read_inputs <- function(input_dir) {
  counts <- readMM(file.path(input_dir, "counts.mtx"))
  cells <- read_table_vector(file.path(input_dir, "cells.tsv"), "cell_id")
  genes <- read_table_vector(file.path(input_dir, "genes.tsv"), "gene")
  dimnames(counts) <- list(cells, genes)

  coords <- read.csv(file.path(input_dir, "coords.csv"), row.names = 1,
                     check.names = FALSE)
  labels <- read.csv(file.path(input_dir, "labels.csv"), stringsAsFactors = FALSE)
  deconv <- read.csv(file.path(input_dir, "deconv.csv"), row.names = 1,
                     check.names = FALSE)

  list(counts = counts, coords = coords, labels = labels, deconv = as.matrix(deconv))
}

write_result <- function(result, output_path) {
  if (is.null(output_path) || is.na(output_path)) {
    return(invisible(NULL))
  }
  if (is.data.frame(result) || is.matrix(result)) {
    write.csv(result, output_path, row.names = TRUE)
  } else {
    saveRDS(result, output_path)
  }
}

if (command == "create_library_matrix") {
  inputs <- read_inputs(params$input_dir)
  cell_type <- cbind(inputs$labels$cell_id, inputs$labels$cell_type)
  library_matrix <- CreateLibraryMatrix(as.matrix(inputs$counts), cell_type)
  write.csv(library_matrix, params$output_path, row.names = TRUE)

} else if (command == "create_object") {
  inputs <- read_inputs(params$input_dir)
  library_matrix <- read.csv(params$library_matrix_path, row.names = 1,
                             check.names = FALSE)
  nde_obj <- CreateNicheDEObject(
    as.matrix(inputs$counts),
    inputs$coords,
    as.matrix(library_matrix),
    inputs$deconv,
    sigma = as.numeric(params$sigma)
  )
  saveRDS(nde_obj, params$output_path)

} else if (command == "calculate_effective_niche") {
  nde_obj <- readRDS(params$object_path)
  if (isTRUE(params$large_scale)) {
    nde_obj <- CalculateEffectiveNicheLargeScale(
      nde_obj,
      batch_size = as.integer(params$batch_size),
      cutoff = as.numeric(params$cutoff)
    )
  } else {
    nde_obj <- CalculateEffectiveNiche(nde_obj)
  }
  saveRDS(nde_obj, params$output_path)

} else if (command == "niche_de") {
  nde_obj <- readRDS(params$object_path)
  nde_obj <- niche_DE(
    nde_obj,
    num_cores = as.integer(params$num_cores),
    outfile = params$outfile,
    C = as.numeric(params$C),
    M = as.numeric(params$M),
    gamma = as.numeric(params$gamma),
    print = isTRUE(params$print),
    Int = isTRUE(params$Int),
    batch = isTRUE(params$batch),
    self_EN = isTRUE(params$self_EN),
    G = as.numeric(params$G)
  )
  saveRDS(nde_obj, params$output_path)

} else if (command == "get_niche_de_genes") {
  nde_obj <- readRDS(params$object_path)
  genes <- get_niche_DE_genes(
    nde_obj,
    test.level = params$test_level,
    index = params$index,
    niche = params$niche,
    pos = isTRUE(params$positive),
    alpha = as.numeric(params$alpha)
  )
  write.csv(genes, params$output_path, row.names = TRUE)

} else if (command == "call_function") {
  fn <- get(params$function_name, envir = asNamespace("NicheDE"))
  positional <- params$args
  keyword <- params$kwargs
  if (!is.null(params$object_path) && !is.na(params$object_path)) {
    positional <- c(list(readRDS(params$object_path)), positional)
  }
  result <- do.call(fn, c(positional, keyword))
  write_result(result, params$output_path)

} else {
  stop(paste("Unknown command:", command))
}
