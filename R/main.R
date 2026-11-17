#' @title Encode DNA k-mers into per-position nucleotide features
#'
#' @description Splits each DNA sequence into single nucleotides and returns
#'   one factor column per nucleotide position. The columns are named nt_pos1,
#'   nt_pos2, ... and use the level order A, T, C, G that the model was
#'   trained with.
#'
#' @param dna_strings A character vector of DNA sequences that all have the
#'   same length, for example c("GGACA", "AAACT"). The length of the first
#'   element determines how many columns are returned.
#'
#' @return A data.frame with one row per input sequence and one factor column
#'   per nucleotide position (nt_pos1 to nt_pos5 for a 5-mer). Every column
#'   has the levels "A", "T", "C" and "G".
#'
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "AAACT"))
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' @title Predict m6A sites for a table of candidate adenosines
#'
#' @description Applies the bundled random forest model to a feature table and
#'   adds two columns: the predicted probability of being an m6A site and the
#'   Positive/Negative call obtained by comparing it with positive_threshold.
#'
#' @param ml_fit A fitted randomForest model that predicts the classes
#'   "Positive" and "Negative", as stored in inst/extdata/rf_fit.rds.
#'
#' @param feature_df A data.frame with one row per candidate site. It must
#'   contain the columns gc_content, RNA_type, RNA_region, exon_length,
#'   distance_to_junction, evolutionary_conservation and DNA_5mer. Note that
#'   exon_length and distance_to_junction must already be log2-transformed,
#'   and DNA_5mer must be a five-nucleotide string with the candidate
#'   adenosine at the central position.
#'
#' @param positive_threshold A single number between 0 and 1. Sites whose
#'   predicted probability is strictly greater than this value are called
#'   "Positive"; sites exactly at the threshold are called "Negative". The
#'   default is 0.5.
#'
#' @return feature_df with two added columns: predicted_m6A_prob (numeric, the
#'   probability of the "Positive" class) and predicted_m6A_status
#'   (character, either "Positive" or "Negative").
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, example_df[1:5, ], positive_threshold = 0.6)
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}
#' @title Predict the m6A status of a single candidate adenosine
#'
#' @description A convenience wrapper around prediction_multiple() that takes
#'   the seven features of one site as separate arguments instead of a data
#'   frame.
#'
#' @param ml_fit A fitted randomForest model that predicts the classes
#'   "Positive" and "Negative", as stored in inst/extdata/rf_fit.rds.
#'
#' @param gc_content A number between 0 and 1 giving the GC content of the
#'   sequence window around the candidate adenosine.
#'
#' @param RNA_type A character string, one of "mRNA", "lincRNA", "lncRNA" or
#'   "pseudogene".
#'
#' @param RNA_region A character string, one of "CDS", "intron", "3'UTR" or
#'   "5'UTR".
#'
#' @param exon_length A number giving the log2-transformed exon length that
#'   contains the candidate site.
#'
#' @param distance_to_junction A number giving the log2-transformed distance
#'   from the candidate site to the nearest exon-intron junction.
#'
#' @param evolutionary_conservation A number between 0 and 1 giving the
#'   conservation score of the candidate site.
#'
#' @param DNA_5mer A five-nucleotide character string such as "GGACA", with
#'   the candidate adenosine at the central position.
#'
#' @param positive_threshold A single number between 0 and 1; sites whose
#'   probability is strictly greater than it are called "Positive". The
#'   default is 0.5.
#'
#' @return A named character vector of length two: predicted_m6A_prob (the
#'   probability, formatted as text because c() coerces the vector to
#'   character) and predicted_m6A_status ("Positive" or "Negative").
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit,
#'                   gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA")
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
