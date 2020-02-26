#!/usr/bin/env Rscript

library(argparse)

## Arguments to script
parser <- ArgumentParser(description = "Filter BLAST output for best hits")

## Defining arguments
parser$add_argument("-b", "--blastDir", 
                    help = "Path to BLAST output directory")
parser$add_argument("-e", "--blastExt", 
                    help = "Extension of BLAST files [default: %(default)s]",
                    default = '.outfmt6')
parser$add_argument("-f", "--fastaDir", 
                    help = "Path to fasta files ")
parser$add_argument("-fe", "--fastaExt",
                    help = "Extension of fasta files [default: %(default)s]",
                    default = '.fasta')
parser$add_argument("-o", "--output",
                    help = "Output directory",
                    default = NULL)
parser$add_argument("-c", "--blastColumnNames", 
                    help = "Vector of column names if a custom blast table is used",
                    nargs = '*',
                    default = NULL)

## Initialising arguments
args <- parser$parse_args()

args$blastDir <- '~/Documents/kate/2002_matt-cds/03_blast-superTranscripts'
args$blastExt <- '.outfmt6'
args$blastColumnNames <- c('qseqid','sseqid','pident','qlen','slen','length','qcovs','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
args$fastaDir <- '~/Documents/kate/2002_matt-cds/02_annotated-superTranscripts'
args$fastaExt <- '.SuperTrans.anno.cds'

## Check if required packages are installed
if(! 'BiocManager' %in% rownames(installed.packages())){install.packages("BiocManager")}
if(! 'Biostrings' %in% rownames(installed.packages())){BiocManager::install("Biostrings")}

suppressPackageStartupMessages(library(tidyverse))

## Check input
assertthat::assert_that(! is.null(args$output), msg = 'No output directory provided... Assign value to `-o`')
if(! dir.exists(args$output)){dir.create(path = args$output, showWarnings = FALSE, recursive = TRUE)}

## Assign default column names if none are given
if(is.null(args$blastColumnNames)){
  
  ## Default column names
  colNames <- c('qseqid','sseqid','pident','length','mismatch',
                'gapopen','qstart','qend','sstart','send','evalue','bitscore')
  msg <- paste0('Setting default BLAST column names: ', paste0(colNames, collapse = ','), '...')
  print(msg)
  args$blastColumnNames <- colNames
  
}

## Normalising paths
args$blastDir <- normalizePath(args$blastDir)
args$fastaDir <- normalizePath(args$fastaDir)

## Listing BLAST files
print('Importing BLAST tables...')
b <- list.files(path = args$blastDir, pattern = args$blastExt, full.names = TRUE, recursive = TRUE)
regx <- paste0('^(.*)', args$blastExt)
names(b) <- sub(pattern = regx, "\\1", basename(b))
b <- map(b, read_tsv, col_types = cols(), col_names = args$blastColumnNames)

## CDS sequences
print('Importing fasta sequences...')
f <- list.files(path = args$fastaDir, pattern = args$fastaExt, full.names = TRUE, recursive = TRUE)
regx <- paste0("(.*)", args$fastaExt)
names(f) <- sub(regx, "\\1", basename(f))
f <- map(f, Biostrings::readDNAStringSet, use.names = TRUE)

## Check that the names of blast tables and fasta files match
logi <- all(names(f) %in% names(b))
assertthat::assert_that(logi, msg = 'Basenames of BLAST tables don\'t mach the basenames of the fasta files... Check your filenames or the extensions you\'ve passed')
 
## Get the best hits
print('Getting best BLAST hits for each transcript...')
bestHits <- map(.x = b, ~{

  o <- .x %>%
    separate(col = sseqid, into = c("org", "symbol"), sep = "_", extra = "merge") %>% # Separate subject into organism and gene symbol
    mutate(symbol = sub("_.*", "", symbol), # Remove everything after the underscores - will be junk
           symbol = toupper(symbol)) %>% # Convert gene symbols to upper case
    separate(col = qseqid, into = c("trinity", 'type'), sep = '::') %>%
    separate(col = trinity, into = c("trinity", "ann"), sep = "\\^", remove = FALSE) %>%
    mutate(symbol = sub('PCDH21', 'CDHR1', symbol),
           len_prop = (qlen/slen) * 100,
           scovs = (length/slen) * 100) %>%
    group_by(trinity, symbol) %>%
    arrange(trinity, symbol, desc(bitscore)) %>%
    slice(1) %>%
    arrange(symbol, desc(bitscore)) %>% # Arrange by gene symbol and bit score
    group_by(symbol) %>%
    slice(1) %>% # Take the first hit - best hit based on bit score
    ungroup()

})

## Subset fasta by best hits
print('Subsetting fasta files...')
seqSubset <- map2(bestHits, f, ~{

  ## ID to subset by
  g <- paste(.x[["trinity"]], .x[["type"]], sep = '::')

  ## ID with gene annotation (new name)
  n <- paste(.x[["trinity"]], .x[["symbol"]], sep = "_")
  
  ## Subset fa
  s <- .y[g]
  names(s) <- n ## Assigning new names with annotation
  
  return(s)

})

## Writing fasta to file
walk(.x = names(seqSubset), ~{
  seqSubset[[.x]] %>%
    Biostrings::writeXStringSet(filepath = paste0(args$output, '/', .x, ".fasta"))
})
