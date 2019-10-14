# Gene Capture Assembly Pipeline

This is a simple pipeline that aligns gene capture data to a transcriptome reference file (Trinity), generating a consensus 'gene sequence' at the other end.

## Usage

The pipeline has the following arguments:

**Key-value file (-kv)**: This is a two columned CSV file that has the columns `sample` and `reference`. For each sample, match it to a **unique** basename of the reference it belongs to. Do not provide the file extension.
    -kv /path/to/key_value.csv

**Reads file directory (-rd)**: A directory path to where the gene capture reads are located.
  
- -rd /path/to/reads_dir

**Reference file directory (-fd)**: A directory path to where the transcriptomic (or any reference) files are located.
  
- -fd /path/to/references_dir

**Output directory (-o)**: Provide the directory path to where you'd like the output of the pipeline to be stored.
  
- -o /path/to/output_dir

**Adapter file (-af)**: Provide the file path to a fasta file containing the adapter sequences used for your sequencing. If you do not have them, I've included a comprehensive adapter file as part of this repo.

- -af /path/to/adapters.fa

**Threads**: If you have a multi-threaded machine, provide the number of threads you'd like to use in multi-threaded processes (e.g. FastQC, Trimmomatic and BWA)

- -j 8
