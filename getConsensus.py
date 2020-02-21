#!/usr/bin/env python3

import os
import sys
import glob
import logging
import argparse
import subprocess
import pandas as pd

def checkOutputExists(outPath, ext, readsList):
    
    ## List sample specific outputs
    b = os.path.basename(readsList[0][0])
    b = b.replace('_L001_R1.fastq.gz', '')
    chk = glob.glob(outPath + '/' + b + '*' + ext)
    
    ## If files exist and have size > 0 --> skip
    if len(chk) > 0 and os.path.exists(chk[0]) and os.path.getsize(chk[0]) > 0:
        return "Skipping: Output exists for" + ' ' + b
    else:
        return 0

def getSeqData(keyValuePath, readsPath, refPath):

    ## Importing key-value
    try:
        df = pd.read_csv(keyValuePath)
    except pd.errors.ParserError as e:
        print("ERROR: Unable to read key-value CSV file", keyValuePath)
        print(e)
    
    ## Iterate over rows listing files
    lst = []
    for idx, row in df.iterrows():
        reads = sorted(glob.glob(readsPath + '/' + row['sample'] + '*')) # list fastq files
        ref = glob.glob(refPath + '/' + row['reference'] + '*' + '.fasta') # list reference file
        
        lst.append([reads, ref, row['sample']])

    return lst


def runFastQC(readsList, outDir, rawTrim, threads):
    
    ## Output directory
    out = outDir + '/FastQC_' + rawTrim
    
    ## Create output directory if it doesn't exist
    if not os.path.exists(out):
        os.makedirs(out, exist_ok=True)

    ## Check if output exists in directory + isn't empty
    chk = checkOutputExists(out, 'html', readsList)
    if chk != 0:
        print(chk)
        return

    ## Raw or trimmed data - conditional
    if rawTrim == 'raw':
        ## Raw data
        R1 = readsList[0][0]
        R2 = readsList[0][1]

        ## Call to FastQC at the command line
        subprocess.run(['fastqc', '-k', '9', '-t', str(threads), '-o',
                         out, R1, R2], stderr=subprocess.DEVNULL)
    else:
        ## Getting trimmed files
        f = sorted(glob.glob(outDir + '/trimmomatic' + '/' + readsList[2] + '*' + '.P.qtrim.fastq.gz'))
        
        R1 = f[0]
        R2 = f[1]

        ## Call to FastQC
        subprocess.run(['fastqc', '-k', '9', '-t', str(threads), '-o',
                         out, R1, R2], stderr=subprocess.DEVNULL)


def runTrimmomatic(readsList, outDir, adapterFile, threads):
    
    ## Output directory
    out = outDir + '/trimmomatic'
    if not os.path.exists(out):
        os.makedirs(out, exist_ok=True)

    ## Check if output exists in directory + isn't empty
    chk = checkOutputExists(out, 'summary', readsList)
    if chk != 0:
        print(chk)
        return

    ## Building read input/outputs
    R1 = readsList[0][0]
    R1_out = out + '/' + readsList[2] + '_R1.P.qtrim.fastq.gz'
    R1_unpaired = out + '/' + readsList[2] + '_R1.U.qtrim.fastq.gz'
    
    R2 = readsList[0][1]
    R2_out = out + '/' + readsList[2] + '_R2.P.qtrim.fastq.gz'
    R2_unpaired = out + '/' + readsList[2] + '_R2.U.qtrim.fastq.gz'

    ## Extra paramters to pass to Trimmomatic
    ic = 'ILLUMINACLIP:' + adapterFile + ':2:30:10'
    log = out + '/' + readsList[2] + '.log'
    summary = out + '/' + readsList[2] + '.summary'

    ## Calling Trimmomatic
    subprocess.run(['trimmomatic', 'PE', '-threads', str(threads), 
                    '-trimlog', log, '-summary', summary, R1, R2, 
                    R1_out, R1_unpaired, R2_out, R2_unpaired, ic, 'SLIDINGWINDOW:4:5', 'LEADING:5', 'TRAILING:5', 'MINLEN:25'], 
                    stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


def runBWAmem(readsList, outDir, threads, scriptDir):

    ## Create output directory
    outRAW = outDir + '/bwa_aligned'
    outFILT = outDir + '/bwa_aligned/filtered'
    if not os.path.exists(outRAW):
        os.makedirs(outRAW, exist_ok=True)

    if not os.path.exists(outFILT):
        os.makedirs(outFILT, exist_ok=True)

    ## Check if output exists in directory + isn't empty
    chk = checkOutputExists(outFILT, 'bam', readsList)
    if chk != 0:
        print(chk)
        return

    ## List trimmed files for alignment
    f = sorted(glob.glob(outDir + '/trimmomatic' + '/' + readsList[2] + '*' + '.P.qtrim.fastq.gz'))
    
    ## Arguments for BWA alignment script
    reference = ''.join(readsList[1]) ## Removing square brackets
    R1 = f[0]
    R2 = f[1]
    outBAM = outRAW + '/' + readsList[2] + '.bam'
    outBAM_stats = outRAW + '/' + readsList[2] + '.flagstat'
    outBAMFILT = outFILT + '/' + readsList[2] + '_filtered_sorted.bam'
    script = scriptDir + "/bwa_align.sh"

    ## Run BWA script
    subprocess.run( [ script, reference, R1, R2, outBAM, outBAM_stats, outBAMFILT, str(threads) ], stderr=subprocess.DEVNULL )


def getConsensus(readsList, outDir):
    
    ## Arguments
    out = outDir + '/consensus'
    if not os.path.exists(out):
        os.makedirs(out, exist_ok=True)

    ## Check if output exists in directory + isn't empty
    chk = checkOutputExists(out, 'fasta', readsList)
    if chk != 0:
        print(chk)
        return

    ## Inputs/outputs
    reference = ''.join(readsList[1])
    bam = ''.join(glob.glob(outDir + '/bwa_aligned/filtered' + '/' +
                           readsList[2] + '*' + '_filtered_sorted.bam'))
    vcf = out + '/' + readsList[2] + '.vcf.gz'
    consensus = out + '/' + readsList[2] + '.fasta'

    ## Calling genotypes
    bcf_call = "bcftools mpileup -Ou -d 5000 -Q 20 -q 20 -f '%s' '%s' | \
        bcftools call -c - | \
        bcftools norm -f '%s' -m +any -Oz -o '%s'" % (reference, bam, reference, vcf)
    os.system(bcf_call) ## Calling genotypes
    subprocess.run(['bcftools', 'index', vcf],
                    stderr=subprocess.DEVNULL)  # Indexing vcf file

    ## Calling consensus using VCF
    subprocess.run(['bcftools', 'consensus', '-f',
                     reference, '-H', '1', '-o', consensus, vcf], stderr=subprocess.DEVNULL)


def main():

    ## Set up argument parser
    desc = str("Pipeline to assemble gene capture data using closely related references")

    epi = str("Developed: Alastair Ludington\n" +
              "Institution: University of Adelaide\n" +
              "Date: 11/10/2019")

    parser = argparse.ArgumentParser(
        description=desc,
        epilog=epi,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-kv", "--keyValuePath",
                        help="File path to key-value table",
                        required=True,
                        type=str)

    parser.add_argument("-rd", "--readsDir",
                        help="Directory path to gene capture data",
                        required=True,
                        type=str)

    parser.add_argument("-fd", "--refDir",
                        help="Directory path to reference fasta files",
                        required=False,
                        type=str)

    parser.add_argument("-o", "--outputDir",
                        help="Directory path for output location",
                        required=True,
                        type=str)

    parser.add_argument("-af", "--adapterFile",
                        help="File path to adapter sequence file for Trimmomatic",
                        required=True,
                        type=str)

    parser.add_argument("-j", "--threads",
                        help="Number of threads for different tasks (default: 2)",
                        required=False,
                        default=2,
                        type=int)

    args = parser.parse_args()

    ## Create ouput directory for logger
    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir, exist_ok=True)

    ## Create logger
    logFile = args.outputDir + '/' + 'geneCapture.log'
    FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\t%(message)s"
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=logFile, format=FORMAT, filemode='w', level=logging.DEBUG)

    ## Handler to print everything to stdout
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    ## Printing arguments
    print("\n" + "*" * 80)
    print("Passed Arguments:\n")
    for arg in vars(args):
        print(arg, getattr(args, arg))

    ## Running the pipeline

    ## Directory this script resides in
    scriptDir = os.path.dirname(os.path.realpath(__file__))

    print("\n" + "*" * 80)
    logger.info("Getting sequence inputs (reads + references)")
    seqData = getSeqData(args.keyValuePath, args.readsDir, args.refDir)
    
    logger.info("FastQC: Raw data")
    # print("\nFastQC: Raw data")
    for i in seqData:
        runFastQC(i, args.outputDir, 'raw', args.threads)

    logger.info("Trimmomatic: Trimming reads (adapters + quality)")
    for i in seqData:
        runTrimmomatic(i, args.outputDir, args.adapterFile, args.threads)

    logger.info("FastQC: Trimmed data")
    for i in seqData:
        runFastQC(i, args.outputDir, 'trim', args.threads)
    
    logger.info("BWA-mem: Aligning reads")
    for i in seqData:
        runBWAmem(i, args.outputDir, args.threads, scriptDir + '/scripts')

    logger.info("BCFtools: Calling consensus sequence for genes")
    for i in seqData:
        getConsensus(i, args.outputDir)

    logger.info("GeneCapture pipeline: Complete")


if __name__ == "__main__":
    main()
