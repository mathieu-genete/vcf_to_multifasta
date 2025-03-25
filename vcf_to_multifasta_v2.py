#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from collections import defaultdict
from Bio import SeqIO
import vcf
import random
import time
import gzip
from io import StringIO

__version__="2.0.0"

def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed arguments object containing user inputs.
    """
    description = "VCF to multifasta"
    parser = argparse.ArgumentParser(prog="vcf_to_multifasta_v2.py", description=description)
    # Define arguments with descriptions
    parser.add_argument("-r", "--reffasta", help="Fasta reference file")
    parser.add_argument("-c", "--contigslist", help="text file with contig or chromosomes name to compute")
    parser.add_argument("-v", "--vcffile", help="Input VCF file (.vcf or .vcf.gz)", required=True)
    parser.add_argument("-o", "--outfolder", help="Output folder for fasta files", required=True)
    parser.add_argument("-n", "--popname", help="Tabulated file containing individual and population names")
    parser.add_argument("-s", "--randomSeed", help="Random seed for reproducibility", type=int)
    parser.add_argument("-p", "--ploidy", help="Set the ploidy of samples (default=2)", type=int, default=2)
    #parser.add_argument("-b", "--biallelic", help="Force biallelic SNPs check", action="store_true")
    return parser.parse_args()

def setup_environment(outfolder, randomSeed):
    """
    Set up the output folder and initialize random seed.

    Args:
        outfolder (str): Path to the output folder for saving results.
        randomSeed (int or None): Random seed for reproducibility.
    """
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    random.seed(randomSeed or int(time.time()))

def load_contigs_file(contigs_file):
    """
    Load contig names from a file.

    Args:
        contigs_file (str): Path to the file containing contig names.

    Returns:
        list: A list of contig names.
    """
    contigs_list=[]
    if contigs_file:
        with open(contigs_file,'r') as inContigs:
            for line in inContigs:
                contigs_list.append(line.strip())
    return contigs_list

def load_population_file(popname_file):
    """
    Load population names from a tab-delimited file.

    Args:
        popname_file (str): Path to the population name file.

    Returns:
        dict: A dictionary mapping individuals to population names.
    """
    popname = {}
    if popname_file:
        with open(popname_file, 'r') as infile:
            for line in infile:
                tmp = line.strip().split('\t')
                if len(tmp) >= 2:
                    popname[tmp[0]] = tmp[1]
    return popname
    
def parse_vcf(vcf_reader,ploidy,outfolder,ref_fasta_datas,contigs_length,popname,contigs_list):
    """
    Parse the VCF file and generate multifasta files.

    Args:
        vcf_reader (vcf.Reader): VCF reader object.
        ploidy (int): Ploidy of samples.
        outfolder (str): Output folder for fasta files.
        ref_fasta_datas (dict): Fasta reference data indexed by contig names.
        contigs_length (dict): Dictionary with contig lengths.
        popname (dict): Population names mapped to individuals.
        contigs_list (list): List of contigs to process.
    """
    vcf_datas=defaultdict(dict)
    chr_name=""
    progress=""
    nbr_contigs=len(contigs_length)
    i=0
    print_progress(i,nbr_contigs)
    
    for record in vcf_reader:
        # Skip contigs not in the list if the list is provided
        if len(contigs_list)>0 and record.CHROM not in contigs_list:
            continue
        change_chr = (record.CHROM != chr_name)        
        if change_chr and len(vcf_datas)>0:
            i+=1
            if (i%10==0):
                print_progress(i,nbr_contigs)
                
            write_multifasta(chr_name,outfolder,vcf_datas,ref_fasta_datas,ploidy,contigs_length,popname)
            vcf_datas=defaultdict(dict)
                            
        chr_name = record.CHROM
        process_vcf_samples(record,vcf_datas)
        
    # Write the last contig
    write_multifasta(chr_name,outfolder,vcf_datas,ref_fasta_datas,ploidy,contigs_length,popname)
    print_progress(i,nbr_contigs)

def print_progress(i,nbr_contigs):
    """
    Print progress to the console.

    Args:
        i (int): Current progress index.
        nbr_contigs (int): Total number of contigs.
    """
    if sys.stdout.isatty():
        progress=round((i/nbr_contigs)*100,1)
        print(f'\r{progress}%',end='')
    
def process_vcf_samples(record,vcf_datas):
    """
    Process samples in a VCF record and extract allele data.

    Args:
        record (vcf.model._Record): VCF record.
        vcf_datas (dict): Data structure for storing sample alleles.
    """
    for sample in record.samples:
        genotype = sample['GT']
        if genotype is not None and (genotype != './.' and genotype != '.'):
            # Extract alleles based on indices
            alleles = genotype.split('|') if '|' in genotype else genotype.split('/')
            alleles = [record.REF if a == '0' else str(record.ALT[int(a) - 1]) for a in alleles]
            alleles = [a.replace('*','N') for a in alleles]
            random.shuffle(alleles)

            if sample.sample not in vcf_datas.keys():
                vcf_datas[sample.sample]={}
            # Adjust position as VCF is 1-based
            vcf_datas[sample.sample][record.POS-1]=alleles
            
def write_multifasta(chr_name,outfolder,vcf_datas,ref_fasta_datas,ploidy,contigs_length,popname):
    """
    Write multifasta files for a specific chromosome or contig.

    Args:
        chr_name (str): Name of the chromosome or contig.
        outfolder (str): Directory to save output fasta files.
        vcf_datas (dict): Sample data extracted from the VCF file.
        ref_fasta_datas (dict): Reference fasta data indexed by contig names.
        ploidy (int): Ploidy level of samples.
        contigs_length (dict): Lengths of contigs or chromosomes.
        popname (dict): Dictionary mapping samples to population names.
    """
    out_fasta_datas=defaultdict(list)
    contig_len=contigs_length[chr_name]

    # Retrieve reference sequence for the contig, if available
    if ref_fasta_datas is not None:
        fasta_seq = ref_fasta_datas[chr_name]
        check_length=contig_len==len(fasta_seq.seq)
    else:
        fasta_seq = None
        check_length=True

    # Validate consistency between reference and contig length
    if check_length:
        for sample,datas in vcf_datas.items():
            for i in range(contig_len):
                if i in vcf_datas[sample].keys():
                    out_fasta_datas[sample].append(vcf_datas[sample][i])
                elif fasta_seq is not None:
                    out_fasta_datas[sample].append(fasta_seq.seq[i])
                else:
                    out_fasta_datas[sample].append('N')
                        
    else:
        print(f"VCF and Fasta length of {chr_name} not corresponding")

    # Write sequences to an output fasta file
    out_fasta_file = os.path.join(outfolder,"{}.fa".format(replace_notAllowed_chars(chr_name)))
    with open(out_fasta_file,'w') as out_fasta:
        for sample,seq_list in out_fasta_datas.items():
            for p in range(ploidy):
                if sample in popname.keys():
                    popname_value=popname[sample]
                else:
                    popname_value="NoPop"

                seq=""
                for b in seq_list:
                    if len(b)>1:
                        seq+=b[p]
                    else:
                        seq+=b
                        
                outseq_ID="{gene}|{pop}|{indiv}|allele{anbr}".format(gene=replace_notAllowed_chars(chr_name),pop=popname_value,indiv=sample,anbr=p+1)
                outseq=">{seqID}\n{sequence}\n".format(seqID=outseq_ID,sequence=seq)
                out_fasta.write(outseq)

def replace_notAllowed_chars(txt):
    """
    Replace characters not allowed in filenames or sequence IDs.

    Args:
        txt (str): Input string.

    Returns:
        str: String with replaced characters.
    """
    NotAllowed={'|':'_'}
    return ''.join(c if c not in NotAllowed.keys() else NotAllowed[c] for c in txt)

def is_indel(rec):
    """
    Check if a VCF record contains an indel (insertion or deletion).

    Args:
        rec (vcf.model._Record): VCF record.

    Returns:
        bool: True if the record is an indel, False otherwise.
    """
    indel_REF= sum([len(str(a).replace('None','N')) for a in list(rec.REF)])!=len(list(rec.REF))
    indel_ALT= sum([len(str(a).replace('None','N')) for a in list(rec.ALT)])!=len(list(rec.ALT))
    return indel_REF or indel_ALT

def is_biallelic(rec):
    """
    Check if a VCF record is biallelic.

    Args:
        rec (vcf.model._Record): VCF record.

    Returns:
        bool: True if the record is biallelic, False otherwise.
    """
    alphabet=["A","T","C","G"]
    alt_bases=[r for r in rec.ALT if (r in alphabet)]
    biallelic=False
    if len(alt_bases)==1:
        biallelic = any(sample['GT'] and '0' in sample['GT'] for sample in rec.samples)
    elif len(alt_bases)==2:
        biallelic = not any(sample['GT'] and '0' in sample['GT'] for sample in rec.samples)
    return biallelic
        
def find_Magic_Bytes(filename,magic_bytes):
    """
    Check if a file starts with specific magic bytes.

    Args:
        filename (str): Path to the file.
        magic_bytes (str): String of expected magic bytes.

    Returns:
        bool: True if the file matches the magic bytes, False otherwise.
    """
    with open(filename,'r',encoding="ISO-8859-1") as infile:
        file_start = infile.read(len(magic_bytes))
    if file_start.startswith(magic_bytes):
        return True
    return False

def is_gzip(filename):
    """
    Check if a file is gzip-compressed.

    Args:
        filename (str): Path to the file.

    Returns:
        bool: True if the file is gzip-compressed, False otherwise.
    """
    magic_bytes = "\x1f\x8b\x08"
    return find_Magic_Bytes(filename,magic_bytes)

def load_vcf(vcf_filename):
    """
    Load a VCF file.

    Args:
        vcf_filename (str): Path to the VCF file.

    Returns:
        vcf.Reader: VCF reader object.
    """
    if is_gzip(vcf_filename):
        with gzip.open(vcf_filename, 'rb') as f:
            text_content = f.read().decode("latin-1")
            vcf_reader = vcf.Reader(StringIO(text_content))
    else:
        vcf_reader = vcf.Reader(open(vcf_filename,'r',encoding="latin-1"))

    return vcf_reader

def load_fasta(fasta_name):
    """
    Load and index a fasta file.

    Args:
        fasta_name (str): Path to the fasta file.

    Returns:
        dict: Indexed fasta sequences.
    """
    return SeqIO.index(fasta_name,"fasta")

def check_vcf_contigs(vcf_reader,fasta_datas):
    """
    Check consistency of contig names between VCF and fasta.

    Args:
        vcf_reader (vcf.Reader): VCF reader object.
        fasta_datas (dict): Indexed fasta data.
    """
    fasta_IDs=set(fasta_datas.keys())
    contigs_IDs=set(vcf_reader.contigs.keys())
    contigs_not_found=contigs_IDs.difference(fasta_IDs)

    if len(contigs_not_found)>0:
        print(contigs_not_found)

def convert_to_int(value):
    """
    Convert a value to an integer, handling errors gracefully.

    Args:
        value: Value to convert.

    Returns:
        int: Converted integer value, or 0 on failure.
    """
    try:
        value = int(value)
    except TypeError:
        value = 0
    return value

def verify_files_exist(*file_paths):
    """
    Verify that all provided file paths exist.

    Args:
        file_paths (str): Paths to files that need verification.

    Raises:
        FileNotFoundError: If any file does not exist.
    """
    for file_path in file_paths:
        if file_path and not os.path.exists(file_path):
            raise FileNotFoundError(f"The file '{file_path}' does not exist. Please check the path.")

def main():
    """
    Main function to orchestrate the script execution.
    """
    args = parse_arguments()

    # Check if files exists
    try:
        verify_files_exist(args.reffasta, args.contigslist, args.vcffile, args.popname)
    except FileNotFoundError as e:
        print(e)
        sys.exit(1)
    
    outfolder = args.outfolder
    ploidy = args.ploidy

    # Set up environment
    setup_environment(outfolder, args.randomSeed)

    # Load auxiliary files
    start=time.time()
    popname=load_population_file(args.popname)
    contigs_list=load_contigs_file(args.contigslist)
    vcf_reader = load_vcf(args.vcffile)
    all_samples = vcf_reader.samples
    contigs_length={contig:int(v.length) for contig,v in vcf_reader.contigs.items()}
    contigs_nbr=len(contigs_length)

    # Load fasta reference if provided
    if args.reffasta:
        ref_fasta_datas = load_fasta(args.reffasta)
        check_vcf_contigs(vcf_reader,ref_fasta_datas)
    else:
        ref_fasta_datas = None

    print("== VCF to multifasta -- version {} ==".format(__version__))
    print("Processing VCF file : {}".format(os.path.basename(args.vcffile)))
    print("ploidy : {}".format(ploidy))
    print("{} contigs\n".format(contigs_nbr))
    print("{} individuals : {}\n".format(len(all_samples),', '.join(all_samples)))
    print("Populations list : {}\n".format(', '.join(sorted(set(popname.values())))))

    # Parse the VCF file and generate multifasta files
    print("*** Start parsing VCF file ***")
    parse_vcf(vcf_reader,ploidy,outfolder,ref_fasta_datas,contigs_length,popname,contigs_list)
    print("\nProcessing complete. Output files generated.")

if __name__=="__main__":
    main()
