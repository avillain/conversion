#!/usr/bin/env python2.7

import os
import pwd
import re
from sys import argv, exit, stdout, stderr, path

import argparse
import SynTView.Formatter as Formatter
import PTT.PTTReader as PTTReader
import SnpEffUtils.Reader as SnpEffReader

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts flanking regions of SNPS/indels.')
    parser.add_argument("ptt_file", help='CDS file in .ptt format (see ftp://ftp.ncbi.nih.gov/genomes/)')
    parser.add_argument("vcf_file", help='SNP/indel file in VCF format, annotated with snpEff (see http://snpeff.sourceforge.net/SnpEff_manual.html#gatk)')
    
    #args = parser.parse_args()    
    
    #parser.add_option("-p", dest="ptt_file", 
    #   help=".ptt file", default=None)
    #parser.add_option("-s", dest="snp_file", 
    #    help="snpEff output file", default=None)
    # grab options and arguments in an usable form
   

    #parser.add_option("-l", dest="log_file", 
    #    help="Optional log file", default=None)
    #parser.add_option("-o", dest="out_prefix", 
    #    help="Prefix for output files", default=None)
    
    parser.add_argument("-l", dest="log_file", 
        help="Optional log file", default=None)
    parser.add_argument("-o", dest="out_prefix", 
        help="Prefix for output files", default=None)
    
    if not argv[1:] :
       argv.append('-h')

    prog_name = os.path.basename(argv[0])
    args = parser.parse_args()
    #print options 
    
               
    if not args.vcf_file :
        parser.error("\nvcf file missing.")
    else :
        if not os.access(args.vcf_file, os.R_OK):
            parser.error("Unable to read vcf file %s " %args.vcf_file)
                
    if not args.out_prefix :
        out_snp   = args.vcf_file[:-4]+".snp"
        out_indel = args.vcf_file[:-4]+".indel"
        if os.access(out_snp, os.F_OK) and not os.access(out_snp, os.W_OK):
            parser.error("\nSNP output file %s exists but is not writable." %out_snp)
        if os.access(out_indel, os.F_OK) and not os.access(out_indel, os.W_OK):
            parser.error("\nINDEL output file %s exists but is not writable." %out_indel)
    else :
        out_snp   = args.out_prefix+".snp"
        out_indel = args.out_prefix+".indel"
        if os.access(out_snp, os.F_OK) and not os.access(out_snp, os.W_OK):
            parser.error("\nSNP output file %s exists but is not writable." %out_snp)
        if os.access(out_indel, os.F_OK) and not os.access(out_indel, os.W_OK):
            parser.error("\nINDEL output file %s exists but is not writable." %out_indel)

    log_fh = stderr
    if args.log_file:
        if os.access(args.log_file, os.F_OK) and not os.access(args.log_file, os.W_OK):
            parser.error("\n log file %s exists but is not writable." %args.log_file)
        else:
            log_fh = open(args.log_file,"w")

    if not args.ptt_file :
        parser.error("\nPTT file missing.\n")
    else :
        if not os.access(args.ptt_file, os.R_OK):
            parser.error("Unable to read PTT file %s " %args.ptt_file)
    
    ptt_reader      = PTTReader.PTTReader()
    snp_reader      = SnpEffReader.SnpEffReader()
    snp_reformatter = Formatter.SynTViewSnpFormatter()
    with open(args.ptt_file, "r") as ptt_fh:
        ptt_dict    = ptt_reader.read_and_check(ptt_fh)
	#from pprint import pprint
	#pprint(ptt_dict)
    
    with open(args.vcf_file, "r") as snp_fh:
        snp_list = snp_reader.read_and_check(snp_fh)
        
    with open(out_snp, "w") as out_snp_fh:
        with open(out_indel, "w") as out_indel_fh:
            for snp in snp_list:
                try:
                    type, reformatted_line = snp_reformatter.format(snp, ptt_dict)
                    if type == "SNP":
                        print >>out_snp_fh, reformatted_line
                    if type == "INS" or type == "DEL":
                        print >>out_indel_fh, reformatted_line
                except Formatter.SynTViewSnpFormatterGeneError, msg:
                    print >>log_fh, str(msg)
                except Formatter.SynTViewSnpFormatterIntragenicError, msg:
                    print >>log_fh, str(msg)+" ignored"
