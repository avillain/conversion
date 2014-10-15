#!/mount/biotools/bin/python2.6

import os
from sys import argv, exit, stdout, stderr

from optparse import OptionParser

from Bio import SeqIO

def write_ptt_and_faa_from_feature(seq_feature, faa_fh, ptt_fh, log_fh,
                                   prefix=None, organism=None, step=70):
    location = ""
    strand   = "+"
    length   = 0
    pid      = "-"         
    gene     = "-"
    synonym  = None  
    code     = "-"
    cog      = "-" 
    product  = "-"
    comment  = ""
    start, stop = seq_feature.location.nofuzzy_start, seq_feature.location.nofuzzy_end 
    location = str(start+1)+".."+str(stop)
    if seq_feature.strand < 0:
        strand = "-"         
    diff   = stop - start
    # for length, suppress stop codon for the length of the protein
    length =  (diff + 1)/3 -1
    if 'db_xref' in seq_feature.qualifiers and len(seq_feature.qualifiers['db_xref']) > 0:
        pid     = seq_feature.qualifiers['db_xref'][0]
        pid     = pid.split(":")[1]  
    if 'gene' in seq_feature.qualifiers and len(seq_feature.qualifiers['gene']) > 0 :
        gene    = seq_feature.qualifiers['gene'][0]
    if prefix :
        synonym = prefix  
    elif 'locus_tag' in seq_feature.qualifiers and len(seq_feature.qualifiers['locus_tag']) == 1 :
        synonym = seq_feature.qualifiers['locus_tag'][0]
    if not synonym :
        print >>log_fh, "Error: CDS without locus_tag and no prefix given %s %s" %(gene, location)
        return -1
    if (not 'translation' in seq_feature.qualifiers) or len(seq_feature.qualifiers['translation']) != 1 :
        print >>log_fh, "Error: CDS without translation %s %s %s" %(synonym, gene, location) 
        return -2
    if 'product' in seq_feature.qualifiers and len(seq_feature.qualifiers['product']) > 0 :
        product = seq_feature.qualifiers['product'][0]
    if organism:
        if product != "hypothetical protein":
            comment = " " + product
        comment += " [" + organism + "]"
    ptt_fh.write( "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s%s" % (
        location, strand, length, pid, gene, synonym,
        code, cog, product, os.linesep) )
    faa_fh.write(">%s%s%s" % ( synonym, comment,os.linesep ))
    translation = seq_feature.qualifiers['translation'][0]
    for i in range(0, len(translation), step):
        faa_fh.write("%s%s"  % ( translation[i:i+step], os.linesep))
    return 0


def write_ffn_from_feature(seq_feature, whole_seq, ffn_fh, log_fh, id=None, step=60 ):
    if not id and  'locus_tag' in seq_feature.qualifiers and len(seq_feature.qualifiers['locus_tag']) == 1 :
        id = seq_feature.qualifiers['locus_tag'][0]
    if not id :
        start, stop = seq_feature.location.nofuzzy_start, seq_feature.location.nofuzzy_end 
        location = str(start+1)+".."+str(stop)
        gene = "-"
        if 'gene' in seq_feature.qualifiers and len(seq_feature.qualifiers['gene']) > 0 :
            gene    = seq_feature.qualifiers['gene'][0]
        print >>log_fh, "Error: CDS without locus_tag and no prefix given %s %s" %(gene, location)
        return -1
    cds_seq = seq_feature.extract(whole_seq)
    ffn_fh.write(">%s%s" % ( id, os.linesep ))
    for i in range(0, len(cds_seq), step):
        ffn_fh.write("%s%s"  % ( cds_seq[i:i+step].lower(), os.linesep))
    return 0


def write_sequin_from_feature(seq_feature, sequin_fh):
    nb_qual_written = 0
    start, stop = seq_feature.location.nofuzzy_start, \
                  seq_feature.location.nofuzzy_end
    start += 1    
    if seq_feature.strand < 0:
        start, stop = stop, start
    sequin_fh.write( "%s\t%s\t%s%s" %(start, stop, seq_feature.type, os.linesep) )
    for qualifier, val_list in seq_feature.qualifiers.iteritems():
        for val in val_list:
            if qualifier != "codon_start" and qualifier != "transl_table":
                val = '"'+val+'"'
            sequin_fh.write("\t\t\t%s\t%s%s" %(qualifier, val, os.linesep))
            nb_qual_written += 1
    return nb_qual_written


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-g", dest="gb_file", 
        help="Genbank input file", default=None)
    parser.add_option("-e", dest="embl_file", 
        help="EMBL input file", default=None)
    parser.add_option("-p", dest="prefix", 
        help="prefix for file names", default=None)
    parser.add_option("-a", dest="ptt_flag",
        help="write the .ptt and .faa files",
        action="store_true",default=False) 
    parser.add_option("-n", dest="nucseq_flag", 
        help="write the .nucseq file",
        action="store_true",default=False)
    parser.add_option("-s", dest="sequin_flag", 
        help="write the sequin file ( .tbl )",
        action="store_true",default=False)
    parser.add_option("-l", dest="log_file", 
        help="Optional log file", default=None)
    parser.add_option("-d", dest="output_dir", 
        help="Output directory", default="./")
    parser.add_option("-c", dest="comment_flag", 
        help="If set, add the product and the organism name "+ \
              os.linesep+"after the ID on the header line of sequences",
        action="store_true", default=False)
    parser.add_option("-i", dest="change_nb", 
        help="If set, replace the locus_tag by the prefix + number",
        action="store_true", default=False)
    
    # grab options and arguments in an usable form
    if not argv[1:] :
       argv.append('-h')

    prog_name = os.path.basename(argv[0])
    (options, args) = parser.parse_args()

    in_format = None
    in_file   = None
    fh_list   = []
    
    if not options.gb_file and not options.embl_file :
        parser.error("\nGenbank or EMBL file missing.\n")
    if options.gb_file and options.embl_file :
        parser.error("\nYou can't use -e and -g at the same time.\n")

    if options.gb_file :
        if not os.access( options.gb_file, os.R_OK):
            parser.error("Unable to read Genbank file %s " %( options.gb_file))
        in_format = "genbank"
        in_file   = options.gb_file
    elif options.embl_file :
        if not os.access( options.embl_file, os.R_OK):
            parser.error("Unable to read EMBL file %s " %( options.embl_file))
        in_format = "embl"
        in_file   = options.embl_file
                
    if not options.prefix :
        parser.error("\nprefix for output file is missing.\n")
    prefix = options.prefix

    if options.ptt_flag:
        faa_file = os.path.join(options.output_dir, prefix + ".faa")
        ptt_file = os.path.join(options.output_dir, prefix + ".ptt")
        if os.access(faa_file, os.F_OK) and not os.access(faa_file, os.W_OK):
            parser.error("\nfaa output file %s exists but is not writable.\n" %faa_file)
            if os.access(ptt_file, os.F_OK) and not os.access(ptt_file, os.W_OK):
                parser.error("\nptt output file %s exists but is not writable.\n" %ptt_file)
        faa_fh = open(faa_file,"w")
        ptt_fh = open(ptt_file,"w")
        fh_list.extend( (faa_fh, ptt_fh) )
        
    if options.nucseq_flag:
        nucseq_file = os.path.join(options.output_dir, prefix + ".nucseq")
        if os.access(nucseq_file, os.F_OK) and not os.access(nucseq_file, os.W_OK):
            parser.error("\nnucseq output file %s exists but is not writable.\n" %nucseq_file)
        nucseq_fh = open(nucseq_file,"w")
        fh_list.append(nucseq_fh)

    if options.sequin_flag:
        tbl_file =  os.path.join(options.output_dir, prefix + ".tbl")
        if os.access(tbl_file, os.F_OK) and not os.access(tbl_file, os.W_OK):
            parser.error("\ntbl sequin output file %s exists but is not writable.\n" %tbl_file)
        tbl_fh = open(tbl_file,"w")
        fh_list.append(tbl_fh)

    log_fh = stderr
    if options.log_file:
        if os.access(options.log_file, os.F_OK) and not os.access(options.log_file, os.W_OK):
            parser.error("\n log file %s exists but is not writable.\n" %options.log_file)
        else:
            log_fh = open(options.log_file,"w")
            fh_list.append(log_fh)
        
        
    with open(in_file,"r") as in_fh:
        nb        = 0
        organism  = None
        accession = None 
        for seq_record in SeqIO.parse(in_fh, in_format) :
            accession = seq_record.annotations['accessions'][0]
            if options.comment_flag :
                organism = seq_record.annotations['organism']
            if options.sequin_flag:
                tbl_fh.write(">Feature %s%s" %( accession, os.linesep ))
            if options.ptt_flag:
                nb_prot = 0
                for seq_feature in seq_record.features :
                    if seq_feature.type == "CDS" :
                        nb_prot += 1
                ptt_fh.write( "%s %s - 1..%d%s" %(seq_record.name, seq_record.description,
                                                  len(seq_record.seq), os.linesep) )
                ptt_fh.write( "%d proteins%s" %(nb_prot, os.linesep) )                
                ptt_fh.write("Location\tStrand\tLength\tPID\tGene")
                ptt_fh.write("\tSynonym\tCode\tCOG\tProduct%s" %os.linesep)
            for seq_feature in seq_record.features :
                # in tbl file, everything is written
                if options.sequin_flag:
                    write_sequin_from_feature(seq_feature, tbl_fh)
                    
                if seq_feature.type == "CDS" :
                    id = None
                    if options.change_nb:
                        nb += 1
                        id = options.prefix + "_%0.4d" %nb
                    if options.ptt_flag:
                        write_ptt_and_faa_from_feature(seq_feature, faa_fh,
                                                       ptt_fh, log_fh, 
                                                       prefix   = id,
                                                       organism = organism)
                    if options.nucseq_flag:
                        write_ffn_from_feature(seq_feature, seq_record.seq,
                                               nucseq_fh, log_fh, id)

    for fh in fh_list:
        fh.close()
