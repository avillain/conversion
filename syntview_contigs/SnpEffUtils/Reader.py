import re
import string
import sys

class SnpEffLineError(Exception):
    pass


class SnpEffReader(object):
   
    def __init__(self):
        pass
    
    def read_and_check(self, snp_fh):
        snp_list  = []
        error_msg = ""
        nb        = 0 
        for line in snp_fh:
            line = line.replace('"','')	
            nb  += 1
            if re.search("^\s*#", line):
                # ignore header
                continue
            try:
               # chr, snp_pos, ref, variant, type, homozyguous, \
                #     qual, cov, warn, gene_id, gene_name, \
                 #    biotype, transcript_id, exon_id, exon_rank, \
                  #   effect, aa_change, codon_change, _ = line.split("\t", 18)
		chr,pos,_,ref,variant,qual,status,info,field,sample=line.split("\t",10)
		#type=="SNP" if len(ref)==len(variant)==1 else "INDEL"
		ty="SNP" if (len(ref)==1 and len(variant)==1) else "INS" if (len(ref)==1 and len(variant)>1) else "DEL"
		try:
			gene_name=re.search("SNPEFF_TRANSCRIPT_ID=(.*);",info).groups()[0]
		except:
			try:
				gene_name=re.search("SNPEFF_TRANSCRIPT_ID=(.*)$",info).groups()[0]
			except:
				gene_name=""
		if gene_name!="":
			try:
				effect=re.search("SNPEFF_EFFECT=([^;]*);",info).groups()[0]
			except:
				effect=""
		else:
			effect=""
		biotype=""
            except ValueError, msg:
                error_msg += "Line %d not well formatted : expecting vcf variant file annotated using snpEff with gatk compatibility mode" %(nb, line)
                continue
            
	    old_aa, new_aa = "", ""

	    if ty == 'SNP':
           	aa_change=re.search("AMINO_ACID_CHANGE=([^;]*);",info)
		oth,all=sample.split(':')[1].split(',')
		#print oth,all
		totalcov=str(int(oth)+int(all))
		#print totalcov
		if effect=="NON_SYNONYMOUS_CODING" or effect=="STOP_GAINED":
                	try:
                    	# we are in a gene
				change=aa_change.groups()[0]
		      	#old_aa, new_aa = aa_change.split("/")
				old_aa,new_aa = re.search("([A-Z\*])\d+([A-Z\*])",change).groups()
            		except:
				pass
		elif effect=="SYNONYMOUS_CODING":
			try:
                        # we are in a gene
                                change=aa_change.groups()[0]
	
                        #old_aa, new_aa = aa_change.split("/")
                                old_aa = change[:1]
				new_aa = old_aa
                        except:
                                pass
	    elif ty=='DEL':
		toto=ref
		ref=variant
		variant=toto
            new_dict = {"chromosome": chr,
			"position"  : string.atoi(pos),
                        "type"      : ty,
                        "reference" : ref,
                        "variant"   : variant,
                        "gene_name" : gene_name,
                        "effect"    : effect,
                        "biotype"   : biotype,
			"cov"	    : all,
			"covtot"    : totalcov}
                
            new_dict["old_aa"]    = old_aa
            new_dict["new_aa"]    = new_aa
            snp_list.append(new_dict)

        if error_msg:
            raise SnpEffLineError, error_msg
        return snp_list
