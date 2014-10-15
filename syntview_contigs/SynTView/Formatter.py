import re
import sys


class SynTViewSnpFormatterGeneError(Exception):
    pass

class SynTViewSnpFormatterIntragenicError(Exception):
    pass

class SynTViewSnpFormatter(object):
    
    def __init__(self):
        self._intragenic = re.compile("intragenic", re.I)
        self._intergenic = re.compile("DOWNSTREAM|UPSTREAM", re.I)
        self._pseudo     = re.compile("pseudo", re.I)
        self._rna        = re.compile("[^m]rna", re.I)
        self._aa_dict =  {'A' : 'Ala',
                         'R' : 'Arg',
                         'N' : 'Asn',
                         'D' : 'Asp',
                         'C' : 'Cys',
                         'Q' : 'Gln',
                         'E' : 'Glu',
                         'G' : 'Gly',
                         'H' : 'His',
                         'I' : 'Ile',
                         'L' : 'Leu',
                         'K' : 'Lys',
                         'M' : 'Met',
                         'F' : 'Phe',
                         'P' : 'Pro',
                         'S' : 'Ser',
                         'T' : 'Thr',
                         'W' : 'Trp',
                         'Y' : 'Tyr',
                         'V' : 'Val',
                         '*' : '***',
                         '-' : '-',
                         '' : ''}

    def format(self, snp_dict, ptt_dict):
        id        = ""
        strand    = None
        gene_name = snp_dict["gene_name"]
        if self._intragenic.search(snp_dict["effect"]):
            # we ignore intragenic SNP
            # because they are redundant with others
            # we don't want to count the same SNP twice
            raise SynTViewSnpFormatterIntragenicError, \
                  "Intragenic element for gene %s" %gene_name
        if self._intergenic.search(snp_dict["effect"]) or \
               self._pseudo.search(snp_dict["biotype"]) or \
               self._rna.search(snp_dict["biotype"]):
            id     = ""
            strand = ""
        elif gene_name in ptt_dict:
	    id     = gene_name
            strand = ptt_dict[gene_name]["strand"]
        else:
            raise SynTViewSnpFormatterGeneError, "Gene name %s unknown in ptt file" %gene_name
                
        if snp_dict["type"] == 'SNP':
            old_aa = self._aa_dict[ snp_dict["old_aa"] ]
            new_aa = self._aa_dict[ snp_dict["new_aa"] ]
	    id=id if id else "-"
	    #from pprint import pprint
	    #pprint(snp_dict)
	    strand=strand if strand else "-"
	    old_aa=old_aa if old_aa else "-"
	    new_aa=new_aa if new_aa else "-"
	    #print strand,"\t",old_aa,"\t",new_aa,"\t",id,"\t",snp_dict["cov"]
            out_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(snp_dict["position"],
                                                      snp_dict["reference"],
                                                      snp_dict["variant"],
                                                      id, strand,
                                                      old_aa, new_aa, snp_dict["cov"], snp_dict["covtot"])
        elif snp_dict["type"] == 'INS':
            variant  = snp_dict["variant"]
            length   = len(variant)
            out_line = "%s\t%s\t%d\t\t%s" %(snp_dict["position"], variant, length, id)
        elif snp_dict["type"] == "DEL" :
            variant  = snp_dict["variant"]
            length   = len(variant)
            out_line = "%s\t%s\t-%d\t\t%s" %(snp_dict["position"], variant, length, id)
       
        return snp_dict["type"], out_line
