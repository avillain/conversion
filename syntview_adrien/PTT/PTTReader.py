import fileinput
import os
import re
import string

class PTTHeaderError(Exception):
    pass

class PTTLineError(Exception):
    pass

class PTTCompletenessError(Exception):
    pass

class PTTReader(object):

    def __init__(self):
        self._header    = [""]*3
        self._header[0] = re.compile("^.* - 1\.\.([1-9][0-9]*)\s*$", re.I)
        self._header[1] = re.compile("^([0-9][0-9]*) proteins\s*$", re.I)
        self._header[2] = re.compile("^Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct$", re.I)
 

    def _check_header(self, in_fh):
        line  = in_fh.readline()
        line  = line.replace('"','')
        match = self._header[0].search(line)
        if not match:
            raise PTTHeaderError, "Header line number 1 not correctly formatted"
        else:
            genome_size = string.atoi( match.group(1) )
        
        line  = in_fh.readline()
        line  = line.replace('"','')
        match = self._header[1].search(line)    
        if not match:
            raise PTTHeaderError, "Header line number 2 not correctly formatted"
        else:
            nb_prot = string.atoi( match.group(1) )
            
        line  = in_fh.readline()
        line  = line.replace('"','')
        match = self._header[2].search(line)    
        if not match:
            raise PTTHeaderError, "Header line number 3 not correctly formatted"
        return (genome_size, nb_prot)


    def read_and_check(self, ptt_fh):
        error_msg = ""
        # let this part throw an exception if necessary
        genome_size, nb_prot = self._check_header(ptt_fh)
            
        ptt_dict = {}
        nb       = 3
        for line in ptt_fh:
            line = line[:-1]
            if not line:
                continue
            try:
                location, strand, length, PID, \
                          gene, synonym, code, COG, product = line.split("\t")
                nb += 1
            except ValueError:
                error_msg += "Line %d not well formatted : 9 columns with tab separator required" %(nb)
                continue
            
	    if gene !="-":
	        nmid=gene
	    else:
		nmid=synonym
            ptt_dict[nmid] = {}
            start, end        = location.split("..")
            ptt_dict[nmid]['start'], ptt_dict[nmid]['end'] = string.atoi(start), string.atoi(end)
            ptt_dict[nmid]['strand'] = strand
            
            ptt_dict[nmid]['PID']     = ""
            ptt_dict[nmid]['gene']    = ""
            ptt_dict[nmid]['COG']     = ""
            ptt_dict[nmid]['product'] = ""
            
            if PID     != "-":
                ptt_dict[nmid]['PID']     = PID
            if gene    != "-" :
                ptt_dict[nmid]['gene']    = gene
            if COG     != "-":
                ptt_dict[nmid]['COG']     = COG
            if product != "-":
                ptt_dict[nmid]['product'] = product

            if error_msg:
                raise PTTLineError, error_msg
            
        if nb-3 != nb_prot:
            raise PTTCompletenessError, "The number of lines %d is not equal to the expected number of proteins %d" %(nb-3, nb_prot)
            
        return ptt_dict
    
            
