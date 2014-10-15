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
	name=line.split()[0]
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
        return (genome_size, nb_prot, name)


    def read_and_check(self, ptt_fh):
        error_msg = ""
        # let this part throw an exception if necessary
        genome_size, nb_prot, name = self._check_header(ptt_fh)
        ptt_dict = {}
        nb       = 3
	for i in range(nb_prot):
	    line = ptt_fh.readline()[:-1]
	    if not line:
                continue
            try:
                location, strand, length, PID, \
                          gene, synonym, code, COG, product = line.split("\t")
                nb += 1
			
            except ValueError:
		error_msg += "Line %d not well formatted : 9 columns with tab separator required" %(nb)
                continue
            
            ptt_dict[synonym] = {}
            start, end        = location.split("..")
            ptt_dict[synonym]['start'], ptt_dict[synonym]['end'] = string.atoi(start), string.atoi(end)
            ptt_dict[synonym]['strand'] = strand
            
            ptt_dict[synonym]['PID']     = ""
            ptt_dict[synonym]['gene']    = ""
            ptt_dict[synonym]['COG']     = ""
            ptt_dict[synonym]['product'] = ""
            
            if PID     != "-":
                ptt_dict[synonym]['PID']     = PID
            if gene    != "-" :
                ptt_dict[synonym]['gene']    = gene
            if COG     != "-":
                ptt_dict[synonym]['COG']     = COG
            if product != "-":
                ptt_dict[synonym]['product'] = product
	    #print nb-3, nb_prot
    	    #if nb-3==nb_prot:
		#break
            if error_msg:
                raise PTTLineError, error_msg
	
        if nb-3 != nb_prot:
            raise PTTCompletenessError, "The number of lines %d is not equal to the expected number of proteins %d" %(nb-3, nb_prot)
        #print line   
        return [name,ptt_dict]
    
    def read_ptt(self, ptt_fh):
	ptt_list=[]
	names_list=[]
	while(1):
	    try:
		name,dic=self.read_and_check(ptt_fh)
		names_list.append(name)
		ptt_list.append(dic)
	    except:
		#raise
		break
	return [ptt_list,names_list]
