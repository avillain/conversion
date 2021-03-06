#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Villain Adrien
# 2011-2014     
# avillain@pasteur.fr


import re # module permettant d'utiliser les expressions régulières
import sys

### FONCTIONS ###
# Complémentaire inverse
def reverse(seq): # séquence à inverser
    rev = ""
    for base in seq:
        if base == "A":
            rev += "T"
        elif base == "T":
            rev += "A"
        elif base == "C":
            rev += "G"
        elif base == "G":
            rev += "C"
    rev = rev[::-1]
    return rev # séquence inversée

# Extraction de l'organisme, des genes et de la sequence à partir d'un fichier genbank
def genbank_extract(fil): # fichier .gbk
    ## Expressions régulières
    renam= re.compile('^ *ORGANISM  (([A-Za-z]+ ?)+)') # ligne contenant le nom de l'organisme
    #regdef= re.compile('^DEFINITION(.*)ACCESSION')
    regen= re.compile('CDS +([0-9]+)\.\.([0-9]+)') # lignes contenant l'emplacement de début et de fin d'un gène
    recomp= re.compile('CDS +complement\(([0-9]+)\.\.([0-9]+)\)') # idem pour gène complémentaire
    reseq= re.compile('^ *[0-9]+(( [ATCGatcg]{0,10})+)') # lignes contenant la séquence nucléotidique
    relen= re.compile('^LOCUS.* +([0-9]+) bp.*') # longueur de la séquence
    reloc= re.compile('locus_tag="(.*)"')
    reso=re.compile('^ *source *([0-9]+\.\.[0-9]+)')
    repid=re.compile('^ */db_xref.*GI:(\d+)".*')
    reprod=re.compile('^ */product="(.*)".*')
    #regene
    
    ## Variables de stockage (sequence, bornes des genes)
    seq=''
    genes=[]

    ## Ouverture du fichier genbank, mode lecture
    filin=open(fil,'r')
    flag_seq=0 # la ligne 'ORIGIN' n'a pas été atteinte
    locus_tag=''
    tag=0
    tagsource=0
    defin=''
    prots=''
    sign='+'
    pid='-'
    gename='-'
    synon='-'
    code='-'
    cog='-'
    prod='-'
    countCDS=0
    ## Extraction des informations du fichier genbank
    for line in filin: # parcours du contenu du fichier
        if renam.findall(line): # repérage du nom de l'organisme
            org=renam.findall(line)
            organism=org[0][0] # 1er element du groupe 0 soit le nom

        elif re.match('ORIGIN',line):
            flag_seq=1 # on se trouve après ORIGIN
	
	elif reprod.findall(line):
	    prod=reprod.findall(line)[0]
	    
	elif repid.findall(line):
	    pid=repid.findall(line)[0]
	    synon=locus_tag
            prots+="\n%d..%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s" %(start,end,sign,length,pid,gename,synon,code,cog,prod)
	# recherche bornes du gène à chaque ligne
        elif regen.findall(line):
            gen=regen.findall(line)
	    start=int(gen[0][0])
	    end=int(gen[0][1])	   
	    length=(end-start)/3
	    gen.append('') # gène non complémentaire
	    gen.append(locus_tag)
            genes.append(gen)
	    countCDS+=1
	    sign='+'
	    #prots+="\n%d..%d\t+\t%d\t%s\t%s\t%s\t%s\t%s\t%s" %(start,end,length,pid,gename,synon,code,cog,prod)
        # idem gene complémentaire
        elif recomp.findall(line):
            comp=recomp.findall(line)
	    start=int(comp[0][0])
	    end=int(comp[0][1])
	    length=(end-start)/3
            comp.append('complement') # gène complémentaire
	    comp.append(locus_tag)
            genes.append(comp)
	    countCDS+=1
	    sign='-'
	    #prots+="\n%d..%d\t-\t%d" %(start,end,length)
        # obtention de la séquence nucléotidique complète
        elif reseq.findall(line) and flag_seq==1:
            sequence=reseq.findall(line)
            seq+=sequence[0][0]

        # obtention de la longueur de la séquence
        elif relen.findall(line):
            length=relen.findall(line)
            leng=length[0]

	elif reloc.findall(line):
	    locus_tag=reloc.findall(line)[0]

	elif re.search('DEFINITION',line):
	    tag=1
	
	elif re.search('ACCESSION',line):
	    tag=0
	
	if tag:
	    defin+=line
        
	elif reso.findall(line):
	    sou=reso.findall(line)
	    defin+=" - %s" %(sou[0])
    defin=' '.join(defin.replace('DEFINITION','').split())
    defin+="\n%d proteins\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct" %(countCDS)
    print "%s%s" %(defin, prots)
    ## Fermeture du fichier source
    filin.close()

    ## Suppression des espaces etc
    seq="".join(seq.split())

#    ## Verification de la taille de la sequence
#    if len(seq) == int(leng):
#        print "Bonne longueur de sequence extraite"

    return [organism,genes,seq,defin] # liste contenant l'organisme, les gènes et la séquence

## Ecriture dans les fichiers .fasta
#def ptt_out(results): # liste contenant l'organisme, les gènes et la séquence
    
   #compt=0 # compteur du numéro de gène
    #for gene in results[1]: # pour tous les gènes
       #compt+=1
        #if gene[1] == 'complement':
            #genseq=reverse((results[2][int(gene[0][0])-1:int(gene[0][1])]).upper()) # complémentaire inverse en majuscules de la séquence du gène
        #else: # non complémentaire
            #genseq=results[2][int(gene[0][0])-1:int(gene[0][1])].upper() # bornes [début-1:fin] car seq est une chaîne de caractères (indices commencent à 0)

        ## Formatage .fasta
        #out='>%s\n'%(gene[2]) # commentaire avec nom organisme, numéro de gène, début et fin de séquence du gène
        #for i in range(0,len(genseq),80): # écriture de 80 bases par ligne
            #out += "".join(genseq[i:i+80]).lower()

        ## Ecriture dans un fichier
        #fout=open('gene%i.fasta' %compt,'w') # nom de la forme gene1.fasta, gene2.fasta ...
        #fout.write("%s" %out)
        #fout.close() # fermeture du fichier
	#print out
    #print results[3]

### PROGRAMME PRINCIPAL ###
## Extraction de l'organisme, des gènes (bornes de début et de fin de chaque gène) et de la séquence nucléotidique de NC_001133.gbk
results = genbank_extract(sys.argv[1])

## Ecriture dans des fichiers .fasta
#@:x
#ptt_out(results)
