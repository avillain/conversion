#! /usr/bin/env python

import unittest    
from test     import test_support
from StringIO import *

import os
from string import *
import re

import SnpEffUtils.Reader as Reader

class TestSNPReader(unittest.TestCase):

    reader = None
   
    def setUp(self):
        TestSNPReader.reader = Reader.SnpEffReader()

    def test_header(self):
        snp_fh = StringIO()
        print >>snp_fh,  "# SnpEff version 2.0.3 (build 2011-10-08), by Pablo Cingolani"
        expected_snp_list = []
        
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)


    def test_header_quote(self):
        snp_fh = StringIO()
        print >>snp_fh,  '"# SnpEff version 2.0.3 (build 2011-10-08)"	" by Pablo Cingolani"	'
        expected_snp_list = []
        
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)

        
    def test_raise_error(self):
        snp_fh = StringIO()
        print >>snp_fh,  "NC_003210\t1745\tC\tA\tSNP\tHom\t222\t765\t\t\t\t\t\t\t\tINTERGENIC"

        snp_fh.seek(0, 0)
        self.assertRaises(Reader.SnpEffLineError, TestSNPReader.reader.read_and_check, snp_fh)


    def test_intra(self): 
        snp_fh = StringIO()
        print >>snp_fh,  "NC_003210\t2519935\tT\tC\tSNP\tHom\t222\t741\t\tlmo2448\tlmo2448\t\t\t\t\tINTRAGENIC: lmo2448\t\t\t\t"
        expected_snp_list = [{"type"      :"SNP",
                              "position"  :2519935,
                              "reference" :"T",
                              "variant"   :"C",
                              "gene_name" :"lmo2448",
                              "effect"    :"INTRAGENIC: lmo2448",
                              "biotype"   :"",
                              "old_aa"    :"",
                              "new_aa"    :""}]
        
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)

        
    def test_pseudo(self):
        snp_fh = StringIO()
        print >>snp_fh,  "NC_003210\t432305\tC\tT\tSNP\tHom\t222\t655\t\tGene_lmo0410\tGene_lmo0410\tpseudogene\tlmo0410\t\t\tINTRON\t\t\t\t\t0\t\t\t"
        expected_snp_list = [{"type"      :"SNP",
                              "position"  :432305,
                              "reference" :"C",
                              "variant"   :"T",
                              "gene_name" :"Gene_lmo0410",
                              "effect"    :"INTRON",
                              "biotype"   :"pseudogene",
                              "old_aa"    :"",
                              "new_aa"    :""}]
        
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)

        
    def test_nonsynonymous(self):
        snp_fh = StringIO()
        print >>snp_fh, "chr1\t523603\t.\tG\tA\t2390\tPASS\tAC=1;AF=1.00;AN=1;DP=62;Dels=0.00;FS=0.000;HaplotypeScore=0.0000;MLEAC=1;MLEAF=1.00;MQ=59.64;MQ0=0;QD=32.41;SNPEFF_AMINO_ACID_CHANGE=A120T;SNPEFF_CODON_CHANGE=Gcg/Acg;SNPEFF_EFFECT=NON_SYNONYMOUS_CODING;SNPEFF_EXON_ID=Exon_1_523245_525866;SNPEFF_FUNCTIONAL_CLASS=MISSENSE;SNPEFF_GENE_NAME=Exon_chromosome1_523246_525867;SNPEFF_IMPACT=MODERATE;SNPEFF_TRANSCRIPT_ID=Transcript_Exon_chromosome1_523246_525867\tGT:AD:DP:GQ:MLPSAC:MLPSAF:PL\t1:0,62:62:99:1:1.00:2420,0"
        expected_snp_list = [{"type"      :"SNP",
                              "position"  :523603,
                              "reference" :"G",
                              "variant"   :"A",
                              "gene_name" :"Exon_chromosome1_523246_525867",
                              "effect"    :"NON_SYNONYMOUS_CODING",
                              "biotype"   :"",
                              "old_aa"    :"A",
                              "new_aa"    :"T"}]
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)
        
        
    def test_snp_in_gene(self):
        snp_fh = StringIO()
        print >>snp_fh, "NC_003210\t17131\tG\tA\tSNP\tHet\t19.1\t98\t\tGene_lmo0013\tlmo0013\t\tTranscript_lmo0013\tExon_NC_003210_16218_17324\t1\tNON_SYNONYMOUS_CODING\tE/K\tGaa/Aaa\t305\t0\t1107\t\t\t"
        expected_snp_list = [{"type"      :"SNP",
                              "position"  :17131,
                              "reference" :"G",
                              "variant"   :"A",
                              "gene_name" :"lmo0013",
                              "effect"    :"NON_SYNONYMOUS_CODING",
                              "biotype"   :"",
                              "old_aa"    :"E",
                              "new_aa"    :"K"}]
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)
                

    def test_snpstop_in_gene(self):
        snp_fh = StringIO()
        print >>snp_fh, "NC_003210\t2830008\tG\tA\tSNP\tHom\t222\t1116\t\tGene_lmo2752\tlmo2752\t\tTranscript_lmo2752\tExon_NC_003210_2828235_2830007\t1\tSYNONYMOUS_STOP\t*/*\ttaG/taA\t591\t1\t1773\t\t\t"
        expected_snp_list = [{"type"      :"SNP",
                              "position"  :2830008,
                              "reference" :"G",
                              "variant"   :"A",
                              "gene_name" :"lmo2752",
                              "effect"    :"SYNONYMOUS_STOP",
                              "biotype"   :"",
                              "old_aa"    :"*",
                              "new_aa"    :"*"}]
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)
        
    def test_ins_in_gene(self):
        snp_fh = StringIO()
        print >>snp_fh, "chr1\t101871\t.\tA\tAG\t2379.97\tPASS\tAC=1;AF=1.00;AN=1;DP=58;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.00;MQ0=0;QD=29.04;RPA=1,2;RU=G;SNPEFF_AMINO_ACID_CHANGE=-558?;SNPEFF_CODON_CHANGE=-/C;SNPEFF_EFFECT=FRAME_SHIFT;SNPEFF_EXON_ID=Exon_1_101871_103543;SNPEFF_FUNCTIONAL_CLASS=NONE;SNPEFF_GENE_NAME=Exon_chromosome1_101872_103544;SNPEFF_IMPACT=HIGH;SNPEFF_TRANSCRIPT_ID=Transcript_Exon_chromosome1_101872_103544;STR\tGT:AD:DP:GQ:MLPSAC:MLPSAF:PL\t1:0,57:58:99:1:1.00:2419,0"
        expected_snp_list = [{"type"      :"INS",
                              "position"  :101871,
                              "reference" :"A",
			      "variant"   :"AG",
			      "gene_name" :"Transcript_Exon_chromosome1_101872_103544",
                              "effect"    :"FRAME_SHIFT",
                              "biotype"   :"",
                              "old_aa"    :"",
                              "new_aa"    :""}]
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)
        
    def test_del_in_inter(self):
        snp_fh = StringIO()
        print >>snp_fh, "NC_007793\t2627119\t.\tTTA\tT\t1406.97\tPASS\tAC=1;AF=1.00;AN=1;DP=40;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=59.52;MQ0=0;QD=17.59;RPA=7,6;RU=TA;SNPEFF_EFFECT=DOWNSTREAM;SNPEFF_FUNCTIONAL_CLASS=NONE;SNPEFF_GENE_NAME=fnbB;SNPEFF_IMPACT=MODIFIER;SNPEFF_TRANSCRIPT_ID=SAUSA300_2440;STR\tGT:AD:DP:GQ:MLPSAC:MLPSAF:PL\t1:0,33:39:99:1:1.00:1446,0"
        expected_snp_list = [{"type"      :"DEL",
                              "position"  :2627119,
                              "reference" :"TTA",
                              "variant"   :"T",
                              "gene_name" :"",
                              "effect"    :"",
                              "biotype"   :"",
                              "old_aa"    :"",
                              "new_aa"    :""}]
        snp_fh.seek(0, 0)
        snp_list = TestSNPReader.reader.read_and_check(snp_fh)
        self.assertEqual(expected_snp_list, snp_list)
        

if __name__ == "__main__":
    test_support.run_unittest(TestSNPReader)
