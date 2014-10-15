#! /usr/bin/env python

import unittest    
from test     import test_support
from StringIO import *

import os
from string import *
import re

import SynTView.Formatter as Formatter

class TestSNP(unittest.TestCase):

    reformatter = Formatter.SynTViewSnpFormatter()
    ptt_dict    = None
    
    def setUp(self):
        TestSNP.ptt_dict = { 'lmo0003' : {"strand" : "+"},
                             'lmo0013' : {"strand" : "+"},
                             'lmo0130' : {"strand" : "+"},
                             'lmo2752' : {"strand" : "+"}}

        
    def test_raise_gene_error(self):
        snp_dict = {"type"      :"SNP",
                    "position"  :17131,
                    "reference" :"G",
                    "variant"   :"A",
                    "gene_name" :"qoxA",
                    "effect"    :"NON_SYNONYMOUS_CODING",
                    "biotype"   :"mRNA",
                    "old_aa"    :"E",
                    "new_aa"    :"K"}

        self.assertRaises(Formatter.SynTViewSnpFormatterGeneError, \
                          TestSNP.reformatter.format, snp_dict, TestSNP.ptt_dict)

        
    def test_raise_intra_error(self):
        snp_dict =  {"type"      :"SNP",
                     "position"  :2519935,
                     "reference" :"T",
                     "variant"   :"C",
                     "gene_name" :"lmo2448",
                     "effect"    :"INTRAGENIC: lmo2448",
                     "biotype"   :"",
                     "old_aa"    :"",
                     "new_aa"    :""}

        self.assertRaises(Formatter.SynTViewSnpFormatterIntragenicError, \
                          TestSNP.reformatter.format, snp_dict, TestSNP.ptt_dict)
        
    def test_pseudo(self):
        snp_dict =  {"type"      :"SNP",
                     "position"  :432305,
                     "reference" :"C",
                     "variant"   :"T",
                     "gene_name" :"Gene_lmo0410",
                     "effect"    :"INTRON",
                     "biotype"   :"pseudogene",
                     "old_aa"    :"",
                     "new_aa"    :""}
        expected_string = "432305\tC\tT\t\t\t\t"
        
        type, reformatted_string = TestSNP.reformatter.format(snp_dict,  TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("SNP", type)
                
        
        
    def test_inter(self):
        snp_dict =  {"type"      :"SNP",
                     "position"  :1745,
                     "reference" :"C",
                     "variant"   :"A",
                     "gene_name" :"",
                     "effect"    :"INTERGENIC",
                     "biotype"   :"",
                     "old_aa"    :"",
                     "new_aa"    :""}
        expected_string = "1745\tC\tA\t\t\t\t"
        
        type, reformatted_string = TestSNP.reformatter.format(snp_dict,  TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("SNP", type)
                

    def test_snp_in_gene(self):
        snp_dict = {"type"      :"SNP",
                    "position"  :17131,
                    "reference" :"G",
                    "variant"   :"A",
                    "gene_name" :"lmo0013",
                    "effect"    :"NON_SYNONYMOUS_CODING",
                    "biotype"   :"mRNA",
                    "old_aa"    :"E",
                    "new_aa"    :"K"}
        expected_string = "17131\tG\tA\tlmo0013\t+\tGlu\tLys"
        
        type, reformatted_string = TestSNP.reformatter.format(snp_dict,  TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("SNP", type)
                
    def test_snp_in_trna_gene(self):
        snp_dict = {"type"      :"SNP",
                    "position"  :82754,
                    "reference" :"T",
                    "variant"   :"G",
                    "gene_name" :"Gene_lmot01",
                    "effect"    :"INTRON",
                    "biotype"   :"tRNA",
                    "old_aa"    :"",
                    "new_aa"    :""}
        expected_string = "82754\tT\tG\t\t\t\t"
        
        type, reformatted_string = TestSNP.reformatter.format(snp_dict,  TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("SNP", type)
                

    def test_snpstop_in_gene(self):
        snp_dict = {"type"      :"SNP",
                    "position"  :2830008,
                    "reference" :"G",
                    "variant"   :"A",
                    "gene_name" :"lmo2752",
                    "effect"    :"SYNONYMOUS_STOP",
                    "biotype"   :"mRNA",
                    "old_aa"    :"*",
                    "new_aa"    :"*"}
        expected_string = "2830008\tG\tA\tlmo2752\t+\t***\t***"
        
        type, reformatted_string = TestSNP.reformatter.format(snp_dict, TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("SNP", type)
        
    def test_ins_in_gene(self):
        snp_dict = {"type"      :"INS",
                    "position"  :132275,
                    "reference" :"*",
                    "variant"   :"+A",
                    "gene_name" :"lmo0130",
                    "effect"    :"FRAME_SHIFT: Transcript_lmo0130",
                    "biotype"   :"",
                    "old_aa"    :"",
                    "new_aa"    :""}
        expected_string = "132275\tA\t1\t\tlmo0130"
        type, reformatted_string = TestSNP.reformatter.format(snp_dict, TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("INS", type)
        
                    
    def test_del_in_inter(self):
        snp_dict = {"type"      :"DEL",
                    "position"  :134723,
                    "reference" :"A",
                    "variant"   :"AG",
                    "gene_name" :"",
                    "effect"    :"INTERGENIC",
                    "biotype"   :"",
                    "old_aa"    :"",
                    "new_aa"    :""}
        expected_string = "134723\tG\t-1\t\t"
        
        type, reformatted_string = TestSNP.reformatter.format(snp_dict,  TestSNP.ptt_dict)
        self.assertEqual(expected_string, reformatted_string)
        self.assertEqual("DEL", type)

        
if __name__ == "__main__":
    test_support.run_unittest(TestSNP)
