#! /usr/bin/env python

import unittest    
from test     import test_support
from StringIO import *


from PTT.PTTReader import *

class TestPTTReader(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_raise_header_exception_genome_size(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome  1..2944528
2 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
"""
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertRaises(PTTHeaderError, reader._check_header, ptt_fh)
        
    def test_raise_header_exception_nb_prot(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome - 1..2944528
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
"""
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertRaises(PTTHeaderError, reader._check_header, ptt_fh)
        
    def test_raise_header_exception_columns(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome - 1..2944528
2 proteins
Location	Strand	Length	PID	Gene	Synonym\sCode	COG	Product
318..1673	+	451	16802049	dnaA	lmo0001 -	COG0593L	chromosomal replication initiation protein
4644..4865	-	73	16802052	-	lmo0004	-	-	-
"""
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertRaises(PTTHeaderError, reader.read_and_check, ptt_fh)

    def test_check_header_ok(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome - 1..2944528
2 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
"""
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertEquals( (2944528, 2), reader._check_header(ptt_fh) )

    def test_raise_exception_line(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome - 1..2944528
2 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
318..1673	+	451	16802049	dnaA	lmo0001\s-	COG0593L	chromosomal replication initiation protein
4644..4865	-	73	16802052	-	lmo0004	-	-	-
"""
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertRaises(PTTLineError, reader.read_and_check, ptt_fh)

    def test_raise_exception_header(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
318..1673	+	451	16802049	dnaA	lmo0001\s-	COG0593L	chromosomal replication initiation protein
4644..4865	-	73	16802052	-	lmo0004	-	-	-
"""
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertRaises(PTTHeaderError, reader.read_and_check, ptt_fh)


    def test_raise_exception_completeness(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome - 1..2944528
2 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
318..1673	+	451	16802049	dnaA	lmo0001	-	COG0593L	chromosomal replication initiation protein
"""
       
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        self.assertRaises(PTTCompletenessError, reader.read_and_check, ptt_fh)


    def test_read_ok(self):
        ptt_fh = StringIO()
        print >>ptt_fh, """Listeria monocytogenes EGD-e, complete genome - 1..2944528
2 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
318..1673	+	451	16802049	dnaA	lmo0001	-	COG0593L	chromosomal replication initiation protein
4644..4865	-	73	16802052	-	lmo0004	-	-	-
"""
        expected_dict = { "lmo0001" :{ "start"     : 318,
                                       "end"       : 1673,
                                       "strand"    : "+",
                                       "PID"       : "16802049",
                                       "gene"      : "dnaA",
                                       "COG"       : "COG0593L",
                                       "product"   : "chromosomal replication initiation protein"},
                          "lmo0004" : {"start"     : 4644,
                                       "end"       : 4865,
                                       "strand"    : "-",
                                       "PID"       : "16802052",
                                       "gene"      : "",
                                       "COG"       : "",
                                       "product"   : ""}
                          }
        
        ptt_fh.seek(0,0)
        reader      = PTTReader()
        result_dict = reader.read_and_check(ptt_fh)
        self.assertEquals(expected_dict, result_dict)

        
 
if __name__ == "__main__":
    test_support.run_unittest(TestPTTReader)
