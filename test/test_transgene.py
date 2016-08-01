#!/usr/bin/env python2.7
# Copyright 2016 Arjun Arkal Rao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : protect/test/test_file_downloads.py
"""
from __future__ import print_function
import os
import argparse
import logging
import unittest
import transgene
import collections

log = logging.getLogger(__name__)

class TransgeneTest(unittest.TestCase):
    """
    Tests functionality of transgene
    """
    def setUp(self):
        self.pep_lens = None
        self.output = None
        self.output_files = []

    def tearDown(self):
        for f in self.output_files:
            os.remove(f)

    def test_transgene(self):
        """
        Tests that transgene runs and checks output files
        """
        params = argparse.Namespace()
        params.prefix = 'unit_test'
        params.pep_lens = '9,10,15'
        params.no_json_dumps = False
        params.input_file = open('test_input/test.pc_translations.fa')
        params.snpeff_file = open('test_input/snpeff_test.vcf')
        transgene.main(params)
        output = {'9mer': {'fasta': 'unit_test_tumor_9_mer_snpeffed.faa',
                           'map': 'unit_test_tumor_9_mer_snpeffed.faa.map'},
                  '10mer': {'fasta': 'unit_test_tumor_10_mer_snpeffed.faa',
                            'map': 'unit_test_tumor_10_mer_snpeffed.faa.map'},
                  '15mer': {'fasta': 'unit_test_tumor_15_mer_snpeffed.faa',
                            'map': 'unit_test_tumor_15_mer_snpeffed.faa.map'}}
        for key,  data in output.iteritems():
            for feature, filename in data.iteritems():
                assert os.path.exists(filename)
                self.output_files.append(filename)
        self.pep_lens = output.keys()
        self.output = output
        self.check_ouptut()


    def check_ouptut(self):
        alpha = 'ARNDCQEGHILKMFPSTWYVBZJUOX'
        expected_fastas = [('9mer', 'test_input/expected_9mer.fa'),
                           ('10mer', 'test_input/expected_10mer.fa'),
                           ('15mer', 'test_input/expected_15mer.fa')]
        expected_peptides = collections.defaultdict(set)
        for kmer, fasta in expected_fastas:
            for _, _, seq in transgene.read_fasta(open(fasta, 'r'), alpha):
                expected_peptides[kmer].add(seq)
        # Compare test output to expected output
        for kmer in self.output.keys():
            observed_fasta = self.output[kmer]['fasta']
            for header, _, seq in transgene.read_fasta(open(observed_fasta, 'r'), alpha):
                assert seq in expected_peptides[kmer], 'Unexpected {}: {}'.format(kmer, seq)


