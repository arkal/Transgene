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

import argparse
import logging
import unittest
import shutil
import os
import tempfile
import transgene

log = logging.getLogger(__name__)


class TransgeneTest(unittest.TestCase):
    """
    Tests functionality of transgene
    """
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.test_dir)
        self.pep_lens = None
        self.output = None
        self.output_files = set()

    def tearDown(self):
        shutil.rmtree(self.test_dir)
        os.chdir(self.cwd)

    def test_transgene_RNA_WXS_FUSION(self):
        print('Testing Transgene with FUSIONS using RNA and DNA.')
        self._test_transgene(use_RNA=True, use_DNA=True, fusions=True)

    def test_transgene_RNA_WXS_xFUSION(self):
        print('Testing Transgene without FUSIONS using RNA and DNA.')
        self._test_transgene(use_RNA=True, use_DNA=True, fusions=False)

    def test_transgene_RNA_xWXS_FUSION(self):
        print('Testing Transgene with FUSIONS using only RNA.')
        self._test_transgene(use_RNA=True, use_DNA=False, fusions=True)

    def test_transgene_RNA_xWXS_xFUSION(self):
        print('Testing Transgene without FUSIONS using only RNA.')
        self._test_transgene(use_RNA=True, use_DNA=False, fusions=False)

    def test_transgene_xRNA_WXS_FUSION(self):
        print('Testing Transgene with FUSIONS using only DNA.')
        self._test_transgene(use_RNA=False, use_DNA=True, fusions=True)

    def test_transgene_xRNA_WXS_xFUSION(self):
        print('Testing Transgene without FUSIONS using only DNA.')
        self._test_transgene(use_RNA=False, use_DNA=True, fusions=False)

    def test_transgene_xRNA_xWXS_FUSION(self):
        print('Testing Transgene with FUSIONS using no sequence-based filtering.')
        self._test_transgene(use_RNA=False, use_DNA=False, fusions=True)

    def test_transgene_xRNA_xWXS_xFUSION(self):
        print('Testing Transgene without FUSIONS using no sequence-based filtering.')
        self._test_transgene(use_RNA=False, use_DNA=False, fusions=False)

    def _test_transgene(self, use_RNA=False, use_DNA=False, fusions=False):
        """
        Tests that transgene runs and checks output files
        """
        params = argparse.Namespace()
        params.prefix = 'unit_test'
        params.pep_lens = '9,10,15'
        params.no_json_dumps = False
        params.log_level = 'DEBUG'
        params.logfile = None

        params.reject_threshold = 5
        params.rna_file = self._get_input_path('test_input/test_rna.bam') if use_RNA else None
        params.rna_min_alt_freq = 0.1 if use_RNA else None
        params.peptide_file = open(self._get_input_path('test_input/test.pc_translations.fa'))
        params.snpeff_file = open(self._get_input_path('test_input/snpeff_test.vcf'))
        params.transcript_file = open(self._get_input_path('test_input/test.pc_transcripts.fa'))

        params.filter_oxog = use_DNA
        params.dna_file = self._get_input_path('test_input/test_dna.bam') if use_DNA else None
        params.oxog_min_alt_freq = 0.1 if use_DNA else None

        params.fusion_file = open(self._get_input_path('test_input/test_fusions.'
                                                       'bedpe')) if fusions else None
        params.annotation_file = open(self._get_input_path('test_input/gencode.v19.chr21.'
                                                           'annotation.gtf')) if fusions else None
        params.genome_file = open(self._get_input_path('test_input/chr21.fa')) if fusions else None

        params.cores = 3
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
                self.output_files.add(filename)
        self.pep_lens = output.keys()
        self.output = output
        self.check_output(use_RNA, use_DNA, fusions)
        params.peptide_file.close()
        params.snpeff_file.close()
        params.transcript_file.close()
        if fusions:
            params.fusion_file.close()

    def check_output(self, test_with_rna_file, test_with_dna_file, test_with_fusions):
        """
        Check the output from transgene

        :param test_with_rna_file: Was this test run with the rna expression?
        :param test_with_dna_file: Was this test run with the DNA bam?
        :param test_with_fusions: Was this test looking at fusions?
        """
        alpha = 'ARNDCQEGHILKMFPSTWYVBZJUOX'
        expected_peptides = {
            '9mer': {'ELAGGGYVPSAPCPGET',  # ENST00000492084.1:L56P
                     'PRLYKIYRGRDSERAPA',  # ENST00000395952.3:E19G
                     'TAVTAPHSNSWDTYHQPRALEKH'},  # ENST00000395952.3:S42NXXXXXY48H
            '10mer': {'GELAGGGYVPSAPCPGETC',  # ENST00000492084.1:L56P
                      'GPRLYKIYRGRDSERAPAS',  # ENST00000395952.3:E19G
                      'PTAVTAPHSNSWDTYHQPRALEKHA'},  # ENST00000395952.3:S42NXXXXXY48H
            '15mer': {'ESLYSGELAGGGYVPSAPCPGETC',  # ENST00000492084.1:L56P
                      'LSCVLGPRLYKIYRGRDSERAPASVPETP',  # ENST00000395952.3:E19G
                      'SVPETPTAVTAPHSNSWDTYHQPRALEKHADSILA'}  # ENST00000395952.3:S42NXXXXXY48H
        }
        if not test_with_rna_file:
            expected_peptides['9mer'].update([
                'SLYSGELADGGYVPSAPCPGET',  # ENST00000492084.1:G51DXXXXL56P
                'SLYSGELADGGYVLSAP',  # ENST00000492084.1:G51D
                'SLYSGELADGGYLSLSK',  # ENST00000440843.2:G51D
                'PRLYKIYRGRDS',  # ENST00000395952.3:E19GXXXY17*
                'TAVTAPHSNSWDTYYQP',  # ENST00000395952.3:S42N
                'HSSSWDTYHQPRALEKH',  # ENST00000395952.3:Y48H
                'PWTVGKNELSQTVGEVF'  # ENST00000229729.10:F>L
            ])
            expected_peptides['10mer'].update([
                'ESLYSGELADGGYVPSAPCPGETC',  # ENST00000492084.1:G51DXXXXL56P
                'ESLYSGELADGGYVLSAPC',  # ENST00000492084.1:G51D
                'ESLYSGELADGGYLSLSKV',  # ENST00000440843.2:G51D
                'GPRLYKIYRGRDS',  # ENST00000395952.3:E19GXXXY17*
                'PTAVTAPHSNSWDTYYQPR',  # ENST00000395952.3:S42N
                'PHSSSWDTYHQPRALEKHA',  # ENST00000395952.3:Y48H
                'DPWTVGKNELSQTVGEVFY'  # ENST00000229729.10:F>L
            ])
            expected_peptides['15mer'].update([
                'LAWRPESLYSGELADGGYVPSAPCPGETC',  # ENST00000492084.1:G51DXXXXL56P
                'LAWRPESLYSGELADGGYVLSAPCPGETC',  # ENST00000492084.1:G51D
                'LAWRPESLYSGELADGGYLSLSKVVPFSH',  # ENST00000440843.2:G51D
                'LSCVLGPRLYKIYRGRDS',  # ENST00000395952.3:E19GXXXY17*
                'SVPETPTAVTAPHSNSWDTYYQPRALEKH',  # ENST00000395952.3:S42N
                'TAVTAPHSSSWDTYHQPRALEKHADSILA',  # ENST00000395952.3:Y48H
                'SSCPEDPWTVGKNELSQTVGEVFYTKNRN'  # ENST00000229729.10:F>L
            ])
            if not test_with_dna_file:
                # Then add the OXOG STUFF
                expected_peptides['9mer'].update([
                    'EFQNDFYRYCIRRSSPQ',  # ENST00000440843.2:S24Y
                    'ATGAPPRRKRVPGRACP'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['10mer'].update([
                    'LEFQNDFYRYCIRRSSPQP',  # ENST00000440843.2:S24Y
                    'VATGAPPRRKRVPGRACPW'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['15mer'].update([
                    'SLYPRLEFQNDFYRYCIRRSSPQPPPNLA',  # ENST00000440843.2:S24Y
                    'EVDTNVATGAPPRRKRVPGRACPWREPIR'  # ENST00000375331.2:Q>K
                ])
        if test_with_fusions:
            # Then add the OXOG STUFF
            expected_peptides['9mer'].update([
                'SQLETYKRQEDPKWEF',    # HOOK3-RET fusion
                'QSSSYGQQTASGDMQT',    # EWSR1-ATF1 fusion
                'VVCTQPKSPSSTPVSP',    # TMPRSS2-ETV1 fusion
                'QSSSYGQQSPPLGGAQ',    # EWSR1-FLI fusion
                'NSKMALNSEALSVVSE'   # TMPRSS2-ERG fusion
            ])
            expected_peptides['10mer'].update([
                'RSQLETYKRQEDPKWEFP',  # HOOK3-RET fusion
                'QQSSSYGQQTASGDMQTY',  # EWSR1-ATF1 fusion
                'PVVCTQPKSPSSTPVSPL',  # TMPRSS2-ETV1 fusion
                'QQSSSYGQQSPPLGGAQT',  # EWSR1-FLI fusion
                'DNSKMALNSEALSVVSED'  # TMPRSS2-ERG fusion
            ])
            expected_peptides['15mer'].update([
                'KANAARSQLETYKRQEDPKWEFPRKNLV',  # HOOK3-RET fusion
                'PSQYSQQSSSYGQQTASGDMQTYQIRTT',  # EWSR1-ATF1 fusion
                'TQASNPVVCTQPKSPSSTPVSPLHHASP',  # TMPRSS2-ETV1 fusion
                'PSQYSQQSSSYGQQSPPLGGAQTISKNT',  # EWSR1-FLI fusion
                'LLDAVDNSKMALNSEALSVVSEDQSLFE'  # TMPRSS2-ERG fusion
            ])

        # Compare test output to expected output
        for kmer in self.output.keys():
            observed_fasta = self.output[kmer]['fasta']
            observed_seqs = set()
            for header, _, seq in transgene.read_fasta(open(observed_fasta, 'r'), alpha):
                observed_seqs.add(seq)
            if observed_seqs != expected_peptides[kmer]:
                if observed_seqs - expected_peptides[kmer]:
                    # There are more observed than expected
                    print('Unexpected {}s called by transgene: {}'.format(kmer, ','.join(
                        observed_seqs - expected_peptides[kmer])))
                if expected_peptides[kmer] - observed_seqs:
                    # There are more expected than observed
                    print('Transgene failed to predict some {}s: {}'.format(kmer, ','.join(
                        expected_peptides[kmer] - observed_seqs)))
                raise RuntimeError

    @staticmethod
    def _get_input_path(file_path):
        """
        Return the absolute path to the input files relative to the project directory

        :param str file_path: relative path to file
        :return: absolute path to file
        """
        project_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(project_path, file_path)
