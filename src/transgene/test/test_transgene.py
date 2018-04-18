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
import os
import shutil
import tempfile
import unittest

from transgene.core import main as transgene_main
from transgene.common import file_type, get_exons, read_fasta
from transgene.indel import get_codon, get_exon_start_pos

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
        self.output_fastas = None
        self.output_files = set()
        self.annotator = '3.6'

    def tearDown(self):
        shutil.rmtree(self.test_dir)
        os.chdir(self.cwd)

    def test_transgene_RNA_WXS_FUSION_GENOME(self):
        print('Testing Transgene with FUSIONS using RNA and DNA with a genome fasta.')
        self._test_transgene(use_RNA=True, use_DNA=True, fusions=True, genome=True)

    # def test_transgene_RNA_WXS_FUSION_xGENOME(self): Fusion True cannot pair with genome false

    def test_transgene_RNA_WXS_xFUSION_GENOME(self):
        print('Testing Transgene without FUSIONS using RNA and DNA with a genome fasta.')
        self._test_transgene(use_RNA=True, use_DNA=True, fusions=False, genome=True)

    def test_transgene_RNA_WXS_xFUSION_xGENOME(self):
        print('Testing Transgene without FUSIONS using RNA and DNA without a genome fasta.')
        self._test_transgene(use_RNA=True, use_DNA=True, fusions=False, genome=False)

    def test_transgene_RNA_xWXS_FUSION_GENOME(self):
        print('Testing Transgene with FUSIONS using only RNA with a genome fasta.')
        self._test_transgene(use_RNA=True, use_DNA=False, fusions=True, genome=True)

    # def test_transgene_RNA_xWXS_FUSION_xGENOME(self): Fusion True cannot pair with genome false

    def test_transgene_RNA_xWXS_xFUSION_GENOME(self):
        print('Testing Transgene without FUSIONS using only RNA with a genome fasta.')
        self._test_transgene(use_RNA=True, use_DNA=False, fusions=False, genome=True)

    def test_transgene_RNA_xWXS_xFUSION_xGENOME(self):
        print('Testing Transgene without FUSIONS using only RNA without a genome fasta.')
        self._test_transgene(use_RNA=True, use_DNA=False, fusions=False, genome=False)

    def test_transgene_xRNA_WXS_FUSION_GENOME(self):
        print('Testing Transgene with FUSIONS using only DNA with a genome fasta.')
        self._test_transgene(use_RNA=False, use_DNA=True, fusions=True, genome=True)

    # def test_transgene_xRNA_WXS_FUSION_xGENOME(self): Fusion True cannot pair with genome false

    def test_transgene_xRNA_WXS_xFUSION_GENOME(self):
        print('Testing Transgene without FUSIONS using only DNA with a genome fasta.')
        self._test_transgene(use_RNA=False, use_DNA=True, fusions=False, genome=True)

    def test_transgene_xRNA_WXS_xFUSION_xGENOME(self):
        print('Testing Transgene without FUSIONS using only DNA without a genome fasta.')
        self._test_transgene(use_RNA=False, use_DNA=True, fusions=False, genome=False)

    def test_transgene_xRNA_xWXS_FUSION_GENOME(self):
        print('Testing Transgene with FUSIONS using no sequence-based filtering with a genome '
              'fasta.')
        self._test_transgene(use_RNA=False, use_DNA=False, fusions=True, genome=True)

    # def test_transgene_xRNA_xWXS_FUSION_xGENOME(self): Fusion True cannot pair with genome false

    def test_transgene_xRNA_xWXS_xFUSION_GENOME(self):
        print('Testing Transgene without FUSIONS using no sequence-based filtering with a genome '
              'fasta.')
        self._test_transgene(use_RNA=False, use_DNA=False, fusions=False, genome=True)

    def test_transgene_xRNA_xWXS_xFUSION_xGENOME(self):
        print('Testing Transgene without FUSIONS using no sequence-based filtering without a '
              'genome fasta.')
        self._test_transgene(use_RNA=False, use_DNA=False, fusions=False, genome=False)

    def _test_transgene(self, use_RNA=False, use_DNA=False, fusions=False, genome=True):
        """
        Tests that transgene runs and checks output files
        """
        if fusions:
            # Genome HAS to be true if fusions has been selected
            genome = True
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
        if self.annotator != 'vep':
            params.snpeff_file = open(
                self._get_input_path('test_input/snpeff_test.vcf')
                if self.annotator == '3.6' else
                self._get_input_path('test_input/snpeff_test_' + self.annotator + '.vcf')
            )
            params.vep_file = None
        else:
            params.snpeff_file = None
            params.vep_file = open(self._get_input_path('test_input/snpeff_test_vep.vcf'))
        params.transcript_file = open(self._get_input_path('test_input/test.pc_transcripts.fa'))
        params.extend_length = 30
        params.filter_oxog = use_DNA
        params.dna_file = self._get_input_path('test_input/test_dna.bam') if use_DNA else None
        params.oxog_min_alt_freq = 0.1 if use_DNA else None

        params.fusion_file = params.annotation_file = params.genome_file = None
        if genome:
            params.annotation_file = file_type(self._get_input_path(
                'test_input/gencode.v19.chr6.chr21.annotation.gtf.gz'))
            params.genome_file = file_type(self._get_input_path('test_input/chr6_chr21.fa.gz'))
            if fusions:
                params.fusion_file = open(self._get_input_path('test_input/test_fusions.'
                                                               'bedpe'))

        params.filter_mt = params.filter_rg = params.filter_ig = params.filter_rt = True
        params.rt_threshold = 100000

        params.cores = 1
        try:
            transgene_main(params)
        except RuntimeError as e:
            if self.annotator == '4.3' and genome is False:
                self.assertTrue(e.message == 'Cannot handle SNPEff4.3 annotated vcfs without a '
                                             'genome fasta and a gencode annotation.')
                return
        output = {'9mer': {'tumor': 'unit_test_tumor_9_mer_peptides.faa',
                           'map': 'unit_test_tumor_9_mer_peptides.faa.map',
                           'normal': 'unit_test_normal_9_mer_peptides.faa'},
                  '10mer': {'tumor': 'unit_test_tumor_10_mer_peptides.faa',
                            'map': 'unit_test_tumor_10_mer_peptides.faa.map',
                            'normal': 'unit_test_normal_10_mer_peptides.faa'},
                  '15mer': {'tumor': 'unit_test_tumor_15_mer_peptides.faa',
                            'map': 'unit_test_tumor_15_mer_peptides.faa.map',
                            'normal': 'unit_test_normal_15_mer_peptides.faa'}}
        for key,  data in output.iteritems():
            for feature, filename in data.iteritems():
                assert os.path.exists(filename)
                self.output_files.add(filename)
        self.pep_lens = output.keys()
        self.output_fastas = output
        self.check_output(use_RNA, use_DNA, fusions, genome, self.annotator)
        params.peptide_file.close()
        if params.snpeff_file:
            params.snpeff_file.close()
        if params.vep_file:
            params.vep_file.close()
        params.transcript_file.close()
        if fusions:
            params.fusion_file.close()

    def check_output(self, test_with_rna_file, test_with_dna_file, test_with_fusions,
                     test_with_genome_files, annotator):
        """
        Check the output from transgene

        :param bool test_with_rna_file: Was this test run with the rna expression?
        :param bool test_with_dna_file: Was this test run with the DNA bam?
        :param bool test_with_fusions: Was this test looking at fusions?
        :param bool test_with_genome_files: Was this test run with a genome fasta?
        :param str annotator: The annotator being tested
        """
        alpha = 'ARNDCQEGHILKMFPSTWYVBZJUOX'
        expected_peptides = {
            '9mer': {'tumor': {'PRLYKIYRGRDSERAPA',  # ENST00000395952.3:E19G
                               'TAVTAPHSNSWDTYHQPRALEKH'},  # ENST00000395952.3:S42NXXXXXY48H
                     'normal': {'PRLYKIYRERDSERAPA',  # ENST00000395952.3:E19G
                                'TAVTAPHSSSWDTYYQPRALEKH'}  # ENST00000395952.3:S42NXXXXXY48H
                     },
            '10mer': {'tumor': {'GPRLYKIYRGRDSERAPAS',  # ENST00000395952.3:E19G
                                'PTAVTAPHSNSWDTYHQPRALEKHA'},  # ENST00000395952.3:S42NXXXXXY48H
                      'normal': {'GPRLYKIYRERDSERAPAS',  # ENST00000395952.3:E19G
                                 'PTAVTAPHSSSWDTYYQPRALEKHA'}  # ENST00000395952.3:S42NXXXXXY48H
                      },
            '15mer': {'tumor': {'LSCVLGPRLYKIYRGRDSERAPASVPETP',  # ENST00000395952.3:E19G
                                'SVPETPTAVTAPHSNSWDTYHQPRALEKHADSILA'},  # ENST00000395952.3:S42NXXXXXY48H
                      'normal': {'LSCVLGPRLYKIYRERDSERAPASVPETP',  # ENST00000395952.3:E19G
                                 'SVPETPTAVTAPHSSSWDTYYQPRALEKHADSILA'}  # ENST00000395952.3:S42NXXXXXY48H
                      }
        }
        if annotator == '3.6':
            # This is not a valid variant in other annotators
            expected_peptides['9mer']['tumor'].update([
                'ELAGGGYVPSAPCPGET'  # ENST00000492084.1:L56P
            ])
            expected_peptides['9mer']['normal'].update([
                'ELAGGGYVLSAPCPGET'  # ENST00000492084.1:L56P
            ])
            expected_peptides['10mer']['tumor'].update([
                'GELAGGGYVPSAPCPGETC'  # ENST00000492084.1:L56P
            ])
            expected_peptides['10mer']['normal'].update([
                'GELAGGGYVLSAPCPGETC'  # ENST00000492084.1:L56P
            ])
            expected_peptides['15mer']['tumor'].update([
                'ESLYSGELAGGGYVPSAPCPGETC'  # ENST00000492084.1:L56P
            ])
            expected_peptides['15mer']['normal'].update([
                'ESLYSGELAGGGYVLSAPCPGETC'  # ENST00000492084.1:L56P
            ])

        if not test_with_rna_file:
            expected_peptides['9mer']['tumor'].update([
                'SLYSGELADGGYLSLSK',  # ENST00000440843.2:G51D
                'PRLYKIYRGRDS',  # ENST00000395952.3:E19GXXXY17*
                'TAVTAPHSNSWDTYYQP',  # ENST00000395952.3:S42N
                'TAVTAPHSTSWDTYYQP',  # ENST00000395952.3:S42T
                'HSSSWDTYHQPRALEKH',  # ENST00000395952.3:Y48H
                'TAVTAPHSTSWDTYHQPRALEKH',  # ENST00000395952.3:S42TXXXXXY48H
                'PWTVGKNELSQTVGEVF',  # ENST00000229729.10:F>L
                'PRSSPYGRRWWGVNAEP',  # ENST00000376560.3:G100R
                'RRQRRERRIRRYLSAGR',  # ENST00000376148.4:F37I
                'WLHNLCDVHLEAVKPVL'  # ENST00000321897.5:Y821H
            ])
            expected_peptides['9mer']['normal'].update([
                'SLYSGELAGGGYLSLSK',  # ENST00000440843.2:G51D
                'PRLYKIYRERDS',  # ENST00000395952.3:E19GXXXY17*
                'TAVTAPHSSSWDTYYQP',  # ENST00000395952.3:S42N
                'TAVTAPHSSSWDTYYQP',  # ENST00000395952.3:S42T
                'HSSSWDTYYQPRALEKH',  # ENST00000395952.3:Y48H
                'TAVTAPHSSSWDTYYQPRALEKH',  # ENST00000395952.3:S42TXXXXXY48H
                'PWTVGKNEFSQTVGEVF',  # ENST00000229729.10:F>L
                'PRSSPYGRGWWGVNAEP',  # ENST00000376560.3:G100R
                'RRQRRERRFRRYLSAGR',  # ENST00000376148.4:F37I
                'WLHNLCDVYLEAVKPVL'  # ENST00000321897.5:Y821H
            ])
            expected_peptides['10mer']['tumor'].update([
                'ESLYSGELADGGYLSLSKV',  # ENST00000440843.2:G51D
                'GPRLYKIYRGRDS',  # ENST00000395952.3:E19GXXXY17*
                'PTAVTAPHSNSWDTYYQPR',  # ENST00000395952.3:S42N
                'PTAVTAPHSTSWDTYYQPR',  # ENST00000395952.3:S42T
                'PHSSSWDTYHQPRALEKHA',  # ENST00000395952.3:Y48H
                'PTAVTAPHSTSWDTYHQPRALEKHA',  # ENST00000395952.3:S42TXXXXXY48H
                'DPWTVGKNELSQTVGEVFY',  # ENST00000229729.10:F>L
                'GPRSSPYGRRWWGVNAEPP',  # ENST00000376560.3:G100R
                'SRRQRRERRIRRYLSAGRL',  # ENST00000376148.4:F37I
                'FWLHNLCDVHLEAVKPVLW'  # ENST00000321897.5:Y821H
            ])
            expected_peptides['10mer']['normal'].update([
                'ESLYSGELAGGGYLSLSKV',  # ENST00000440843.2:G51D
                'GPRLYKIYRERDS',  # ENST00000395952.3:E19GXXXY17*
                'PTAVTAPHSSSWDTYYQPR',  # ENST00000395952.3:S42N
                'PTAVTAPHSSSWDTYYQPR',  # ENST00000395952.3:S42T
                'PHSSSWDTYYQPRALEKHA',  # ENST00000395952.3:Y48H
                'PTAVTAPHSSSWDTYYQPRALEKHA',  # ENST00000395952.3:S42TXXXXXY48H
                'DPWTVGKNEFSQTVGEVFY',  # ENST00000229729.10:F>L
                'GPRSSPYGRGWWGVNAEPP',  # ENST00000376560.3:G100R
                'SRRQRRERRFRRYLSAGRL',  # ENST00000376148.4:F37I
                'FWLHNLCDVYLEAVKPVLW'  # ENST00000321897.5:Y821H
            ])
            expected_peptides['15mer']['tumor'].update([
                'LAWRPESLYSGELADGGYLSLSKVVPFSH',  # ENST00000440843.2:G51D
                'LSCVLGPRLYKIYRGRDS',  # ENST00000395952.3:E19GXXXY17*
                'SVPETPTAVTAPHSNSWDTYYQPRALEKH',  # ENST00000395952.3:S42N
                'SVPETPTAVTAPHSTSWDTYYQPRALEKH',  # ENST00000395952.3:S42T
                'TAVTAPHSSSWDTYHQPRALEKHADSILA',  # ENST00000395952.3:Y48H
                'SVPETPTAVTAPHSTSWDTYHQPRALEKHADSILA',  # ENST00000395952.3:S4INXXXXXY48H
                'SSCPEDPWTVGKNELSQTVGEVFYTKNRN',  # ENST00000229729.10:F>L
                'IRRGLGPRSSPYGRRWWGVNAEPPFPGPG',  # ENST00000376560.3:G100R
                'SMASTSRRQRRERRIRRYLSAGRLVRAQA',  # ENST00000376148.4:F37I
                'HALHHFWLHNLCDVHLEAVKPVLWHSPRP'  # ENST00000321897.5:Y821H
            ])
            expected_peptides['15mer']['normal'].update([
                'LAWRPESLYSGELAGGGYLSLSKVVPFSH',  # ENST00000440843.2:G51D
                'LSCVLGPRLYKIYRERDS',  # ENST00000395952.3:E19GXXXY17*
                'SVPETPTAVTAPHSSSWDTYYQPRALEKH',  # ENST00000395952.3:S42N
                'SVPETPTAVTAPHSSSWDTYYQPRALEKH',  # ENST00000395952.3:S42T
                'TAVTAPHSSSWDTYYQPRALEKHADSILA',  # ENST00000395952.3:Y48H
                'SVPETPTAVTAPHSSSWDTYYQPRALEKHADSILA',  # ENST00000395952.3:S4INXXXXXY48H
                'SSCPEDPWTVGKNEFSQTVGEVFYTKNRN',  # ENST00000229729.10:F>L
                'IRRGLGPRSSPYGRGWWGVNAEPPFPGPG',  # ENST00000376560.3:G100R
                'SMASTSRRQRRERRFRRYLSAGRLVRAQA',  # ENST00000376148.4:F37I
                'HALHHFWLHNLCDVYLEAVKPVLWHSPRP'  # ENST00000321897.5:Y821H
            ])
            if annotator == '3.6':
                # This is not a valid variant in other annotators
                expected_peptides['9mer']['tumor'].update([
                    'SLYSGELADGGYVPSAPCPGET',  # ENST00000492084.1:G51DXXXXL56P
                    'SLYSGELADGGYVLSAP',  # ENST00000492084.1:G51D
                ])
                expected_peptides['9mer']['normal'].update([
                    'SLYSGELAGGGYVLSAPCPGET',  # ENST00000492084.1:G51DXXXXL56P
                    'SLYSGELAGGGYVLSAP',  # ENST00000492084.1:G51D
                ])
                expected_peptides['10mer']['tumor'].update([
                    'ESLYSGELADGGYVPSAPCPGETC',  # ENST00000492084.1:G51DXXXXL56P
                    'ESLYSGELADGGYVLSAPC',  # ENST00000492084.1:G51D
                ])
                expected_peptides['10mer']['normal'].update([
                    'ESLYSGELAGGGYVLSAPCPGETC',  # ENST00000492084.1:G51DXXXXL56P
                    'ESLYSGELAGGGYVLSAPC',  # ENST00000492084.1:G51D
                ])
                expected_peptides['15mer']['tumor'].update([
                    'LAWRPESLYSGELADGGYVPSAPCPGETC',  # ENST00000492084.1:G51DXXXXL56P
                    'LAWRPESLYSGELADGGYVLSAPCPGETC',  # ENST00000492084.1:G51D
                ])
                expected_peptides['15mer']['normal'].update([
                    'LAWRPESLYSGELAGGGYVLSAPCPGETC',  # ENST00000492084.1:G51DXXXXL56P
                    'LAWRPESLYSGELAGGGYVLSAPCPGETC',  # ENST00000492084.1:G51D
                ])
            if not test_with_dna_file:
                # Then add the OXOG STUFF
                expected_peptides['9mer']['tumor'].update([
                    'EFQNDFYRYCIRRSSPQ',  # ENST00000440843.2:S24Y
                    'ATGAPPRRKRVPGRACP'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['9mer']['normal'].update([
                    'EFQNDFYRSCIRRSSPQ',  # ENST00000440843.2:S24Y
                    'ATGAPPRRQRVPGRACP'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['10mer']['tumor'].update([
                    'LEFQNDFYRYCIRRSSPQP',  # ENST00000440843.2:S24Y
                    'VATGAPPRRKRVPGRACPW'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['10mer']['normal'].update([
                    'LEFQNDFYRSCIRRSSPQP',  # ENST00000440843.2:S24Y
                    'VATGAPPRRQRVPGRACPW'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['15mer']['tumor'].update([
                    'SLYPRLEFQNDFYRYCIRRSSPQPPPNLA',  # ENST00000440843.2:S24Y
                    'EVDTNVATGAPPRRKRVPGRACPWREPIR'  # ENST00000375331.2:Q>K
                ])
                expected_peptides['15mer']['normal'].update([
                    'SLYPRLEFQNDFYRSCIRRSSPQPPPNLA',  # ENST00000440843.2:S24Y
                    'EVDTNVATGAPPRRQRVPGRACPWREPIR'  # ENST00000375331.2:Q>K
                ])
        else:
            # Indels will only show up if there was an rna file AND a genome file
            if test_with_genome_files:
                expected_peptides['9mer']['tumor'].update(
                    ['PRSSPYGRRW',  # ENST00000376560.3:G100RXWG102*
                     'WLHNLCDVHLEAVKPVL',  # ENST00000321897.5:Y821H
                     'ASTSRRQRQRTSHSSLLVCRTAGPGPGPPPATPRPRCR',  # ENST00000376148.4:R33Q?XXX_SNV_
                     ])
                expected_peptides['9mer']['normal'].update(
                    ['PRSSPYGRGW',  # ENST00000376560.3:G100RXWG102*
                     'WLHNLCDVYLEAVKPVL',  # ENST00000321897.5:Y821H
                     'ASTSRRQRR',  # ENST00000376148.4:R33Q?XXX_SNV_
                     ])
                expected_peptides['10mer']['tumor'].update(
                    ['GPRSSPYGRRW',  # ENST00000376560.3:G100RXWG102*
                     'FWLHNLCDVHLEAVKPVLW',  # ENST00000321897.5:Y821H
                     'MASTSRRQRQRTSHSSLLVCRTAGPGPGPPPATPRPRCR',  # ENST00000376148.4:R33Q?XXX_SNV_
                     ])
                expected_peptides['10mer']['normal'].update(
                    ['GPRSSPYGRGW',  # ENST00000376560.3:G100RXWG102*
                     'FWLHNLCDVYLEAVKPVLW',  # ENST00000321897.5:Y821H
                     'MASTSRRQRR',  # ENST00000376148.4:R33Q?XXX_SNV_
                     ])
                expected_peptides['15mer']['tumor'].update(
                    ['IRRGLGPRSSPYGRRW',  # ENST00000376560.3:G100RXWG102*
                     'HALHHFWLHNLCDVHLEAVKPVLWHSPRP',  # ENST00000321897.5:Y821H
                     'RPKSSMASTSRRQRQRTSHSSLLVCRTAGPGPGPPPATPRPRCR', # ENST00000376148.4:R33Q?XXX_SNV_
                     ])
                expected_peptides['15mer']['normal'].update(
                    ['IRRGLGPRSSPYGRGW',  # ENST00000376560.3:G100RXWG102*
                     'HALHHFWLHNLCDVYLEAVKPVLWHSPRP',  # ENST00000321897.5:Y821H
                     'RPKSSMASTSRRQRR',  # ENST00000376148.4:R33Q?XXX_SNV_
                     ])
                if annotator == '4.0':
                    if test_with_genome_files:
                        expected_peptides['9mer']['tumor'].update(
                            ['AARILDISTVTQLRSLVIDLERTAKELDTNTLKVHTLR',  # ENST00000307859.4:RA168T
                             'HFWLHNLCQDVYLVSEAGGGLVFPCLLLIPLEISKAES',  # ENST00000321897.5:D819QD
                             'QVLFKGQGPSTHVLLTHTISRIAVSYQTKVNLLSAIK',  # ENST00000449264.2:GC144G
                             'MGEEINAAQARHLVEQRRGSPMEPVHYVLVPFHDPGKC',  # ENST00000375688.4:AKI334A
                             'ATLPRPPHQLAPPGPAAGAPPPGCATLPRMPPDPYLQE'  # ENST00000211413.5:H87QL
                             ])
                        expected_peptides['9mer']['normal'].update(
                            ['AARILDISR',  # ENST00000307859.4:RA168T
                             'HFWLHNLCD',  # ENST00000321897.5:D819QD
                             'QVLFKGQGC',  # ENST00000449264.2:GC144G
                             'MGEEINAAK',  # ENST00000375688.4:AKI334A
                             'ATLPRPPHH'  # ENST00000211413.5:H87QL
                             ])
                        expected_peptides['10mer']['tumor'].update(
                            ['EAARILDISTVTQLRSLVIDLERTAKELDTNTLKVHTLR',  # ENST00000307859.4:RA168T
                             'HHFWLHNLCQDVYLVSEAGGGLVFPCLLLIPLEISKAES',  # ENST00000321897.5:D819QD
                             'SQVLFKGQGPSTHVLLTHTISRIAVSYQTKVNLLSAIK',  # ENST00000449264.2:GC144G
                             'SMGEEINAAQARHLVEQRRGSPMEPVHYVLVPFHDPGKC',  # ENST00000375688.4:AKI334A
                             'SATLPRPPHQLAPPGPAAGAPPPGCATLPRMPPDPYLQE'  # ENST00000211413.5:H87QL
                             ])
                        expected_peptides['10mer']['normal'].update(
                            ['EAARILDISR',  # ENST00000307859.4:RA168T
                             'HHFWLHNLCD',  # ENST00000321897.5:D819QD
                             'SQVLFKGQGC',  # ENST00000449264.2:GC144G
                             'SMGEEINAAK',  # ENST00000375688.4:AKI334A
                             'SATLPRPPHH'  # ENST00000211413.5:H87QL
                             ])
                        expected_peptides['15mer']['tumor'].update(
                            ['EHMPAEAARILDISTVTQLRSLVIDLERTAKELDTNTLKVHTLR',
                             # ENST00000307859.4:RA168T
                             'VTHALHHFWLHNLCQDVYLVSEAGGGLVFPCLLLIPLEISKAES',
                             # ENST00000321897.5:D819QD
                             'LYLIYSQVLFKGQGPSTHVLLTHTISRIAVSYQTKVNLLSAIK',
                             # ENST00000449264.2:GC144G
                             'LDTTGSMGEEINAAQARHLVEQRRGSPMEPVHYVLVPFHDPGKC',
                             # ENST00000375688.4:AKI334A
                             'RGPSSSATLPRPPHQLAPPGPAAGAPPPGCATLPRMPPDPYLQE'
                             # ENST00000211413.5:H87QL
                             ])
                        expected_peptides['15mer']['normal'].update(
                            ['EHMPAEAARILDISR',  # ENST00000307859.4:RA168T
                             'VTHALHHFWLHNLCD',  # ENST00000321897.5:D819QD
                             'LYLIYSQVLFKGQGC',  # ENST00000449264.2:GC144G
                             'LDTTGSMGEEINAAK',  # ENST00000375688.4:AKI334A
                             'RGPSSSATLPRPPHH',  # ENST00000211413.5:H87QL
                             ])
                else:
                    expected_peptides['9mer']['tumor'].update(
                        ['AARILDISTVTQLRSLV',  # ENST00000307859.4:RA168T
                         'HFWLHNLCQDVYLEAVK',  # ENST00000321897.5:D819QD
                         'QVLFKGQGPSTHVLLT',  # ENST00000449264.2:GC144G
                         'MGEEINAAQARHLVEQ',  # ENST00000375688.4:AKI334A
                         'ATLPRPPHQLAPPGPAAG'  # ENST00000211413.5:H87QL
                         ])
                    expected_peptides['9mer']['normal'].update(
                        ['AARILDISRAVTQLRSLV',  # ENST00000307859.4:RA168T
                         'HFWLHNLCDVYLEAVK',  # ENST00000321897.5:D819QD
                         'QVLFKGQGCPSTHVLLT',  # ENST00000449264.2:GC144G
                         'MGEEINAAKIQARHLVEQ',  # ENST00000375688.4:AKI334A
                         'ATLPRPPHHAPPGPAAG'  # ENST00000211413.5:H87QL
                         ])
                    expected_peptides['10mer']['tumor'].update(
                        ['EAARILDISTVTQLRSLVI',  # ENST00000307859.4:RA168T
                         'HHFWLHNLCQDVYLEAVKP',  # ENST00000321897.5:D819QD
                         'SQVLFKGQGPSTHVLLTH',  # ENST00000449264.2:GC144G
                         'SMGEEINAAQARHLVEQR',  # ENST00000375688.4:AKI334A
                         'SATLPRPPHQLAPPGPAAGA'  # ENST00000211413.5:H87QL
                         ])
                    expected_peptides['10mer']['normal'].update(
                        ['EAARILDISRAVTQLRSLVI',  # ENST00000307859.4:RA168T
                         'HHFWLHNLCDVYLEAVKP',  # ENST00000321897.5:D819QD
                         'SQVLFKGQGCPSTHVLLTH',  # ENST00000449264.2:GC144G
                         'SMGEEINAAKIQARHLVEQR',  # ENST00000375688.4:AKI334A
                         'SATLPRPPHHAPPGPAAGA'  # ENST00000211413.5:H87QL
                         ])
                    expected_peptides['15mer']['tumor'].update(
                        ['EHMPAEAARILDISTVTQLRSLVIDLERT',  # ENST00000307859.4:RA168T
                         'VTHALHHFWLHNLCQDVYLEAVKPVLWHS',  # ENST00000321897.5:D819QD
                         'LYLIYSQVLFKGQGPSTHVLLTHTISRI',  # ENST00000449264.2:GC144G
                         'LDTTGSMGEEINAAQARHLVEQRRGSPM',  # ENST00000375688.4:AKI334A
                         'RGPSSSATLPRPPHQLAPPGPAAGAPPPGC'  # ENST00000211413.5:H87QL
                         ])
                    expected_peptides['15mer']['normal'].update(
                        ['EHMPAEAARILDISRAVTQLRSLVIDLERT',  # ENST00000307859.4:RA168T
                         'VTHALHHFWLHNLCDVYLEAVKPVLWHS',  # ENST00000321897.5:D819QD
                         'LYLIYSQVLFKGQGCPSTHVLLTHTISRI',  # ENST00000449264.2:GC144G
                         'LDTTGSMGEEINAAKIQARHLVEQRRGSPM',  # ENST00000375688.4:AKI334A
                         'RGPSSSATLPRPPHHAPPGPAAGAPPPGC',  # ENST00000211413.5:H87QL
                         ])
            else:
                # The chained INDEL test SNVS will show up in the output as just SNVs
                expected_peptides['9mer']['tumor'].update(
                    ['PRSSPYGRRWWGVNAEP',  # ENST00000376560.3:G100R
                     'RRQRRERRIRRYLSAGR',  # ENST00000376148.4:F37I
                     'WLHNLCDVHLEAVKPVL'  # ENST00000321897.5:Y821H
                     ])
                expected_peptides['9mer']['normal'].update(
                    ['PRSSPYGRGWWGVNAEP',  # ENST00000376560.3:G100R
                     'RRQRRERRFRRYLSAGR',  # ENST00000376148.4:F37I
                     'WLHNLCDVYLEAVKPVL'  # ENST00000321897.5:Y821H
                     ])
                expected_peptides['10mer']['tumor'].update(
                    ['GPRSSPYGRRWWGVNAEPP',  # ENST00000376560.3:G100R
                     'SRRQRRERRIRRYLSAGRL',  # ENST00000376148.4:F37I
                     'FWLHNLCDVHLEAVKPVLW'  # ENST00000321897.5:Y821H
                     ])
                expected_peptides['10mer']['normal'].update(
                    ['GPRSSPYGRGWWGVNAEPP',  # ENST00000376560.3:G100R
                     'SRRQRRERRFRRYLSAGRL',  # ENST00000376148.4:F37I
                     'FWLHNLCDVYLEAVKPVLW'  # ENST00000321897.5:Y821H
                     ])
                expected_peptides['15mer']['tumor'].update(
                    ['IRRGLGPRSSPYGRRWWGVNAEPPFPGPG',  # ENST00000376560.3:G100R
                     'SMASTSRRQRRERRIRRYLSAGRLVRAQA',  # ENST00000376148.4:F37I
                     'HALHHFWLHNLCDVHLEAVKPVLWHSPRP'  # ENST00000321897.5:Y821H
                     ])
                expected_peptides['15mer']['normal'].update(
                    ['IRRGLGPRSSPYGRGWWGVNAEPPFPGPG',  # ENST00000376560.3:G100R
                     'SMASTSRRQRRERRFRRYLSAGRLVRAQA',  # ENST00000376148.4:F37I
                     'HALHHFWLHNLCDVYLEAVKPVLWHSPRP'  # ENST00000321897.5:Y821H
                     ])
        if test_with_fusions:
            # Then add the OXOG STUFF
            expected_peptides['9mer']['tumor'].update([
                'SQLETYKRQEDPKWEF',    # HOOK3-RET fusion
                'QSSSYGQQTASGDMQT',    # EWSR1-ATF1 fusion
                'VVCTQPKSPSSTPVSP',    # TMPRSS2-ETV1 fusion
                'QSSSYGQQSPPLGGAQ',    # EWSR1-FLI fusion
                'NSKMALNSEALSVVSE'   # TMPRSS2-ERG fusion
            ])
            expected_peptides['10mer']['tumor'].update([
                'RSQLETYKRQEDPKWEFP',  # HOOK3-RET fusion
                'QQSSSYGQQTASGDMQTY',  # EWSR1-ATF1 fusion
                'PVVCTQPKSPSSTPVSPL',  # TMPRSS2-ETV1 fusion
                'QQSSSYGQQSPPLGGAQT',  # EWSR1-FLI fusion
                'DNSKMALNSEALSVVSED'  # TMPRSS2-ERG fusion
            ])
            expected_peptides['15mer']['tumor'].update([
                'KANAARSQLETYKRQEDPKWEFPRKNLV',  # HOOK3-RET fusion
                'PSQYSQQSSSYGQQTASGDMQTYQIRTT',  # EWSR1-ATF1 fusion
                'TQASNPVVCTQPKSPSSTPVSPLHHASP',  # TMPRSS2-ETV1 fusion
                'PSQYSQQSSSYGQQSPPLGGAQTISKNT',  # EWSR1-FLI fusion
                'LLDAVDNSKMALNSEALSVVSEDQSLFE'  # TMPRSS2-ERG fusion
            ])

        # Compare test output to expected output
        for ttype in 'tumor', 'normal':
            for kmer in self.output_fastas.keys():
                observed_fasta = self.output_fastas[kmer][ttype]
                observed_seqs = set()
                for header, _, seq in core.read_fasta(open(observed_fasta, 'r'), alpha):
                    observed_seqs.add(seq)
                if observed_seqs != expected_peptides[kmer][ttype]:
                    if observed_seqs - expected_peptides[kmer][ttype]:
                        # There are more observed than expected
                        print('Unexpected {} {}s called by transgene: {}'.format(
                            ttype, kmer, ','.join(observed_seqs - expected_peptides[kmer][ttype])))
                    if expected_peptides[kmer][ttype] - observed_seqs:
                        # There are more expected than observed
                        print('Transgene failed to predict some {} {}s: {}'.format(
                            ttype, kmer, ','.join(expected_peptides[kmer][ttype] - observed_seqs)))
                    raise RuntimeError

    def test_get_ref_pos_alt_aa(self):
        want = [('A123B', ['A', 123, 'B']),
                ('A123', ['A', 123, '']),
                ('ABC132BAT', ['ABC', 132, 'BAT']),
                ('BCA132B?', ['BCA', 132, 'B?']),
                ('B132*', ['B', 132, '*'])
                ]
        for string, parsed_string in want:
            get = core.get_ref_pos_alt_aa(string)
            self.assertListEqual(parsed_string, get)

    def test_get_exon_start_pos(self):
        genome_file = file_type(self._get_input_path('test_input/chr6_chr21.fa.gz'))
        annotation_file = file_type(self._get_input_path(
            'test_input/gencode.v19.chr6.chr21.annotation.gtf.gz'))
        genes_to_translate = {'PRRT1',  # negative strand
                              'TRIM40'  # positive strand
                              }
        exons, cds_starts = get_exons(genome_file, annotation_file,
                                      genes_to_translate)

        want = [(32119598, False),
                (32119586, False),
                (32118240, False),
                (32117499, False),
                (30105105, True),
                (30113771, True),
                (30114984, True),
                ]
        get = [get_exon_start_pos('ENST00000211413.5', 3, exons, cds_starts),
               get_exon_start_pos('ENST00000211413.5', 7, exons, cds_starts),
               get_exon_start_pos('ENST00000211413.5', 155, exons, cds_starts),
               get_exon_start_pos('ENST00000211413.5', 187, exons, cds_starts),
               get_exon_start_pos('ENST00000307859.4', 98, exons, cds_starts),
               get_exon_start_pos('ENST00000307859.4', 116, exons, cds_starts),
               get_exon_start_pos('ENST00000307859.4', 193, exons, cds_starts),
        ]
        self.assertListEqual(want, get)

    def test_get_codon(self):
        genome_file = file_type(self._get_input_path('test_input/chr6_chr21.fa.gz'))
        annotation_file = file_type(self._get_input_path(
            'test_input/gencode.v19.chr6.chr21.annotation.gtf.gz'))
        genes_to_translate = {'PRRT1',  # negative strand
                              'TRIM40'  # positive strand
                              }
        exons, cds_starts = get_exons(genome_file, annotation_file,
                                      genes_to_translate)

        genome_file.seek(0)
        chrom = ''
        for header, comment, seq in read_fasta(genome_file, 'ACGTN'):
            if header == 'chr6':
                chrom = seq

        want = ['GCT',
                'GGC',
                'GTG',
                'GGA',  # splice
                'CAC',
                'AGG',  # splice
                'ACA',
                ]

        get = [get_codon('ENST00000211413.5', 264, exons, cds_starts, chrom),
               get_codon('ENST00000211413.5', 239, exons, cds_starts, chrom),
               get_codon('ENST00000211413.5', 179, exons, cds_starts, chrom),
               get_codon('ENST00000211413.5', 7, exons, cds_starts, chrom),
               get_codon('ENST00000307859.4', 113, exons, cds_starts, chrom),
               get_codon('ENST00000307859.4', 201, exons, cds_starts, chrom),
               get_codon('ENST00000307859.4', 191, exons, cds_starts, chrom),
        ]
        self.assertListEqual(want, get)

    @staticmethod
    def _get_input_path(file_path):
        """
        Return the absolute path to the input files relative to the project directory

        :param str file_path: relative path to file
        :return: absolute path to file
        """
        project_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(project_path, file_path)
