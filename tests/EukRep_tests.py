#!/usr/bin/env python3

"""
Tests to check if package is installed and working correctly.
patrickwest@berkeley.edu
"""

import argparse
import os
import sys
from Bio import SeqIO

import EukRep.EukRep as EukRep

def load_reference_set(kmer_length):
    args.model = resource_stream(__name__, 'models/linsvm_160_%smer_1.0.pickle' % args.kmer_len)
    fh = open(kmer_length + '_reference.fa')

def compare_outputs(reference_db, test_db):
    reference_ids = []
    for record in SeqIO.parse(reference_db, "fasta"):
        reference_ids.append(record.id)

    for record in SeqIO.parse(test_db, "fasta"):
        if record.id not in reference_ids:
            return True

    return False

class PredictionTests():
    '''
    Ensure prediction results are as expected on a few test scaffolds
    '''

    def setUp(self):
        self.kmer_lengths = [3,4,5,6]
        self.dir_path = os.path.dirname(os.path.realpath(__file__))

    def run(self):
        self.setUp()
        for kmer in self.kmer_lengths:
            self.unit_test_01(kmer)
        self.tearDown()

    def unit_test_01(self,kmer_length):
        args = EukRep.Parse_Args(['-i', self.dir_path + '/test_sequences/test_scaffolds.fa','-o', '%smer_out.fa' % str(kmer_length),'-k', str(kmer_length)])
        EukRep.main(args)

        assert not compare_outputs(open(self.dir_path + '/test_sequences/test_scaffolds.fa.pred%s' % kmer_length), open('%smer_out.fa' % kmer_length)), \
            'Prediction with kmer length %s returned unexpected result' % kmer_length

    def tearDown(self):
        for kmer in self.kmer_lengths:
            if os.path.isfile('%smer_out.fa' % kmer):
                os.remove('%smer_out.fa' % kmer)


def prediction_test():
    ''' run simple prediction tests'''
    PredictionTests().run()


if __name__ == '__main__':
    prediction_test()

    print('EukRep prediction appears to be in working order')


#NEED TO TEST predictions
    #-All kmers
    #-test input/output at same time
    #-preferably make it easy to update test comparison datasets
