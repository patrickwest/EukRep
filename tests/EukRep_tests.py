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

def load_test_sequences(file_path):
    records = []
    fh = open(file_path)
    for record in SeqIO.parse(fh, "fasta"):
        records.append(record)
    return records

def compare_outputs(reference_db, test_db, out_name):
    reference_ids = {}
    print('Expect the following sequences to be predicted to be eukaryotic:')
    for record in SeqIO.parse(reference_db, "fasta"):
        reference_ids[record.id] = 0
        print(record.id)

    #Check for enexptected output
    print('The following sequences were predicted to be eukaryotic:')
    for record in SeqIO.parse(test_db, "fasta"):
        print(record.id)
        if record.id not in reference_ids:
            print('\nUnexpected scaffold %s in %s' % (record.id, out_name))
            return True
        else:
            reference_ids[record.id] = 1

    #Check for missing output
    for scaffold in reference_ids:
        if reference_ids[scaffold] is 0:
            print('\nMissing scaffold %s in %s' % (scaffold, out_name))
            return True
    return False

class PredictionTests():
    '''
    Ensure prediction results are as expected on a few test scaffolds
    '''

    def setUp(self):
        self.kmer_lengths = [3,4,5,6]
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.test_seqs = load_test_sequences(self.dir_path + '/test_sequences/test_scaffolds.fa')

    def run(self):
        self.setUp()

        print('\n------------------------------------------------------------')
        print('Prediction tests to be performed on %s total sequences' % len(self.test_seqs) )
        print('Test sequences located in ' + self.dir_path + '/test_sequences/test_scaffolds.fa')
        print('------------------------------------------------------------')

        for kmer in self.kmer_lengths:
            self.unit_test_01(kmer)
        self.tearDown()

    def unit_test_01(self,kmer_length):
        print('\nRunning test predictions with %smers...' % kmer_length)

        args = EukRep.Parse_Args(['-i', self.dir_path + '/test_sequences/test_scaffolds.fa','-o', '%smer_out.fa' % str(kmer_length),'-k', str(kmer_length)])
        EukRep.main(args)

        assert not compare_outputs(open(self.dir_path + '/test_sequences/test_scaffolds.fa.pred%s' % kmer_length), open('%smer_out.fa' % kmer_length), '%smer_out.fa' % str(kmer_length)), \
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

    print('\nEukRep prediction appears to be in working order\n')


#NEED TO TEST predictions
    #-All kmers
    #-test input/output at same time
    #-preferably make it easy to update test comparison datasets
