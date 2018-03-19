#!/usr/bin/env python3

#Patrick West 2017

import os
import sys
import numpy
import pickle
import argparse
import warnings
from argparse import RawTextHelpFormatter
from kpal.klib import Profile
from numpy import array
from sklearn import svm
from Bio import SeqIO
from io import StringIO
from random import randint
from pkg_resources import resource_stream

def main(args):

    #check arguments are valid
    args = check_args(args)

    #Open files
    outfile = open(args.o,'w')
    #open prokarya output file
    if args.prokarya is None:
        prok_fh = None
    else:
        prok_fh = open(args.prokarya,'w')

    #Open model
    #model_fh = open(args.model,'rb')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        model = pickle.load(args.model)

    #Perform predictions
    euk_scafs, prok_scafs = Make_Predictions(args.i, args.min, 5000, args.kmer_len, model, args.tie)

    #Write output
    if args.seq_names == True:
        print_seq_names(outfile, euk_scafs, prok_scafs, prok_fh)
    else:
        print_contigs_as_fa(args.i, outfile, euk_scafs, prok_fh)

    #cleanup
    args.model.close()
    outfile.close()
    if prok_fh is not None:
        prok_fh.close()

def check_args(args):
    '''
    Ensure user provided arguments won't cause issues
    '''

    #Outfile
    if os.path.isfile(args.o) and not args.ff:
        print('Outfile: %s already exists.' \
            % args.o, file=sys.stderr)
        exit()

    #Prokaryote Outfile
    if args.prokarya is not None:
        if os.path.isfile(args.prokarya) and not args.ff:
            print('Outfile: %s already exists.' \
                % args.prokarya, file=sys.stderr)
            exit()

    #Kmer length
    if args.kmer_len is None:
        args.kmer_len = 5
    elif int(args.kmer_len) > 7 or int(args.kmer_len) < 3 and args.model is None:
        print('Specified kmer length: %s is invalid. Please choose a length between 3-7 unless using a custom trained model.' \
            % args.kmer_len, file=sys.stderr)
        exit()
    else:
        args.kmer_len = int(args.kmer_len)

    #Model
    if args.model is None:
        #Mode
        if args.mode is None:
            args.model = resource_stream(__name__, 'models/linsvm_160_%smer_1.1.balanced.pickle' % args.kmer_len)
        elif args.mode is 'lenient':
            args.model = resource_stream(__name__, 'models/linsvm_160_%smer_1.1.pickle' % args.kmer_len)
        elif args.mode is 'strict':
            args.model = resource_stream(__name__, 'models/linsvm_160_%smer_1.1.strict.pickle' % args.kmer_len)
        elif args.mode is 'balanced':
            args.model = resource_stream(__name__, 'models/linsvm_160_%smer_1.1.balanced.pickle' % args.kmer_len)
        #args.model = os.path.join(sys.path[0], "models/linsvm_160_%smer_1.0.pickle" % args.kmer_len)
    else:
        args.model = open(args.model,'rb')

    #Minimum sequence length
    if args.min is None:
        args.min = 3000
    else:
        args.min = int(args.min)

    #Tie
    if args.tie is None:
        args.tie = "euk"
    elif args.tie != "euk"\
        and args.tie != "prok"\
        and args.tie != "rand"\
        and args.tie != "skip":

        print('%s is an invalid --tie parameter. Please choose euk, prok, rand, or skip.'\
            % args.tie, file=sys.stderr)
        exit()

    return args

def print_contigs_as_fa(fa_file_name, out_file, euk_ids, prokarya):
    '''
    Write predicted euk and predicted prok scaffolds in fa format to their respective output files
    '''

    fa_fh = open(fa_file_name)

    for record in SeqIO.parse(fa_fh, "fasta"):
        if record.id in euk_ids:
            out_file.write('>' + record.description + '\n')
            out_file.write(str(record.seq) + '\n')

        elif prokarya is not None:
            prokarya.write('>' + record.description + '\n')
            prokarya.write(str(record.seq) + '\n')

def print_seq_names(outfile, euk_ids, prok_ids, prokarya):
    '''
    Write predicted euk and predicted prok scaffolds names to their respective output files
    '''

    for line in euk_ids:
        outfile.write(line + "\n")

    if prokarya is not None:
        for line in prok_ids:
            prokarya.write(line + "\n")


def Make_Predictions(fa_file_name, min_size, max_size, kmer_size, model, tie):
    '''
    Read in fasta file, chop into 5kb parts, calculate kmer frequencies,
    and make predictions using provided trained machine learning model
    '''

    #Classified Sequences
    euk_seqs = []
    prok_seqs = []

    #open file
    fh = open(fa_file_name)

    #Iterate through fasta
    for record in SeqIO.parse(fh, "fasta"):
        s = StringIO(str(record.seq))
        split_seqs = []
        kmer_freqs = []
        seq_name = record.id

        #split sequence into 5kb max_size chunks
        split_seqs = chunk_sequence(s, min_size, max_size)

        #Calculate kmer frequences for each chunk
        kmer_freqs = calc_kmer_freqs(split_seqs, kmer_size)

        #Make prediction on each chunk
        #scikit learn requires .reshape(1,-1) for 1d arrays so if only one seq exists for a contig transform it
        if len(kmer_freqs) == 1:
            predictions = model.predict(numpy.asarray(kmer_freqs).reshape(1,-1))
        elif len(kmer_freqs) >= 1:
            predictions = model.predict(kmer_freqs)

        if len(kmer_freqs) >= 1:
            euk_seqs, prok_seqs = classify_by_majority_rule(predictions, seq_name, euk_seqs, prok_seqs, tie)

    fh.close()

    return euk_seqs, prok_seqs

def chunk_sequence(sequence, min_size, max_size):
    '''
    Cut sequences longer than 5kb into 5kb chunks and exclude trailing sequences
    if shorter than user specified min_length
    '''

    split_seqs = []
    while True:
        chunk = sequence.read(max_size)
        if len(chunk) >= min_size:
            split_seqs.append(chunk)
        else:
            break

    return split_seqs

def calc_kmer_freqs(split_seqs, kmer_size):
    '''
    Use kpal to calculate kmer frequencies for split sequences
    '''

    kmer_freqs = []
    for seq in split_seqs:
        temp_list = []

        #for some reason this kmer counter function only works on iterable(str) type objects.
        temp_list.append(str(seq))
        ktable = Profile.from_sequences(temp_list, kmer_size, name=None)

        #skip sequences with a lot of Ns/characters besides A|T|C|G
        if len(str(seq)) < 3000:
            if ktable.total >= len(str(seq))/2:
                ktable.counts = [count/ktable.total for count in ktable.counts]
                kmer_freqs.append(ktable.counts)
        else:
            if ktable.total >= 1500:
                ktable.counts = [count/ktable.total for count in ktable.counts]
                kmer_freqs.append(ktable.counts)

    return kmer_freqs

def classify_by_majority_rule(predictions, seq_name, euk_seqs, prok_seqs, tie):
    '''
    Tally predictions for each 5kb chunk comprising a fasta sequence and determine
    classification by majority rule of the 5kb chunks
    '''

    #Classify sequence by majority rule
    prok_total = 0
    euk_total = 0
    chunk_total = 0
    for pred in predictions:
        if pred == "bact" or pred == "arch" or "bact" in pred or "arch" in pred:
            prok_total += 1
        else:
            euk_total += 1
        chunk_total += 1
    if euk_total > prok_total:
        euk_seqs.append(seq_name)
    elif euk_total == prok_total:
        if tie == "euk":
            euk_seqs.append(seq_name)
        elif tie == "prok":
            prok_seqs.append(seq_name)
        elif tie == "rand":
            if (randint(0, 100) >= 50):
                euk_seqs.append(seq_name)
            else:
                prok_seqs.append(seq_name)
        #if tie == skip, skip it
    else:
        prok_seqs.append(seq_name)

    return euk_seqs, prok_seqs

class Log(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("%s_log.txt", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass

def Parse_Args(args):
    # argument parser
    parser = argparse.ArgumentParser(description='Identify sequences of predicted eukaryotic origin from a nucleotide fasta file. Individual sequences are split into 5kb chunks. Prediction is performed on each 5kb chunk and sequence origin is determined by majority rule of the chunks.', formatter_class=RawTextHelpFormatter)

    # input/output
    parser.add_argument(\
        '-i', required = True,\
        help = 'input fasta file')
    parser.add_argument(\
        '-o', required = True,\
        help = 'output file name')

    # options
    parser.add_argument(\
        '-ff',\
        action = 'store_true', help = 'Force overwrite of existing output files')
    parser.add_argument(\
        '--min',\
        help = 'Minimum sequence length cutoff for sequences to be included in prediction. Default is 3kb')
    parser.add_argument(\
        '--model',\
        help = 'Path to an alternate trained linear SVM model. Default is lin_svm_160_3.0.pickle')
    parser.add_argument(\
        '-k', '--kmer_len',\
        help = 'Kmer length to use for making predictions. Lengths between 3-7bp are available by default. If using a custom trained model, specify kmer length here.')
    parser.add_argument(\
        '--prokarya',\
        help = 'Name of file to output predicted prokaryotic sequences to. Default is to not output prokaryotic sequences.')
    parser.add_argument(\
        '--seq_names',\
        action = 'store_true', help = 'Only output fasta headers of identified sequences. Default is full fasta entry')
    parser.add_argument('--mode',\
        choices=['strict', 'balanced', 'lenient'],\
        help = 'Not compatable with --model.\n\
        How stringent the algorithm is in identifying eukaryotic scaffolds. Strict has a lower false positive rate and true positive rate; vice verso for leneient. Default is balanced.')
    #parser.add_argument(\
        #'--pre_merged', action = 'store_true',\
        #help = 'Output predictions for each sub-chunk of input sequences. Default is predictions for full merged sequences')
    parser.add_argument(\
        '--tie',\
        help = 'Specify how to handle cases where an equal number of a sequences chunks are predicted to be of eukaryotic and prokaryotic origin (Generally occurs infrequently).\n\
        euk = classify as euk\n\
        prok = classify as prok\n\
        rand = assign randomly\n\
        skip = do not classify\n\
        Default is to classify as eukaryotic.')

    return parser.parse_args(args)

if __name__ == '__main__':

    args = Parse_Args(sys.argv[1:])
    main(args)
