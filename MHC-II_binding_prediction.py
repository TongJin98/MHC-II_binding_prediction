'''
Parse command line to get
(1) Name of sequence file (seqfile)
(2) Name of prediction method (pmethod)
(3) Name of MHC molecule (geneName)

Loop over sequences
     Loop over 15-mers within sequence
          Test that 15-mer contains only valid amino acids
               Call function that takes 15-mer and geneName and returns affinity, 9-mer
'''

import argparse
import sys
from Bio import SeqIO
import csv
import math
import re
import os

def main(args):

    # Read in command line arguments to variables
    input_dir = args.i
    output_dir = args.o

    if output_dir[-1] is not "/":
        output_dir += "/"

    try:
        os.makedirs(output_dir)
    except:
        pass

    filename = input_dir.split('.')[-2]
    output_filename = f"{output_dir}{filename}"
    output_data = parse_input(input_dir)
    generate_output_file(output_filename, output_data)



def parse_input(infile):
    '''Parse the original input file of peptide sequence to store all valid
       fifteen mers within each sequence.

    Attributes
    ----------
        infile (str): full path to input file
    '''

    # dictionary to store resulting data
    data = {'seqfile': [],
            'pmethod': [],
            'geneName' : [],
            'geneName': []
            }

    for seq_record in SeqIO.parse(infile, "fasta"):
        num_fifteen_mer = len(seq_record.seq)-14
        nonstd_aas =  'BJOUXZ'

        for i in range(num_fifteen_mer):
            fifteen_mer = seq_record.seq[i:i+15]
            if not any(c in fifteen_mer for c in nonstd_aas):
                data['seqfile'].append(fifteen_mer)
                data['geneName'].append(seq_record.description)

    return data



def generate_output_file(fileName, output_data):
    '''generate the out file

    Attributes
    ----------
        fileName (str): file name for the anchor file generated
        data (dictionary): list of parsed fifteen_mer along with their gene_name
    '''
    fileName = fileName + '.csv'
    with open(fileName, 'w') as csv_file:
        fieldnames = ['gene','fifteen_mer']
        csv_writer = csv.DictWriter(csv_file,fieldnames=fieldnames,delimiter=';')
        csv_writer.writeheader()

        for geneName, seqfile in zip( output_data['geneName'], output_data['seqfile']):
            csv_writer.writerow({'gene': geneName, 'fifteen_mer': seqfile})



if __name__ == '__main__':

    # Set commend line arugments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'path to the input file')
    parser.add_argument('-o', help = 'path to the output file')
    args = parser.parse_args()

    if (args.i == None or args.o == None):
        print("Command line arugment error\nCorrect Usage:\npython MHC-II_binding_prediction -i <full path of input file> -o <full path to output file>")
        sys.exit()
    main(args)
