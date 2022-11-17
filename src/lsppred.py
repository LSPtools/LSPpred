#!/usr/bin/env python
'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) ANDREW LONSDALE, 17 Sep 2020
License     : MIT
Maintainer  : ANDREW.LONSDALE@LONSBIO.COM.AU
Portability : POSIX

Predict LSP status for plant proteins, provided in a single FASTA file.
'''

from argparse import ArgumentParser, FileType
import sys
import logging
import pandas as pd
import sys
import os
import pickle
from Bio import SeqIO
from FeatureGen import Get_Protein_Feat,normalize_zscore
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaIterator, FastaWriter
import numpy as np
from tabulate import tabulate
import joblib
from collections import defaultdict
import pkg_resources

all_feature_dict = {}

add_control_seq = False
control_seq = """MTKNYPTVSEDYKKAVEKCRRKLRGLIAEKNCAPIMVRLAWHSAGTFDCQSRTGGPFGTM
RFDAEQAHGANSGIHIALRLLDPIREQFPTISFADFHQLAGVVAVEVTGGPDIPFHPGRE
DKPQPPPEGRLPDATKGCDHLRDVFAKQMGLSDKDIVALSGAHTLGRCHKDRSGFEGAWT
SNPLIFDNSYFKELLSGEKEGLLQLVSDKALLDDPVFRPLVEKYAADEDAFFADYAEAHM
KLSELGFADA"""
from pathlib import Path
script_dir=Path(os.path.realpath(__file__)).parent 



# Settings
mode_dict = defaultdict(dict)
mode_dict["LSPpred"]["normfile"]=(script_dir / "../model/trainingSetNormParams_LSPpred.csv").resolve()
mode_dict["LSPpred"]["featfile"]=(script_dir / "../model/LSPpred2_v20190730features.csv").resolve()
mode_dict["LSPpred"]["libfile"]=(script_dir / "../model/LSPpred2_v20190730.joblib").resolve()
mode_dict["LSPpred"]["threshold_val"]=0.68
mode_dict["SPLpred"]["normfile"]=(script_dir /"../model/trainingSetNormParams_SPLpred.csv").resolve()
mode_dict["SPLpred"]["featfile"]=(script_dir /"../model/SPLpred3_v20190730features.csv").resolve()
mode_dict["SPLpred"]["libfile"]=(script_dir /"../model/SPLpred3_v20190730.joblib").resolve()
mode_dict["SPLpred"]["threshold_val"]=0.7


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
HEADER = "LSPpred version: 1.0.0"
PROGRAM_NAME = "lsppred"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Predict LSPs in a FASTA file'
    parser = ArgumentParser(description=description)

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--low', action='store_true', help='Add low confidence predictions')
    parser.add_argument('--output',
                        metavar='OUT_FILE',
                        type=FileType('w'), default= sys.stdout,
                        help='save output in CSV format to OUT_FILE')
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
#    parser.add_argument('--format',
#                            metavar='LOG_FILE',
#                            type=str,
#                            default="csv",
#                            help='output format as csv or text table ')
    parser.add_argument('fasta_file',
                        metavar='FASTA_FILE',
                        type=FileType('r'),
                        help='Input FASTA protein file')
    return parser.parse_args()

def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z")
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def process_files(options):


    #with open(options.fasta_file) as fasta:
    with options.fasta_file as fasta:
              for record in SimpleFastaParser(fasta): #SeqIO.SimpleFastaParser(fasta):
                  sequence = record[1]
                  title = record[0]
                  seq_id = title.split(None, 1)[0]
    #for record in SeqIO.parse(fasta, "fasta"):
    #              sequence = record.seq
    #              title = record.id
    #              seq_id = title.split("|")[0]
                  #print(seq_id ,title)
                  if len(sequence) > 45:
                      if sequence.find('X') < 0 and sequence.find('B') < 0  and sequence.find('Z') < 0: # no ambiguous
                          all_feature_dict[seq_id]=Get_Protein_Feat(sequence)
    #print(all_feature_dict)

    all_feature_dict["control-seq-lsppred"] = Get_Protein_Feat(control_seq.replace("\n",""))

    df = pd.DataFrame(all_feature_dict).transpose()

    combined_results =[]

    for pred_mod in ["LSPpred","SPLpred"]:

        clf = joblib.load(mode_dict[pred_mod]["libfile"])
        # normalsiazton data needs to be loaded and the new features normalsiexed
        normParamDf = pd.read_csv(mode_dict[pred_mod]["normfile"], index_col=[0])
        z_df, normParams = normalize_zscore(df.copy(),saveCSV = False,filename='NA.csv',normParams = normParamDf)

        # only want to run on the seletced features
        features_df = pd.read_csv(mode_dict[pred_mod]["featfile"])

        results = pd.DataFrame()
        # trey to remov missing 0 but htey may not be ther:
        try:
            pred = clf.predict_proba(z_df.loc[:,features_df['Features']])
        except Exception as e:
            print(e)
            print("unable to predict")
            results = z_df.copy()
            results[pred_mod] = "Unable"
            results[pred_mod+"_probability"] = -1

        if results.empty:
            results = z_df.copy()
            results[pred_mod+"_probability"] = pred[:,1]
            results[pred_mod] =  results[pred_mod+'_probability'] >= mode_dict[pred_mod]["threshold_val"]
            results[pred_mod+"_LowConf"] =  results[pred_mod+'_probability'] >= 0.5

        results.index.name="Sequence"
        #print(results)
        combined_results.append(results.copy())
        #combined_results.append(results.loc[results.index != "control-seq-lsppred",['LSP','LSP_probability']])

    results = pd.merge(combined_results[0],combined_results[1],left_index=True, right_index=True)
    #print(results.loc[results.index != "control-seq-lsppred",['LSPpred','LSPpred_probability','SPLpred','SPLpred_probability']])
    #print(results.loc[results.index != "control-seq-lsppred",['LSP','LSP_probability']].to_string(header=True))
    #print(results.loc[results.index != "control-seq-lsppred",['LSP','LSP_probability']].to_csv(header=True))
    #print("# threshold:",threshold_val)
    results['Consensus'] = results.LSPpred & results.SPLpred
    results['Consensus_LowConf'] = results.LSPpred_LowConf & results.SPLpred_LowConf
    results['Either'] = results.LSPpred | results.SPLpred
    results['Either_LowConf'] = results.LSPpred_LowConf | results.SPLpred_LowConf

#    if (options.format == "text"):
    if (False):
        print(tabulate(results.loc[results.index != "control-seq-lsppred",['LSPpred_probability','LSPpred_LowConf','LSPpred','SPLpred_probability','SPLpred_LowConf','SPLpred','Either_LowConf','Either','Consensus_LowConf','Consensus']],headers='keys', tablefmt='psql'))
        #print(tabulate(results.loc[results.index != "control-seq-lsppred",['LSPpred','LSPpred_LowConf','LSPpred_probability','SPLpred','SPLpred_LowConf','SPLpred_probability','Consensus','Consensus_LowConf','Either','Either_LowConf']],headers='keys', tablefmt='psql'))


    else:
        if options.low:
            print("LSPpred_LowConf: "+str((results['LSPpred_LowConf']).sum()-1))
        print("LSPpred: "+str((results['LSPpred']).sum()-1))
        if options.low:
            print("SPLpred_LowConf: "+str((results['SPLpred_LowConf']).sum()))
        print("SPLpred: "+str((results['SPLpred']).sum()))
        if options.low:
            print("Either_LowConf: "+str((results['Either_LowConf']).sum()-1))
        print("Either: "+str((results['Either']).sum()-1))
        if options.low:
            print("Consensus_LowConf: "+str((results['Consensus_LowConf']).sum()))
        print("Consensus: "+str((results['Consensus']).sum()))
        print("Total sequences: "+str(len(results.index)-1)) # remove control seq
        print("Output written to: "+ options.output.name)
        if options.low:

            print(results.loc[results.index != "control-seq-lsppred",['LSPpred_probability','LSPpred_LowConf','LSPpred','SPLpred_probability','SPLpred_LowConf','SPLpred','Either_LowConf','Either','Consensus_LowConf','Consensus']].to_csv(),file=options.output)


        else:
            print(results.loc[results.index != "control-seq-lsppred",['LSPpred_probability','LSPpred','SPLpred_probability','SPLpred','Either','Consensus']].to_csv(),file=options.output)
def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    print(HEADER)
    process_files(options)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
