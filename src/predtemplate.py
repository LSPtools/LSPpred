import pandas as pd
import sys
import pickle
from Bio import SeqIO
from FeatureGen import Get_Protein_Feat,normalize_zscore
from Bio.SeqIO.FastaIO import SimpleFastaParser, FastaIterator, FastaWriter
import numpy as np
from tabulate import tabulate
import joblib
from collections import defaultdict


predmode = sys.argv[1]
filename = sys.argv[2]

all_feature_dict = {}

add_control_seq = False
control_seq = """MTKNYPTVSEDYKKAVEKCRRKLRGLIAEKNCAPIMVRLAWHSAGTFDCQSRTGGPFGTM
RFDAEQAHGANSGIHIALRLLDPIREQFPTISFADFHQLAGVVAVEVTGGPDIPFHPGRE
DKPQPPPEGRLPDATKGCDHLRDVFAKQMGLSDKDIVALSGAHTLGRCHKDRSGFEGAWT
SNPLIFDNSYFKELLSGEKEGLLQLVSDKALLDDPVFRPLVEKYAADEDAFFADYAEAHM
KLSELGFADA"""


mode_dict = defaultdict(dict)
mode_dict["LSPpred"]["normfile"]="trainingSetNormParams_LSPpred.csv"
mode_dict["LSPpred"]["featfile"]="LSPpred2_v20190730features.csv"
mode_dict["LSPpred"]["libfile"]="LSPpred2_v20190730.joblib"
mode_dict["LSPpred"]["threshold_val"]=0.68
mode_dict["SPLpred"]["normfile"]="trainingSetNormParams_SPLpred.csv"
mode_dict["SPLpred"]["featfile"]="SPLpred3_v20190730features.csv"
mode_dict["SPLpred"]["libfile"]="SPLpred3_v20190730.joblib"
mode_dict["SPLpred"]["threshold_val"]=0.7

with open(filename) as fasta:
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
                  if sequence.find('X') < 0: # no ambiguous
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

print(tabulate(results.loc[results.index != "control-seq-lsppred",['LSPpred_probability','LSPpred_LowConf','LSPpred','SPLpred_probability','SPLpred_LowConf','SPLpred','Either_LowConf','Either','Consensus_LowConf','Consensus']],headers='keys', tablefmt='psql'))
#print(tabulate(results.loc[results.index != "control-seq-lsppred",['LSPpred','LSPpred_LowConf','LSPpred_probability','SPLpred','SPLpred_LowConf','SPLpred_probability','Consensus','Consensus_LowConf','Either','Either_LowConf']],headers='keys', tablefmt='psql'))


print("# Summary")
print("# LSPpred_LowConf: "+str((results['LSPpred_LowConf']).sum()))
print("# LSPpred: "+str((results['LSPpred']).sum()))
print("# SPLpred_LowConf: "+str((results['SPLpred_LowConf']).sum()))
print("# SPLpred: "+str((results['SPLpred']).sum()))
print("# Either_LowConf: "+str((results['Either_LowConf']).sum()))
print("# Either: "+str((results['Either']).sum()))
print("# Consensus_LowConf: "+str((results['Consensus_LowConf']).sum()))
print("# Consensus: "+str((results['Consensus']).sum()))
print("# Total sequences: "+str(len(results.index)))
