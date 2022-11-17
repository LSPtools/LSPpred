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
mode_dict["LSPpred"][normfile]="trainingSetNormParams_LSPpred.csv"
mode_dict["LSPpred"][featfile]="LSPpred2_v20190730features.csv"
mode_dict["LSPpred"][libfile]="LSPpred2_v20190730.joblib"
mode_dict["LSPpred"][threshold_val]=0.68
mode_dict["SPLpred"][normfile]="trainingSetNormParams_SPLpred.csv"
mode_dict["SPLpred"][featfile]="SPLpred3_v20190730features.csv"
mode_dict["SPLpred"][libfile]="SPLpred3_v20190730.joblib"
mode_dict["SPLpred"][threshold_val]=0.7

if (predmode == "lsp"):

    normfile="trainingSetNormParams_LSPpred.csv"
    featfile="LSPpred2_v20190730features.csv"
    libfile="LSPpred2_v20190730.joblib"
    threshold_val=0.68
else:

    libfile="SPLpred3_v20190730.joblib"
    featfile="SPLpred3_v20190730features.csv"
    normfile="trainingSetNormParams_SPLpred.csv"
    threshold_val=0.7

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

clf = joblib.load(libfile)
df = pd.DataFrame(all_feature_dict).transpose()


# normalsiazton data needs to be loaded and the new features normalsiexed
normParamDf = pd.read_csv(normfile, index_col=[0])
z_df, normParams = normalize_zscore(df,saveCSV = False,filename='NA.csv',normParams = normParamDf)
df = z_df


#featfile


#print(df.head())
#print(df[(np.isnan(df))])
#print(np.all(np.isfinite(df)))
#df = df.fillna(0)
#df = df.dropna(axis=1)


# only want to run on the seletced features
# featfile
features_df = pd.read_csv(featfile)


results = pd.DataFrame()
# trey to remov missing 0 but htey may not be ther:
try:
    pred = clf.predict_proba(df.loc[:,features_df['Features']])
    #pred = clf.predict_proba(df.drop(missing, axis=1, inplace=False))
except Exception as e:
    #try:
    #    pred = clf.predict_proba(df)
    #except ValueError:
    print(e)
    print("unable to predict")
    print(df.shape)
    print(df.columns)
    results = df
    results["LSP"] = "Unable"
    results["LSP_probability"] = -1
#print(pred[:,1])

if results.empty:
    results = df.copy()
    results['LSP_probability'] = pred[:,1]
    results['LSP'] =  results['LSP_probability'] >= threshold_val
results.index.name="Sequence"
#print(results.loc[results.index != "control-seq-lsppred",['LSP','LSP_probability']].to_string(header=True))
#print(results.loc[results.index != "control-seq-lsppred",['LSP','LSP_probability']].to_csv(header=True))
print("# threshold:",threshold_val)
print(tabulate(results.loc[results.index != "control-seq-lsppred",['LSP','LSP_probability']],headers='keys', tablefmt='psql'))

