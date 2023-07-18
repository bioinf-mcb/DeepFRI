import tempfile
import numpy as np

from deepfrier.ONNXPredictor import Predictor

from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO

import argparse

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

from logging import getLogger

logger = getLogger(__name__)
THRESHOLD = 10.0

def load_predicted_PDB(pdbfile):
    # Generate (diagonalized) C_alpha distance matrix from a pdbfile
    parser = PDBParser()
    structure = parser.get_structure(pdbfile.split('/')[-1].split('.')[0], pdbfile)
    residues = [r for r in structure.get_residues()]

    # sequence from atom lines
    records = SeqIO.parse(pdbfile, 'pdb-atom')
    seqs = [str(r.seq) for r in records]

    distances = np.empty((len(residues), len(residues)))
    for x in range(len(residues)):
        for y in range(len(residues)):
            one = residues[x]["CA"].get_coord()
            two = residues[y]["CA"].get_coord()
            distances[x, y] = np.linalg.norm(one-two)

    return distances, seqs[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Predict contact maps from PDB files')

    # add two arguments - `inputs` and `inputs_file`, where `inputs_file` would be a file with paths to PDB files. One of the arguments is required
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input', type=str, help='Input PDB file')
    group.add_argument('-f', '--input_file', type=str, help='Input file with paths to PDB files')

    parser.add_argument('-m', '--models', nargs='+', type=str, required=True, help='Models to use')
    args = parser.parse_args()

    assert len(args.models) > 0, "No models specified"
    assert all([m in ['cc', 'mf', 'bp'] for m in args.models]), "Invalid model name"  

    # load PDB file
    model_template = "trained_models/DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_{0}_cpu"

    models = {
        "cc": Predictor(model_template.format('cc'), gcn=True, threads=1),
        "mf": Predictor(model_template.format('mf'), gcn=True, threads=1),
        "bp": Predictor(model_template.format('bp'), gcn=True, threads=1),
    }
    
    if args.input:
        fnames = [args.input]
    else:
        fnames = [line.strip() for line in open(args.input_file)]
    for fname in fnames:
        c_alpha, seq = load_predicted_PDB(fname)
        cmap = (c_alpha < THRESHOLD).astype(np.int32)

        for model in args.models:
            predictor = models[model]
            predictor.predict_function(seq, cmap, fname)
            
            with tempfile.NamedTemporaryFile() as f:
                predictor.export_csv(f.name)
                for idx, line in enumerate(open(f.name)):
                    if idx == 0:
                        continue
                    print(model+","+line.strip())