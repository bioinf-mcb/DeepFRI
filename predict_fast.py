import tempfile
import numpy as np

from deepfrier.ONNXPredictor import Predictor
from deepfrier.fold import fold_protein

from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO

import argparse

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

from logging import getLogger
import obonet
import networkx

graph = obonet.read_obo("go-basic.obo")

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

def explode_goterm(goterm):
    # get name of the GO term
    return list(networkx.predecessor(graph, goterm).keys())

def explode_all_goterms(goterms):
    res = {}
    for goterm, score in goterms:
        for parent in explode_goterm(goterm):
            if parent not in res:
                res[parent] = score
            else:
                res[parent] = max(res[parent], score)
    return list(res.items())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Predict contact maps from PDB files')

    # add two arguments - `inputs` and `inputs_file`, where `inputs_file` would be a file with paths to PDB files. One of the arguments is required
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input', type=str, help='Input PDB file')
    group.add_argument('-f', '--input_file', type=str, help='Input file with paths to PDB files')
    group.add_argument('-s', '--seq', type=str, help='Input sequence. It will be folded using ESMFold (requires internet access)')
    group.add_argument('--fasta_fn', type=str, help='Input fasta file')

    parser.add_argument('-m', '--models', nargs='+', type=str, required=True, help='Models to use')
    parser.add_argument('-p', '--propagate', action='store_true', help='Propagate predictions to neighboring residues')
    
    args = parser.parse_args()

    assert len(args.models) > 0, "No models specified"
    assert all([m in ['cc', 'mf', 'bp'] for m in args.models]), "Invalid model name"  

    # load PDB file
    model_template = "onnx_models/DeepFRI-UNIPROT_GraphConv_gcd_512-512-512_fcd_1024_ca_10.0_ext_desc_{0}_cpu"

    models = {
        "cc": Predictor(model_template.format('cc'), gcn=True, threads=1),
        "mf": Predictor(model_template.format('mf'), gcn=True, threads=1),
        "bp": Predictor(model_template.format('bp'), gcn=True, threads=1),
    }
    
    fnames = []

    if args.input:
        fnames = [args.input]
    elif args.seq:
        fname = fold_protein(args.seq)
        fnames = [fname]
    elif args.input_file:
        fnames = [line.strip() for line in open(args.input_file)]
    elif args.fasta_fn:
        fnames = []
        for record in SeqIO.parse(args.fasta_fn, "fasta"):
            fname = fold_protein(str(record.seq))
            fnames.append(fname)
    for fname in fnames:
        try:
            c_alpha, seq = load_predicted_PDB(fname)
        except KeyError:
            print("Failed to load PDB file "+fname)
            continue
        cmap = (c_alpha < THRESHOLD).astype(np.int32)

        for model in args.models:
            predictor = models[model]
            predictor.predict_function(seq, cmap, fname)
            predicted_goterms = [(i[0], float(i[2])) for i in predictor.prot2goterms[fname]]
            original_goterms = [i[0] for i in predicted_goterms]
            predictor.flush()

            if args.propagate:
                predicted_goterms = explode_all_goterms(predicted_goterms)

            print("model,go_id,score,go_name,propagated,depth")
            for goterm, score in predicted_goterms:
                goname = graph.nodes[goterm]['name']
                propagated = "F"
                if goterm not in original_goterms:
                    propagated = "T"
                num_parents = len(explode_goterm(goterm))
                print(model+","+goterm+","+str(score)+","+goname+","+propagated+","+str(num_parents))