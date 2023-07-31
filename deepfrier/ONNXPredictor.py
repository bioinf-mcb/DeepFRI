import os
import csv
import glob
import json
import gzip
import secrets

import numpy as np
import tensorflow as tf

from .utils import load_catalogue, load_FASTA, load_predicted_PDB, seq2onehot
from .layers import MultiGraphConv, GraphConv, FuncPredictor, SumPooling

import onnxruntime as rt

class GradCAM(object):
    """
    GradCAM for protein sequences.
    [Adjusted for GCNs based on https://arxiv.org/abs/1610.02391]
    """
    def __init__(self, model, layer_name="GCNN_concatenate"):
        self.grad_model = tf.keras.models.Model([model.inputs], [model.get_layer(layer_name).output, model.output])

    def _get_gradients_and_filters(self, inputs, class_idx, use_guided_grads=False):
        with tf.GradientTape() as tape:
            conv_outputs, predictions = self.grad_model(inputs)
            loss = predictions[:, class_idx, 0]
        grads = tape.gradient(loss, conv_outputs)

        if use_guided_grads:
            grads = tf.cast(conv_outputs > 0, "float32")*tf.cast(grads > 0, "float32")*grads

        return conv_outputs, grads

    def _compute_cam(self, output, grad):
        weights = tf.reduce_mean(grad, axis=1)
        # perform weighted sum
        cam = tf.reduce_sum(tf.multiply(weights, output), axis=-1).numpy()

        return cam

    def heatmap(self, inputs, class_idx, use_guided_grads=False):
        output, grad = self._get_gradients_and_filters(inputs, class_idx, use_guided_grads=use_guided_grads)
        cam = self._compute_cam(output, grad)
        heatmap = (cam - cam.min())/(cam.max() - cam.min())

        return heatmap.reshape(-1)


class Predictor(object):
    """
    Class for loading trained models and computing GO/EC predictions and class activation maps (CAMs).
    """
    def __init__(self, model_prefix, gcn=True, threads=1):
        self.model_prefix = model_prefix
        self.gcn = gcn

        if rt.get_device() == 'CPU':
            self.threads = threads
        elif rt.get_device() == 'GPU':
            self.threads = 0

        self._load_model()
        self.prot2goterms = {}
        self.goidx2chains = {}

    def _load_model(self):
        session_options = rt.SessionOptions()
        session_options.intra_op_num_threads = self.threads
        self.session = rt.InferenceSession(
            self.model_prefix + '.onnx',
            providers=['CUDAExecutionProvider', 'CPUExecutionProvider'],
            sess_options=session_options,
        )

        # load parameters
        with open(self.model_prefix + "_model_params.json") as json_file:
            metadata = json.load(json_file)

        self.gonames = np.asarray(metadata['gonames'])
        self.goterms = np.asarray(metadata['goterms'])
        self.thresh = 0.1*np.ones(len(self.goterms))

    def _load_cmap(self, filename, cmap_thresh=10.0):
        if filename.endswith('.pdb'):
            D, seq = load_predicted_PDB(filename)
            A = np.double(D < cmap_thresh)
        elif filename.endswith('.npz'):
            cmap = np.load(filename)
            if 'C_alpha' not in cmap:
                raise ValueError("C_alpha not in *.npz dict.")
            D = cmap['C_alpha']
            A = np.double(D < cmap_thresh)
            seq = str(cmap['seqres'])
        elif filename.endswith('.pdb.gz'):
            rnd_fn = "".join([secrets.token_hex(10), '.pdb'])
            with gzip.open(filename, 'rb') as f, open(rnd_fn, 'w') as out:
                out.write(f.read().decode())
            D, seq = load_predicted_PDB(rnd_fn)
            A = np.double(D < cmap_thresh)
            os.remove(rnd_fn)
        else:
            raise ValueError("File must be given in *.npz or *.pdb format.")
        # ##
        S = seq2onehot(seq)
        S = S.reshape(1, *S.shape)
        A = A.reshape(1, *A.shape)

        return A, S, seq

    def predict_function(
        self,
        seqres: str,
        cmap = None,
        chain: str = "",
    ):
        """
        Computes GO/EC predictions for a single protein chain from sequence and contact map.

        Args:
            seqres (str): protein sequence
            cmap (np.array): contact map
            chain (str): chain ID

        Returns:
            None
        """

        self.Y_hat = np.zeros((1, len(self.goterms)), dtype=float)
        self.data = {}
        self.test_prot_list = [chain]
        S = seq2onehot(seqres)
        S = S.reshape(1, *S.shape)
        inputDetails = self.session.get_inputs()
        if cmap is not None and self.gcn:
            A = cmap.reshape(1, *cmap.shape)
            prediction = self.session.run(
                None, {
                    inputDetails[0].name: A.astype(np.float32),
                    inputDetails[1].name: S.astype(np.float32)
                })[0]
            self.data[chain] = [[A, S], seqres]
        else:
            prediction = self.session.run(
                None, {inputDetails[0].name: S.astype(np.float32)})[0]
            self.data[chain] = [[S], seqres]

        y = prediction[:, :, 0].reshape(-1)
        self.Y_hat[0] = y
        self.prot2goterms[chain] = []
        go_idx = np.where(y >= self.thresh)[0]
        for idx in go_idx:
            if idx not in self.goidx2chains:
                self.goidx2chains[idx] = set()
            self.goidx2chains[idx].add(chain)
            self.prot2goterms[chain].append(
                (self.goterms[idx], self.gonames[idx], float(y[idx])))
            

    def save_predictions(self, output_fn):
        print ("### Saving predictions to *.json file...")
        # pickle.dump({'pdb_chains': self.test_prot_list, 'Y_hat': self.Y_hat, 'goterms': self.goterms, 'gonames': self.gonames}, open(output_fn, 'wb'))
        with open(output_fn, 'w') as fw:
            out_data = {'pdb_chains': self.test_prot_list,
                        'Y_hat': self.Y_hat.tolist(),
                        'goterms': self.goterms.tolist(),
                        'gonames': self.gonames.tolist()}
            json.dump(out_data, fw, indent=1)

    def save_file(self, output_fn: str, delimiter: str, quotechar: str = '"'):
        """
        Exports predictions to .csv or .tsv format

        Args:
            output_fn (str): output file name
            delimiter (str): delimiter for .csv or .tsv file
            quotechar (str): quotechar for .csv file

        Returns:
            None
        """

        with open(output_fn, 'w', encoding="utf-8") as f:
            writer = csv.writer(f, delimiter=delimiter, quotechar=quotechar)
            writer.writerow([
                'Protein', 'GO_term/EC_number', 'Score',
                'GO_term/EC_number name'
            ])
            for prot, goterms in self.prot2goterms.items():
                sorted_rows = sorted(goterms, key=lambda x: x[2], reverse=True)
                for row in sorted_rows:
                    writer.writerow(
                        [prot, row[0], '{:.5f}'.format(row[2]), row[1]])
        self.flush()

    def export_csv(self, output_fn: str):
        self.save_file(output_fn, ",")

    def export_tsv(self, output_fn: str):
        self.save_file(output_fn, "\t")

    def flush(self):
        self.prot2goterms = {}
        self.goidx2chains = {}