import requests
import tempfile

from logging import getLogger

logger = getLogger(__name__)

def fold_protein(sequence, fname=""):
    """
    Fold a protein sequence using ESMFOLD.
    """

    logger.debug("Folding protein sequence")

    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    r = requests.post(url, data=sequence)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", prefix=fname+"_") as fp:
        fp.write(r.content)
        fp.flush()
        return fp.name