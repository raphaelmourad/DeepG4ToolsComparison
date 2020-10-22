import numpy as np
from Bio import SeqIO
def sequence_to_ohe(
        sequence,
        channel={
            'A': 0,
            'T': 1,
            'U': 1,
            'C': 2,
            'G': 3
        }
):
    """
  fun builds the one hot encoding numpy array for a sequence.
  :param sequence
  :param channel: the coding of nucleotides.
  """

    sequence_size = len(sequence)
    ohe_dataset = np.zeros((1, sequence_size, 4))
    for pos, nucleotide in enumerate(sequence):
        nucleotide = nucleotide.upper()
        if nucleotide == 'N':
            continue
        try:
            ohe_dataset[0, pos, channel[nucleotide]] = 1
        except KeyError:
            print()
            print("Unknown nucleotide", nucleotide)
            print("Counting as N.")
            continue
    return ohe_dataset
