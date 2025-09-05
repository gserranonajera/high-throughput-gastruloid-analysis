import os
import sys
import numpy as np
from PIL import Image


def extractProbabilitiesCellposeSeg(filename, in_folder, out_folder):
    # loads a _seg.npy file created by cellpose and save the probabilities (0-1) as an 8bit image (0-255)

    infile_path = os.path.join(in_folder, filename+'_seg.npy')
    outfile_path = os.path.join(out_folder, 'prob_'+filename+'.tif')

    data = np.load(infile_path, allow_pickle=True).item()
    probs = np.squeeze(data['flows'][1]) #it has 3 dimensions. the first one is a singleton
    prob_imag = Image.fromarray(probs)
    prob_imag.save(outfile_path)

if __name__ == "__main__":
    filename = sys.argv[1]
    in_folder = sys.argv[2]
    out_folder = sys.argv[3]
    extractProbabilitiesCellposeSeg(filename, in_folder, out_folder)