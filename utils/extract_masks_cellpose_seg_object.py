import os
import sys
import numpy as np
from PIL import Image


def extractProbabilitiesCellposeSeg(filename, in_folder, out_folder):
    # loads a _seg.npy file created by cellpose and save the masks

    infile_path = os.path.join(in_folder, filename+'_seg.npy')
    outfile_path = os.path.join(out_folder, filename+'_cp_masks.tif')

    data = np.load(infile_path, allow_pickle=True).item()
    mask = np.squeeze(data['masks'])
    mask_imag = Image.fromarray(mask)
    mask_imag.save(outfile_path)

if __name__ == "__main__":
    filename = sys.argv[1]
    in_folder = sys.argv[2]
    out_folder = sys.argv[3]
    extractProbabilitiesCellposeSeg(filename, in_folder, out_folder)