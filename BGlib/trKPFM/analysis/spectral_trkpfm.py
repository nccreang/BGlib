"""
Analysis code for spectral tr-KPFM
"""

class SpectralTRKPFM():
    def __init__(self,h5_main):
        self.h5_main = h5_main

    def unload_attributes(self):
        import numpy as np
        from scipy.io import loadmat

        self.dc_amp_vec = np.squeeze(usid.hdf_utils.find_dataset(h5_main,'Spectroscopic_Values')[0][:])
        self.pnts_per_pix = self.dc_amp_vec.shape[1]
        self.h5_f = loadmat(parms_path)
        self.num_rows = int(h5_f['num_rows'])
        self.num_cols = int(h5_f['num_cols'])
        self.num_pix = self.num_cols*self.num_rows

        bias = h5_f['dc_amp_vec']
        self.bias_vec = bias[0][1::2]
        self.bias_positive = bias[0][1::2]
        self.bias_negative = np.flip(bias[0][3::4].T,0)
        self.num_dc_step = np.size(self.bias_vec)



