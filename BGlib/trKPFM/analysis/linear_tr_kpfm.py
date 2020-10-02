"""
test file
"""


class LinearKPFM():
    def __init__(self,h5_main,scan_rate,scan_size):
        self.h5_main = h5_main
        self.scan_rate = scan_rate
        self.scan_size = scan_size


    def split_voltages(self):
        import numpy as np
        t = np.array([ii for ii in range(self.ndim_form[0])]) * (1 / self.scan_rate) * 2
        y = np.linspace(0, self.scan_size, self.ndim_form[0])

        cnttot = 0
        count = 1
        cntmax = 0
        indx = []
        vs = []
        for ii in range(1, len(t)):
            if np.rint(self.volt[ii, 0]) != np.rint(self.volt[ii - 1, 0]):
                cnttot = cnttot + 1  # counts the total number of voltage switches in the data set
                if np.rint(self.volt[ii, 0]) != np.rint(self.volt[ii, -1]):  # checks to see if the voltage is changed in the middle fo the scan line
                    indx.append(ii + 1)
                else:
                    indx.append(ii)
                    vs.append(self.volt[ii, 0])
                if count > cntmax:
                    cntmax = count
                    count = 1
                else:
                    count = 1
            else:
                count = count + 1
        self.vs = vs
        self.cntmax = cntmax
        self.cnttot = cnttot
        self.indx = indx
        self.y = y
        self.t = t
        # return vs, cntmax, cnttot, indx, t, y

    def unpack_data(self):

        """
          Not sure how to get the scan rate and scan size out of the igor file...

          scan_rate: in Hz
          scan_size: in um
        """
        import numpy as np

        desr = [ii.data_descriptor for ii in self.h5_main]
        subs = ['Amplitude', 'Potential', 'UserIn', 'HeightRetrace', 'Phase']
        indx0 = [desr.index(list(filter(lambda x: ii in x, desr))[0]) for ii in subs]
        amp_data = self.h5_main[indx0[0]]
        pot_data = self.h5_main[indx0[1]]
        volt_data = self.h5_main[indx0[2]]
        height_data = self.h5_main[indx0[3]]
        phase_data = self.h5_main[indx0[4]]
        self.ndim_form = volt_data.get_n_dim_form().shape
        self.volt = (np.flipud(np.reshape(volt_data, self.ndim_form[:2])))
        self.pot = (np.flipud(np.reshape(pot_data, self.ndim_form[:2])))
        self.height = (np.flipud(np.reshape(height_data, self.ndim_form[:2])))
        self.phase = (np.flipud(np.reshape(phase_data, self.ndim_form[:2])))
        self.amp = (np.flipud(np.reshape(amp_data,self.ndim_form[:2])))


        # data_dict = {'Voltage': volt, 'CPD': pot, 'Height': height, 'Phase': phase, 'Amplitude': amp_data, 'ndim_form': ndim_form,
        #              'scan_rate': self.scan_rate, 'scan_size': self.scan_size}

        # vs, cntmax, cnttot, indx, t, y = self.split_voltages(data_dict)
        self.split_voltages()
        # data_dict.update({'indx': indx, 'vs': vs, 'count max': cntmax,'count total':cnttot,'t':t,'y':y})
        # return data_dict

    def compute_voltage_averages(self):
        import numpy as np

        avgs = []
        zeroavg = np.mean(self.pot[:self.indx[0] - 1, :], axis=0)
        for ii in range(len(self.indx) - 1):
            avgs.append(np.mean(self.pot[self.indx[ii]:self.indx[ii + 1] + 1, :], axis=0) - zeroavg)

        # data_dict.update({'zeroavg':zeroavg,'avgs':avgs})
        self.zeroavg = zeroavg
        self.avgs = self.avgs


    def plot_CPD_voltages(self,method='Raw',window=13,poly=3):
        import matplotlib.pyplot as plt
        import numpy as np

        cmap = plt.cm.get_cmap('plasma', self.cntmax)
        jj = 0
        fig, axs = plt.subplots(nrows=int(self.cnttot / 2 + 1), ncols=2, sharex='col', figsize=(15, 10))
        axs[0, 0].axis('off')
        axs[1, 0].set_title('Biasing')
        axs[0, 1].set_title('Zero Voltage After Bias')
        axs[0, 1].text(0.3, 0.8, '0 V', transform=axs[0, 1].transAxes)
        axs[int(self.cnttot / 2), 0].set_xlabel('$\mu$m')
        axs[int(self.cnttot / 2), 1].set_xlabel('$\mu$m')
        cbaxs = fig.add_axes([0.91, 0.095, 0.02, 0.71])
        fig.subplots_adjust(hspace=0.25, wspace=0.25)
        axs[0, 1].set_ylabel('CPD (V)', rotation=90, labelpad=2)
        axs[0, 1].axvspan(self.y[0], 0, facecolor='0.5', alpha=0.5)

        col = 1
        row = 0
        lab = 0
        for ii in range(len(self.t)):
            if np.rint(self.volt[ii, -1]) != np.rint(self.volt[ii, 0]):
                jj = jj + 1
                lab = 1
                continue
            if ii != 0:
                if lab == 1:
                    prev = np.rint(self.volt[ii - 2, 0])
                else:
                    prev = np.rint(self.volt[ii - 1, 0])
                if np.rint(self.volt[ii, 0]) != prev:
                    jj = 0
                    if np.rint(self.volt[ii, 0]) != 0:
                        col = 0
                        row = row + 1
                    else:
                        col = 1
                        zero_pot = self.zeroavg
                    S = np.array2string(np.rint(self.volt[ii, 0])) + ' V'
                    axs[row, col].text(0.3, 0.8, S, transform=axs[row, col].transAxes)
                    axs[row, col].set_ylabel('CPD (V)', rotation=90, labelpad=2)
                if ii == 0:
                    zero_pot = self.zeroavg
            if method == 'Raw':
                yy = self.pot[ii, :]
                axs[row, col].plot(data_dict['y'], yy, c=cmap(jj))
            elif method == 'Static_rm':
                yy = self.pot[ii, :] - self.zeroavg
                axs[row, col].plot(self.y, yy, c=cmap(jj))
            elif method == 'Efield':
                smooth = si.savgol_filter(self.pot[ii,:]-zero_pot,window,poly)
                yy = np.diff(smooth)/np.diff(data_dict['y'])*1e4
                axs[row, col].plot(self.y[:,-1],yy,c=cmap(jj))

            lab = 0
            jj = jj + 1
        scbar = plt.cm.ScalarMappable(cmap=plt.cm.plasma,
                                      norm=plt.Normalize(vmin=0, vmax=data_dict['t'][data_dict['count max']]))
        scbar._A = []
        cbar = plt.colorbar(scbar, cax=cbaxs)
        cbar.ax.set_ylabel('Relaxation Time (s)', rotation=270, labelpad=15, size=12)
        cbar.ax.tick_params(labelsize=12)
        fig.align_ylabels(axs)
        fig.subplots_adjust(wspace=.3)
        return fig, axs