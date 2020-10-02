"""
test file
"""


class LinearKPFM():
    def __init__(self,h5_main):
        self.h5_main = h5_main


    def split_voltages(self,data_dict):
        import numpy as np
        t = np.array([ii for ii in range(data_dict['ndim_form'][0])]) * (1 / data_dict['scan_rate']) * 2
        y = np.linspace(0, data_dict['scan_size'], data_dict['ndim_form'][0])

        cnttot = 0
        count = 1
        cntmax = 0
        indx = []
        vs = []
        for ii in range(1, len(t)):
            if np.rint(data_dict['Voltage'][ii, 0]) != np.rint(data_dict['Voltage'][ii - 1, 0]):
                cnttot = cnttot + 1  # counts the total number of voltage switches in the data set
                if np.rint(data_dict['Voltage'][ii, 0]) != np.rint(data_dict['Voltage'][
                                                                       ii, -1]):  # checks to see if the voltage is changed in the middle fo the scan line
                    indx.append(ii + 1)
                else:
                    indx.append(ii)
                    vs.append(data_dict['Voltage'][ii, 0])
                if count > cntmax:
                    cntmax = count
                    count = 1
                else:
                    count = 1
            else:
                count = count + 1
        return vs, cntmax, cnttot, indx, t, y

    def unpack_data(self,h5_main,scan_rate,scan_size):

        """
          Not sure how to get the scan rate and scan size out of the igor file...

          scan_rate: in Hz
          scan_size: in um
        """
        import numpy as np

        scan_rate = scan_rate
        desr = [ii.data_descriptor for ii in h5_main]
        subs = ['Amplitude', 'Potential', 'UserIn', 'HeightRetrace', 'Phase']
        indx0 = [desr.index(list(filter(lambda x: ii in x, desr))[0]) for ii in subs]
        amp_data = h5_main[indx0[0]]
        pot_data = h5_main[indx0[1]]
        volt_data = h5_main[indx0[2]]
        height_data = h5_main[indx0[3]]
        phase_data = h5_main[indx0[4]]
        ndim_form = volt_data.get_n_dim_form().shape
        volt = (np.flipud(np.reshape(volt_data, ndim_form[:2])))
        pot = (np.flipud(np.reshape(pot_data, ndim_form[:2])))
        height = (np.flipud(np.reshape(height_data, ndim_form[:2])))
        phase = (np.flipud(np.reshape(phase_data, ndim_form[:2])))

        data_dict = {'Voltage': volt, 'CPD': pot, 'Height': height, 'Phase': phase, 'Amplitude': amp_data, 'ndim_form': ndim_form,
                     'scan_rate': scan_rate, 'scan_size': scan_size}

        vs, cntmax, cnttot, indx, t, y = self.split_voltages(data_dict)

        data_dict.update({'indx': indx, 'vs': vs, 'count max': cntmax,'count total':cnttot,'t':t,'y':y})
        return data_dict

    def compute_voltage_averages(self,data_dict):
        import numpy as np

        avgs = []
        zeroavg = np.mean(data_dict['CPD'][:data_dict['indx'][0] - 1, :], axis=0)
        for ii in range(len(data_dict['indx']) - 1):
            avgs.append(np.mean(data_dict['CPD'][data_dict['indx'][ii]:data_dict['indx'][ii + 1] + 1, :], axis=0) - zeroavg)

        data_dict.update({'zeroavg':zeroavg,'avgs':avgs})

    def plot_CPD_voltages(data_dict,method='Raw',window=13,poly=3):
        import matplotlib.pyplot as plt
        import numpy as np

        cmap = plt.cm.get_cmap('plasma', data_dict['count max'])
        jj = 0
        fig, axs = plt.subplots(nrows=int(data_dict['count total'] / 2 + 1), ncols=2, sharex='col', figsize=(15, 10))
        axs[0, 0].axis('off')
        axs[1, 0].set_title('Biasing')
        axs[0, 1].set_title('Zero Voltage After Bias')
        axs[0, 1].text(0.3, 0.8, '0 V', transform=axs[0, 1].transAxes)
        axs[int(data_dict['count total'] / 2), 0].set_xlabel('$\mu$m')
        axs[int(data_dict['count total'] / 2), 1].set_xlabel('$\mu$m')
        cbaxs = fig.add_axes([0.91, 0.095, 0.02, 0.71])
        fig.subplots_adjust(hspace=0.25, wspace=0.25)
        axs[0, 1].set_ylabel('CPD (V)', rotation=90, labelpad=2)
        axs[0, 1].axvspan(data_dict['y'][0], 0, facecolor='0.5', alpha=0.5)

        col = 1
        row = 0
        lab = 0
        for ii in range(len(data_dict['t'])):
            if np.rint(data_dict['Voltage'][ii, -1]) != np.rint(data_dict['Voltage'][ii, 0]):
                jj = jj + 1
                lab = 1
                continue
            if ii != 0:
                if lab == 1:
                    prev = np.rint(data_dict['Voltage'][ii - 2, 0])
                else:
                    prev = np.rint(data_dict['Voltage'][ii - 1, 0])
                if np.rint(data_dict['Voltage'][ii, 0]) != prev:
                    jj = 0
                    if np.rint(data_dict['Voltage'][ii, 0]) != 0:
                        col = 0
                        row = row + 1
                    else:
                        col = 1
                        zero_pot = data_dict['zeroavg']
                    S = np.array2string(np.rint(data_dict['Voltage'][ii, 0])) + ' V'
                    axs[row, col].text(0.3, 0.8, S, transform=axs[row, col].transAxes)
                    axs[row, col].set_ylabel('CPD (V)', rotation=90, labelpad=2)
                if ii == 0:
                    zero_pot = data_dict['zeropot']
            if method == 'Raw':
                yy = data_dict['CPD'][ii, :]
                axs[row, col].plot(data_dict['y'], yy, c=cmap(jj))
            elif method == 'Static_rm':
                yy = data_dict['CPD'][ii, :] - data_dict['zeroavg']
                axs[row, col].plot(data_dict['y'], yy, c=cmap(jj))
            elif method == 'Efield':
                smooth = si.savgol_filter(data_dict['CPD'][ii,:]-zero_pot,window,poly)
                yy = np.diff(smooth)/np.diff(data_dict['y'])*1e4
                axs[row, col].plot(data_dict['y'][:,-1],yy,c=cmap(jj))

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