import mne
import numpy as np
import matplotlib.pyplot as plt
from eegTools import Channel 

raw = mne.io.read_raw_edf("../data/tutorial_eeg.edf", preload=True)
#raw.plot()
raw.info['bads'] = ['EOG EOG1', 'EOG EOG2', 'EMG EMG', 'ECG ECG', 'STI 014'] 

#montage = mne.channels.read_montage('GSN-HydroCel-129')
#montage.ch_names[3:132] = raw.ch_names[:129]
#raw.plot_sensors()
#raw.plot_sensors('3d')
#picks = mne.pick_types(raw.info, meg=False, eeg=True, eog=False)

iir_params = dict(order=2, ftype='butter')
raw.filter(l_freq=0.3, h_freq=50., method='iir', iir_params=iir_params)
#raw.notch_filter(phase='zero', fir_window='hamming',filter_length='auto', freqs=np.arange(50, 100, 50), n_jobs=1)
#raw.plot()
raw.set_eeg_reference([]) #averaging using the default electrodes
#raw.plot()
#raw.filter(0.5,80)
dataArray = raw._data
channels, samples = dataArray.shape

# Isolating a sample channel
picked = 5 # channel index to pick
sampleCh=dataArray[picked,:]
print("PICKED CHANNEL: %s" %(raw.info["ch_names"][picked]))
c1 = Channel(sampleCh,raw.info['sfreq'])
c1.neoOp(30,'halftime')
