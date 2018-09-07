
# coding: utf-8

# In[1]:


import mne
import numpy as np
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
from eegTools import Channel 


# In[2]:


raw = mne.io.read_raw_edf("../data/tutorial_eeg.edf", preload=True)
#raw.info['bads'] = ['EOG EOG1', 'EOG EOG2', 'EMG EMG', 'ECG ECG', 'STI 014']
bads  = ['EOG EOG1', 'EOG EOG2', 'EMG EMG', 'ECG ECG', 'STI 014']


# In[3]:


iir_params = dict(order=2, ftype='butter')
raw.filter(l_freq=0.3, h_freq=50., method='iir', iir_params=iir_params)
raw.set_eeg_reference([])


# In[4]:


dataArray = raw._data
channels, samples = dataArray.shape


# In[5]:


c=[0]*6
picked = 0 # channel index to pick
for picked in range(0,6):
    sampleCh=dataArray[picked,:]
    print("PICKED CHANNEL: %s" %(raw.info["ch_names"][picked]))
    c[picked]= Channel(sampleCh,raw.info['sfreq'],raw.info["ch_names"][picked])
    c[picked].neoOp(50,'halftime')


# In[6]:


for picked in range(0,6):
    c[picked].neoOp(70,'halftime')   


# In[8]:


c[0].plotSpyke()

