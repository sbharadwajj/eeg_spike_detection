import numpy as np
import matplotlib.pyplot as plt
import random
from random import randint
from IPython.display import clear_output
import matplotlib as mpl
from scipy.signal import argrelextrema

###################	class, representing a single EEG channel	###################

class Channel: 

	def __init__(self,ch,sampFreq,name):		#initializing data_members

		self.data = ch
		self.samples = len(ch)
		self.sampFreq = sampFreq
		self.setVariables()
		self.name=name

	def setVariables(self):
		ch = self.data
		sampFreq = self.sampFreq
		self.duration = len(ch)/sampFreq
		self.epochLen = int(np.ceil(0.35*sampFreq)) 
		#Assuming the spyke period to be 20ms to 70ms and background data of 4*max(spyke_period)

	
	def operations(self):
                ch = self.data
                sampFreq = self.sampFreq
                self.lastMarking = 0
                self.featureVec = ['0']*len(ch)
                self.der1 = np.gradient(ch)
                self.der2 = np.gradient(self.der1)
                self.spykeMark = [0]*self.samples
                self.spykes=0
                self.visited=[0]*self.samples
	
                d1=self.der1
                d2=self.der2
                l=self.samples 
                lastState=0
                nowState=0
                for i in argrelextrema(ch,np.greater)[0]:
                    self.featureVec[i] = 'M'
                for i in argrelextrema(ch,np.less)[0]:
                    self.featureVec[i] = 'm'

	#	for i in range(0,l):
	#		if d1[i] < 0:
	#			nowState=-1
	#		else:
	#			nowState=1
	#		if i==0:
	#			lastState = nowState
	#		else:
	#			if nowState != lastState:			###### state transition wrt monotonicity ######

	#				if d2[i] < 0:				###### type of edge (-ve raising or +ve falling)
	#					self.featureVec[i] = 'M'
	#				else:
	#					self.featureVec[i] = 'm'
	#				lastState = nowState

	def plotSampleEpoch(self):
		data=self.data
		dur=self.epochLen/2
		maximas=0
		for i in range(0,self.samples):
			if self.featureVec[i]=='M':
				maximas = maximas+1
		print(maximas)
		n=randint(1,maximas+1)
		cnt=0
		for i in range(0,self.samples): ################### to stop at nth Maxima ###################
			if cnt == n:
				plt.plot(data[int(i-dur):int(i+dur)])
				plt.show()
				break	
			if self.featureVec[i]=='M':
				cnt = cnt+1

	def plotEpoch(self,inst):
		data=self.data
		dur=self.epochLen/2
		plt.plot(data[int(inst-dur):int(inst+dur)])
		plt.show()

	def neoOp(self,const,method='derivative'):

		self.resetAll()
		l=self.samples
		self.neo = [0.00]*l
		x = self.data

###################		    Algorithm to calculate NEO			###################
###################		NEO(x[n]) = (x[n]^2)-(x[n-1]*x[n+1])		###################

		for i in range (0,l):
			if i == 0:
				self.neo[i] = 0

			elif i==l-1:
				self.neo[i] = 0
			else:
				self.neo[i] = (np.square(x[i]))-(x[i-1]*x[i+1])

###################			Tresholding				###################

		neo=self.neo
		self.mark = [0]*l
		self.tresh = (const/l)*(sum(self.neo))

		for i in range(0,l):
			if neo[i] > self.tresh:
				self.mark[i] =1
		print("Executing getspykes()")
		self.getSpykes(method)

	def getSpykes(self,method):

		#### Get a continuous length of boolean 1 ####

		l=self.samples
		cont = 0
		for i in range (0,l):
			if self.mark[i] == 1 and cont ==0:
				cont = 1
				beg = i
			elif self.mark[i] == 0 and cont ==1:
				cont = 0
				end = i
				if method == 'derivative':
					self.checkSpykeDer(beg,end)
				elif method == 'halftime':
					self.checkSpykeHft(beg,end)
		c = 0
		for i in range(0,l):
			if self.spykeMark[i] == 1:
				c = c+1
		self.spykes=c
		print("Total number of suspected spykes = %d" %c)
	####   Iterate over (begin-5) to (end+5)  ####

	def checkSpykeDer(self,beg,end):

		beg= beg-5
		end = end+5
		#print(self.featureVec[beg:end])

		for i in range(beg,end):
			if self.featureVec[i] == 'M':
				#print(i)
				self.spykeMark[i] = 1

	def checkSpykeHft(self,beg,end):

		self.spykeMark[int((beg+end)/2)] = 1

	def plotSpyke(self,dur=15):
		dur=dur*self.sampFreq
		dr=self.epochLen/2
		l=self.samples
		data=self.data
		spykes = self.spykes
		n=randint(1,spykes+1)
		cnt = 0
		for i in range(int(dur+1),l):
			if self.spykeMark[i] == 1:
				cnt = cnt+1
			if cnt == n:
				plt.plot(data[int(i-dur):int(i+dur)])
				plt.axhline(y=0, xmin=100-int(dr), xmax=100+int(dr), color='red', zorder=1)
				plt.vlines(x=dur, ymin=-0.00004, ymax=data[i], color='red', zorder=2)
				plt.show()
				break

	def resetAll(self):
		self.setVariables()
		

def formDataset(channels,cno,i,dur=15):
	dur=dur*channels[0].sampFreq
	dr=channels[0].epochLen/2
	data=channels[cno].data
	for chno in range(0,3):
		data=channels[chno].data
		plt.subplot(len(channels),1,chno+1)
		axes = plt.gca()
		axes.set_ylim([-0.00005,0.00005]) 
		t = "Plot for channel %d" %(cno)
		plt.rc('figure', figsize=(20.0, 5.0))
		plt.title(t)
		plt.plot(data[int(i-dur):int(i+dur)],linewidth=0.5)
		#plt.axhline(y=0, xmin=100-int(dr), xmax=100+int(dr), color='red', zorder=1)
		#plt.vlines(x=dur, ymin=-0.00004, ymax=data[i], color='red', zorder=2)
		#plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
		plt.show()
	resp=input("Is this a spike in channel %d? (y/n): " %(cno))
	data=channels[cno].data
	if resp == 'y':
		clear_output()
		return data[int(i-dr):int(i+dr)].tolist(),[1]

	elif resp == 'n':
		clear_output()
		return data[int(i-dr):int(i+dr)].tolist(),[0]

def getMontage():
        time=np.linspace(25,70)
        offset=[0.652,0.761,0.473,0.3532,0.6923,0.355]
        amplitude=np.linspace(35,70)
        sign=[-1,1]
        period=random.choice(time)
        os=random.choice(offset)
        amp=random.choice(amplitude)
        s=random.choice(sign)
        sinsamples=int(256*period*0.001)
        x = np.linspace(-np.pi, np.pi, sinsamples)
        k=np.sin(x)
        k +=os
        noise = np.random.normal(0,1,sinsamples)
        k = np.add(k,0.5*noise)
        k *=amp*0.0000005
        k *=s
        u = np.argmax(k)
        return k,u

def importMontage(channels):
	marks = [0]*200
	length = len(channels[0].data)
	chs = len(channels)
	for m in range (0,200):
		placer=randint(4000,(length-100))
		spk,ofst = getMontage()
		marks[m]=placer+ofst
		for rng in range(0,len(spk)):
			channels[0].data[placer+rng] += spk[rng]
			channels[1].data[placer+rng] += 0.65*spk[rng]
			channels[2].data[placer+rng] += 0.5*spk[rng]
			channels[3].data[placer+rng] += 0.35*spk[rng]
	return channels,marks
