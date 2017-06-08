#import fileinput
import glob
import ROOT as R
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
import math as m
R.gROOT.SetBatch(1)


pileup = str(sys.argv[1])

samplefiles = '/user_data/ktaylor/GridTree/JetGridTrees_'+pileup+'PU/*'
output='/user_data/ktaylor/test/images_SDmass_'+pileup+'PU'

entries=[]


def plot_mean_jet(rec, field = 'images', title = 'Average Jet Image'):
	fig = plt.figure(figsize=(8, 8), dpi=100)
	ax = fig.add_subplot(111)
	im = ax.imshow(np.mean(rec, axis = 0),  norm=LogNorm(vmin=0.00001, vmax=1), interpolation='nearest')
	plt.title(r''+title)
	return fig


def plot_jet(rec, title = 'Jet Image', log=True):
	fig = plt.figure(figsize=(8, 8), dpi=100)
	ax = fig.add_subplot(111)
	if log:
		im = ax.imshow(rec,  norm=LogNorm(vmin=0.00001, vmax=1), interpolation='nearest')
	else:
		im = ax.imshow(rec, interpolation='nearest')
	plt.title(r''+title)
	return fig



jet_pt=[]
jet_eta=[]
jet_phi=[]
jet_mass=[]
jet_sdmass=[]
tau21=[]
Sublead_eta=[]
Sublead_phi=[]
PCeta=[]
PCphi=[]
signal=[]

Delroot=glob.glob(samplefiles)
#print Delroot
total_count=0
for files in Delroot:
	print files
	tFile = R.TFile.Open(files, "read")
	tTree = tFile.Get("jetgrids")
	count=0
	for event in tTree:
		if (event.JetAK8_SDMASS < 65 or event.JetAK8_SDMASS > 95):
			continue
		if (event.JetAK8_PT > 325):
			continue
		total_count+=1
		count+=1
		
		Jetarray=np.array(np.zeros((25,25)))
		
		etacenter=event.JetAK8_LsubjetEta
		phicenter=event.JetAK8_LsubjetPhi
		
		Sublead_eta.append(np.float32(event.JetAK8_SubLeadingEta))
		Sublead_phi.append(np.float32(event.JetAK8_SubLeadingPhi))
		PCeta.append(np.float32(event.JetAK8_PCEta))
		PCphi.append(np.float32(event.JetAK8_PCPhi))
		jet_pt.append(np.float32(event.JetAK8_PT))
		jet_eta.append(np.float32(etacenter))
		jet_phi.append(np.float32(phicenter))
		jet_mass.append(np.float32(event.JetAK8_MASS))
		jet_sdmass.append(np.float32(event.JetAK8_SDMASS))
		tau21.append(np.float32(event.JetAK8_Tau21))
		if 'Wprime' in files : 
			signal.append(1)
		else:
			signal.append(0)
		for calo in range(len(event.CaloTower_ET)):
			etadistance=int(m.floor((event.CaloTower_ETA[calo]-etacenter)/.0714 +.5))
			phidistance= int(m.floor((event.CaloTower_PHI[calo]-phicenter)/(m.pi/44)+.5))
			if abs(etadistance)<=12 and abs(phidistance)<=12:
				Jetarray[12+etadistance][12+phidistance]+=event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo])
		
		entries.append(Jetarray)
 
	print 'Events: ',count
	print 'Total Events: ', total_count
            
#plot_mean_jet(entries).savefig(output+'.pdf')
#plot_jet(entries[0]).savefig('Wprime_few.pdf')
np.savez(output,image=entries,signal=signal,PCEta=PCeta,PCPhi=PCphi,SubLeadingEta=Sublead_eta,SubLeadingPhi=Sublead_phi,jet_pt=jet_pt,jet_eta=jet_eta,jet_phi=jet_phi,jet_mass=jet_mass,tau21=tau21)


#                
                       

