#import fileinput
import glob
import ROOT as R
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import math as m
from jettools import plot_mean_jet,flip_jet,rotate_jet
R.gROOT.SetBatch(1)

#ljmetDir = '/mnt/hadoop/users/nelson/TreeGrids/'
samplefiles = '../WQCDGrids_35PU/*'
output='Jetimages_35PU'

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
tau21=[]
#Sublead_eta=[]
#Sublead_phi=[]
#PCeta=[]
#PCphi=[]
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
        total_count+=1
        count+=1

        Jetarray=np.array(np.zeros((25,25)))
        etacenter=event.JetAK8_LsubjetEta
        phicenter=event.JetAK8_LsubjetPhi
        
            #print 'jet#:',ijet,'etacenter:',event.JetAK8_ETA[ijet],'phicenter',event.JetAK8_PHI[ijet]'
        for calo in range(len(event.CaloTower_ET)):
            etadistance=int(abs(event.CaloTower_ETA[calo]-etacenter)/.0714+.5)
            phidistance= int(abs(event.CaloTower_PHI[calo]-phicenter)/(m.pi/44)+.5)

            if (etacenter >= event.CaloTower_ETA[calo] and phicenter >= event.CaloTower_PHI[calo]) and (etadistance<=12 and phidistance<=12):
                Jetarray[12-etadistance][12-phidistance]=event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo])
                            #print event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo]-etacenter), event.CaloTower_ETA[calo]-etacenter
                            #print event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo]-etacenter), np.cosh(event.CaloTower_ETA[calo]-etacenter)
            elif (etacenter >= event.CaloTower_ETA[calo] and phicenter < event.CaloTower_PHI[calo]) and (etadistance<=12 and phidistance<=12):
                Jetarray[12-etadistance][12+phidistance]=event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo])
                           #continue
            elif (etacenter < event.CaloTower_ETA[calo] and phicenter >= event.CaloTower_PHI[calo]) and (etadistance<=12 and phidistance<=12):
                Jetarray[12+etadistance][12-phidistance]=event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo])
                             #print etacenter, phicenter ,event.CaloTower_ETA[calo], event.CaloTower_PHI[calo]
                             #print etadistance, phidistance
                           #continue
            elif (etacenter < event.CaloTower_ETA[calo] and phicenter < event.CaloTower_PHI[calo]) and (etadistance<=12 and phidistance<=12):
                Jetarray[12+etadistance][12+phidistance]=event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo])
                       #print event.CaloTower_ET[calo]/np.cosh(event.CaloTower_ETA[calo]-etadistance)
                       #print 'LowerBinETA:', event.CaloEdgeEtaMin[calo], ' HigherBinEta:',event.CaloEdgeEtaMax[calo], ' LowerBinPhi:',event.CaloEdgePhiMin[calo], ' CaloET:',event.CaloTower_ET[calo]
                           #continue
                
        if (event.JetAK8_SubLeadingEta == -99) | (event.JetAK8_SubLeadingPhi ==-99):
            e, p = (event.JetAK8_PCEta, event.JetAK8_PCPhi)
        else:
            e, p = (event.JetAK8_SubLeadingEta, event.JetAK8_SubLeadingPhi)
            
        if e==0:continue
    
        angle = np.arctan(p / e) + 2.0 * np.arctan(1.0)
            
        if (-np.sin(angle) * e + np.cos(angle) * p) > 0:
            angle += -4.0 * np.arctan(1.0)

        image = flip_jet(rotate_jet(np.array(Jetarray), -angle, dim=25),'r')
                #e_norm = np.linalg.norm(image)   
        jet_pt.append(np.float32(event.JetAK8_PT))
        jet_eta.append(np.float32(etacenter))
        jet_phi.append(np.float32(phicenter))
        jet_mass.append(np.float32(event.JetAK8_MASS))
        tau21.append(np.float32(event.JetAK8_Tau21))
        if 'Wprime' in files :
            signal.append(1)
        else:            
            signal.append(0)
                    
        entries.append(image)
 
    print 'Events: ',count
    print 'Total Events: ', total_count
            
#plot_mean_jet(entries).savefig(output+'.pdf')
#plot_jet(entries[0]).savefig('Wprime_few.pdf')
np.savez(output,image=entries,signal=signal,jet_pt=jet_pt,jet_eta=jet_eta,jet_phi=jet_phi,jet_mass=jet_mass,tau21=tau21)


#                
                       

