
import os,sys,shutil,datetime,time
import getpass
from ROOT import *

start_time = time.time()

processes = ['QCD_HT200to300', 'QCD_HT300to500', 'QCD_HT500to700', 'QCD_HT700to1000', 'WprimeToWZToWhadZinv_M-600', 'WprimeToWZToWhadZinv_M-800'] 
pileups = ['_0PU', '_35PU', '_140PU', '_200PU']

count=0

for process in processes:
	for pileup in pileups:

		shift = process + pileup
		
		#IO directories must be full paths
		
		#relbase = '/user_data/nelson/CMSSW_8_0_4/'
		inputDir='/mnt/hadoop/users/nelson/DelphesJets/'+shift+'/'
		outputDir= '/user_data/ktaylor/GridTree/'+shift+'/'
		
		runDir=os.getcwd()
		
		#gROOT.ProcessLine('.x compileStep3.C')
		
		cTime=datetime.datetime.now()
		date='%i_%i_%i_%i_%i_%i'%(cTime.year,cTime.month,cTime.day,cTime.hour,cTime.minute,cTime.second)
		
		condorDir=runDir+'/condorLogs/'+shift+'/'
		
		print 'Getting proxy'
		#os.system('voms-proxy-init -valid 168:00')
		proxyPath="blabla"#os.popen('voms-proxy-info -path')
		proxyPath="blabla"#proxyPath.readline().strip()
		
		print 'Starting submission'
		#count=0
		
		rootfiles = os.popen('ls '+inputDir)
		os.system('mkdir -p '+outputDir)
		os.system('mkdir -p '+condorDir)
		
		for file in rootfiles:
			print file    
			rawname = file[:-6]
			print rawname    
		
			count+=1
			dict={'RUNDIR':runDir, 'CONDORDIR':condorDir, 'INPUTDIR':inputDir, 'FILENAME':rawname, 'PROXY':proxyPath, 'OUTPUTDIR':outputDir}
			jdfName=condorDir+'/%(FILENAME)s.job'%dict
			print jdfName
			jdf=open(jdfName,'w')
			jdf.write(
		        
"""universe = vanilla
Executable = %(RUNDIR)s/makeJetGrids.sh 
Should_Transfer_Files = NO
Transfer_Input_Files = 
Output = %(FILENAME)s.out
Error = %(FILENAME)s.err
Log = %(FILENAME)s.log
Notification = Never
Arguments = %(FILENAME)s.root %(FILENAME)s.root %(INPUTDIR)s %(OUTPUTDIR)s %(RUNDIR)s
Queue 1""" %dict)
			
			jdf.close()
			os.chdir('%s/'%(condorDir))
			os.system('condor_submit %(FILENAME)s.job'%dict)
			os.system('sleep 0.5')                                
			os.chdir('%s'%(runDir))
			print count, "jobs submitted..."
		
		
print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))
