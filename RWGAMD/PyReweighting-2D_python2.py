#! /usr/bin/env python

import math
import scipy
import scipy.stats as stats
import numpy as np
import sys
import matplotlib.pyplot as plt
import csv
from argparse import ArgumentParser
from scipy.optimize import curve_fit
## from scipy.optimize import *

print "============================================================"
print "PyReweighting: Python scripts used to reweight accelerated \n\
               and scaled molecular dynamics simulations."
print "  "
print "Authors: Yinglong Miao <yinglong.miao@gmail.com>"
print "         Bill Sinko <wsinko@gmail.com>"
print "  "
print "Last Update: Dec 30, 2014"
print "\n\
Citations:  \n\
1. Sinko W, Miao Y, de Oliveira CAF, McCammon JA (2013) Population Based Reweighting of Scaled Molecular Dynamics. The Journal of Physical Chemistry B 117(42):12759-12768. \n\
2. Miao Y, Sinko W, Pierce L, Bucher D, Walker RC, McCammon JA (2014) Improved reweighting of accelerated molecular dynamics simulations for free energy calculation. J Chemical Theory and Computation. 10(7): 2677-2689."
print " "

###########MAIN
def main():
## Set control parameters
    plt_figs = 1

    args = cmdlineparse()   
    
    inputfile=loadfiletoarray(args.input)
    data=inputfile[:,0]
    
    rows = len(data)
    weights,dV = weightparse(rows, args)
    if args.Xdim:
        binsX, discX= assignbinsX(args.Xdim, args)
    else:
        max_data = np.ceil(np.amax(inputfile[:,0]))
        min_data = np.floor(np.amin(inputfile[:,0]))
        binsX, discX= assignbinsX([min_data,max_data], args)  ## Default bin size
        # print "(max_data,min_data) = ", max_data,min_data
        # binsX, discX= assignbinsX([-180,180], args)  ## Default bin size
    if args.Ydim:
        binsY, discY= assignbinsY(args.Ydim, args)
    else:
        max_data = np.ceil(np.amax(inputfile[:,1]))
        min_data = np.floor(np.amin(inputfile[:,1]))
        binsY, discY= assignbinsY([min_data,max_data], args)  ## Default bin size
        # binsY, discY= assignbinsY([-180,180], args)  ## Default bin size

##  SET MAX ENERGY FOR ALL INFINITY VALUES
    if args.Emax:
        cb_max=float(args.Emax)
    else :
        cb_max = 8

##  SET HISTOGRAM CUTOFF
    if args.cutoff:
        hist_min=int(args.cutoff)
    else :
        hist_min = 10	# minimum number of configurations in one bin

##  SET ORDER of McLaurin series expansion
    if args.order:
        order=int(args.order)
    else :
        order = 10	# default

##  SET TEMPERATURE
    if args.T:
        T=float(args.T)
    else :
        T = 300	# simulation temperature
    beta = 1.0/(0.001987*T)

##REWEIGHTING
##  SET flag for Gaussian fitting of deltaV
    if args.fit:
        fit=args.fit
    else :
        fit=False	# simulation temperature
##    print "gaussian fitting:", fit

##REWEIGHTING
    if args.job == "amdweight_CE":
    	hist2,newedgesX,newedgesY,c1,c2,c3 = reweight_CE(inputfile,hist_min,binsX,discX,binsY,discY,dV,T,fit)
	pmf = hist2pmf2D(hist2,hist_min,T)
	c1 = -np.multiply(1.0/beta,c1)
	c2 = -np.multiply(1.0/beta,c2)
	c3 = -np.multiply(1.0/beta,c3)

	c12 = np.add(c1,c2)
	c123 = np.add(c12,c3)
	pmf_c1 = np.add(pmf, c1)
	print "pmf_min-c1 = ", np.min(pmf_c1)
	pmf_c1 = normalize2D(pmf_c1,cb_max)
	pmf_c2 = np.add(pmf, c12)
	print "pmf_min-c2 = ", np.min(pmf_c2)
	pmf_c2 = normalize2D(pmf_c2,cb_max)
	pmf_c3 = np.add(pmf, c123)
	print "pmf_min-c3 = ", np.min(pmf_c3)
	pmf_c3 = normalize2D(pmf_c3,cb_max)
    elif args.job == "amdweight_MC":
        n=order
	MCweight=np.zeros(len(dV))
        beta_dV=np.multiply(dV,beta)
	for x in range(0,n+1):
	  MCweight=np.add(MCweight,(np.divide(np.power(beta_dV, x), float(scipy.misc.factorial(x)))))
	weights=MCweight
    	hist2,newedgesX,newedgesY = np.histogram2d(inputfile[:,0], inputfile[:,1], bins = (binsX, binsY), weights=weights)
        hist2=prephist(hist2,T,cb_max)
    elif args.job == "amdweight":
    	hist2,newedgesX,newedgesY = np.histogram2d(inputfile[:,0], inputfile[:,1], bins = (binsX, binsY), weights=weights)
        hist2=prephist(hist2,T,cb_max)
    else :
    	hist2,newedgesX,newedgesY = np.histogram2d(inputfile[:,0], inputfile[:,1], bins = (binsX, binsY), weights=None)
        hist2=prephist(hist2,T,cb_max)

##SAVE FREE ENERGY DATA INTO A FILE
    if args.job == "amdweight_MC" or args.job == "amdweight" or args.job == "noweight" :
    	pmffile = 'pmf-'+str(args.input)+'.xvg'
	output_pmf2D(pmffile,hist2,binsX,binsY)
    if args.job == "amdweight_CE" :
	hist2 = pmf_c1
    	pmffile = 'pmf-c1-'+str(args.input)+'.xvg'
	output_pmf2D(pmffile,hist2,binsX,binsY)

	hist2 = pmf_c3
    	pmffile = 'pmf-c3-'+str(args.input)+'.xvg'
	output_pmf2D(pmffile,hist2,binsX,binsY)

	hist2 = pmf_c2
    	pmffile = 'pmf-c2-'+str(args.input)+'.xvg'
	output_pmf2D(pmffile,hist2,binsX,binsY)

    if args.job == "amd_dV":
        plt_figs = 0
    	hist2,newedgesX,newedgesY,binfX,binfY,dV_avg,dV_std,dV_anharm,dV_mat = reweight_dV(inputfile,hist_min,binsX,binsY,discX,discY,dV,T)

    	pmffile = 'dV-hist-2D-'+str(args.input) + '.xvg'
	output_dV(pmffile,dV)

        alpha = anharm(dV)
        print "Anharmonicity of all dV = " + str(alpha)

    	pmffile = 'dV-anharm-2D-'+str(args.input)+'.xvg'
	output_dV_anharm2D(pmffile,binsX,binsY,dV_anharm)

    	pmffile = 'dV-stat-2D-'+str(args.input)+'.xvg'
	output_dV_stat2D(pmffile,binsX,binsY,dV_avg,dV_std,dV_anharm)

    	pmffile = 'dV-mat-2D-'+str(args.input)+'.xvg'
	output_dV_mat2D(pmffile,binsX,binsY,hist2,dV_avg,dV_std,dV_anharm,dV_mat)

###PLOTTING FUNCTION FOR FREE ENERGY FIGURE
    if plt_figs :
    	cbar_ticks=[0, cb_max*.25, cb_max*.5, cb_max*.75, DUMMY_EMAX]
      	plt.figure(2, figsize=(11,8.5))
    	extent = [newedgesX[0], newedgesX[-1], newedgesY[-1], newedgesY[0]]
    	print "Range of values in 2D Free Energy Map ",extent
      	plt.imshow(hist2.transpose(), extent=extent, interpolation='gaussian')
    	cb = plt.colorbar(ticks=cbar_ticks, format=('% .1f'), aspect=10) # grab the Colorbar instance
    	imaxes = plt.gca()
    	plt.axes(cb.ax)
    	plt.clim(vmin=0,vmax=DUMMY_EMAX)
    	plt.yticks(fontsize=18)
    	plt.axes(imaxes)
    	axis=(min(binsX), max(binsX), min(binsY), max(binsY))
    	plt.axis(axis)
    	plt.xticks(size='18')
    	plt.yticks(size='18')
    	plt.xlabel('RC1',fontsize=18)
    	plt.ylabel('RC2',fontsize=18)
##    	plt.xlabel(r'$\phi$',fontsize=18)
##    	plt.ylabel(r'$\psi$',fontsize=18)
##    	plt.xlabel(r'$\chi$1',fontsize=18)
##    	plt.ylabel(r'$\chi$2',fontsize=18)
    	plt.savefig('2D_Free_energy_surface.png',bbox_inches=0)
    	print "FIGURE SAVED 2D_Free_energy_surface.png"
    
###PLOTTING FUNCTION FOR WEIGHTS histogram
    	[hist, edges] = np.histogram(weights, bins=100)
    	width=np.absolute(np.subtract(edges[0], edges[1]))
    	plt.figure(1, figsize=(11,8.5))
    	plt.bar(edges[:100], hist, width=width, log=True)
    	plt.yscale('log')   ###if typerror is thrown delete .matplotlib/fontList.cache  file
    	plt.xticks(fontsize='18')
    	plt.yticks(fontsize='18')
	plt.savefig('weights.png',bbox_inches=0)
    	print "FIGURE SAVED weights.png"

    print " "
    print "END"

########READ datafiles and print weights

def cmdlineparse():
    parser = ArgumentParser(description="command line arguments")
    parser.add_argument("-input", dest="input", required=True, help="2D input file", metavar="<2D input file>")
    parser.add_argument("-job", dest="job", required=True, help="Reweighting method to use: <noweight>, <weighthist>, <amd_time>, <amd_dV>, <amdweight>, <amdweight_MC>, <amdweight_CE>", metavar="<Job type reweighting method>")
    parser.add_argument("-weight", dest="weight", required=False, help="weight file", metavar="<weight file>")
    parser.add_argument("-Xdim", dest="Xdim", required=False, nargs="+", help="Xdimensions", metavar="<Xmin Xmax >")
    parser.add_argument("-Ydim", dest="Ydim", required=False, nargs="+", help="Ydimension", metavar="<Ymin Ymax >")
    parser.add_argument("-discX", dest="discX", required=False,  help="Discretization size in X dimension", metavar="<discretization-X>")
    parser.add_argument("-discY", dest="discY", required=False,  help="Discretization size in Y dimension", metavar="<discretization-Y>")
    parser.add_argument("-cutoff", dest="cutoff", required=False,  help="histogram cutoff", metavar="<cutoff>")
    parser.add_argument("-T", dest="T", required=False,  help="Temperature", metavar="<Temperature>")
    parser.add_argument("-Emax", dest="Emax", required=False,  help="Maximum free energy", metavar="<Emax>")
    parser.add_argument("-fit", dest="fit", required=False, help="Fit deltaV distribution", metavar="<fit>")
    parser.add_argument("-order", dest="order", required=False, help="Order of Maclaurin series", metavar="<order>")
    args=parser.parse_args()
    return args
    
    
def loadfiletoarray(file):
    loaded=np.loadtxt(file, usecols=[0,1])
    print "DATA LOADED:    "+file
    return loaded

def weightparse(rows, args):
    if args.job == "weighthist":
        data=np.loadtxt(args.weight)
        weights=data[:,0]
        dV = np.zeros(rows)
    elif args.job == "amd_time" or args.job == "amd_dV" or args.job == "amdweight" or args.job == "amdweight_MC" or args.job == "amdweight_CE" :
        data=np.loadtxt(args.weight)
        weights = np.exp(data[:,0])
        dV = data[:,2]
    elif args.job == "noweight":
        weights = np.zeros(rows)
        weights = weights + 1
        dV = np.zeros(rows)
    else:
        print "ERROR JOBTYPE"+ args.job+ " NOT RECOGNIZED"
        del data
        del weights
        del dV
    return weights,dV

def assignbinsX(dim, args):
    minimum=float(dim[0])
    maximum=float(dim[1])
    if args.discX:
        discX=float(args.discX)
    else :
        discX = 6
    binsX =np.arange(minimum,(maximum+discX),discX)
    return binsX,discX

def assignbinsY(dim, args):
    minimum=float(dim[0])
    maximum=float(dim[1])
    if args.discY:
        discY=float(args.discY)
    else :
        discY = 6
    binsY =np.arange(minimum,(maximum+discY),discY)
    return binsY,discY

def normalize2D(pmf,cb_max):
    pmf=pmf-np.min(pmf)  ## zero value to lowest energy state
    temphist=pmf
    # print "rows = ", len(temphist[:,0])
    # print "cols = ", len(temphist[0,:])
    #set infinity free energy values to is cb_max
    for jy in range(len(temphist[0,:])):
      for jx in range(len(temphist[:,0])):
        if np.isinf(temphist[jx,jy]):
                temphist[jx,jy]=cb_max
    return temphist

def prephist(hist2,T,cb_max):
    hist2=np.add(hist2,0.000000000000000001)  ###so that distrib
    hist2=(0.001987*T)*np.log(hist2) ####Convert to free energy in Kcal/mol
    hist2=np.max(hist2)-hist2  ## zero value to lowest energy state
##    print np.max(hist2)
    temphist2=hist2
    #set infinity free energy values to is cb_max
    for jy in range(len(temphist2[0,:])):
        for jx in range(len(temphist2[:,0])):
            if np.isinf(temphist2[jx,jy]):
                temphist2[jx,jy]=cb_max
    return temphist2

def reweight_CE(data,hist_min,binsX,discX,binsY,discY,dV,T,fit):
    hist2, newedgesX, newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=None)

    hist_max = np.max(hist2)
##    print np.max(hist2)

    beta = 1.0/(0.001987*T)
    nf = len(data[:,0])
    nbinsX = len(hist2[:,0])
    nbinsY = len(hist2[0,:])

    c1 = np.zeros((nbinsX,nbinsY)) 
    c2 = np.zeros((nbinsX,nbinsY)) 
    c3 = np.zeros((nbinsX,nbinsY)) 

    binfX = np.zeros(nf) # array for storing assigned bin of each frame
    binfY = np.zeros(nf) # array for storing assigned bin of each frame
    nA = np.zeros((nbinsX,nbinsY)) # nA is equivalent to hist here
    dV_avg = np.zeros((nbinsX,nbinsY)) 
    dV_avg2 = np.zeros((nbinsX,nbinsY)) 
    dV_avg3 = np.zeros((nbinsX,nbinsY)) 
    dV_std = np.zeros((nbinsX,nbinsY)) 
    dV_anharm = np.zeros((nbinsX,nbinsY)) 
    dV_mat = np.zeros((nbinsX,nbinsY,hist_max)) # matrix for storing dV of each assigned 

    dV_avg_all=np.average(dV)
    dV_std_all=np.std(dV)
    print 'dV all: avg = ', dV_avg_all, 'std = ', dV_std_all

    diff_tol_avg = 10
    diff_tol_std = 1
    dV_binsize = 50

    for i in range(len(data[:,0])):
	jx = int((data[i,0]-binsX[0])/discX)
	jy = int((data[i,1]-binsY[0])/discY)
	if jx < nbinsX and jy < nbinsY :
	  binfX[i] = jx
	  binfY[i] = jy
	  dV_mat[jx,jy,nA[jx,jy]] = dV[i]
	  nA[jx,jy] = nA[jx,jy]+1

    for jx in range(nbinsX):
      for jy in range(nbinsY):
        dV_anharm[jx,jy] = 100
        if nA[jx,jy]>=hist_min :
	  num = int(nA[jx,jy])
          atemp = np.zeros(num)
          atemp2 = np.zeros(num)
          atemp3 = np.zeros(num)
          for k in range(num):
	     atemp[k] = dV_mat[jx,jy,k]
	     atemp2[k] = dV_mat[jx,jy,k]**2
	     atemp3[k] = dV_mat[jx,jy,k]**3
 	  dV_avg[jx,jy] = np.average(atemp)
	  dV_std[jx,jy]=np.std(atemp)
          dV_anharm[jx,jy] = anharm(atemp)

	  if np.absolute(dV_avg[jx,jy]-dV_avg_all)>diff_tol_avg or np.absolute(dV_std[jx,jy]-dV_std_all)>diff_tol_std :
	     dV_avg[jx,jy]=0
	     dV_std[jx,jy]=0
##	       print binsX[j], ': corrected mean = ', dV_avg[j], 'sigma = ', dV_std[j]

 	  dV_avg2[jx,jy]=np.average(atemp2)
 	  dV_avg3[jx,jy]=np.average(atemp3)
	  del atemp
	  del atemp2
	  del atemp3
	  c1[jx,jy] = beta*dV_avg[jx,jy]
	  c2[jx,jy] = 0.5*beta**2*dV_std[jx,jy]**2
	  c3[jx,jy] = (1.0/6.0)*beta**3*(dV_avg3[jx,jy]-3.0*dV_avg2[jx,jy]*dV_avg[jx,jy]+2.0*dV_avg[jx,jy]**3)
    return hist2,newedgesX,newedgesY,c1,c2,c3

def reweight_dV(data,hist_min,binsX,binsY,discX,discY,dV,T):
    hist2, newedgesX, newedgesY = np.histogram2d(data[:,0], data[:,1], bins = (binsX, binsY), weights=None)
    hist_max = np.max(hist2)
##    print np.max(hist2)

    nf = len(data[:,0])
    nbinsX = len(hist2[:,0])
    nbinsY = len(hist2[0,:])

    binfX = np.zeros(nf) # array for storing assigned bin of each frame
    binfY = np.zeros(nf) # array for storing assigned bin of each frame
    nA = np.zeros((nbinsX,nbinsY)) # nA is equivalent to hist here
    dV_avg = np.zeros((nbinsX,nbinsY)) 
    dV_std = np.zeros((nbinsX,nbinsY)) 
    dV_anharm = np.zeros((nbinsX,nbinsY)) 
    dV_mat = np.zeros((nbinsX,nbinsY,hist_max)) # matrix for storing dV of each assigned 

    for i in range(len(data[:,0])):
	jx = int((data[i,0]-binsX[0])/discX)
	jy = int((data[i,1]-binsY[0])/discY)
	if jx < nbinsX and jy < nbinsY :
	  binfX[i] = jx
	  binfY[i] = jy
	  dV_mat[jx,jy,nA[jx,jy]] = dV[i]
	  nA[jx,jy] = nA[jx,jy]+1

    for jx in range(nbinsX):
      for jy in range(nbinsY):
        dV_anharm[jx,jy] = 100
        if nA[jx,jy]>=hist_min :
	  num = int(nA[jx,jy])
          atemp = np.zeros(num)
          for k in range(num):
	     atemp[k] = dV_mat[jx,jy,k]
 	  dV_avg[jx,jy] = np.average(atemp)
	  dV_std[jx,jy]=np.std(atemp)
          dV_anharm[jx,jy] = anharm(atemp)
	  del atemp
    return hist2,newedgesX,newedgesY,binfX,binfY,dV_avg,dV_std,dV_anharm,dV_mat

##  Convert histogram to free energy in Kcal/mol
def hist2pmf2D(hist,hist_min,T):
        nbinsX = len(hist[:,0])
        nbinsY = len(hist[0,:])
    	pmf = np.zeros((nbinsX,nbinsY))
	pmf_min = 100
    	for jx in range(len(hist[:,0])):
    	  for jy in range(len(hist[0,:])):
      	   if hist[jx,jy]>=hist_min :
	      pmf[jx,jy]=-(0.001987*T)*np.log(hist[jx,jy])
     	   if pmf_min > pmf[jx,jy] :
	      pmf_min=pmf[jx,jy]
##        pmf=pmf-pmf_min  ## zero value to lowest energy state
    	return pmf

def output_pmf2D(pmffile,hist,binsX,binsY):
    	fpmf = open(pmffile, 'w')
    	strpmf='#RC1\tRC2\tPMF(kcal/mol)\n\n@    xaxis  label \"RC1\"\n@    yaxis  label \"RC2\"\n@TYPE xy\n'
    	fpmf.write(strpmf)
    	for jx in range(len(hist[:,0])):
    	  for jy in range(len(hist[0,:])):
		strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(hist[jx,jy]) + '\n'
		fpmf.write(strpmf)
    	fpmf.closed
    	return fpmf

def output_dV(pmffile,dV):
    	fpmf = open(pmffile, 'w')
    	strpmf='#dV \tp(dV) \n\n@    xaxis  label \"dV\"\n@    yaxis  label \"p(dV)\"\n@TYPE xy\n'
	hist_dV, bin_dV = np.histogram(dV, bins=50)
        for k in range(len(hist_dV)):
		strpmf=strpmf + str(bin_dV[k]) + ' \t' + str(hist_dV[k]) + ' \n'
	fpmf.write(strpmf)
    	fpmf.closed
    	return fpmf

def output_dV_anharm2D(pmffile,binsX,binsY,dV_anharm):
    	fpmf = open(pmffile, 'w')
    	strpmf='#RC \tdV_anharm \tError\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV_anmarm\"\n@TYPE xy\n'
    	fpmf.write(strpmf)
    	for jx in range(len(dV_anharm[:,0])):
    	  for jy in range(len(dV_anharm[0,:])):
		strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(dV_anharm[jx,jy]) + '\n'
		fpmf.write(strpmf)
    	fpmf.closed
    	return fpmf

def output_dV_stat2D(pmffile,binsX,binsY,dV_avg,dV_std,dV_anharm):
    	fpmf = open(pmffile, 'w')
    	strpmf='#RC \tdV_avg(kcal/mol) \tError\n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xydy\n'
    	fpmf.write(strpmf)
    	for jx in range(len(dV_anharm[:,0])):
    	  for jy in range(len(dV_anharm[0,:])):
		strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(dV_avg[jx,jy]) + ' \t' + str(dV_std[jx,jy]) + ' \t' + str(dV_anharm[jx,jy]) + '\n'
		fpmf.write(strpmf)
    	fpmf.closed
    	return fpmf

def output_dV_mat2D(pmffile,binsX,binsY,hist,dV_avg,dV_std,dV_anharm,dV_mat):
    	fpmf = open(pmffile, 'w')
    	strpmf='#RC \tNf \tdV_avg \tdV_std \tdV_ij \n\n@    xaxis  label \"RC\"\n@    yaxis  label \"dV(kcal/mol)\"\n@TYPE xy\n'
    	fpmf.write(strpmf)
    	for jx in range(len(hist[:,0])):
    	  for jy in range(len(hist[0,:])):
		nf_j = int(hist[jx,jy])
		strpmf=str(binsX[jx]) + ' \t' + str(binsY[jy]) + ' \t' + str(hist[jx,jy]) + ' \t' + str(dV_avg[jx,jy]) + ' \t' + str(dV_std[jx,jy]) + ' \t' + str(dV_anharm[jx,jy])
        	for k in range(int(nf_j)):
			strpmf=strpmf + ' \t' + str(dV_mat[jx,jy,k])
		strpmf=strpmf + '\n'
		fpmf.write(strpmf)
    	fpmf.closed
    	return fpmf

def anharm(data):
#    print "Compute anharmonicity"
    var=np.var(data)
    hist, edges=np.histogram(data, 50, normed=True)
    hist=np.add(hist,0.000000000000000001)  ###so that distrib
    dx=edges[1]-edges[0]
    S1=-1*np.trapz(np.multiply(hist, np.log(hist)),dx=dx)
    S2=0.5*np.log(2.00*np.pi*np.exp(1.0)*var+0.000000000000000001)
    alpha=S2-S1
    if np.isinf(alpha):
       alpha = 100
    return alpha
 
if __name__ == '__main__':
    main()
    
