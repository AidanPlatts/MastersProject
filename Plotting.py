import numpy as np
import matplotlib.pyplot as plt
#import pyfits as F
from astropy.table import Table
from astropy.io import fits
from scipy import stats as st
import math

rosssample= fits.open(r'full_sample_.fits')
data=rosssample[1].data

bhsample = fits.open(r'blackholes.fits')
bhdata=bhsample[1].data


vollim = fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\fits files\fullvollim.fits')
vollimdata=vollim[1].data
notvollim=fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\fits files\notvollim.fits')
notvoldata=notvollim[1].data

bhlim=fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\fits files\blackholevollim.fits')
bhlimdata=bhlim[1].data
bhnot=fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\fits files\bhnot.fits')
bhnotdata=bhnot[1].data

samplepsi=fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\fits files\fullsamplepsi.fits')
datapsi=samplepsi[1].data
bhsamplepsi=fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\fits files\blackholespsi.fits')
bhdatapsi=bhsamplepsi[1].data

bhsamp = fits.open(r'C:\Users\aidpl\Documents\Uni\Masters Project\blackholess.fits')
bhsampdata=bhsamp[1].data

#redshift vs stellar mass plot, with black holes marked
def zmass():
    plt.figure(figsize=(5,5))
    plt.scatter(notvoldata['z_cone'], notvoldata['logM'],marker='o',facecolors='none',edgecolors='blue',alpha=0.6)
    plt.scatter(vollimdata['z_cone'],vollimdata['logM'],marker='s',facecolors='none', edgecolors='midnightblue',alpha=0.6)
    plt.scatter(bhlimdata['z_cone'],bhlimdata['logM'],c='red',marker='^')
    plt.scatter(bhnotdata['z_cone'],bhnotdata['logM'],c='red',marker='^')
    plt.plot((0.015,0.06),(10,10),c='black',linestyle='dashed',linewidth=2) #volume limited sample line
    plt.xlabel(r'z')
    plt.ylabel(r'$\log (M_{*}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('redshift-massplot.png')
    rvalue(data['z_cone'],data['logM'])
    rvalue(bhdata['z_cone'],bhdata['logM'])

#black hole mass vs stellar mass plot, with splitted styles for the volume limited samples
def bhmass():
    plt.figure(figsize=(5,5))
    plt.scatter(bhlimdata['logM'],bhlimdata['log(MBH)'],marker='^',facecolors='none',edgecolors='darkred',alpha=0.6)
    plt.scatter(bhnotdata['logM'],bhnotdata['log(MBH)'],marker='^',facecolors='none',edgecolors='red',alpha=0.6)
    plt.plot((10,10),(5.3,8.1),c='black',linestyle='dashed',linewidth=2)
    
    #calculate lines of best fit and plot them
    x=np.linspace(9.4,10.9,100)
    
    xlim=np.linspace(10,10.9,100)
    m,c=np.polyfit(bhlimdata['logM'],bhlimdata['log(MBH)'],1)
    plt.plot(xlim, m*xlim+c,c='darkred',linestyle=('dotted'),label='Volume limited sample')
    
    xnot=np.linspace(9.4,10,100)
    m,c=np.polyfit(bhnotdata['logM'],bhnotdata['log(MBH)'],1)
    plt.plot(xnot, m*xnot+c,c='red',linestyle=('dotted'),label='Outside volume limited sample')
    
    plt.plot(x,8.20+1.12*x-11,c='grey',linestyle='dashdot',label='Haring & Rix 2004')
    plt.plot(x,7.45+1.05*x-11,c='darkgrey',linestyle='dashed',label='Reines & Volonteri 2015')
    
    plt.ylabel(r'$\log (M_{BH}/\rm{M}_{\odot})$')
    plt.xlabel(r'$\log (M_{*}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.legend(fontsize='small',framealpha=0.0)
    plt.savefig('mass-massplot.png')
    
    #find correlation values
    rvalue(bhdata['logM'],bhdata['log(MBH)'])
    rvalue(bhlimdata['logM'],bhlimdata['log(MBH)'])
    rvalue(bhnotdata['logM'],bhnotdata['log(MBH)'])


#plot windedness vs pitch angle plot
def psiw():
    volimpsi=[]
    notvolimpsi=[]
    wavg=[]
    wavglim=[]
    wavgnotlim=[]
    psisparc=[]
    
    #calculate w_avg for the different samples and append to empty lists for ease of plotting
    for i in range(0,len(datapsi['psi_gz2'])):
        if datapsi['logM'][i] >= 10:
            volimpsi.append(datapsi['psi_gz2'][i])
            wavglim.append(datapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*datapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*datapsi['t10_arms_winding_a30_loose_debiased'][i]))
        else:
            wavgnotlim.append(datapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*datapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*datapsi['t10_arms_winding_a30_loose_debiased'][i]))
            notvolimpsi.append(datapsi['psi_gz2'][i])
    for i in range(0,len(datapsi['psi_sparcfire'])):
        if datapsi['psi_sparcfire'][i]>0:
            wavg.append(datapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*datapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*datapsi['t10_arms_winding_a30_loose_debiased'][i]))
            psisparc.append(datapsi['psi_sparcfire'][i])
    
    plt.figure(figsize=(5,5))
    plt.scatter(wavglim,volimpsi,marker='s',facecolors='none',edgecolors='midnightblue',alpha=0.6)
    plt.scatter(wavgnotlim,notvolimpsi,marker='o',facecolors='none',edgecolors='b',alpha=0.6)
    plt.scatter(wavg,psisparc,marker='s',facecolors='none',edgecolors='red',alpha=0.3)
    plt.xlabel(r'$w_{avg}$')
    plt.ylabel(r'$\psi\,[\,^o\,]$')
    plt.tight_layout()
    plt.savefig('windedness.png')
    rvalue(wavg,psisparc)

#1d histogram for the redshift
def histogram():
    plt.figure(figsize=(5,5))
    plt.hist(data['z_cone'])
    plt.xlabel(r'z')
    plt.tight_layout()
    plt.savefig('hist1d.png')
   
#2d histogram for redshift and stellar mass
def hist2d():
    x=data['z_cone']
    y=data['logM']
    x_min = np.min(x) 
    x_max = np.max(x) 
  
    y_min = np.min(y) 
    y_max = np.max(y) 
  
    #calculate bins with binsize = 20
    x_bins = np.linspace(x_min, x_max, 20) 
    y_bins = np.linspace(y_min, y_max, 20) 
    plt.figure(figsize=(5,5))
    plt.hist2d(data['z_cone'],data['logM'],bins=(x_bins,y_bins), cmap=plt.cm.Greys)
    plt.colorbar()
    plt.xlabel(r'z')
    plt.ylabel(r'$\log (M_{*}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('hist2d.png')
    rvalue(data['z_cone'],data['logM'])

#black hole mass vs windedness/pitch angle subplot (not used in report)    
def bhwsub():
    wavg=[]
    for i in range(0,len(bhdatapsi['psi_gz2'])):
        wavg.append(bhdatapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*bhdatapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*bhdatapsi['t10_arms_winding_a30_loose_debiased'][i]))
    fig=plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    
    ax1.scatter(wavg,bhdatapsi['log(MBH)'],marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(bhdatapsi['log(MBH)'],wavg,1)
    ax1.plot(m*bhdatapsi['log(MBH)']+c,bhdatapsi['log(MBH)'],c='black',linestyle='dashed')
    
    ax2.scatter(bhdatapsi['psi_gz2'],bhdatapsi['log(MBH)'],marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'],1)
    ax2.plot(m*bhdatapsi['log(MBH)']+c,bhdatapsi['log(MBH)'],c='black',linestyle='dashed')
    
    ax1.set_xlabel(r'$w_{avg}$')
    ax2.set_xlabel(r'$\psi\,[\,^o\,]$')
    ax1.set_ylabel(r'$\log (M_{BH}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('bhmasswsub.png')
    rvalue(wavg,bhdatapsi['log(MBH)'])
    rvalue(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'])

#stellar mass vs windedness/pitch angle plot (not used in report)  
def smwsub():
    wavg=[]
    for i in range(0,len(datapsi['psi_gz2'])):
        wavg.append(datapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*datapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*datapsi['t10_arms_winding_a30_loose_debiased'][i]))
    fig=plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    #f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.scatter(wavg,datapsi['logM'],marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(datapsi['logM'],wavg,1)
    ax1.plot(m*datapsi['logM']+c,datapsi['logM'],c='black',linestyle='dashed')
    
    ax2.scatter(datapsi['psi_gz2'],datapsi['logM'],marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(datapsi['logM'],datapsi['psi_gz2'],1)
    ax2.plot(m*datapsi['logM']+c,datapsi['logM'],c='black',linestyle='dashed')
    
    ax1.set_xlabel(r'$w_{avg}$')
    ax2.set_xlabel(r'$\psi\,[\,^o\,]$')
    ax1.set_ylabel(r'$\log (M_{*}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('stellarmasswsub.png')
    rvalue(datapsi['logM'],wavg)
    rvalue(datapsi['logM'],datapsi['psi_gz2'])

#black hole mass vs windedness plot (not used in report)    
def bhw():
    wavg=[]
    for i in range(0,len(bhdatapsi['psi_gz2'])):
        wavg.append(bhdatapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*bhdatapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*bhdatapsi['t10_arms_winding_a30_loose_debiased'][i]))
    plt.figure(figsize=(5,5))
    plt.scatter(bhdatapsi['log(MBH)'],wavg,marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(bhdatapsi['log(MBH)'],wavg,1)
    plt.plot(bhdatapsi['log(MBH)'], m*bhdatapsi['log(MBH)']+c,c='black',linestyle='dashed')
    plt.ylabel(r'$w_{avg}$')
    plt.xlabel(r'$\log (M_{BH}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('bhmassw.png')
    rvalue(bhdatapsi['log(MBH)'],wavg)

#black hole mass vs pitch angle plot
def bhpsi():
    x=np.array([22.5,14.3,7.1,17.3,13.3,28.8,20.1,8.6,39.8,17.0,13.8,18.6,10.6,34.6,17.8,14.4,9.3,30.1,10.7,24.4,14.4])
    y=np.array([3.7e6,1.6e7,1.7e8,2e7,3.8e7,1.4e6,1.6e7,3.9e7,5e4,8.1e6,2.2e7,1.7e7,4.2e7,8.4e5,9.1e6,1.2e7,7.1e7,1.3e6,6.5e7,3e6,1.7e7])
    
    ylimdown=np.array([1.5e3,5e5])
    xlimdown=np.array([42.2,37.1])
    
    ylimup=np.array([1.5e3,4e3,1.8e3,4.8e3])
    xlimup=np.array([33.7,41.2,43.4,41.4])
    
    y=np.log10(y)
    ylimup=np.log10(ylimup)
    ylimdown=np.log10(ylimdown)
    
    xline=[]
    yline=[]
    for i in range(0,len(x)):
        if x[i]<38:
            xline.append(x[i])
            yline.append(y[i])
            
    for i in range(0,len(xlimup)):
        if xlimup[i]<38:
            xline.append(xlimup[i])
            yline.append(ylimup[i])
            
    for i in range(0,len(xlimdown)): 
        if xlimdown[i]<38:
            xline.append(xlimdown[i])
            yline.append(ylimdown[i])
    
    plt.figure(figsize=(5,5))
    mbh=[]
    psi=[]
    for i in range(0,len(bhdatapsi['log(MBH)'])):
        if bhdatapsi['logM'][i] > 10:
            mbh.append(bhdatapsi['log(MBH)'][i])
            psi.append(bhdatapsi['psi_gz2'][i])
    x=np.linspace(5.5,8.2,100)
    plt.scatter(psi,mbh,marker='^',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(mbh,psi,1)
    plt.plot(m*x+c,x, c='black',linestyle='dashed',label='Volume limited sample')
    m,c=np.polyfit(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'],1)
    plt.plot(m*x+c,x, c='grey',linestyle='dotted',label='Full sample')
    
    xl=np.linspace(16,25,100)
    m,c=np.polyfit(xline,yline,1)
    plt.plot(xl,m*xl+c ,c='blue',linestyle='dashdot',label='Seigar et al. 2008')
    
    plt.xlabel(r'$\psi\,[\,^o\,]$')
    plt.ylabel(r'$\log (\rm{M}_{BH}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.legend()
    plt.savefig('bhmasspsi.png')
    rvalue(mbh,psi)
    rvalue(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'])

def smw():
    wavg=[]
    for i in range(0,len(datapsi['psi_gz2'])):
        wavg.append(datapsi['t10_arms_winding_a28_tight_debiased'][i]+(2*datapsi['t10_arms_winding_a29_medium_debiased'][i])+(3*datapsi['t10_arms_winding_a30_loose_debiased'][i]))
    plt.figure(figsize=(5,5))
    plt.scatter(wavg,datapsi['logM'],marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(wavg,datapsi['logM'],1)
    plt.plot(np.array(wavg), m*np.array(wavg)+c,c='black',linestyle='dashed')
    plt.xlabel(r'$w_{avg}$')
    plt.ylabel(r'$\log (M_{*}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('stellarmassw.png')
    rvalue(wavg,datapsi['logM'])

def smpsi():
    plt.figure(figsize=(5,5))
    volimM=[]
    notvolimM=[]
    volimpsi=[]
    notvolimpsi=[]
    for i in range(0,len(datapsi)):
        if datapsi['logM'][i] >= 10:
            volimM.append(datapsi['logM'][i])
            volimpsi.append(datapsi['psi_gz2'][i])
        else:
            notvolimM.append(datapsi['logM'][i])
            notvolimpsi.append(datapsi['psi_gz2'][i])
    #plt.scatter(datapsi['psi_gz2'],datapsi['logM'],marker='o',facecolors='none',edgecolors='blue',alpha=0.6)
    
    
    plt.scatter(volimpsi,volimM,marker='s',facecolors='none',edgecolors='midnightblue',alpha=0.6)
    plt.scatter(notvolimpsi,notvolimM,marker='o',facecolors='none',edgecolors='blue',alpha=0.6)
    plt.plot((15,27),(10,10),c='black',linestyle='dashed',linewidth=2)
    
    m,c=np.polyfit(datapsi['logM'],datapsi['psi_gz2'],1)
    plt.plot(m*datapsi['logM']+c,datapsi['logM'],c='black',linestyle='dotted')
    
    plt.xlabel(r'$\psi\,[\,^o\,]$')
    plt.ylabel(r'$\log (M_{*}/\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('stellarmasspsi.png')
    rvalue(datapsi['logM'],datapsi['psi_gz2'])

def rvalue(x,y):
    #calculate the spearman correlation coefficient and corresponding p value and output them
    print(st.spearmanr(x,y)[0])
    print(st.spearmanr(x,y)[1])

def seigarplotlog():
    #reproduce the plot from Seigar et. al. 2008
    x=[22.5,14.3,7.1,17.3,13.3,28.8,20.1,8.6,39.8,17.0,13.8,18.6,10.6,34.6,17.8,14.4,9.3,30.1,10.7,24.4,14.4]
    y=[3.7e6,1.6e7,1.7e8,2e7,3.8e7,1.4e6,1.6e7,3.9e7,5e4,8.1e6,2.2e7,1.7e7,4.2e7,8.4e5,9.1e6,1.2e7,7.1e7,1.3e6,6.5e7,3e6,1.7e7]
    
    xerr=[2.5,3,0.4,2.2,0.9,3.6,0.8,2.4,1.6,2.4,1.6,2.7,1.4,3.8,0.9,2.5,1.1,2.6,2,2,0.7]
    yerr=[0.2e6,0.9e7,0.6e8,1e7,3.2e7,0.8e6,0.9e7,0.3e7,0.5e4,4.6e6,0.4e7,0.4e7,0.5e7,4.3e5,1.7e6,0.4e7,0.5e7,0.4e6,1e7,0.4e6,0.4e7]
    
    ylimdown=[1.5e3,5e5]
    xlimdown=[42.2,37.1]
    ylimup=[1.5e3,4e3,1.8e3,4.8e3]
    xlimup=[33.7,41.2,43.4,41.4]
    xlimerrdown=[3,1.3]
    ylimerrdown=[0.5e3,1.6e5]
    xlimerrup=[1.7,4.2,1.4,3.9]
    ylimerrup=[0.5e3,1.3e3,0.6e3,1.6e3]
    fig=plt.figure(figsize=(5,5))
    ax=fig.add_subplot(1,1,1)
    ax.errorbar(x,y,xerr=xerr,yerr=yerr,fmt='s',capsize=3,c='black')
    ax.errorbar(xlimdown,ylimdown,xerr=xlimerrdown,yerr=ylimerrdown,fmt='s',uplims=False, lolims=True,capsize=3,c='black')
    ax.errorbar(xlimup,ylimup,xerr=xlimerrup,yerr=ylimerrup,fmt='s',uplims=True, lolims=False,capsize=3,c='black')
    xline=[]
    yline=[]
    xrange=[]
    yrange=[]
    for i in range(0,len(x)):
        if x[i]<38:
            xline.append(x[i])
            yline.append(math.log10(y[i]))
        if 16<x[i]<25:
            xrange.append(x[i])
            yrange.append(math.log10(y[i]))
    for i in range(0,len(xlimup)):
        if xlimup[i]<38:
            xline.append(xlimup[i])
            yline.append(math.log10(ylimup[i]))
        if 16<xlimup[i]<25:
            xrange.append(xlimup[i])
            yrange.append(math.log10(ylimup[i]))
    for i in range(0,len(xlimdown)):
        if xlimdown[i]<38:
            xline.append(xlimdown[i])
            yline.append(math.log10(ylimdown[i]))
        if 16<xlimdown[i]<25:
            xrange.append(xlimdown[i])
            yrange.append(math.log10(ylimdown[i]))
    xl=np.linspace(5,38,100)
    m,c=np.polyfit(xline,yline,1)
    plt.plot(xl,10**(m*xl+c),c='black',linestyle='dashed')
    m,c=np.polyfit(xrange,yrange,1)
    plt.plot(xl,10**(m*xl+c),c='red',linestyle='dashdot')
    m,c=np.polyfit(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'],1)
    #plt.plot((m*xl+c),xl, c='grey',linestyle='dotted',label='Full sample')
    ax.set_yscale("log")
    ax.set_xlabel(r'$\psi\,[\,^o\,]$')
    ax.set_ylabel(r'$\rm{M}_{\rm{BH}} \,(\rm{M}_{\odot})$')
    plt.tight_layout()
    plt.savefig('seigar.png')
    rvalue(xrange,yrange)
    #print(xline)
    
def seigarplot():
    #same as above but a slightly different method
    x=[22.5,14.3,7.1,42.2,37.1,17.3,13.3,28.8,20.1,8.6,39.8,17.0,13.8,18.6,10.6,34.6,17.8,14.4,9.3,30.1,10.7,24.4,14.4,33.7,41.2,43.3,41.4]
    y=[3.7e6,1.6e7,1.7e8,1.5e3,5e5,2e7,3.8e7,1.4e6,1.6e7,3.9e7,5e4,8.1e6,2.2e7,1.7e7,4.2e7,8.4e5,9.1e6,1.2e7,7.1e7,1.3e6,6.5e7,3e6,1.7e7,1.5e5,4e3,1.8e3,4.8e3]
    #xerr=[2.5,3,0.4,3,1.3,2.2,0.9,3.6,0.8,2.4,1.6,2.4,1.6,2.7,1.4,3.8,0.9,2.5,1.1,2.6,2,2,0.7,1.7,4.2,1.4,3.9]
    #yerr=[0.2e6,0.9e7,0.6e8,1e3,1e5,1e7,3.2e7,0.8e6,0.9e7,0.3e7,0.5e4,4.6e6,0.4e7,0.4e7,0.5e7,4.3e5,1.7e6,0.4e7,0.5e7,0.4e6,1e7,0.4e6,0.4e7,1e3,1e3,1e3,1e3]
    xline=[]
    yline=[]
    for i in range(0,len(y)):
        y[i]=math.log10(y[i])
        if x[i]<38:
            xline.append(x[i])
            yline.append(y[i])
    plt.figure(figsize=(5,5))
    plt.scatter(x,y,marker='o',facecolors='none',edgecolors='red',alpha=0.6)
    xl=np.linspace(5,38,100)
    m,c=np.polyfit(xline,yline,1)
    plt.plot(xl,m*xl+c,c='black',linestyle='dashed')
    xn=np.linspace(3,8.2,100)
    m,c=np.polyfit(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'],1)
    plt.plot((m*xn+c),xn,c='grey',linestyle='dotted',label='Full sample')
    plt.ylabel(r'$(M_{BH}/\rm{M}_{\odot})$')
    plt.xlabel(r'$\psi\,[\,^o\,]$')
    rvalue(x,y)
    
def seigarsub():
    #produce a subplot of our bhpsi and the Seigar et. al. plot
    x=[22.5,14.3,7.1,42.2,37.1,17.3,13.3,28.8,20.1,8.6,39.8,17.0,13.8,18.6,10.6,34.6,17.8,14.4,9.3,30.1,10.7,24.4,14.4,33.7,41.2,43.3,41.4]
    y=[3.7e6,1.6e7,1.7e8,1.5e3,5e5,2e7,3.8e7,1.4e6,1.6e7,3.9e7,5e4,8.1e6,2.2e7,1.7e7,4.2e7,8.4e5,9.1e6,1.2e7,7.1e7,1.3e6,6.5e7,3e6,1.7e7,1.5e3,4e3,1.8e3,4.8e3]
    #xerr=[2.5,3,0.4,3,1.3,2.2,0.9,3.6,0.8,2.4,1.6,2.4,1.6,2.7,1.4,3.8,0.9,2.5,1.1,2.6,2,2,0.7,1.7,4.2,1.4,3.9]
    #yerr=[0.2e6,0.9e7,0.6e8,1e3,1e5,1e7,3.2e7,0.8e6,0.9e7,0.3e7,0.5e4,4.6e6,0.4e7,0.4e7,0.5e7,4.3e5,1.7e6,0.4e7,0.5e7,0.4e6,1e7,0.4e6,0.4e7,1e3,1e3,1e3,1e3]
    xline=[]
    yline=[]
    for i in range(0,len(y)):
        y[i]=math.log10(y[i])
        if x[i]<38:
            xline.append(x[i])
            yline.append(y[i])
    fig=plt.figure(figsize=(10,5))
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)
    ax1.scatter(bhdatapsi['psi_gz2'],bhdatapsi['log(MBH)'],marker='^',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(bhdatapsi['log(MBH)'],bhdatapsi['psi_gz2'],1)
    ax1.plot(m*bhdatapsi['log(MBH)']+c,bhdatapsi['log(MBH)'],c='black',linestyle='dashed')

    ax2.scatter(x,y,marker='s',facecolors='none',edgecolors='black',alpha=0.6)
    m,c=np.polyfit(xline,yline,1)
    ax2.plot(np.array(xline),m*np.array(xline)+c,c='black',linestyle='dashed')
    ax1.set_ylabel(r'$\log (\rm{M}_{\rm{BH}}\,/\,\rm{M}_{\odot})$')
    ax1.set_xlabel(r'$\psi\,[\,^o\,]$')
    ax2.set_xlabel(r'$\psi\,[\,^o\,]$')
    plt.tight_layout()
    
    rvalue(datapsi['logM'],datapsi['psi_gz2'])
    rvalue(xline,yline)
    plt.savefig('seigarsubplot.png')

    
def subseigar():
    #as above but utilising the error bars this time
    x=np.array([22.5,14.3,7.1,17.3,13.3,28.8,20.1,8.6,39.8,17.0,13.8,18.6,10.6,34.6,17.8,14.4,9.3,30.1,10.7,24.4,14.4])
    y=np.array([3.7e6,1.6e7,1.7e8,2e7,3.8e7,1.4e6,1.6e7,3.9e7,5e4,8.1e6,2.2e7,1.7e7,4.2e7,8.4e5,9.1e6,1.2e7,7.1e7,1.3e6,6.5e7,3e6,1.7e7])
    
    xerr=np.array([2.5,3,0.4,2.2,0.9,3.6,0.8,2.4,1.6,2.4,1.6,2.7,1.4,3.8,0.9,2.5,1.1,2.6,2,2,0.7])
    yerr=np.array([0.2e6,0.9e7,0.6e8,1e7,3.2e7,0.8e6,0.9e7,0.3e7,0.5e4,4.6e6,0.4e7,0.4e7,0.5e7,4.3e5,1.7e6,0.4e7,0.5e7,0.4e6,1e7,0.4e6,0.4e7])
    
    ylimdown=np.array([1.5e3,5e5])
    xlimdown=np.array([42.2,37.1])
    ylimup=np.array([1.5e3,4e3,1.8e3,4.8e3])
    xlimup=np.array([33.7,41.2,43.4,41.4])
    xlimerrdown=np.array([3,1.3])
    ylimerrdown=np.array([0.5e3,1.6e5])
    xlimerrup=np.array([1.7,4.2,1.4,3.9])
    ylimerrup=np.array([0.5e3,1.3e3,0.6e3,1.6e3])
    
    yerr=yerr/y
    ylimerrup=ylimerrup/ylimup
    ylimerrdown=ylimerrdown/ylimdown
    y=np.log10(y)
    ylimup=np.log10(ylimup)
    ylimdown=np.log10(ylimdown)
    
    GDm=[]
    GDp=[]
    sparcx=[]
    sparcy=[]
    for i in range(0,len(bhsampdata['Grand Design?'])):
        if bhsampdata['Grand Design?'][i] == True and bhsampdata['psi_gz2'][i] > 0:
            GDm.append(bhsampdata['log(MBH)'][i])
            GDp.append(bhsampdata['psi_gz2'][i])
        if bhdatapsi['psi_gz2'][i]>0:
            sparcy.append(bhdatapsi['log(MBH)'][i])
            sparcx.append(bhdatapsi['psi_gz2'][i])
           
    fig=plt.figure(figsize=(10,5))
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)
    
    y1=np.linspace(5.3,8.2,100)
    ax1.scatter(sparcx,sparcy,marker='^',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(sparcy,sparcx,1)
    ax1.plot(m*y1+c,y1,c='black',linestyle='dashed',label='Full sample')
    ax1.scatter(GDp,GDm,marker='^',facecolors='red',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(GDm,GDp,1)
    ax1.plot(m*y1+c,y1,c='red',linestyle='dotted',label='Grand Design galaxies')
    
    #calculate the lines of best fit for the Seigar plot
    xline=[]
    yline=[]
    for i in range(0,len(x)):
        if x[i]<38:
            xline.append(x[i])
            yline.append(y[i])
            
    for i in range(0,len(xlimup)):
        if xlimup[i]<38:
            xline.append(xlimup[i])
            yline.append(ylimup[i])
            
    for i in range(0,len(xlimdown)): 
        if xlimdown[i]<38:
            xline.append(xlimdown[i])
            yline.append(ylimdown[i])
            
    ax2.errorbar(x,y,xerr=xerr,yerr=yerr,fmt='s',capsize=3,c='black')
    ax2.errorbar(xlimdown,ylimdown,xerr=xlimerrdown,yerr=ylimerrdown,fmt='s',uplims=False, lolims=True,capsize=3,c='black')
    ax2.errorbar(xlimup,ylimup,xerr=xlimerrup,yerr=ylimerrup,fmt='s',uplims=True, lolims=False,capsize=3,c='black')
    xl=np.linspace(5,38,100)
    m,c=np.polyfit(xline,yline,1)
    ax2.plot(xl,(m*xl+c),c='black',linestyle='dashed') 
    xl=np.linspace(16,25,100)
    m,c=np.polyfit(xline,yline,1)
    ax1.plot(xl,m*xl+c ,c='blue',linestyle='dashdot',label='Seigar et al. 2008')
    ax1.set_ylabel(r'$\log (\rm{M}_{\rm{BH}}\,/\,\rm{M}_{\odot})$')
    ax1.set_xlabel(r'$\psi\,[\,^o\,]$')
    ax2.set_xlabel(r'$\psi\,[\,^o\,]$')
    ax1.legend()
    plt.tight_layout()
    plt.savefig('seigarsubplot.png')
    
def bhsparc():
    #plot bh mass vs pitch angles measured from sparcfire (not used in report)
    sparcx=[]
    sparcy=[]
    for i in range(0,len(bhdatapsi['psi_sparcfire'])):
        if bhdatapsi['psi_sparcfire'][i] > -99:
            sparcx.append(bhdatapsi['psi_sparcfire'][i])
            sparcy.append(bhdatapsi['log(MBH)'][i])
    
    x=np.array([22.5,14.3,7.1,17.3,13.3,28.8,20.1,8.6,39.8,17.0,13.8,18.6,10.6,34.6,17.8,14.4,9.3,30.1,10.7,24.4,14.4])
    y=np.array([3.7e6,1.6e7,1.7e8,2e7,3.8e7,1.4e6,1.6e7,3.9e7,5e4,8.1e6,2.2e7,1.7e7,4.2e7,8.4e5,9.1e6,1.2e7,7.1e7,1.3e6,6.5e7,3e6,1.7e7])
    
    ylimdown=np.array([1.5e3,5e5])
    xlimdown=np.array([42.2,37.1])
    
    ylimup=np.array([1.5e3,4e3,1.8e3,4.8e3])
    xlimup=np.array([33.7,41.2,43.4,41.4])
    
    y=np.log10(y)
    ylimup=np.log10(ylimup)
    ylimdown=np.log10(ylimdown)
    
    xline=[]
    yline=[]
    for i in range(0,len(x)):
        if x[i]<38:
            xline.append(x[i])
            yline.append(y[i])
            
    for i in range(0,len(xlimup)):
        if xlimup[i]<38:
            xline.append(xlimup[i])
            yline.append(ylimup[i])
            
    for i in range(0,len(xlimdown)): 
        if xlimdown[i]<38:
            xline.append(xlimdown[i])
            yline.append(ylimdown[i])
            
    plt.figure(figsize=(5,5))        
    xl=np.linspace(5,35,100)
    m,c=np.polyfit(xline,yline,1)
    plt.plot(xl,m*xl+c ,c='blue',linestyle='dashdot',label='Seigar et al. 2008')
    #x=np.linspace(5.5,8.2,100)
        
    plt.scatter(sparcx,sparcy,marker='^',facecolors='none',edgecolors='red',alpha=0.6)
    m,c=np.polyfit(sparcx,sparcy,1)
    print(m,c)
    plt.plot(xl,m*xl+c, c='black',linestyle='dashed')
    plt.xlabel(r'$\psi\,[\,^o\,]$')
    plt.ylabel(r'$\log (\rm{M}_{BH}/\rm{M}_{\odot})$')
    plt.tight_layout()
    
    rvalue(sparcx,sparcy)
 