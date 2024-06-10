import numpy as np
import pylab as plt
import matplotlib
import scipy.stats as stats
from scipy.stats import norm
from astroML.stats import binned_statistic_2d


def plotXYZ(df):

    X = df['X']
    Y = df['Y']
    Z = df['Z']
    
    fig=plt.figure(2,figsize=(12,8))
    fig.subplots_adjust(wspace=0.2,hspace=0.2,top=0.97,bottom=0.1,left=0.09,right=0.95)
    plt.rcParams['font.size'] = 14

    ax2=fig.add_subplot(221)
    bla = ax2.hexbin(X, Y, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Y [kpc]')
    ax2.set_xlabel('X [kpc]')
    ax2.set_xlim(3,13)
    ax2.set_ylim(-6,6)

    
    ax2=fig.add_subplot(222)
    bla = ax2.hexbin(X, Z, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_xlabel('X [kpc]')
    ax2.set_xlim(3,13)
    ax2.set_ylim(-4,4)

    ax2=fig.add_subplot(223)
    bla = ax2.hexbin(Y, Z, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_xlabel('Y [kpc]')
    ax2.set_xlim(-6,6)
    ax2.set_ylim(-4,4)

    ax2=fig.add_subplot(224)
    bla = ax2.hexbin(np.sqrt(X**2+Y**2), Z, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_xlabel('R [kpc]')
    ax2.set_xlim(3,13)
    ax2.set_ylim(-4,4)

    plt.savefig('../plots/plotXYZ.png')
    print('made plot: ../plots/plotXYZ.png')
    plt.show() 
    plt.close("all")
    return 



def plotXYZgiants(df):

    X = df['X']
    Y = df['Y']
    Z = df['Z']
    
    fig=plt.figure(2,figsize=(12,8))
    fig.subplots_adjust(wspace=0.2,hspace=0.2,top=0.97,bottom=0.1,left=0.09,right=0.95)
    plt.rcParams['font.size'] = 14

    ax2=fig.add_subplot(221)
    bla = ax2.hexbin(X, Y, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Y [kpc]')
    ax2.set_xlabel('X [kpc]')
    ax2.set_xlim(-5,20)
    ax2.set_ylim(-15,15)

    
    ax2=fig.add_subplot(222)
    bla = ax2.hexbin(X, Z, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_xlabel('X [kpc]')
    ax2.set_xlim(-5,20)
    ax2.set_ylim(-15,15)

    ax2=fig.add_subplot(223)
    bla = ax2.hexbin(Y, Z, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_xlabel('Y [kpc]')
    ax2.set_xlim(-15, 15)
    ax2.set_ylim(-15,15)

    ax2=fig.add_subplot(224)
    bla = ax2.hexbin(np.sqrt(X**2+Y**2), Z, gridsize=200, bins='log', cmap='viridis')
    cb = fig.colorbar(bla)
    cb.set_label('N')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_xlabel('R [kpc]')
    ax2.set_xlim(0, 25)
    ax2.set_ylim(-15,15)

    plt.savefig('../plots/plotXYZ.png')
    print('made plot: ../plots/plotXYZ.png')
    plt.show() 
    plt.close("all")
    return 



def plot3VvsZ(df, title=""):

    X = df['X']
    Y = df['Y']
    Z = np.abs(df['Z'])
    v_phi = df['v_phi']
    v_R = df['v_R']
    v_Z = df['v_Z']

    def getStats(X, V, Xmin=0, Xmax=1500, Npts=250):
        mean = stats.binned_statistic(X, V, statistic='mean', bins=np.arange(Xmin, Xmax, Npts))
        std = stats.binned_statistic(X, V, statistic="std", bins=np.arange(Xmin, Xmax, Npts))
        return mean, std

    def plotPanel(ax, X, V, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax):
        mean, std = getStats(X, V)
        ax.scatter(X, V, c="red", s=.1)
        ax.scatter(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic, marker='^', color="yellow", edgecolor='black', s=60)
        ax.plot(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic+2*std.statistic, ls='dashed', color="black", lw=3)
        ax.plot(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic-2*std.statistic, ls='dashed', color="black", lw=3)
        ax.set_xlabel(Xlabel)
        ax.set_ylabel(Ylabel)
        ax.set_xlim(Xmin, Xmax) 
        ax.set_ylim(Ymin, Ymax)
        return
    
    fig=plt.figure(1,figsize=(12,4))
    fig.subplots_adjust(wspace=0.2,hspace=0.2,top=0.97,bottom=0.1,left=0.09,right=0.95)
    plt.rcParams['font.size'] = 14

    ## vPhi
    ax=fig.add_subplot(131)
    plotPanel(ax, Z*1000, v_phi, '|Z| [pc]', '$v_{\\phi}$ [km/s]', 0, 2000, -300,0)
    ax.set_title(title)
 
    # eq. 17 from Bond+(2010)
    Zg = np.linspace(0, 2000, 10)
    vPhiM = -205 + 19.2*(Zg/1000)**1.25
    sigmaM = 30 + 3.0*(Zg/1000)**2.0
    ax.plot(Zg, 0*vPhiM, ls='dashdot', color="blue", lw=2)
    ax.plot(Zg, vPhiM, ls='dashdot', color="blue", lw=2)
    ax.plot(Zg, vPhiM-2*sigmaM, ls='dashdot', color="blue", lw=1)
    ax.plot(Zg, vPhiM+2*sigmaM, ls='dashdot', color="blue", lw=1)

    # new model
    Zg2 = np.linspace(0, 2.0, 100)
    vPhiMean, vPhiDisp = getvPhiMeanSig(Zg2)
    #ax.plot(Zg2*1000, vPhiMean, ls='dashdot', color="green", lw=2)
    ax.plot(Zg2*1000, vPhiMean-2*vPhiDisp, ls='dashed', color="green", lw=1)
    ax.plot(Zg2*1000, vPhiMean+2*vPhiDisp, ls='dashed', color="green", lw=1)
    print('new model:', np.median(vPhiMean), np.median(vPhiDisp)) 
    
    ## vR 
    ax=fig.add_subplot(132)
    plotPanel(ax, Z*1000, v_R, '|Z| [pc]', '$v_R$ [km/s]', 0, 2000, -200,200)

    vM = 0*(Zg/1000)**1.25
    sigmaM = 40 + 5.0*(Zg/1000)**1.5
    ax.plot(Zg, vM, ls='dashdot', color="blue", lw=2)
    ax.plot(Zg, vM-2*sigmaM, ls='dashdot', color="blue", lw=1)
    ax.plot(Zg, vM+2*sigmaM, ls='dashdot', color="blue", lw=1)

    ## vZ 
    ax=fig.add_subplot(133)
    plotPanel(ax, Z*1000, v_Z, '|Z| [pc]', '$v_Z$ [km/s]', 0, 2000, -200,200)

    vM = 0*(Zg/1000)**1.25
    sigmaM = 25 + 4.0*(Zg/1000)**1.5
    # sigmaM = 18 + 4.0*(Zg/1000)**1.5
    ax.plot(Zg, vM, ls='dashdot', color="blue", lw=2)
    ax.plot(Zg, vM-2*sigmaM, ls='dashdot', color="blue", lw=1)
    ax.plot(Zg, vM+2*sigmaM, ls='dashdot', color="blue", lw=1)
    
    plt.tight_layout()
    plt.savefig('../plots/plot3VvsZ.png')
    print('made plot: ../plots/plot3VvsZ.png')
    plt.show() 
    plt.close("all")
    return 



def plot3VvsZgiants(df, title="", disk=True):

    X = df['X']
    Y = df['Y']
    Z = np.abs(df['Z'])
    v_phi = df['v_phi']
    v_R = df['v_R']
    v_Z = df['v_Z']
    if disk!=True:
        v_r = df['v_r']
        v_th = df['v_th']

    def getStats(X, V, Xmin=0, Xmax=1500, Npts=250):
        Xok = X[np.abs(V)<350]
        Vok = V[np.abs(V)<350]
        mean = stats.binned_statistic(Xok, Vok, statistic='mean', bins=np.arange(Xmin, Xmax, Npts))
        std = stats.binned_statistic(Xok, Vok, statistic="std", bins=np.arange(Xmin, Xmax, Npts))
        return mean, std

    def plotPanel(ax, X, V, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, disk=True, gridsize=100):
        mean, std = getStats(X, V, Xmax=Xmax-500, Npts=450)
        if (disk):
            bla = ax.hexbin(X, V, gridsize=gridsize, bins='log', cmap='viridis')
            cb = fig.colorbar(bla)
            cb.set_label('N')
            ax.scatter(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic, marker='^', color="yellow", edgecolor='black', s=60)
        else:
            ax.scatter(X, V, c="blue", s=.1)
            ax.scatter(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic, marker='^', color="red", edgecolor='white', s=80)

        ax.plot(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic+2*std.statistic, ls='dashed', color="red", lw=3)
        ax.plot(mean.bin_edges[:-1]+(mean.bin_edges[1]- mean.bin_edges[0])/2, mean.statistic-2*std.statistic, ls='dashed', color="red", lw=3)
        ax.hlines(y=0.01, xmin=0.01, xmax=max(X), linestyles='dashdot', colors='yellow', lw=2) 
        ax.set_xlabel(Xlabel)
        ax.set_ylabel(Ylabel)
        ax.set_xlim(Xmin, Xmax) 
        ax.set_ylim(Ymin, Ymax)

    
    fig=plt.figure(1,figsize=(12,4))
    fig.subplots_adjust(wspace=0.2,hspace=0.2,top=0.97,bottom=0.1,left=0.09,right=0.95)
    plt.rcParams['font.size'] = 14

    
    ### vPhi ###
    ax=fig.add_subplot(131)
    ax.set_title(title)
    if (disk):
        Zmax = 5001
        plotPanel(ax, Z*1000, v_phi, '|Z| [pc]', '$v_{\\phi}$ [km/s]', 0, Zmax, -300, 100)
    else: 
        Zmax = 7001
        plotPanel(ax, Z*1000, v_phi, '|Z| [pc]', '$v_{\\phi}$ [km/s]', 0, Zmax, -300, 300, disk=False)

    Zg = np.linspace(0, Zmax-500, 10)

    if disk: 
        # eq. 17 from Bond+(2010)
        vPhiM = -205 + 19.2*(Zg/1000)**1.25
        sigmaM = 30 + 3.0*(Zg/1000)**2.0
        ax.plot(Zg, vPhiM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM-2*sigmaM, ls='dashdot', color="cyan", lw=1)
        ax.plot(Zg, vPhiM+2*sigmaM, ls='dashdot', color="cyan", lw=1)

        # new model
        Zg2 = np.linspace(0, 4.5, 100)
        vPhiMean, vPhiDisp = getvPhiMeanSig(Zg2)
        #ax.plot(Zg2*1000, vPhiMean-2*vPhiDisp, ls='dashed', color="green", lw=1)
        #ax.plot(Zg2*1000, vPhiMean+2*vPhiDisp, ls='dashed', color="green", lw=1)
        # print('new disk model:', np.median(vPhiMean), np.median(vPhiDisp)) 
    else:
        vPhiM = 0 + 0*(Zg/1000) 
        sigmaM = 85 + 0*(Zg/1000)
        ax.plot(Zg, vPhiM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM-2*sigmaM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM+2*sigmaM, ls='dashdot', color="cyan", lw=2)
 
        
    ### vR ###
    ax=fig.add_subplot(132)
    if (disk):
        plotPanel(ax, Z*1000, v_R, '|Z| [pc]', '$v_R$ [km/s]', 0, Zmax, -200,200)
    else:
        plotPanel(ax, Z*1000, v_r, '|Z| [pc]', '$v_r$ [km/s]', 0, Zmax, -350,350, disk=False)

    if disk:
        vM = 0*(Zg/1000)**1.25
        sigmaM = 40 + 5.0*(Zg/1000)**1.5
        ax.plot(Zg, vM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vM-2*sigmaM, ls='dashdot', color="cyan", lw=1)
        ax.plot(Zg, vM+2*sigmaM, ls='dashdot', color="cyan", lw=1)
    else:
        vPhiM = 0 + 0*(Zg/1000) 
        sigmaM = 141 + 0*(Zg/1000)
        ax.plot(Zg, vPhiM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM-2*sigmaM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM+2*sigmaM, ls='dashdot', color="cyan", lw=2)

    ### vZ ###
    ax=fig.add_subplot(133)
    if (disk):
        plotPanel(ax, Z*1000, v_Z, '|Z| [pc]', '$v_Z$ [km/s]', 0, Zmax, -200,200)
    else:
        plotPanel(ax, Z*1000, v_th, '|Z| [pc]', '$v_{\\theta}$ [km/s]', 0, Zmax, -300, 300, disk=False)

    if disk:
        vM = 0*(Zg/1000)**1.25
        sigmaM = 25 + 4.0*(Zg/1000)**1.5
        # sigmaM = 18 + 4.0*(Zg/1000)**1.5
        ax.plot(Zg, vM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vM-2*sigmaM, ls='dashdot', color="cyan", lw=1)
        ax.plot(Zg, vM+2*sigmaM, ls='dashdot', color="cyan", lw=1)
    else:
        vPhiM = 0 + 0*(Zg/1000) 
        sigmaM = 75 + 0*(Zg/1000)
        ax.plot(Zg, vPhiM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM-2*sigmaM, ls='dashdot', color="cyan", lw=2)
        ax.plot(Zg, vPhiM+2*sigmaM, ls='dashdot', color="cyan", lw=2)

    plt.tight_layout()
    plt.savefig('../plots/plot3VvsZ.png')
    print('made plot: ../plots/plot3VvsZ.png')
    plt.show() 
    plt.close("all")
    return




def plotVvsV2(vX1, vY1, tilt1, vX2, vY2, tilt2, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, title1, title2, gridsize=1000):

    fig2=plt.figure(2,figsize=(12,8))
    fig2.subplots_adjust(wspace=0.2,hspace=0.2,top=0.97,bottom=0.1,left=0.09,right=0.95)
    plt.rcParams['font.size'] = 14

    ax2=fig2.add_subplot(221)
    ax2.set_title(title1)
    bla = ax2.hexbin(vX1, vY1, gridsize=gridsize, bins='log', cmap='viridis')
    ax2.plot([Xmin, Xmax], [tilt1*Xmin, tilt1*Xmax], ls='--', c='red')
    
    ax2.set_xlabel(Xlabel)
    ax2.set_ylabel(Ylabel)
    ax2.set_xlim(Xmin, Xmax)
    ax2.set_ylim(Ymin, Ymax)

    ax2=fig2.add_subplot(222)
    ax2.set_title(title2)
    bla = ax2.hexbin(vX2, vY2, gridsize=gridsize, bins='log', cmap='viridis')
    ax2.plot([Xmin, Xmax], [tilt2*Xmin, tilt2*Xmax], ls='--', c='red')

    cb = fig2.colorbar(bla)
    cb.set_label('N')
    ax2.set_xlabel(Xlabel)
    ax2.set_ylabel(Ylabel)
    ax2.set_xlim(Xmin, Xmax)
    ax2.set_ylim(Ymin, Ymax)


    plt.tight_layout()
    plt.savefig('../plots/plotVvsV2.png')
    print('made plot: ../plots/plotVvsV2.png')
    plt.show() 
    plt.close("all")
    return 





def plotVvsVtriple(vX, vY, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, title="", gridsize=1000, giants=False):

    def plotPanel(fig, ax, vX, vY, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, title, gridsize):
        #bla = ax.hexbin(vX, vY, gridsize=gridsize, bins='log', cmap='viridis')
        bla = ax.hexbin(vX, vY, gridsize=gridsize, bins='log', cmap='viridis')
        cb = fig.colorbar(bla)
        cb.set_label('N')
        ax.set_xlabel(Xlabel)
        ax.set_ylabel(Ylabel)
        ax.set_xlim(Xmin, Xmax)
        ax.set_ylim(Ymin, Ymax)
        ax.set_title(title)
        
    fig2=plt.figure(2,figsize=(18,5.))
    fig2.subplots_adjust(wspace=0.02,hspace=0.02,top=0.97,bottom=0.1,left=0.0,right=0.99)
    plt.rcParams['font.size'] = 16

    ax2=fig2.add_subplot(131)
    plotPanel(fig2, ax2, vX[0], vY[0], Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, title[0], gridsize)
    if giants: ax2.plot([Xmin, Xmax], [0, 0], ls='--', lw=1, c='red')

    ax2=fig2.add_subplot(132)
    if giants:
        gS = gridsize
    else:
        gS = int(gridsize/10.0)
    plotPanel(fig2, ax2, vX[1], vY[1], Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, title[1], gS)
    if giants: ax2.plot([Xmin, Xmax], [0, 0], ls='--', lw=1, c='red')

    ax2=fig2.add_subplot(133)
    if giants:
        gS = int(gridsize/2)
    else:
        gS = gridsize 
    plotPanel(fig2, ax2, vX[2], vY[2], Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax, title[2], gS)
    if giants: ax2.plot([Xmin, Xmax], [0, 0], ls='--', lw=1, c='red')

    plt.tight_layout()
    plt.savefig('../plots/plotVvsVtriple.png')
    print('made plot: ../plots/plotVvsVtriple.png')
    plt.show() 
    plt.close("all")
    return 



def plotHRdiagram(df):

    # Setting up the data columns for use. 
    lum = df['lum_flame']
    logL = np.log10(lum)
    Teff = df['teff_gspphot']
    FeH = df['mh_gspphot']
    age = df['age_flame']
    red = df['ebpminrp_gspphot']

    # astroML binned statistics to set up color coding for N.
    N, xedges, yedges = binned_statistic_2d(Teff[~Teff.isna() & ~logL.isna() & ~FeH.isna()], 
                                        logL[~Teff.isna() & ~logL.isna() & ~FeH.isna()], 
                                        FeH[~Teff.isna() & ~logL.isna() & ~FeH.isna()],
                                        'count', bins=200)

    # Color coding for FeH. Gives a warning but still works in the plots.
    FeH_mean, xedges, yedges = binned_statistic_2d(Teff[~Teff.isna() & ~logL.isna() & ~FeH.isna()], 
                                               logL[~Teff.isna() & ~logL.isna() & ~FeH.isna()], 
                                               FeH[~Teff.isna() & ~logL.isna() & ~FeH.isna()],
                                               'mean', bins=200)

    # Color coding for reddening.
    red_mean, xedges, yedges = binned_statistic_2d(Teff[~Teff.isna() & ~logL.isna() & ~red.isna()], 
                                               logL[~Teff.isna() & ~logL.isna() & ~red.isna()], 
                                               red[~Teff.isna() & ~logL.isna() & ~red.isna()],
                                               'mean', bins=200)

    # Color coding for age.
    age_mean, xedges, yedges = binned_statistic_2d(Teff[~Teff.isna() & ~logL.isna() & ~age.isna()], 
                                          logL[~Teff.isna() & ~logL.isna() & ~age.isna()], 
                                          age[~Teff.isna() & ~logL.isna() & ~age.isna()],
                                               'mean', bins=200)
    # color maps  
    cmap = plt.get_cmap('viridis')
    cmap2 = plt.get_cmap('seismic')
    cmap3 = plt.get_cmap('Reds')
    cmap4 = plt.get_cmap('gist_rainbow')

    # Setting up plots for Figure 9, Panel 4 recreation.
    fig = plt.figure(figsize=(15, 12))
    fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95,
                    bottom=0.07, top=0.95)
    plt.rcParams['font.size'] = 14


    #------------------------------------------------------------
    # Create first plot color coded by N. Gives an error, but the plot still works and looks like what we want.
    plt.subplot(221)
    plt.imshow(N.T, origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', norm=matplotlib.colors.LogNorm(), interpolation='nearest', cmap=cmap)

    plt.xlim(xedges[-1]*1.05, xedges[0]*0.95)
    plt.ylim(yedges[0]*1.1, yedges[-1]*1.1)
    plt.xlabel('Effective Temperature (K)')
    plt.ylabel('log$_{10}$(Luminosity) (L$_\\odot$)')

    cb = plt.colorbar()  #ticks=[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.0, 4.5], orientation='vertical')
    cb.set_label('N')


    #------------------------------------------------------------
    # Create second plot color coded by FeH.
    plt.subplot(222)
    plt.imshow(FeH_mean.T, origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', interpolation='nearest', cmap=cmap2)

    plt.xlim(xedges[-1]*1.05, xedges[0]*0.95)
    plt.ylim(yedges[0]*1.1, yedges[-1]*1.1)
    plt.xlabel('Effective Temperature (K)')
    plt.ylabel('log$_{10}$(Luminosity) (L$_\\odot$)')

    cb = plt.colorbar()  #ticks=[-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8], format=r'$%.1f$', orientation='vertical')
    cb.set_label('Mean [Fe/H] in pixel')


    #------------------------------------------------------------
    # Create third plot color coded by reddening.
    plt.subplot(223)
    plt.imshow(red_mean.T, origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', interpolation='nearest', cmap=cmap3)

    plt.xlim(xedges[-1]*1.05, xedges[0]*0.95)
    plt.ylim(yedges[0]*1.1, yedges[-1]*1.1)
    plt.xlabel('Effective Temperature (K)')
    plt.ylabel('log$_{10}$(Luminosity) (L$_\\odot$)')

    cb = plt.colorbar()   #ticks=[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5], format=r'$%.1f$', orientation='vertical')
    cb.set_label('Mean Reddening in pixel')


    #------------------------------------------------------------
    # Create fourth plot color coded by age.
    plt.subplot(224)
    plt.imshow(age_mean.T, origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', interpolation='nearest', cmap=cmap4)

    plt.xlim(xedges[-1]*1.05, xedges[0]*0.95)
    plt.ylim(yedges[0]*1.1, yedges[-1]*1.1)
    plt.xlabel('Effective Temperature (K)')
    plt.ylabel('log$_{10}$(Luminosity) (L$_\\odot$)')

    cb = plt.colorbar()  #ticks=[0, 2, 4, 6, 8, 10, 12, 14], orientation='vertical')
    cb.set_label('Mean Age (Gyr) in pixel')

    plt.tight_layout()
    plt.savefig('../plots/HRdiagram.png')
    print('made plot: ../plots/HRdiagram.png')
    plt.show() 
    plt.close("all")
    return 




def plotVphiDistr4(Z, vPhi, vPhiMin, vPhiMax, titles):

    fig2=plt.figure(2,figsize=(12,8))
    fig2.subplots_adjust(wspace=0.2,hspace=0.2,top=0.97,bottom=0.1,left=0.09,right=0.95)
    plt.rcParams['font.size'] = 14

    def plotPanel(ax, Z, vPhi, vPhiMin, vPhiMax, title):
        vPhiGrid = np.linspace(vPhiMin, vPhiMax, 270)
        Zok = Z[(vPhi>vPhiMin)&(vPhi<vPhiMax)]
        medZ = np.median(Zok)
        p1, p2, pModel, fD1, fD2 = vPhiDistrModelNew(medZ, vPhiGrid)
        renorm = 0.92 # does ax.hist give a bias for density=True? bin edge issue?
        pModel = renorm*(fD1*p1+fD2*p2)
        ax.plot(vPhiGrid, pModel, c='b')
        ax.plot(vPhiGrid, fD1*p1, ls='dashed', c='red')
        ax.plot(vPhiGrid, fD2*p2, ls='dashed', c='green')
        ax.hist(vPhi, density=True, bins=30, histtype='stepfilled', alpha=0.2)
        ax.set_xlabel('$v_\\phi$')
        ax.set_ylabel('n (km/s)$^{-1}$')
        ax.set_xlim(vPhiMin, vPhiMax)
        plt.title(title, x=0.30, y=1.008)
        
   
    ax2=fig2.add_subplot(221)
    k=0
    plotPanel(ax2, Z[k], vPhi[k], vPhiMin, vPhiMax, titles[k])
    ax2=fig2.add_subplot(222)
    k=1
    plotPanel(ax2, Z[k], vPhi[k], vPhiMin, vPhiMax, titles[k])
    ax2=fig2.add_subplot(223)
    k=2
    plotPanel(ax2, Z[k], vPhi[k], vPhiMin, vPhiMax, titles[k])
    ax2=fig2.add_subplot(224)
    k=3
    plotPanel(ax2, Z[k], vPhi[k], vPhiMin, vPhiMax, titles[k])
        
    plt.tight_layout()
    plt.savefig('../plots/plotVphiDistr4.png')
    print('made plot: ../plots/plotVphiDistr4.png')
    plt.show() 
    plt.close("all")
    return 


def vPhiDistrModelOld(Z, vPhiGrid):

    ## component normalization
    fD2 = 0.25
    fD1 = 1 - fD2

    ## abc parameters from Table 1 in Bond+2012
    # sig1:
    a1 = 12; b1 = 1.8; c1 = 2.0
    # typo in Bond, as implied by SM code
    a1 = 22
    # sig2:
    a2 = 34; b2 = 1.2; c2 = 2.0
    # model:
    sig1 = a1 + b1*np.abs(Z)**c1
    sig2 = a2 + b2*np.abs(Z)**c2
    vn = -194.0 + 19.2 * np.abs(Z)**1.25
    mu1 = vn 
    mu2 = vn - 34.0 
    # call pdf for Gaussians
    p1 = gaussPDF(vPhiGrid, mu1, sig1) 
    p2 = gaussPDF(vPhiGrid, mu2, sig2)
    pModel = fD1*p1 + fD2*p2 
    return p1, p2, pModel, fD1, fD2



def vPhiDistrModelNew(Z, vPhiGrid):

    ## component normalization
    fD2 = 0.60
    fD1 = 1 - fD2

    ## abc parameters, with new a values 
    # sig1:
    a1 = 22; b1 = 1.8; c1 = 2.0
    # sig2:
    a2 = 17; b2 = 1.2; c2 = 2.0
    # model:
    sig1 = a1 + b1*np.abs(Z)**c1
    sig2 = a2 + b2*np.abs(Z)**c2
    vn = -186.0 + 19.2 * np.abs(Z)**1.25
    mu1 = vn 
    mu2 = vn - 38.0
    # call pdf for Gaussians
    p1 = gaussPDF(vPhiGrid, mu1, sig1) 
    p2 = gaussPDF(vPhiGrid, mu2, sig2) 
    pModel = fD1*p1 + fD2*p2 
    return p1, p2, pModel, fD1, fD2



def gaussPDF(x, mu, sig):
    return 1.0/np.sqrt(6.28)/sig*np.exp(-(x-mu)**2/2/sig**2)
    

#        mu1off   mu2off   a1     a2    vn01   vn02
# OLD     0.0     -34.0    12     34   -194    -228
# NEW     8.0     -30.0    22     17   -186    -216 
#  normalization 0.4:0.6 instead of old 0.75:0.25 

# vn = -194.0 + 19.2 * np.abs(Z)**1.25


def getStats(x,pdf):
    mean = np.sum(x*pdf)/np.sum(pdf)
    V = np.sum((x-mean)**2*pdf)/np.sum(pdf)
    return mean, np.sqrt(V)



def getvPhiMeanSig(Z):
    vPhiGrid = np.linspace(-400, 0, 400)
    vPhiMean = 0*Z
    vPhiDisp = 0*Z
    for i in range(0,np.size(Z)):
        p1, p2, pModel, fD1, fD2 = vPhiDistrModelNew(Z[i], vPhiGrid)
        vPhiMean[i], vPhiDisp[i] = getStats(vPhiGrid, pModel)
    return vPhiMean, vPhiDisp

