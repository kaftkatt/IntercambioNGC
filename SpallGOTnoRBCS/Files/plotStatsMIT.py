import pandas as pd # pandas 1.5.2
import numpy as np  # numpy 1.23.5
import matplotlib.pyplot as plt

# ----------- Functions ----------------------------
def make_dFrame(statsfile, col_names):
    '''Generate a pandas data frame from 
    MITgcm statistics text file statsfile, 
    with column names col_names (list)'''

    with open(statsfile,'r') as file:
        data=file.read().split('\n\n')
    data=[i.split('\n') for i in data]
    data_numeric = [np.array(sublist, dtype=np.float64) for sublist in data[:-1]]
    df = pd.DataFrame(data_numeric).T
    df.columns = col_names
    return(df)

def plot_stats(ax,x,ymin,ymax,ymean,ystd,varname,xlab='',ylab='', color1='darkred', color2='k', color3='midnightblue'):
    '''ax:: axes to plot in.
    x:: nparray: 1d array to plot in x axis
    ymin:: nparray: 1d array of min values of y 
    ymax:: nparray: 1d array of max values of y 
    ymean:: nparray: 1d array of mean values of y 
    ystd:: nparray: 1d array of std values of y 
    varname:: str: variable name for y
    xlab:: str: xlabel
    ylab:: str: ylabel
    '''
    #ax.plot(x,ymin,'.', color=color3)
    #ax.plot(x,ymax,'.', color=color1)
    ax.plot(x,ymin,'-', color=color3, alpha=0.3)
    ax.plot(x,ymax,'-', color=color1, alpha=0.3)
    ax.errorbar(x,ymean,yerr=ystd, capsize=4, marker='.', color=color2, alpha=0.5)
    ax.plot(x,ymean,'.', color=color2)
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.text(1.02,0.5,varname,transform=ax.transAxes, fontweight='bold', fontsize=14)
    return(ax)

def plotts(ax,time,variables,colors,ylab, errorbars=False, std=None):
        if errorbars==False:
            for var, col in zip(variables,colors):
                ax.plot(time, df[var], '-', color=col, alpha=0.5)
                ax.plot(time, df[var], '^', color=col, label=var, alpha=1)
        else:
            for var, err, col in zip(variables, std, colors):
                ax.errorbar(time,df[var],yerr=df[err], capsize=4, marker='o', color=col,label=var)
                ax.plot(time, df[var], '-', color=col, alpha=0.5)
        ax.grid(color='0.9')
        ax.set_xlim(np.min(time), np.max(time))
        ax.set_ylabel(ylab)
        return ax
    
def plot_derived(fig, ax1, ax2, ax3, ax4, ax5, df):
    axs = [ax1, ax2, ax3, ax4, ax5]
    kevars = ['pe_b_mean', 'ke_mean']
    relvars = ['vort_r_min','vort_r_max']
    vortvars = ['vort_a_mean','vort_p_mean']
    errvars = ['vort_a_sd', 'vort_p_sd']
    surfvars = ['surfExpan_theta_mean','surfExpan_salt_mean']
    
    colors = ['#1f78b4','#33a02c','#a6cee3','#b2df8a']
    time = df['time_secondsf']/3600
    
    ax1 = plotts(ax1,time, kevars,colors,'energy')
    ax2 = plotts(ax2,time, ['ke_vol'],[colors[1]],'energy')
    ax3 = plotts(ax3,time, relvars,colors,'vorticity')
    ax4 = plotts(ax4,time, vortvars,colors,'vorticity', errorbars=True, std=errvars)
    ax5 = plotts(ax5,time, surfvars,colors,' meters')
    
    ax1.set_title('MITgcm energy and vorticity monitoring')
    ax5.set_xlabel('time (hr)')
    ax1.legend(handletextpad=0)
    ax2.legend(handletextpad=0)
    ax3.legend(handletextpad=0)
    ax4.legend(handletextpad=0)
    ax5.legend(handletextpad=0)
    
    return(ax1,ax2,ax3,ax4, ax5)

def plot_cfl(fig, ax1, ax2, df):
    axs = [ax1,ax2]
    cflvar = ['u', 'v', 'w', 'W_hf']
    colors = ['coral', 'olivedrab','goldenrod','brown']
    trcfls = ['trAdv_CFL_u_max','trAdv_CFL_v_max','trAdv_CFL_w_max']
    advcfls = ['advcfl_u_max','advcfl_v_max','advcfl_w_max','advcfl_W_hf_max']
    time = df['time_secondsf']/3600
    for var, col in zip(cflvar[:-1],colors):
        ax1.plot(time, df[f'trAdv_CFL_{var}_max'], '-', color=col, alpha=0.5)
        ax1.plot(time, df[f'trAdv_CFL_{var}_max'], '^', color=col, label=var)
        ax1.grid(color='0.9')
        ax1.set_xlim(np.min(time), np.max(time))
        ax1.set_ylabel('trAdv_CFL')

    for var, col in zip(cflvar[:],colors):
        ax2.plot(time, df[f'advcfl_{var}_max'], '-', color=col, alpha=0.5)
        ax2.plot(time, df[f'advcfl_{var}_max'], '^', color=col, label=var)
        ax2.grid(color='0.9')
        ax2.set_xlim(np.min(time), np.max(time))
        ax2.set_ylabel('advcfl_CFL')
        ax2.set_xlabel('time (hr)')

    ax1.set_title('MITgcm CFL monitoring')
    ax2.set_xlabel('time (hr)')
    ax1.legend(handletextpad=0)
    ax2.legend(handletextpad=0)
    return(ax1,ax2)


# -------------- Main -----------

statsfile = 'statsMIT.txt'
col_names=['time_tsnumber','time_secondsf','eta_max', 'eta_min', 'eta_mean','eta_sd',
           'uvel_max', 'uvel_min', 'uvel_mean','uvel_sd',
           'vvel_max', 'vvel_min', 'vvel_mean','vvel_sd',
           'wvel_max', 'wvel_min', 'wvel_mean','wvel_sd',
           'theta_max', 'theta_min', 'theta_mean','theta_sd',
           'salt_max', 'salt_min', 'salt_mean','salt_sd',
           'trAdv_CFL_u_max','trAdv_CFL_v_max','trAdv_CFL_w_max',
           'advcfl_u_max','advcfl_v_max','advcfl_w_max','advcfl_W_hf_max',
           'pe_b_mean','ke_max', 'ke_mean','ke_vol','vort_r_min','vort_r_max',
           'vort_a_mean','vort_a_sd','vort_p_mean','vort_p_sd',
           'surfExpan_theta_mean','surfExpan_salt_mean',      
           ]

# Get dataframe from stats text file
df = make_dFrame(statsfile, col_names)

# Make plots dyn vars
fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1, sharex=True, figsize=(7,7))
axs = [ax1,ax2,ax3,ax4,ax5,ax6]
dynvars = ['eta','uvel','vvel','wvel','theta','salt']
ylabels = [r'$\eta$ (m)',r'$u$ (m/s)',r'$v$ (m/s)',r'$w$ (m/s)',r'$\theta$ ($^{\circ}C$)','salinity']
for ax, var, ylab in zip(axs, dynvars, ylabels):
    ax = plot_stats(ax,df['time_secondsf']/3600,df[f'{var}_min'],df[f'{var}_max'],df[f'{var}_mean'],df[f'{var}_sd'],var,ylab=ylab)
    ax.grid(color='0.9')
ax1.set_title('MITgcm dynstat monitoring')
ax6.set_xlabel('time (hr)')
plt.tight_layout()
plt.savefig('dyn_mit_stats.png', dpi=400)
plt.show(block=False)

# make plots CFL
fig2, (ax1b,ax2b) = plt.subplots(2,1, sharex=True, figsize=(7,4))
ax1b, ax2b = plot_cfl(fig2, ax1b, ax2b, df)
plt.tight_layout()
plt.savefig('cfl_mit_stats.png', dpi=400)
plt.show(block=False)

# Make plots energy, vorticity and surface expansion from T and S 
fig3, (ax1c,ax2c,ax3c,ax4c,ax5c) = plt.subplots(5,1, sharex=True, figsize=(7,7))
ax1c, ax2c, ax3c, ax4c, ax5c = plot_derived(fig3, ax1c, ax2c, ax3c, ax4c,ax5c, df)
plt.tight_layout()
plt.savefig('ek_vort_mit_stats.png', dpi=400)
plt.show()
