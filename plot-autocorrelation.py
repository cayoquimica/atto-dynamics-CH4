from tqdm import tqdm
from joblib import Parallel, delayed
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
from mpl_toolkits.mplot3d import axes3d
import matplotlib.style
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
from argparse import ArgumentParser

def packet(E0,t,nu0,T):
      return E0*np.cos(2*np.pi*nu0*t)*np.exp(-np.pi*t**2/T**2)

def gaussian(const,x,x0,sigma):
    gauss = const * np.exp( -(x-x0)**2/(2*sigma**2) )
    return gauss

def ft(b, a, c:int ):
# This functions performs the fourier transform of b(a). 
# It returns the fourier transform fou(fre) and the respective frequencies in cycles/a_unit
#define d(a) and d(w)
    da = a[-1] - a[-2]
    dnu = 1 / ( a[-1] - a[0] )
    fou = np.zeros(b.shape[0], dtype=complex)
#     if (a.shape[0] % 2) == 1:
    fre = np.arange(-(a.shape[0])/2,(a.shape[0])/2) * dnu
#      else:
#         fre = arange(-(a.shape[0]/2,a.shape[0]/2) * dw
# Do the integral
    for nu in (np.arange(a.shape[0])):
#         ww = nu - (a.shape[0]-1)/2
        for t in np.arange(a.shape[0]):
            fou[nu] = fou[nu] + b[t] * np.exp(-complex(0,1) * 2* np.pi * fre[nu] * a[t])*da
    return fou, fre

# Define a function that creates the subgrids
def def_subgrid(p1i,p1f,p2i,p2f):
    np1 = p1f - p1i+1
    np2 = p2f - p2i+1
    xr = np.zeros(np2*np1)
    yr = np.zeros(np2*np1)
    for j in np.arange(np2):
        for i in np.arange(np1):
            xr[j*np1+i] = (nq1/2 - np.arange(p1i,p1f+1)[i]) * dq1 + dq1/2 #Creating position vectors in the subgrid in bohr
            yr[j*np1+i] = (np.arange(p2i,p2f+1)[j] - 114 -1) * dq2 #Creating position vectors in the subgrid in bohr
    return xr,yr,np1,np2

start = 1
xx = np.loadtxt('simu_parameters')
tf, nt, nf, nq1, nq2 = xx
stop = int(nf)
nq1 = int(nq1)
nq2 = int(nq2)
tstep = tf / nt * (nt / nf)      #time step between each file saved
tf = int(tf)
t_au = 24.18884326505 / 1000 #conversion of time from a.u to femtoseconds
nf = int(nf)
s = nq1 * nq2
nst = 3
passoq1 = 0.08
passoq2 = 0.07
dq1 = passoq1
dq2 = passoq2
r1 = ( 46 - np.linspace(1,92,92) ) * passoq1 + passoq1/2
r2 = (np.linspace(1,144,144)-94-1)*passoq2
n1 = 27
n2 = 20
r1n = np.arange(r1[-1]-n1*passoq1, r1[0]+n1*passoq1, passoq1) # np.arange excludes the stop point!!!!!!!!!
r1n = np.flip(r1n) # have to flip the vector
r2n = np.arange(r2[0]-n2*passoq2, r2[-1]+n2*passoq2, passoq2) # with passoq2=0.07 and n1=20, it is going to be from -7.98 to 4.83
x, y = np.meshgrid(r1n, r2n)
#num_cores = 10
#inputs = tqdm(np.arange(start,stop,1))
inputs = np.arange(start, stop+2,1)
alldata = np.loadtxt('alldata')
time = alldata[:,0]

#qq1 = np.array([0.005468383941194,-0.005539667067822,0.099940066929625,0.352088100072972,0.315861943182830,-0.495856225504668,-0.382798634963684,-0.284702415219885,-0.495881834018862,-0.056657307983618,0.123960758550165,-0.097711469008833,0.022204480454270,-0.089107486989717,-0.101474677166113])
#qq2 = np.array([0.005503122218018,-0.005574858175503,0.000798651732665,-0.173145514483078,-0.209601800471665,0.414670087187852,0.142239888853687,0.240959271431981,0.414644315993849,-0.174783204984944,0.242514205248471,-0.417522164824840,0.140111513551201,-0.207439525935984,-0.421309279015076])
#d0 = 2.088356031 # equilibrium distance on neutral CH4 in bohr
#r0 = d0/np.sqrt(3.)
#geo_0 = np.asarray([0.,0.,0., r0,r0,r0, -r0,-r0,r0, r0,-r0,-r0, -r0,r0,-r0])
#geo_0_r = geo_0.reshape(5,3)
#h1_0 = [geo_0_r[1,0]-geo_0_r[0,0], geo_0_r[1,1]-geo_0_r[0,1], geo_0_r[1,2]-geo_0_r[0,2]]
#h2_0 = [geo_0_r[2,0]-geo_0_r[0,0], geo_0_r[2,1]-geo_0_r[0,1], geo_0_r[2,2]-geo_0_r[0,2]]
#h3_0 = [geo_0_r[3,0]-geo_0_r[0,0], geo_0_r[3,1]-geo_0_r[0,1], geo_0_r[3,2]-geo_0_r[0,2]]
#h4_0 = [geo_0_r[4,0]-geo_0_r[0,0], geo_0_r[4,1]-geo_0_r[0,1], geo_0_r[4,2]-geo_0_r[0,2]]
#H1CH2_0 = 180. / np.pi * np.arccos( ( h1_0[0] * h2_0[0] + h1_0[1] * h2_0[1] + h1_0[2] * h2_0[2] ) / ( np.sqrt( h1_0[0] ** 2 + h1_0[1] ** 2 + h1_0[2] ** 2 ) * np.sqrt ( h2_0[0] ** 2 + h2_0[1] ** 2 + h2_0[2] ** 2 ) ) )
#H1CH3_0 = 180. / np.pi * np.arccos( ( h1_0[0] * h3_0[0] + h1_0[1] * h3_0[1] + h1_0[2] * h3_0[2] ) / ( np.sqrt( h1_0[0] ** 2 + h1_0[1] ** 2 + h1_0[2] ** 2 ) * np.sqrt ( h3_0[0] ** 2 + h3_0[1] ** 2 + h3_0[2] ** 2 ) ) )
#H1CH4_0 = 180. / np.pi * np.arccos( ( h1_0[0] * h4_0[0] + h1_0[1] * h4_0[1] + h1_0[2] * h4_0[2] ) / ( np.sqrt( h1_0[0] ** 2 + h1_0[1] ** 2 + h1_0[2] ** 2 ) * np.sqrt ( h4_0[0] ** 2 + h4_0[1] ** 2 + h4_0[2] ** 2 ) ) )
#H2CH3_0 = 180. / np.pi * np.arccos( ( h2_0[0] * h3_0[0] + h2_0[1] * h3_0[1] + h2_0[2] * h3_0[2] ) / ( np.sqrt( h2_0[0] ** 2 + h2_0[1] ** 2 + h2_0[2] ** 2 ) * np.sqrt ( h3_0[0] ** 2 + h3_0[1] ** 2 + h3_0[2] ** 2 ) ) )
#H2CH4_0 = 180. / np.pi * np.arccos( ( h2_0[0] * h4_0[0] + h2_0[1] * h4_0[1] + h2_0[2] * h4_0[2] ) / ( np.sqrt( h2_0[0] ** 2 + h2_0[1] ** 2 + h2_0[2] ** 2 ) * np.sqrt ( h4_0[0] ** 2 + h4_0[1] ** 2 + h4_0[2] ** 2 ) ) )
#H3CH4_0 = 180. / np.pi * np.arccos( ( h3_0[0] * h4_0[0] + h3_0[1] * h4_0[1] + h3_0[2] * h4_0[2] ) / ( np.sqrt( h3_0[0] ** 2 + h3_0[1] ** 2 + h3_0[2] ** 2 ) * np.sqrt ( h4_0[0] ** 2 + h4_0[1] ** 2 + h4_0[2] ** 2 ) ) )

ng = int(s)
n_p = int(nf)
tt = np.arange(n_p)
time_final = int(stop)
n_evaluations = int(nt)
t_step = time_final / n_evaluations
tstep_files = n_evaluations/n_p * t_step #in atomic units
tt = tt * tstep_files * 24.18884326505 / 1000. # time axis in femtoseconds
#
#c2v0 = np.zeros((n_p,11))
#d2d0 = np.zeros((n_p,11))
#td0 = np.zeros((n_p,11))
#cs0 = np.zeros((n_p,11))
#c2v1 = np.zeros((n_p,11))
#d2d1 = np.zeros((n_p,11))
#td1 = np.zeros((n_p,11))
#cs1 = np.zeros((n_p,11))
#c2v2 = np.zeros((n_p,11))
#d2d2 = np.zeros((n_p,11))
#td2 = np.zeros((n_p,11))
#cs2 = np.zeros((n_p,11))
#fg0 = np.zeros((n_p,11))
#fg1 = np.zeros((n_p,11))
#fg2 = np.zeros((n_p,11))
#c2v0c = np.zeros((n_p,11))
#d2d0c = np.zeros((n_p,11))
#c2v1c = np.zeros((n_p,11))
#d2d1c = np.zeros((n_p,11))
#c2v2c = np.zeros((n_p,11))
#d2d2c = np.zeros((n_p,11))
#nac210 = np.zeros((n_p,11))
#nac211 = np.zeros((n_p,11))
#nac212 = np.zeros((n_p,11))
#c2v0[:,0] = tt
#d2d0[:,0] = tt
#td0[:,0] = tt
#cs0[:,0] = tt
#c2v1[:,0] = tt
#d2d1[:,0] = tt
#td1[:,0] = tt
#cs1[:,0] = tt
#c2v2[:,0] = tt
#d2d2[:,0] = tt
#td2[:,0] = tt
#cs2[:,0] = tt
#fg0[:,0] = tt
#fg1[:,0] = tt
#fg2[:,0] = tt
#c2v1c[:,0] = tt
#c2v1nc = c2v1c
#d2d1c[:,0] = tt
#d2d1nc = d2d1c
#nac210[:,0] = tt
#nac211[:,0] = tt
#nac212[:,0] = tt
#
#popall0 = np.zeros(stop)
#popall1 = np.zeros(stop)
#popall2 = np.zeros(stop)
#popc2v0 = np.zeros(stop)
#popc2v1 = np.zeros(stop)
#popc2v2 = np.zeros(stop)
#popd2d0 = np.zeros(stop)
#popd2d1 = np.zeros(stop)
#popd2d2 = np.zeros(stop)
#poptd0 = np.zeros(stop)
#poptd1 = np.zeros(stop)
#poptd2 = np.zeros(stop)
#popcs0 = np.zeros(stop)
#popcs1 = np.zeros(stop)
#popcs2 = np.zeros(stop)
#popc2v0c = np.zeros(stop)
#popc2v1c = np.zeros(stop)
#popc2v2c = np.zeros(stop)
#popd2d0c = np.zeros(stop)
#popd2d1c = np.zeros(stop)
#popd2d2c = np.zeros(stop)
#posg0 = np.zeros((stop,2))
#posg1 = np.zeros((stop,2))
#posg2 = np.zeros((stop,2))
#posc2v0 = np.zeros((stop,2))
#posc2v1 = np.zeros((stop,2))
#posc2v2 = np.zeros((stop,2))
#posd2d0 = np.zeros((stop,2))
#posd2d1 = np.zeros((stop,2))
#posd2d2 = np.zeros((stop,2))
#postd0 = np.zeros((stop,2))
#postd1 = np.zeros((stop,2))
#postd2 = np.zeros((stop,2))
#poscs0 = np.zeros((stop,2))
#poscs1 = np.zeros((stop,2))
#poscs2 = np.zeros((stop,2))
#posg0 = np.zeros((stop,2))
#posg1 = np.zeros((stop,2))
#posg2 = np.zeros((stop,2))
#posc2v0c = np.zeros((stop,2))
#posc2v1c = np.zeros((stop,2))
#posc2v2c = np.zeros((stop,2))
#posd2d0c = np.zeros((stop,2))
#posd2d1c = np.zeros((stop,2))
#posd2d2c = np.zeros((stop,2))
#popnac210 = np.zeros(stop)
#popnac211 = np.zeros(stop)
#popnac212 = np.zeros(stop)
#posnac210 = np.zeros((stop,2))
#posnac211 = np.zeros((stop,2))
#posnac212 = np.zeros((stop,2))
#
#
##Minimum C2v is in v1(64,95) or v1(91,115) in the bigger grid
##TD is in v1(46,95) or v1(73,115) in the bigger grid
##Minimum D2d is in v1(46,79) or (73,99) in the bigger grid
##Minimum Cs is in v1(33,113) or (60,133) in the bigger grid
#
##########################################################################
##c2v region
#q1ic2v = 80
#q1fc2v = 100
#q2ic2v = 108
#q2fc2v = 125
#xc2v,yc2v,n1c2v,n2c2v = def_subgrid(q1ic2v,q1fc2v,q2ic2v,q2fc2v)
##########################################################################
##D2d region
#q1id2d = 61
#q1fd2d = 85
#q2id2d = 90
#q2fd2d = 107
#xd2d,yd2d,n1d2d,n2d2d = def_subgrid(q1id2d,q1fd2d,q2id2d,q2fd2d)
##########################################################################
##Td region
#q1itd = 67 
#q1ftd = 77 
#q2itd = 110 
#q2ftd = 120 
#xtd,ytd,n1td,n2td = def_subgrid(q1itd,q1ftd,q2itd,q2ftd)
##########################################################################
##Cs region
#q1ics = 52 
#q1fcs = 70 
#q2ics = 123 
#q2fcs = 142 
#xcs,ycs,n1cs,n2cs = def_subgrid(q1ics,q1fcs,q2ics,q2fcs)
##########################################################################
##c2v closed region
#q1ic2vc = 89 
#q1fc2vc = 93 
#q2ic2vc = 113 
#q2fc2vc = 117 
#xc2vc,yc2vc,n1c2vc,n2c2vc = def_subgrid(q1ic2vc,q1fc2vc,q2ic2vc,q2fc2vc)
##########################################################################
##D2d region
#q1id2dc = 71 
#q1fd2dc = 75 
#q2id2dc = 97 
#q2fd2dc = 101 
#xd2dc,yd2dc,n1d2dc,n2d2dc = def_subgrid(q1id2dc,q1fd2dc,q2id2dc,q2fd2dc)
##########################################################################
##High NAC21 region apart from Td point
#q1inac21 = 96 
#q1fnac21 = 105 
#q2inac21 = 128 
#q2fnac21 = 135 
#xnac21,ynac21,n1nac21,n2nac21 = def_subgrid(q1inac21,q1fnac21,q2inac21,q2fnac21)
##########################################################################
#q1x = np.arange(nq1)
#q2y = np.arange(nq2)
#z = 0
#xf = np.zeros(nq1*nq2)
#yf = np.zeros(nq1*nq2)
#for j in np.arange(nq2):
#    for i in np.arange(nq1):
#        xf[z] = (nq1/2 - q1x[i]) * dq1 + dq1/2 #Creating position vectors in the main grid in bohr
#        yf[z] = (q2y[j] - 114 -1) * dq2 #Creating position vectors in the main grid in bohr 
#        z = z + 1
##########################################################################
#
auto_h   = np.loadtxt('autoc_h.txt')
auto_d   = np.loadtxt('autoc_d.txt')
diag_h   = np.loadtxt('diago_h.txt')
diag_d   = np.loadtxt('diago_d.txt')
auto_h0  = np.loadtxt('cross_h_0.txt')
auto_d0  = np.loadtxt('cross_d_0.txt')
auto_h1  = np.loadtxt('cross_h_1.txt')
auto_d1  = np.loadtxt('cross_d_1.txt')
auto_h2  = np.loadtxt('cross_h_2.txt')
auto_d2  = np.loadtxt('cross_d_2.txt')
auto_h01 = np.loadtxt('cross_h_01.txt')
auto_d01 = np.loadtxt('cross_d_01.txt')
auto_h02 = np.loadtxt('cross_h_02.txt')
auto_d02 = np.loadtxt('cross_d_02.txt')
auto_h12 = np.loadtxt('cross_h_12.txt')
auto_d12 = np.loadtxt('cross_d_12.txt')

normh = auto_h[0]
normd = auto_d[0]
auto_h   = auto_h   / normh
auto_d   = auto_d   / normd
diag_h   = diag_h   / normh
diag_d   = diag_d   / normd
auto_h0  = auto_h0  / normh
auto_d0  = auto_d0  / normd
auto_h1  = auto_h1  / normh
auto_d1  = auto_d1  / normd
auto_h2  = auto_h2  / normh
auto_d2  = auto_d2  / normd
auto_h01 = auto_h01 / normh
auto_d01 = auto_d01 / normd
auto_h02 = auto_h02 / normh
auto_d02 = auto_d02 / normd
auto_h12 = auto_h12 / normh
auto_d12 = auto_d12 / normd

np.savetxt('autoc_h.txt',   auto_h  )
np.savetxt('autoc_d.txt',   auto_d  )
np.savetxt('diago_h.txt',   diag_h  )
np.savetxt('diago_d.txt',   diag_d  )
np.savetxt('cross_h_0.txt', auto_h0 )
np.savetxt('cross_d_0.txt', auto_d0 )
np.savetxt('cross_h_1.txt', auto_h1 )
np.savetxt('cross_d_1.txt', auto_d1 )
np.savetxt('cross_h_2.txt', auto_h2 )
np.savetxt('cross_d_2.txt', auto_d2 )
np.savetxt('cross_h_01.txt',auto_h01)
np.savetxt('cross_d_01.txt',auto_d01)
np.savetxt('cross_h_02.txt',auto_h02)
np.savetxt('cross_d_02.txt',auto_d02)
np.savetxt('cross_h_12.txt',auto_h12)
np.savetxt('cross_d_12.txt',auto_d12)


outputname = 'autocorrelation_ratio.png'
fig, ax = plt.subplots(3, figsize=(12,12))
ax[2].plot(time*t_au, auto_d/auto_h, 'k', lw=1.5)
ax[2].set_ylabel(r'$\vert C(t)_{CD_4}\vert^2 / \vert C(t)_{CH_4}\vert^2$', fontsize=15)
ax[2].set_xlabel('Time / fs', fontsize=15)
ax[2].set_xlim([0,time[nf]*t_au])
ax[2].locator_params(axis='x', nbins=14)
#ax.axes.get_xaxis().set_visible(False)

ax[1].plot(time*t_au, auto_d/auto_h, 'k', lw=1.5)
ax[1].set_ylabel(r'$\vert C(t)_{CD_4}\vert^2 / \vert C(t)_{CH_4}\vert^2$', fontsize=15)
ax[1].set_xlabel('Time / fs', fontsize=15)
ax[1].set_xlim([0,4])
ax[1].set_ylim([0,8])
ax[1].locator_params(axis='x', nbins=14)
#ax.axes.get_xaxis().set_visible(False)

ax[0].plot(time*t_au, auto_d, 'r', label=r'CD$_4$', lw=1.5)
ax[0].plot(time*t_au, auto_h, 'k', label=r'CH$_4$', lw=1.5)
ax[0].set_ylabel(r'$\vert C(t)\vert^2$', fontsize=15)
ax[0].set_xlabel('Time / fs', fontsize=15)
ax[0].set_xlim([0,4])
ax[0].legend()
ax[0].locator_params(axis='x', nbins=14)
fig.tight_layout()
fig.savefig(outputname)
plt.close('all')


outputname = 'autoc_ratio_diagonal.png'
fig, ax = plt.subplots(1, figsize=(6,4))
#ax.plot(time*t_au, (diag_d+auto_d0+auto_d1+auto_d2+auto_d01+auto_d02+auto_d12-auto_d))
ax.plot(time*t_au, diag_d/diag_h,     'k', lw=1.5, label=r'$\frac{\sum_{i,g}\vert c_{i,g}^D(0) c_{i,g}^D(t)\vert^2}{\sum_{i,g}\vert c_{i,g}^H(0) c_{i,g}^H(t)\vert^2}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
#ax.set_ylim([0,10])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')


outputname = 'autoc_ratio_cross_0.png'
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(time*t_au, auto_d0/auto_h0,   'r', lw=1.5, label=r'$\frac{\sum_g\sum_{g\prime\neq g}c_{0,g}^D(0)* c_{0,g}^D(t) c_{0,g\prime }^D(0) c_{0,g\prime }^D(t)*}{\sum_g\sum_{g\prime\neq g}c_{0,g}^D(0)* c_{0,g}^D(t) c_{0,g\prime }^D(0) c_{0,g\prime }^D(t)*}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
#ax.set_ylim([0,10])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')

outputname = 'autoc_ratio_cross_1.png'
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(time*t_au, auto_d1/auto_h1,   'g', lw=1.5, label=r'$\frac{\sum_g\sum_{g\prime\neq g}c_{1,g}^D(0)* c_{1,g}^D(t) c_{1,g\prime}^D(0) c_{1,g\prime}^D(t)*}{\sum_g\sum_{g\prime\neq g}c_{1,g}^D(0)* c_{1,g}^D(t) c_{1,g\prime}^D(0) c_{1,g\prime}^D(t)*}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
#ax.set_ylim([0,10])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')

outputname = 'autoc_ratio_cross_2.png'
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(time*t_au, auto_d2/auto_h2,   'b', lw=1.5, label=r'$\frac{\sum_g\sum_{g\prime\neq g}c_{2,g}^D(0)* c_{2,g}^D(t) c_{2,g\prime}^D(0) c_{2,g\prime}^D(t)*}{\sum_g\sum_{g\prime\neq g}c_{2,g}^D(0)* c_{2,g}^D(t) c_{2,g\prime}^D(0) c_{2,g\prime}^D(t)*}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
ax.set_ylim([0,20])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')

outputname = 'autoc_ratio_cross_01.png'
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(time*t_au, auto_d01/auto_h01, 'c', lw=1.5, label=r'$\frac{\sum_g\sum_{g\prime}c_{0,g}^D(0)* c_{0,g}^D(t) c_{1,g\prime}^D(0) c_{1,g\prime}^D(t)*}{\sum_g\sum_{g\prime}c_{0,g}^D(0)* c_{0,g}^D(t) c_{1,g\prime}^D(0) c_{1,g\prime}^D(t)*}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
ax.set_ylim([0,20])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')

outputname = 'autoc_ratio_cross_02.png'
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(time*t_au, auto_d02/auto_h02, 'm', lw=1.5, label=r'$\frac{\sum_g\sum_{g\prime}c_{0,g}^D(0)* c_{0,g}^D(t) c_{2,g\prime}^D(0) c_{2,g\prime}^D(t)*}{\sum_g\sum_{g\prime}c_{0,g}^D(0)* c_{0,g}^D(t) c_{2,g\prime}^D(0) c_{2,g\prime}^D(t)*}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
ax.set_ylim([0,20])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')

outputname = 'autoc_ratio_cross_12.png'
fig, ax = plt.subplots(1, figsize=(6,4))
ax.plot(time*t_au, auto_d12/auto_h12, color='purple', lw=1.5, label=r'$\frac{\sum_g\sum_{g\prime}c_{1,g}^D(0)* c_{1,g}^D(t) c_{2,g\prime}^D(0) c_{2,g\prime}^D(t)*}{\sum_g\sum_{g\prime}c_{1,g}^D(0)* c_{1,g}^D(t) c_{2,g\prime}^D(0) c_{2,g\prime}^D(t)*}$')
ax.set_ylabel('Autocorrelation', fontsize=15)
ax.set_xlabel('Time / fs', fontsize=15)
ax.legend()
ax.set_xlim([0,time[nf]*t_au])
ax.set_ylim([0,20])
ax.locator_params(axis='x', nbins=14)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')





a = np.zeros((nf+1,7))
a[:,0] = time*t_au
a[:,1] = auto_h0
a[:,2] = auto_h1
a[:,3] = auto_h2
a[:,4] = auto_d0
a[:,5] = auto_d1
a[:,6] = auto_d2
np.savetxt('autocorrelation_states.txt',a)
