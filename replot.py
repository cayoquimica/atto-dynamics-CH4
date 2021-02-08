from tqdm import tqdm
from joblib import Parallel, delayed
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import matplotlib.style
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
from argparse import ArgumentParser

print('check if time in alldata has tstep substract, if not, comment line 174 ')

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

def graph_pop(t, nq2, nq1, x , y, z, r1n, r2n, tstep, t_au):
    for i in [0,1,2]:
        # graph
        outputname = "pop-st{}-{:06}.png".format(i,t)
        fig = plt.figure(figsize=(16,10))
        ax = fig.gca(projection='3d')
        #print(x.shape, y.shape, z[1,:,:].shape)
        surf = ax.plot_surface(x, y, z[i,:,:], cmap=cm.jet)
        axes = plt.gca()
        axes.set_xlim([r1n[-1],r1n[0]])
        axes.set_ylim([r2n[0],r2n[-1]])
        axes.set_zlim([0,0.004])
        axes.set_xlabel('q1')
        axes.set_ylabel('q2')
        axes.set_zlabel('Population')
        title_name = "Population through time. T = {:7.2f} fs".format(t*tstep*t_au)
        plt.title(title_name)
        ax.view_init(35, 30)
        fig.tight_layout()

        fig.savefig(outputname)
        plt.close()
def graph_pop_diff(t, nq2, nq1, x , y, z, r1n, r2n, tstep, t_au):
    for i in [0,1,2]:
        # graph
        outputname = "diff-st{}-{:06}.png".format(i,t)
        fig = plt.figure(figsize=(16,10))
        ax = fig.gca(projection='3d')
        #print(x.shape, y.shape, z[1,:,:].shape)
        surf = ax.plot_surface(x, y, z[i,:,:], cmap=cm.jet)
        axes = plt.gca()
        axes.set_xlim([r1n[-1],r1n[0]])
        axes.set_ylim([r2n[0],r2n[-1]])
        #axes.set_zlim([0,0.004])
        axes.set_xlabel('q1')
        axes.set_ylabel('q2')
        axes.set_zlabel('Population')
        title_name = "CH4+ 9th harmonic. T = {:7.2f} fs".format(t*tstep*t_au)
        plt.title(title_name)
        ax.view_init(35, 30)
        fig.tight_layout()

        fig.savefig(outputname)
        plt.close()
def dist_ang(geo):
    b1 = [bond(0,x) for x in range(1,5)]
    b2 = [ angles(i,j) for i in range(1,4) for j in range(i+1,5)]
    return ( np.concatenate((np.array(b1),np.array(b2))) )

def bond(a,b):
    return np.sqrt((geo[a,0]-geo[b,0])**2+(geo[a,1]-geo[b,1])**2+(geo[a,2]-geo[b,2])**2)

def pos4a(a):
    return [geo[a,0]-0., geo[a,1]-0., geo[a,2]-0.]

def angles(a,b):
    return 180. / np.pi * np.arccos( ( pos4a(a)[0] * pos4a(b)[0] + pos4a(a)[1] * pos4a(b)[1] + pos4a(a)[2] * pos4a(b)[2] ) / ( np.sqrt( pos4a(a)[0] ** 2 + pos4a(a)[1] ** 2 + pos4a(a)[2] ** 2 ) * np.sqrt ( pos4a(b)[0] ** 2 + pos4a(b)[1] ** 2 + pos4a(b)[2] ** 2 ) ) )

start = 1
xx = np.loadtxt('simu_parameters')
tf, nt, nf, nq1, nq2 = xx
stop = int(nf)
nq1 = int(nq1)
nq2 = int(nq2)
tstep = tf / nt * (nt / nf)      #time step between each file saved
t_au = 24.18884326505 / 1000 #conversion of time from a.u to femtoseconds
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

#e0 = np.zeros(int(nf)+1); e00 = 0.; e01 = 0.; e02=0.
#e1 = e0; e2 = e0; time = e0; pulse = e0; norm0 = e0; norm1 = e0; norm2 = e0; am0 = e0; am1 = e0; am2 = e0; lmq10 = e0; lmq20 = e0; lmq11 = e0; lmq21 = e0; lmq12 = e0; lmq22 = e0;
#for i in np.arange(1,800):#301,1):
#    file="alldata.{:04d}".format(i)
#    with open(file, 'r') as f :
#        line1 = f.readline()
#    line1 = line1.strip()
#    columns = line1.split()
#    e00 = e00 + float(columns[1])
#    e01 = e01 + float(columns[2])
#    e02 = e02 + float(columns[3])
#    alldata = np.loadtxt(file)
#    time = alldata[:,0]
#    pulse = pulse + alldata[:,1]
#    norm0 = norm0 + alldata[:,2]
#    norm1 = norm1 + alldata[:,3]
#    norm2 = norm2 + alldata[:,4]
#    e0 = e0 + alldata[:,5]
#    e1 = e1 + alldata[:,6]
#    e2 = e2 + alldata[:,7]
#    am0 = am0 + alldata[:,8]
#    am1 = am1 + alldata[:,9]
#    am2 = am2 + alldata[:,10]
#    lmq10 = lmq10 + alldata[:,11]
#    lmq20 = lmq20 + alldata[:,12]
#    lmq11 = lmq11 + alldata[:,13]
#    lmq21 = lmq21 + alldata[:,14]
#    lmq12 = lmq12 + alldata[:,15]
#    lmq22 = lmq22 + alldata[:,16]
with open('alldata', 'r') as f :
    line1 = f.readline()

line1 = line1.strip()
columns = line1.split()
e00 = float(columns[1])
e01 = float(columns[2])
e02 = float(columns[3])
alldata = np.loadtxt('alldata')
time = alldata[:,0]
pulse = alldata[:,1]
norm0 = alldata[:,2]
norm1 = alldata[:,3]
norm2 = alldata[:,4]
e0 = alldata[:,5]
e1 = alldata[:,6]
e2 = alldata[:,7]
am0 = alldata[:,8]
am1 = alldata[:,9]
am2 = alldata[:,10]
lmq10 = alldata[:,11]
lmq20 = alldata[:,12]
lmq11 = alldata[:,13]
lmq21 = alldata[:,14]
lmq12 = alldata[:,15]
lmq22 = alldata[:,16]

normalization = norm0[0]+norm1[0]+norm2[0]
norm0[:] = norm0[:] / normalization
norm1[:] = norm1[:] / normalization
norm2[:] = norm2[:] / normalization
e0[:]    =    e0[:] / normalization
e1[:]    =    e1[:] / normalization
e2[:]    =    e2[:] / normalization
am0[:]   =   am0[:] / normalization
am1[:]   =   am1[:] / normalization
am2[:]   =   am2[:] / normalization
lmq10[:] = lmq10[:] / normalization
lmq20[:] = lmq20[:] / normalization
lmq11[:] = lmq11[:] / normalization
lmq21[:] = lmq21[:] / normalization
lmq12[:] = lmq12[:] / normalization
lmq22[:] = lmq22[:] / normalization


#outputname = "norms-conservation.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, norm0+norm1+norm2-1, '-', lw=2.0)
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Norm conservation', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#fig.tight_layout()
#fig.savefig(outputname)
##fig.close()
#
#outputname = "pop-all-states.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, norm0,'k', time*t_au, norm1, 'b', time*t_au, norm2, 'r', lw=2.0)
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Population', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
##ax.legend(['st0','st1','st2'],bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) # this is to choose where the legend box is located
#ax.legend(['st0','st1','st2'],loc='center right')
#fig.tight_layout()
#fig.savefig(outputname)
##plt.close()
#
##outputname = "pop-all-states-zoom.png"
##fig, ax = plt.subplots(1,1, figsize=(6,4))
##ax.plot(time*t_au, norm0,'k', time*t_au, norm1, 'b', time*t_au, norm2, 'r', lw=2.0)
##ax.set_xlim(0,time[-1]*t_au)
##ax.set_ylabel('Population', fontsize=15)
##ax.set_xlabel('Time in femtoseconds', fontsize=15)
##ax.legend(['st0','st1','st2'])
##ax.set_xlim([0,15])
##fig.tight_layout()
##fig.savefig(outputname)
#
#outputname = "energy-conservation.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, e0+e1+e2-(e00+e01+e02), '-', lw=2.0)
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Energy in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#fig.tight_layout()
#fig.savefig(outputname)
#
#outputname = "energy-all-states.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, e0,'k', time*t_au, e1, 'b', time*t_au, e2, 'r', lw=2.0)
#ax = plt.gca()
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Energy in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#ax.legend(['st0','st1','st2'])
#fig.tight_layout()
#fig.savefig(outputname)
#
##outputname = "energy-all-states-zoom.png"
##fig, ax = plt.subplots(1,1, figsize=(6,4))
##ax.plot(time*t_au, e0,'k', time*t_au, e1, 'b', time*t_au, e2, 'r', lw=2.0)
##ax = plt.gca()
##ax.set_xlim(0,time[-1]*t_au)
##ax.set_ylabel('Energy in au.', fontsize=15)
##ax.set_xlabel('Time in femtoseconds', fontsize=15)
##ax.legend(['st0','st1','st2'])
##ax.set_xlim([0,15])
##fig.tight_layout()
##fig.savefig(outputname) 
#
#outputname = "angular-momentum.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, am0+am1+am2, '-', lw=2.0)
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Angular momentum in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#fig.tight_layout()
#fig.savefig(outputname)
#
#outputname = "angular-momentum-all-states.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, am0,'k', time*t_au, am1, 'b', time*t_au, am2, 'r', lw=2.0)
#ax = plt.gca()
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Energy in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#ax.legend(['st0','st1','st2'])
#fig.tight_layout()
#fig.savefig(outputname)
#
#outputname = "linear-momentum-q1.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, lmq10+lmq11+lmq12, '-', lw=2.0)
#ax = plt.gca()
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Linear momentum in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#ax.legend(['Along q1'])
#fig.tight_layout()
#fig.savefig(outputname)
#
#outputname = "linear-momentum-q2.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, lmq20+lmq21+lmq22, '-', lw=2.0)
#ax = plt.gca()
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Linear momentum in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#ax.legend(['Along q2'])
#fig.tight_layout()
#fig.savefig(outputname)
#
#outputname = "linear-momentum-q1-all-states.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, lmq10, 'k', time*t_au, lmq11, 'b', time*t_au, lmq12, 'r', lw=2.0)
#ax = plt.gca()
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Linear momentum in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#ax.legend(['Along q1 - st1','Along q1 - st1','Along q1 - st2'])
#fig.tight_layout()
#fig.savefig(outputname)
#
#outputname = "linear-momentum-q2-all-states.png"
#fig, ax = plt.subplots(1,1, figsize=(6,4))
#ax.plot(time*t_au, lmq20, 'k', time*t_au, lmq21, 'b', time*t_au, lmq22, 'r', lw=2.0)
#ax = plt.gca()
#ax.set_xlim(0,time[-1]*t_au)
#ax.set_ylabel('Linear momentum in au.', fontsize=15)
#ax.set_xlabel('Time in femtoseconds', fontsize=15)
#ax.legend(['Along q2 - st1','Along q2 - st1','Along q2 - st2'])
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
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
#
#ng = int(s)
#n_p = int(nf)+1
#tt = np.arange(n_p)
#time_final = int(stop)
#n_evaluations = int(nt)
#t_step = time_final / n_evaluations
#tstep_files = n_evaluations/n_p * t_step #in atomic units
#tt = tt * tstep_files * 24.18884326505 / 1000. # time axis in femtoseconds
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
#cs0c = np.zeros((n_p,11))
#cs1c = np.zeros((n_p,11))
#cs2c = np.zeros((n_p,11))
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
#cs0c[:,0] = tt
#cs1c[:,0] = tt
#cs2c[:,0] = tt
#
#popall0 = np.zeros(stop+1)
#popall1 = np.zeros(stop+1)
#popall2 = np.zeros(stop+1)
#popc2v0 = np.zeros(stop+1)
#popc2v1 = np.zeros(stop+1)
#popc2v2 = np.zeros(stop+1)
#popd2d0 = np.zeros(stop+1)
#popd2d1 = np.zeros(stop+1)
#popd2d2 = np.zeros(stop+1)
#poptd0 = np.zeros(stop+1)
#poptd1 = np.zeros(stop+1)
#poptd2 = np.zeros(stop+1)
#popcs0 = np.zeros(stop+1)
#popcs1 = np.zeros(stop+1)
#popcs2 = np.zeros(stop+1)
#popc2v0c = np.zeros(stop+1)
#popc2v1c = np.zeros(stop+1)
#popc2v2c = np.zeros(stop+1)
#popd2d0c = np.zeros(stop+1)
#popd2d1c = np.zeros(stop+1)
#popd2d2c = np.zeros(stop+1)
#posg0 = np.zeros((stop+1,2))
#posg1 = np.zeros((stop+1,2))
#posg2 = np.zeros((stop+1,2))
#posc2v0 = np.zeros((stop+1,2))
#posc2v1 = np.zeros((stop+1,2))
#posc2v2 = np.zeros((stop+1,2))
#posd2d0 = np.zeros((stop+1,2))
#posd2d1 = np.zeros((stop+1,2))
#posd2d2 = np.zeros((stop+1,2))
#postd0 = np.zeros((stop+1,2))
#postd1 = np.zeros((stop+1,2))
#postd2 = np.zeros((stop+1,2))
#poscs0 = np.zeros((stop+1,2))
#poscs1 = np.zeros((stop+1,2))
#poscs2 = np.zeros((stop+1,2))
#posg0 = np.zeros((stop+1,2))
#posg1 = np.zeros((stop+1,2))
#posg2 = np.zeros((stop+1,2))
#posc2v0c = np.zeros((stop+1,2))
#posc2v1c = np.zeros((stop+1,2))
#posc2v2c = np.zeros((stop+1,2))
#posd2d0c = np.zeros((stop+1,2))
#posd2d1c = np.zeros((stop+1,2))
#posd2d2c = np.zeros((stop+1,2))
#popnac210 = np.zeros(stop+1)
#popnac211 = np.zeros(stop+1)
#popnac212 = np.zeros(stop+1)
#posnac210 = np.zeros((stop+1,2))
#posnac211 = np.zeros((stop+1,2))
#posnac212 = np.zeros((stop+1,2))
#popcs0c = np.zeros(stop+1)
#popcs1c = np.zeros(stop+1)
#popcs2c = np.zeros(stop+1)
#poscs0c = np.zeros((stop+1,2))
#poscs1c = np.zeros((stop+1,2))
#poscs2c = np.zeros((stop+1,2))
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
##D2d closed  region
#q1id2dc = 71 
#q1fd2dc = 75 
#q2id2dc = 97 
#q2fd2dc = 101 
#xd2dc,yd2dc,n1d2dc,n2d2dc = def_subgrid(q1id2dc,q1fd2dc,q2id2dc,q2fd2dc)
##########################################################################
##Cs closed region
#q1icsc = 58
#q1fcsc = 62
#q2icsc = 131
#q2fcsc = 135
#xcsc,ycsc,n1csc,n2csc = def_subgrid(q1icsc,q1fcsc,q2icsc,q2fcsc)
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
#        xf[z] = (nq1/2 - (1+q1x[i])) * dq1 + dq1/2 #Creating position vectors in the main grid in bohr
#        yf[z] = ((1+q2y[j]) - 114 -1) * dq2 #Creating position vectors in the main grid in bohr 
#        z = z + 1
##########################################################################
#print(xf)
#print(yf)
#for t in tqdm(inputs):
# # #     parse files
#    dsname = "time-pop-{:06}".format(t)
#    fname = '{}.h5'.format(dsname)
#    file = h5.File(fname, 'r')
#    aux2 = np.asarray(file[dsname])
#    aux2 = aux2 / np.linalg.norm(aux2, ord=1)
#    z = aux2.reshape(3, int(nq2), int(nq1))
#
#    s = 1
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg01 = cg0
#    cg0 = cg0 / (np.linalg.norm(cg0, ord=1))
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#
#    popall0[t-1] = 0.0
#    posg = np.zeros(2)
#
#    for i in np.arange(nq1*nq2):
#        popall0[t-1] = popall0[t-1] + cg01[i]
#        posg = posg + cg0[i] * np.array([xf[i], yf[i]])
#
#    posg0[t-1,:] = posg
#    geo = (geo_0+posg[0]*qq1+posg[1]*qq2).reshape(5, 3)
#    fg0[t-1,1:11]=dist_ang(geo)
#    if t == 500:
#        print(posg)
#        print(fg0[t-1,:])
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popc2v0[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1c2v*n2c2v):
#        popc2v0[t-1] = popc2v0[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xc2v[i] , yc2v[i] ]) #Calculating mean position in the subgrid
#    posc2v0[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    c2v0[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popd2d0[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1d2d*n2d2d):
#        popd2d0[t-1] = popd2d0[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xd2d[i] , yd2d[i] ]) #Calculating mean position in the subgrid
#    posd2d0[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    d2d0[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #cs region
#    sub = cg[q1ics-1:q1fcs,q2ics-1:q2fcs] #Making a subgrid
#    sv = np.zeros(n1cs*n2cs)
#    sv = sub.reshape(n1cs*n2cs, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popcs0[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1cs*n2cs):
#        popcs0[t-1] = popcs0[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xcs[i] , ycs[i] ]) #Calculating mean position in the subgrid
#    poscs0[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    cs0[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #td region
#    sub = cg[q1itd-1:q1ftd,q2itd-1:q2ftd] #Making a subgrid
#    sv = np.zeros(n1td*n2td)
#    sv = sub.reshape(n1td*n2td, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    poptd0[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1td*n2td):
#        poptd0[t-1] = poptd0[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xtd[i] , ytd[i] ]) #Calculating mean position in the subgrid
#    postd0[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    td0[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #c2v closed region
#    sub = cg[q1ic2vc-1:q1fc2vc,q2ic2vc-1:q2fc2vc] #Making a subgrid
#    sv = np.zeros(n1c2vc*n2c2vc)
#    sv = sub.reshape(n1c2vc*n2c2vc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popc2v0c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1c2vc*n2c2vc):
#        popc2v0c[t-1] = popc2v0c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xc2vc[i] , yc2vc[i] ]) #Calculating mean position in the subgrid
#    posc2v0c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    c2v0c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #D2d closed region
#    sub = cg[q1id2dc-1:q1fd2dc,q2id2dc-1:q2fd2dc] #Making a subgrid
#    sv = np.zeros(n1d2dc*n2d2dc)
#    sv = sub.reshape(n1d2dc*n2d2dc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popd2d0c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1d2dc*n2d2dc):
#        popd2d0c[t-1] = popd2d0c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xd2dc[i] , yd2dc[i] ]) #Calculating mean position in the subgrid
#    posd2d0c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    d2d0c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #High NAC21 region region
#    sub = cg[q1inac21-1:q1fnac21,q2inac21-1:q2fnac21] #Making a subgrid
#    sv = np.zeros(n1nac21*n2nac21)
#    sv = sub.reshape(n1nac21*n2nac21, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popnac210[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1nac21*n2nac21):
#        popnac210[t-1] = popnac210[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xnac21[i] , ynac21[i] ]) #Calculating mean position in the subgrid
#    posnac210[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    nac210[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #cs closed region
#    sub = cg[q1icsc-1:q1fcsc,q2icsc-1:q2fcsc] #Making a subgrid
#    sv = np.zeros(n1csc*n2csc)
#    sv = sub.reshape(n1csc*n2csc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popcs0c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1csc*n2csc):
#        popcs0c[t-1] = popcs0c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xcsc[i] , ycsc[i] ]) #Calculating mean position in the subgrid
#    poscs0c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    cs0c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    s = 2
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg01 = cg0
#    cg0 = cg0 / (np.linalg.norm(cg0, ord=1))
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#
#    popall1[t-1] = 0.0
#    posg = np.zeros(2)
#
#    for i in np.arange(nq1*nq2):
#        popall1[t-1] = popall1[t-1] + cg01[i]
#        posg = posg + cg0[i] * np.array([xf[i], yf[i]])
#
#    posg1[t-1,:] = posg
#    geo = (geo_0+posg[0]*qq1+posg[1]*qq2).reshape(5, 3)
#    fg1[t-1,1:11]=dist_ang(geo)
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popc2v1[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1c2v*n2c2v):
#        popc2v1[t-1] = popc2v1[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xc2v[i] , yc2v[i] ]) #Calculating mean position in the subgrid
#    posc2v1[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    c2v1[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popd2d1[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1d2d*n2d2d):
#        popd2d1[t-1] = popd2d1[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xd2d[i] , yd2d[i] ]) #Calculating mean position in the subgrid
#    posd2d1[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    d2d1[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #cs region
#    sub = cg[q1ics-1:q1fcs,q2ics-1:q2fcs] #Making a subgrid
#    sv = np.zeros(n1cs*n2cs)
#    sv = sub.reshape(n1cs*n2cs, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popcs1[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1cs*n2cs):
#        popcs1[t-1] = popcs1[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xcs[i] , ycs[i] ]) #Calculating mean position in the subgrid
#    poscs1[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    cs1[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #td region
#    sub = cg[q1itd-1:q1ftd,q2itd-1:q2ftd] #Making a subgrid
#    sv = np.zeros(n1td*n2td)
#    sv = sub.reshape(n1td*n2td, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    poptd1[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1td*n2td):
#        poptd1[t-1] = poptd1[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xtd[i] , ytd[i] ]) #Calculating mean position in the subgrid
#    postd1[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    td1[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #c2v closed region
#    sub = cg[q1ic2vc-1:q1fc2vc,q2ic2vc-1:q2fc2vc] #Making a subgrid
#    sv = np.zeros(n1c2vc*n2c2vc)
#    sv = sub.reshape(n1c2vc*n2c2vc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popc2v1c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1c2vc*n2c2vc):
#        popc2v1c[t-1] = popc2v1c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xc2vc[i] , yc2vc[i] ]) #Calculating mean position in the subgrid
#    posc2v1c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    c2v1c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #D2d closed region
#    sub = cg[q1id2dc-1:q1fd2dc,q2id2dc-1:q2fd2dc] #Making a subgrid
#    sv = np.zeros(n1d2dc*n2d2dc)
#    sv = sub.reshape(n1d2dc*n2d2dc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popd2d1c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1d2dc*n2d2dc):
#        popd2d1c[t-1] = popd2d1c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xd2dc[i] , yd2dc[i] ]) #Calculating mean position in the subgrid
#    posd2d1c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    d2d1c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #High NAC21 region region
#    sub = cg[q1inac21-1:q1fnac21,q2inac21-1:q2fnac21] #Making a subgrid
#    sv = np.zeros(n1nac21*n2nac21)
#    sv = sub.reshape(n1nac21*n2nac21, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popnac211[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1nac21*n2nac21):
#        popnac211[t-1] = popnac211[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xnac21[i] , ynac21[i] ]) #Calculating mean position in the subgrid
#    posnac211[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    nac211[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #cs closed region
#    sub = cg[q1icsc-1:q1fcsc,q2icsc-1:q2fcsc] #Making a subgrid
#    sv = np.zeros(n1csc*n2csc)
#    sv = sub.reshape(n1csc*n2csc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popcs1c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1csc*n2csc):
#        popcs1c[t-1] = popcs1c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xcsc[i] , ycsc[i] ]) #Calculating mean position in the subgrid
#    poscs1c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    cs1c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    s = 3
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg01 = cg0
#    cg0 = cg0 / (np.linalg.norm(cg0, ord=1))
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#
#    popall2[t-1] = 0.0
#    posg = np.zeros(2)
#
#    for i in np.arange(nq1*nq2):
#        popall2[t-1] = popall2[t-1] + cg01[i]
#        posg = posg + cg0[i] * np.array([xf[i], yf[i]])
#
#    posg2[t-1,:] = posg
#    geo = (geo_0+posg[0]*qq1+posg[1]*qq2).reshape(5, 3)
#    fg2[t-1,1:11]=dist_ang(geo)
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popc2v2[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1c2v*n2c2v):
#        popc2v2[t-1] = popc2v2[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xc2v[i] , yc2v[i] ]) #Calculating mean position in the subgrid
#    posc2v2[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    c2v2[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popd2d2[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1d2d*n2d2d):
#        popd2d2[t-1] = popd2d2[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xd2d[i] , yd2d[i] ]) #Calculating mean position in the subgrid
#    posd2d2[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    d2d2[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #cs region
#    sub = cg[q1ics-1:q1fcs,q2ics-1:q2fcs] #Making a subgrid
#    sv = np.zeros(n1cs*n2cs)
#    sv = sub.reshape(n1cs*n2cs, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popcs2[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1cs*n2cs):
#        popcs2[t-1] = popcs2[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xcs[i] , ycs[i] ]) #Calculating mean position in the subgrid
#    poscs2[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    cs2[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #td region
#    sub = cg[q1itd-1:q1ftd,q2itd-1:q2ftd] #Making a subgrid
#    sv = np.zeros(n1td*n2td)
#    sv = sub.reshape(n1td*n2td, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    poptd2[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1td*n2td):
#        poptd2[t-1] = poptd2[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xtd[i] , ytd[i] ]) #Calculating mean position in the subgrid
#    postd2[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    td2[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #c2v closed region
#    sub = cg[q1ic2vc-1:q1fc2vc,q2ic2vc-1:q2fc2vc] #Making a subgrid
#    sv = np.zeros(n1c2vc*n2c2vc)
#    sv = sub.reshape(n1c2vc*n2c2vc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popc2v2c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1c2vc*n2c2vc):
#        popc2v2c[t-1] = popc2v2c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xc2vc[i] , yc2vc[i] ]) #Calculating mean position in the subgrid
#    posc2v2c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    c2v2c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #D2d closed region
#    sub = cg[q1id2dc-1:q1fd2dc,q2id2dc-1:q2fd2dc] #Making a subgrid
#    sv = np.zeros(n1d2dc*n2d2dc)
#    sv = sub.reshape(n1d2dc*n2d2dc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popd2d2c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1d2dc*n2d2dc):
#        popd2d2c[t-1] = popd2d2c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xd2dc[i] , yd2dc[i] ]) #Calculating mean position in the subgrid
#    posd2d2c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    d2d2c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #High NAC21 region region
#    sub = cg[q1inac21-1:q1fnac21,q2inac21-1:q2fnac21] #Making a subgrid
#    sv = np.zeros(n1nac21*n2nac21)
#    sv = sub.reshape(n1nac21*n2nac21, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popnac212[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1nac21*n2nac21):
#        popnac212[t-1] = popnac212[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xnac21[i] , ynac21[i] ]) #Calculating mean position in the subgrid
#    posnac212[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    nac212[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
#    #cs closed region
#    sub = cg[q1icsc-1:q1fcsc,q2icsc-1:q2fcsc] #Making a subgrid
#    sv = np.zeros(n1csc*n2csc)
#    sv = sub.reshape(n1csc*n2csc, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#
#    popcs2c[t-1] = 0.
#    pos = np.zeros(2)
#
#    for i in np.arange(n1csc*n2csc):
#        popcs2c[t-1] = popcs2c[t-1] + sv[i]
#        pos = pos + sv1[i] * np.array([xcsc[i] , ycsc[i] ]) #Calculating mean position in the subgrid
#    poscs2c[t-1,:] = pos
#    geo = (geo_0+pos[0]*qq1+pos[1]*qq2).reshape(5, 3)
#    cs2c[t-1,1:11]=dist_ang(geo)
#    del(sv,sub,sv1)
#    #########################################################################
##    dsname = "diff-pop-{:06}".format(t)
##    fname = '{}.h5'.format(dsname)
##    file = h5.File(fname, 'r')
##    aux2 = np.asarray(file[dsname])
##    z_diff = aux2.reshape(3, int(nq2), int(nq1))
##    graph_pop_diff(t, nq2, nq1, x, y, z_diff, r1n, r2n, tstep, t_au)
##    #########################################################################
#    graph_pop(t, nq2, nq1, x, y, z, r1n, r2n, tstep, t_au)
################################################################################################################################################
#
#np.savetxt('c2v-st0.txt',c2v0)
#np.savetxt('c2v-st1.txt',c2v1)
#np.savetxt('c2v-st2.txt',c2v2)
#np.savetxt('d2d-st0.txt',d2d0)
#np.savetxt('d2d-st1.txt',d2d1)
#np.savetxt('d2d-st2.txt',d2d2)
#np.savetxt('cs-st0.txt',cs0)
#np.savetxt('cs-st1.txt',cs1)
#np.savetxt('cs-st2.txt',cs2)
#np.savetxt('td-st0.txt',td0)
#np.savetxt('td-st1.txt',td1)
#np.savetxt('td-st2.txt',td2)
#np.savetxt('c2vc-st0.txt',c2v0c)
#np.savetxt('c2vc-st1.txt',c2v1c)
#np.savetxt('c2vc-st2.txt',c2v2c)
#np.savetxt('d2dc-st0.txt',d2d0c)
#np.savetxt('d2dc-st1.txt',d2d1c)
#np.savetxt('d2dc-st2.txt',d2d2c)
#np.savetxt('full-grid-st0.txt',fg0)
#np.savetxt('full-grid-st1.txt',fg1)
#np.savetxt('full-grid-st2.txt',fg2)
#np.savetxt('csc-st0.txt',cs0c)
#np.savetxt('csc-st1.txt',cs1c)
#np.savetxt('csc-st2.txt',cs2c)
#
#np.savetxt('pop-c2v-st0.txt',popc2v0)
#np.savetxt('pop-c2v-st1.txt',popc2v1)
#np.savetxt('pop-c2v-st2.txt',popc2v2)
#np.savetxt('pop-d2d-st0.txt',popd2d0)
#np.savetxt('pop-d2d-st1.txt',popd2d1)
#np.savetxt('pop-d2d-st2.txt',popd2d2)
#np.savetxt('pop-cs-st0.txt',popcs0)
#np.savetxt('pop-cs-st1.txt',popcs1)
#np.savetxt('pop-cs-st2.txt',popcs2)
#np.savetxt('pop-td-st0.txt',poptd0)
#np.savetxt('pop-td-st1.txt',poptd1)
#np.savetxt('pop-td-st2.txt',poptd2)
#np.savetxt('pop-c2vc-st0.txt',popc2v0c)
#np.savetxt('pop-c2vc-st1.txt',popc2v1c)
#np.savetxt('pop-c2vc-st2.txt',popc2v2c)
#np.savetxt('pop-d2dc-st0.txt',popd2d0c)
#np.savetxt('pop-d2dc-st1.txt',popd2d1c)
#np.savetxt('pop-d2dc-st2.txt',popd2d2c)
#np.savetxt('pop-full-grid-st0.txt',popall0)
#np.savetxt('pop-full-grid-st1.txt',popall1)
#np.savetxt('pop-full-grid-st2.txt',popall2)
#np.savetxt('pop-cs-st0.txt',popcs0)
#np.savetxt('pop-cs-st1.txt',popcs1)
#np.savetxt('pop-cs-st2.txt',popcs2)
#
#np.savetxt('pos-c2v-st0.txt',posc2v0)
#np.savetxt('pos-c2v-st1.txt',posc2v1)
#np.savetxt('pos-c2v-st2.txt',posc2v2)
#np.savetxt('pos-d2d-st0.txt',posd2d0)
#np.savetxt('pos-d2d-st1.txt',posd2d1)
#np.savetxt('pos-d2d-st2.txt',posd2d2)
#np.savetxt('pos-cs-st0.txt',poscs0)
#np.savetxt('pos-cs-st1.txt',poscs1)
#np.savetxt('pos-cs-st2.txt',poscs2)
#np.savetxt('pos-td-st0.txt',postd0)
#np.savetxt('pos-td-st1.txt',postd1)
#np.savetxt('pos-td-st2.txt',postd2)
#np.savetxt('pos-c2vc-st0.txt',posc2v0c)
#np.savetxt('pos-c2vc-st1.txt',posc2v1c)
#np.savetxt('pos-c2vc-st2.txt',posc2v2c)
#np.savetxt('pos-d2dc-st0.txt',posd2d0c)
#np.savetxt('pos-d2dc-st1.txt',posd2d1c)
#np.savetxt('pos-d2dc-st2.txt',posd2d2c)
#np.savetxt('pos-full-grid-st0.txt',posg0)
#np.savetxt('pos-full-grid-st1.txt',posg1)
#np.savetxt('pos-full-grid-st2.txt',posg2)
#np.savetxt('pos-cs-st0.txt',poscs0)
#np.savetxt('pos-cs-st1.txt',poscs1)
#np.savetxt('pos-cs-st2.txt',poscs2)
#
#np.savetxt('popnac21-st0.txt',popnac210)
#np.savetxt('popnac21-st1.txt',popnac211)
#np.savetxt('popnac21-st2.txt',popnac212)
#np.savetxt('posnac21-st0.txt',posnac210)
#np.savetxt('posnac21-st1.txt',posnac211)
#np.savetxt('posnac21-st2.txt',posnac212)

c2v0 = np.loadtxt('c2v-st0.txt')
c2v1 = np.loadtxt('c2v-st1.txt')
c2v2 = np.loadtxt('c2v-st2.txt')
d2d0 = np.loadtxt('d2d-st0.txt')
d2d1 = np.loadtxt('d2d-st1.txt')
d2d2 = np.loadtxt('d2d-st2.txt')
cs0 = np.loadtxt('cs-st0.txt')
cs1 = np.loadtxt('cs-st1.txt')
cs2 = np.loadtxt('cs-st2.txt')
td0 = np.loadtxt('td-st0.txt')
td1 = np.loadtxt('td-st1.txt')
td2 = np.loadtxt('td-st2.txt')
c2v0c = np.loadtxt('c2vc-st0.txt')
c2v1c = np.loadtxt('c2vc-st1.txt')
c2v2c = np.loadtxt('c2vc-st2.txt')
d2d0c = np.loadtxt('d2dc-st0.txt')
d2d1c = np.loadtxt('d2dc-st1.txt')
d2d2c = np.loadtxt('d2dc-st2.txt')
fg0 = np.loadtxt('full-grid-st0.txt')
fg1 = np.loadtxt('full-grid-st1.txt')
fg2 = np.loadtxt('full-grid-st2.txt')
cs0 = np.loadtxt('cs-st0.txt')
cs1 = np.loadtxt('cs-st1.txt')
cs2 = np.loadtxt('cs-st2.txt')

popc2v0 = np.loadtxt('pop-c2v-st0.txt')
popc2v1 = np.loadtxt('pop-c2v-st1.txt')
popc2v2 = np.loadtxt('pop-c2v-st2.txt')
popd2d0 = np.loadtxt('pop-d2d-st0.txt')
popd2d1 = np.loadtxt('pop-d2d-st1.txt')
popd2d2 = np.loadtxt('pop-d2d-st2.txt')
popcs0 = np.loadtxt('pop-cs-st0.txt')
popcs1 = np.loadtxt('pop-cs-st1.txt')
popcs2 = np.loadtxt('pop-cs-st2.txt')
poptd0 = np.loadtxt('pop-td-st0.txt')
poptd1 = np.loadtxt('pop-td-st1.txt')
poptd2 = np.loadtxt('pop-td-st2.txt')
popc2v0c = np.loadtxt('pop-c2vc-st0.txt')
popc2v1c = np.loadtxt('pop-c2vc-st1.txt')
popc2v2c = np.loadtxt('pop-c2vc-st2.txt')
popd2d0c = np.loadtxt('pop-d2dc-st0.txt')
popd2d1c = np.loadtxt('pop-d2dc-st1.txt')
popd2d2c = np.loadtxt('pop-d2dc-st2.txt')
popall0 = np.loadtxt('pop-full-grid-st0.txt')
popall1 = np.loadtxt('pop-full-grid-st1.txt')
popall2 = np.loadtxt('pop-full-grid-st2.txt')
popcs0 = np.loadtxt('pop-cs-st0.txt')
popcs1 = np.loadtxt('pop-cs-st1.txt')
popcs2 = np.loadtxt('pop-cs-st2.txt')

posc2v0 = np.loadtxt('pos-c2v-st0.txt')
posc2v1 = np.loadtxt('pos-c2v-st1.txt')
posc2v2 = np.loadtxt('pos-c2v-st2.txt')
posd2d0 = np.loadtxt('pos-d2d-st0.txt')
posd2d1 = np.loadtxt('pos-d2d-st1.txt')
posd2d2 = np.loadtxt('pos-d2d-st2.txt')
poscs0 = np.loadtxt('pos-cs-st0.txt')
poscs1 = np.loadtxt('pos-cs-st1.txt')
poscs2 = np.loadtxt('pos-cs-st2.txt')
postd0 = np.loadtxt('pos-td-st0.txt')
postd1 = np.loadtxt('pos-td-st1.txt')
postd2 = np.loadtxt('pos-td-st2.txt')
posc2v0c = np.loadtxt('pos-c2vc-st0.txt')
posc2v1c = np.loadtxt('pos-c2vc-st1.txt')
posc2v2c = np.loadtxt('pos-c2vc-st2.txt')
posd2d0c = np.loadtxt('pos-d2dc-st0.txt')
posd2d1c = np.loadtxt('pos-d2dc-st1.txt')
posd2d2c = np.loadtxt('pos-d2dc-st2.txt')
posg0 = np.loadtxt('pos-full-grid-st0.txt')
posg1 = np.loadtxt('pos-full-grid-st1.txt')
posg2 = np.loadtxt('pos-full-grid-st2.txt')
poscs0 = np.loadtxt('pos-cs-st0.txt')
poscs1 = np.loadtxt('pos-cs-st1.txt')
poscs2 = np.loadtxt('pos-cs-st2.txt')

popnac210 = np.loadtxt('popnac21-st0.txt')
popnac211 = np.loadtxt('popnac21-st1.txt')
popnac212 = np.loadtxt('popnac21-st2.txt')
posnac210 = np.loadtxt('posnac21-st0.txt')
posnac211 = np.loadtxt('posnac21-st1.txt')
posnac212 = np.loadtxt('posnac21-st2.txt')


outputname = "dCH1.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,1], 'k', time*t_au, d2d0[:,1], 'b', time*t_au, cs0[:,1], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Distance C-H1 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "dCH2.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,2], 'k', time*t_au, d2d0[:,2], 'b', time*t_au, cs0[:,2], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Distance C-H2 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "dCH3.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,3], 'k', time*t_au, d2d0[:,3], 'b', time*t_au, cs0[:,3], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Distance C-H3 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "dCH4.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,4], 'k', time*t_au, d2d0[:,4], 'b', time*t_au, cs0[:,4], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Distance C-H4 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H1CH2.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,5], 'k', time*t_au, d2d0[:,5], 'b', time*t_au, cs0[:,5], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Angle H1-C-H2 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H1CH3.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,6], 'k', time*t_au, d2d0[:,6], 'b', time*t_au, cs0[:,6], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Angle H1-C-H3 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H1CH4.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,7], 'k', time*t_au, d2d0[:,7], 'b', time*t_au, cs0[:,7], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Angle H1-C-H4 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H2CH3.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,8], 'k', time*t_au, d2d0[:,8], 'b', time*t_au, cs0[:,8], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Angle H2-C-H3 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H2CH4.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,9], 'k', time*t_au, d2d0[:,9], 'b', time*t_au, cs0[:,9], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Angle H2-C-H4 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H3CH4.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, c2v0[:,10], 'k', time*t_au, d2d0[:,10], 'b', time*t_au, cs0[:,10], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Angle H3-C-H4 through time for every region', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

plt.close('all')

outputname = "dCH1-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,1], 'k', time*t_au, fg1[:,1], 'b', time*t_au, fg2[:,1], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Distance C-H1 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "dCH2-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,2], 'k', time*t_au, fg1[:,2], 'b', time*t_au, fg2[:,2], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Distance C-H2 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "dCH3-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,3], 'k', time*t_au, fg1[:,3], 'b', time*t_au, fg2[:,3], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Distance C-H3 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "dCH4-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,4], 'k', time*t_au, fg1[:,4], 'b', time*t_au, fg2[:,4], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Distance C-H4 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H1CH2-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,5], 'k', time*t_au, fg1[:,5], 'b', time*t_au, fg2[:,5], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Angle H1-C-H2 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H1CH3-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,6], 'k', time*t_au, fg1[:,6], 'b', time*t_au, fg2[:,6], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Angle H1-C-H3 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H1CH4-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,7], 'k', time*t_au, fg1[:,7], 'b', time*t_au, fg2[:,7], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Angle H1-C-H4 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H2CH3-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,8], 'k', time*t_au, fg1[:,8], 'b', time*t_au, fg2[:,8], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Angle H2-C-H3 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H2CH4-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,9], 'k', time*t_au, fg1[:,9], 'b', time*t_au, fg2[:,9], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Angle H2-C-H4 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "H3CH4-fullgrid.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,10], 'k', time*t_au, fg1[:,10], 'b', time*t_au, fg2[:,10], 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degrees', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['st0','st1','st2'])
ax.set_title('Angle H3-C-H4 through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-C2v-D2d-points.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0c, 'k', time*t_au, popd2d0c, 'b', time*t_au, popcs0c, 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
#ax.set_title('Population in the C2v and D2d points through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "full-grid-distances-st0.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,1], 'k', time*t_au, fg0[:,2], 'b', time*t_au, fg0[:,3], 'r', time*t_au, fg0[:,4], 'g', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Distance in au.', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['CH1 st0','CH2 st0','CH3 st0','CH4 st0'])
#ax.set_title('Bond distance through time for ground state', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "full-grid-angles-st0.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, fg0[:,5], 'k', time*t_au, fg0[:,6], 'b', time*t_au, fg0[:,7], 'r', time*t_au, fg0[:,8], 'g', time*t_au, fg0[:,9], 'm', time*t_au, fg0[:,10], 'c', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Angle in degree', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['H1CH2 st0','H1CH3 st0','H1CH4 st0','H2CH3 st0','H2CH4 st0','H3CH4 st0'])#,bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.) # this is to choose where the legend box is located
#ax.set_title('Bond angles through time for ground state', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-C2v-D2d-Cs-regions.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0, 'k', time*t_au, popd2d0, 'b', time*t_au, popcs0, 'r', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0'])
ax.set_title('Population in the C2v, D2d and Cs regions through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-C2v-D2d-Nac21-points.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0c, 'k',   time*t_au, popd2d0c, 'b',   lw=2.0)
ax.plot(time*t_au, popc2v1c, 'k:',  time*t_au, popd2d1c, 'b:',  lw=2.0)
ax.plot(time*t_au, popc2v2c, 'k--', time*t_au, popd2d2c, 'b--', lw=2.0)
ax.plot(time*t_au, popnac210, 'r',  time*t_au, popnac211, 'r:', time*t_au, popnac212, 'r--', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','C2v st1','D2d st1','C2v st2','D2d st2','NAC21 region st0','NAC21 region st1','NAC21 region st2'])
ax.set_title('Population in the C2v and D2d points \n and the NAC21 region through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-all-regions.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0, 'k',   time*t_au, popd2d0, 'b',   time*t_au, popcs0, 'r',   time*t_au, poptd0, 'g',   lw=2.0)
ax.plot(time*t_au, popc2v1, 'k:',  time*t_au, popd2d1, 'b:',  time*t_au, popcs1, 'r:',  time*t_au, poptd1, 'g:',  lw=2.0)
ax.plot(time*t_au, popc2v2, 'k--', time*t_au, popd2d2, 'b--', time*t_au, popcs2, 'r--', time*t_au, poptd2, 'g--', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Cs st0','Td st0','C2v st1','D2d st1','Cs st1','Td st1','C2v st2','D2d st2','Cs st2','Td st2'])
ax.set_title('Population in the C2v, D2d, Td  and Cs regions through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-C2v-D2d-td-regions-no-Cs.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0, 'k',   time*t_au, popd2d0, 'b',   time*t_au, poptd0, 'g',   lw=2.0)
ax.plot(time*t_au, popc2v1, 'k:',  time*t_au, popd2d1, 'b:',  time*t_au, poptd1, 'g:',  lw=2.0)
ax.plot(time*t_au, popc2v2, 'k--', time*t_au, popd2d2, 'b--', time*t_au, poptd2, 'g--', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','Td st0','C2v st1','D2d st1','Td st1','C2v st2','D2d st2','Td st2'])
ax.set_title('Population in the C2v, D2d and Td regions through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-Td-C2v-D2d-Nac21-points.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0c, 'k',   time*t_au, popd2d0c, 'b',   lw=2.0)
ax.plot(time*t_au, popc2v1c, 'k:',  time*t_au, popd2d1c, 'b:',  lw=2.0)
ax.plot(time*t_au, popc2v2c, 'k--', time*t_au, popd2d2c, 'b--', lw=2.0)
ax.plot(time*t_au, popnac210, 'r',  time*t_au, popnac211, 'r:', time*t_au, popnac212, 'r--', lw=2.0)
ax.plot(time*t_au, poptd0, 'g',  time*t_au, poptd1, 'g:', time*t_au, poptd2, 'g--', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','C2v st1','D2d st1','C2v st2','D2d st2','NAC21 st0','NAC21 st1','NAC21 st2','Td st0','Td st1','Td st2'])
ax.set_title('Population in the C2v and D2d points \n and the NAC21 region through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)

outputname = "pop-C2v-D2d-regions.png"
fig, ax = plt.subplots(1,1, figsize=(6,4))
ax.plot(time*t_au, popc2v0, 'k',   time*t_au, popd2d0, 'b',   lw=2.0)
ax.plot(time*t_au, popc2v1, 'k:',  time*t_au, popd2d1, 'b:',  lw=2.0)
ax.plot(time*t_au, popc2v2, 'k--', time*t_au, popd2d2, 'b--', lw=2.0)
ax = plt.gca()
ax.set_xlim(0,time[-1]*t_au)
ax.set_ylabel('Population', fontsize=15)
ax.set_xlabel('Time in femtoseconds', fontsize=15)
ax.legend(['C2v st0','D2d st0','C2v st1','D2d st1','C2v st2','D2d st2'],loc='center right')
#ax.set_title('Population in the C2v and D2d points \n and the NAC21 region through time', fontsize=15)
fig.tight_layout()
fig.savefig(outputname)


plt.close('all')
