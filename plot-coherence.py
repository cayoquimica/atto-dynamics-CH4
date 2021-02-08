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
with open('alldata', 'r') as f :
    line1 = f.readline()

line1 = line1.strip()
columns = line1.split()
e00 = float(columns[1])
e01 = float(columns[2])
e02 = float(columns[3])
alldata = np.loadtxt('alldata')
time = alldata[:,0]

qq1 = np.array([0.005468383941194,-0.005539667067822,0.099940066929625,0.352088100072972,0.315861943182830,-0.495856225504668,-0.382798634963684,-0.284702415219885,-0.495881834018862,-0.056657307983618,0.123960758550165,-0.097711469008833,0.022204480454270,-0.089107486989717,-0.101474677166113])
qq2 = np.array([0.005503122218018,-0.005574858175503,0.000798651732665,-0.173145514483078,-0.209601800471665,0.414670087187852,0.142239888853687,0.240959271431981,0.414644315993849,-0.174783204984944,0.242514205248471,-0.417522164824840,0.140111513551201,-0.207439525935984,-0.421309279015076])
d0 = 2.088356031 # equilibrium distance on neutral CH4 in bohr
r0 = d0/np.sqrt(3.)
geo_0 = np.asarray([0.,0.,0., r0,r0,r0, -r0,-r0,r0, r0,-r0,-r0, -r0,r0,-r0])
geo_0_r = geo_0.reshape(5,3)
h1_0 = [geo_0_r[1,0]-geo_0_r[0,0], geo_0_r[1,1]-geo_0_r[0,1], geo_0_r[1,2]-geo_0_r[0,2]]
h2_0 = [geo_0_r[2,0]-geo_0_r[0,0], geo_0_r[2,1]-geo_0_r[0,1], geo_0_r[2,2]-geo_0_r[0,2]]
h3_0 = [geo_0_r[3,0]-geo_0_r[0,0], geo_0_r[3,1]-geo_0_r[0,1], geo_0_r[3,2]-geo_0_r[0,2]]
h4_0 = [geo_0_r[4,0]-geo_0_r[0,0], geo_0_r[4,1]-geo_0_r[0,1], geo_0_r[4,2]-geo_0_r[0,2]]
H1CH2_0 = 180. / np.pi * np.arccos( ( h1_0[0] * h2_0[0] + h1_0[1] * h2_0[1] + h1_0[2] * h2_0[2] ) / ( np.sqrt( h1_0[0] ** 2 + h1_0[1] ** 2 + h1_0[2] ** 2 ) * np.sqrt ( h2_0[0] ** 2 + h2_0[1] ** 2 + h2_0[2] ** 2 ) ) )
H1CH3_0 = 180. / np.pi * np.arccos( ( h1_0[0] * h3_0[0] + h1_0[1] * h3_0[1] + h1_0[2] * h3_0[2] ) / ( np.sqrt( h1_0[0] ** 2 + h1_0[1] ** 2 + h1_0[2] ** 2 ) * np.sqrt ( h3_0[0] ** 2 + h3_0[1] ** 2 + h3_0[2] ** 2 ) ) )
H1CH4_0 = 180. / np.pi * np.arccos( ( h1_0[0] * h4_0[0] + h1_0[1] * h4_0[1] + h1_0[2] * h4_0[2] ) / ( np.sqrt( h1_0[0] ** 2 + h1_0[1] ** 2 + h1_0[2] ** 2 ) * np.sqrt ( h4_0[0] ** 2 + h4_0[1] ** 2 + h4_0[2] ** 2 ) ) )
H2CH3_0 = 180. / np.pi * np.arccos( ( h2_0[0] * h3_0[0] + h2_0[1] * h3_0[1] + h2_0[2] * h3_0[2] ) / ( np.sqrt( h2_0[0] ** 2 + h2_0[1] ** 2 + h2_0[2] ** 2 ) * np.sqrt ( h3_0[0] ** 2 + h3_0[1] ** 2 + h3_0[2] ** 2 ) ) )
H2CH4_0 = 180. / np.pi * np.arccos( ( h2_0[0] * h4_0[0] + h2_0[1] * h4_0[1] + h2_0[2] * h4_0[2] ) / ( np.sqrt( h2_0[0] ** 2 + h2_0[1] ** 2 + h2_0[2] ** 2 ) * np.sqrt ( h4_0[0] ** 2 + h4_0[1] ** 2 + h4_0[2] ** 2 ) ) )
H3CH4_0 = 180. / np.pi * np.arccos( ( h3_0[0] * h4_0[0] + h3_0[1] * h4_0[1] + h3_0[2] * h4_0[2] ) / ( np.sqrt( h3_0[0] ** 2 + h3_0[1] ** 2 + h3_0[2] ** 2 ) * np.sqrt ( h4_0[0] ** 2 + h4_0[1] ** 2 + h4_0[2] ** 2 ) ) )

ng = int(s)
n_p = int(nf)
tt = np.arange(n_p)
time_final = int(stop)
n_evaluations = int(nt)
t_step = time_final / n_evaluations
tstep_files = n_evaluations/n_p * t_step #in atomic units
tt = tt * tstep_files * 24.18884326505 / 1000. # time axis in femtoseconds

c2v0 = np.zeros((n_p,11))
d2d0 = np.zeros((n_p,11))
td0 = np.zeros((n_p,11))
cs0 = np.zeros((n_p,11))
c2v1 = np.zeros((n_p,11))
d2d1 = np.zeros((n_p,11))
td1 = np.zeros((n_p,11))
cs1 = np.zeros((n_p,11))
c2v2 = np.zeros((n_p,11))
d2d2 = np.zeros((n_p,11))
td2 = np.zeros((n_p,11))
cs2 = np.zeros((n_p,11))
fg0 = np.zeros((n_p,11))
fg1 = np.zeros((n_p,11))
fg2 = np.zeros((n_p,11))
c2v0c = np.zeros((n_p,11))
d2d0c = np.zeros((n_p,11))
c2v1c = np.zeros((n_p,11))
d2d1c = np.zeros((n_p,11))
c2v2c = np.zeros((n_p,11))
d2d2c = np.zeros((n_p,11))
nac210 = np.zeros((n_p,11))
nac211 = np.zeros((n_p,11))
nac212 = np.zeros((n_p,11))
c2v0[:,0] = tt
d2d0[:,0] = tt
td0[:,0] = tt
cs0[:,0] = tt
c2v1[:,0] = tt
d2d1[:,0] = tt
td1[:,0] = tt
cs1[:,0] = tt
c2v2[:,0] = tt
d2d2[:,0] = tt
td2[:,0] = tt
cs2[:,0] = tt
fg0[:,0] = tt
fg1[:,0] = tt
fg2[:,0] = tt
c2v1c[:,0] = tt
c2v1nc = c2v1c
d2d1c[:,0] = tt
d2d1nc = d2d1c
nac210[:,0] = tt
nac211[:,0] = tt
nac212[:,0] = tt

popall0 = np.zeros(stop)
popall1 = np.zeros(stop)
popall2 = np.zeros(stop)
popc2v0 = np.zeros(stop)
popc2v1 = np.zeros(stop)
popc2v2 = np.zeros(stop)
popd2d0 = np.zeros(stop)
popd2d1 = np.zeros(stop)
popd2d2 = np.zeros(stop)
poptd0 = np.zeros(stop)
poptd1 = np.zeros(stop)
poptd2 = np.zeros(stop)
popcs0 = np.zeros(stop)
popcs1 = np.zeros(stop)
popcs2 = np.zeros(stop)
popc2v0c = np.zeros(stop)
popc2v1c = np.zeros(stop)
popc2v2c = np.zeros(stop)
popd2d0c = np.zeros(stop)
popd2d1c = np.zeros(stop)
popd2d2c = np.zeros(stop)
posg0 = np.zeros((stop,2))
posg1 = np.zeros((stop,2))
posg2 = np.zeros((stop,2))
posc2v0 = np.zeros((stop,2))
posc2v1 = np.zeros((stop,2))
posc2v2 = np.zeros((stop,2))
posd2d0 = np.zeros((stop,2))
posd2d1 = np.zeros((stop,2))
posd2d2 = np.zeros((stop,2))
postd0 = np.zeros((stop,2))
postd1 = np.zeros((stop,2))
postd2 = np.zeros((stop,2))
poscs0 = np.zeros((stop,2))
poscs1 = np.zeros((stop,2))
poscs2 = np.zeros((stop,2))
posg0 = np.zeros((stop,2))
posg1 = np.zeros((stop,2))
posg2 = np.zeros((stop,2))
posc2v0c = np.zeros((stop,2))
posc2v1c = np.zeros((stop,2))
posc2v2c = np.zeros((stop,2))
posd2d0c = np.zeros((stop,2))
posd2d1c = np.zeros((stop,2))
posd2d2c = np.zeros((stop,2))
popnac210 = np.zeros(stop)
popnac211 = np.zeros(stop)
popnac212 = np.zeros(stop)
posnac210 = np.zeros((stop,2))
posnac211 = np.zeros((stop,2))
posnac212 = np.zeros((stop,2))


#Minimum C2v is in v1(64,95) or v1(91,115) in the bigger grid
#TD is in v1(46,95) or v1(73,115) in the bigger grid
#Minimum D2d is in v1(46,79) or (73,99) in the bigger grid
#Minimum Cs is in v1(33,113) or (60,133) in the bigger grid

#########################################################################
#c2v region
q1ic2v = 80
q1fc2v = 100
q2ic2v = 108
q2fc2v = 125
xc2v,yc2v,n1c2v,n2c2v = def_subgrid(q1ic2v,q1fc2v,q2ic2v,q2fc2v)
#########################################################################
#D2d region
q1id2d = 61
q1fd2d = 85
q2id2d = 90
q2fd2d = 107
xd2d,yd2d,n1d2d,n2d2d = def_subgrid(q1id2d,q1fd2d,q2id2d,q2fd2d)
#########################################################################
#Td region
q1itd = 67 
q1ftd = 77 
q2itd = 110 
q2ftd = 120 
xtd,ytd,n1td,n2td = def_subgrid(q1itd,q1ftd,q2itd,q2ftd)
#########################################################################
#Cs region
q1ics = 52 
q1fcs = 70 
q2ics = 123 
q2fcs = 142 
xcs,ycs,n1cs,n2cs = def_subgrid(q1ics,q1fcs,q2ics,q2fcs)
#########################################################################
#c2v closed region
q1ic2vc = 89 
q1fc2vc = 93 
q2ic2vc = 113 
q2fc2vc = 117 
xc2vc,yc2vc,n1c2vc,n2c2vc = def_subgrid(q1ic2vc,q1fc2vc,q2ic2vc,q2fc2vc)
#########################################################################
#D2d closed region
q1id2dc = 71 
q1fd2dc = 75 
q2id2dc = 97 
q2fd2dc = 101 
xd2dc,yd2dc,n1d2dc,n2d2dc = def_subgrid(q1id2dc,q1fd2dc,q2id2dc,q2fd2dc)
#########################################################################
#High NAC21 region apart from Td point
q1inac21 = 96 
q1fnac21 = 105 
q2inac21 = 128 
q2fnac21 = 135 
xnac21,ynac21,n1nac21,n2nac21 = def_subgrid(q1inac21,q1fnac21,q2inac21,q2fnac21)
#########################################################################
q1x = np.arange(nq1)
q2y = np.arange(nq2)
z = 0
xf = np.zeros(nq1*nq2)
yf = np.zeros(nq1*nq2)
for j in np.arange(nq2):
    for i in np.arange(nq1):
        xf[z] = (nq1/2 - q1x[i]) * dq1 + dq1/2 #Creating position vectors in the main grid in bohr
        yf[z] = (q2y[j] - 114 -1) * dq2 #Creating position vectors in the main grid in bohr 
        z = z + 1
#########################################################################
#cs closed region
q1icsc = 58 
q1fcsc = 62 
q2icsc = 131 
q2fcsc = 135 
xcsc,ycsc,n1csc,n2csc = def_subgrid(q1icsc,q1fcsc,q2icsc,q2fcsc)
#########################################################################
#Cs region 9 points 
q1ics9 = 59 
q1fcs9 = 61 
q2ics9 = 132 
q2fcs9 = 134 
xcs9,ycs9,n1cs9,n2cs9 = def_subgrid(q1ics9,q1fcs9,q2ics9,q2fcs9)
#########################################################################
#c2v region 9 points
q1ic2v9 = 90 
q1fc2v9 = 92 
q2ic2v9 = 114 
q2fc2v9 = 116 
xc2v9,yc2v9,n1c2v9,n2c2v9 = def_subgrid(q1ic2v9,q1fc2v9,q2ic2v9,q2fc2v9)
#########################################################################
#D2d region 9 points
q1id2d9 = 72 
q1fd2d9 = 74 
q2id2d9 = 98 
q2fd2d9 = 101 
xd2d9,yd2d9,n1d2d9,n2d2d9 = def_subgrid(q1id2d9,q1fd2d9,q2id2d9,q2fd2d9)
#########################################################################

coh01 = np.zeros(stop+1)
coh02 = np.zeros(stop+1)
coh12 = np.zeros(stop+1)
coh01_dx = np.zeros(stop+1)
coh02_dx = np.zeros(stop+1)
coh12_dx = np.zeros(stop+1)
coh01_dy = np.zeros(stop+1)
coh02_dy = np.zeros(stop+1)
coh12_dy = np.zeros(stop+1)
coh01_dz = np.zeros(stop+1)
coh02_dz = np.zeros(stop+1)
coh12_dz = np.zeros(stop+1)
coh01_fg_c2v = np.zeros(stop+1)
coh02_fg_c2v = np.zeros(stop+1)
coh12_fg_c2v = np.zeros(stop+1)
coh01_fg_c2v_dx = np.zeros(stop+1)
coh02_fg_c2v_dx = np.zeros(stop+1)
coh12_fg_c2v_dx = np.zeros(stop+1)
coh01_fg_c2v_dy = np.zeros(stop+1)
coh02_fg_c2v_dy = np.zeros(stop+1)
coh12_fg_c2v_dy = np.zeros(stop+1)
coh01_fg_c2v_dz = np.zeros(stop+1)
coh02_fg_c2v_dz = np.zeros(stop+1)
coh12_fg_c2v_dz = np.zeros(stop+1)
coh01_fg_d2d = np.zeros(stop+1)
coh02_fg_d2d = np.zeros(stop+1)
coh12_fg_d2d = np.zeros(stop+1)
coh01_fg_d2d_dx = np.zeros(stop+1)
coh02_fg_d2d_dx = np.zeros(stop+1)
coh12_fg_d2d_dx = np.zeros(stop+1)
coh01_fg_d2d_dy = np.zeros(stop+1)
coh02_fg_d2d_dy = np.zeros(stop+1)
coh12_fg_d2d_dy = np.zeros(stop+1)
coh01_fg_d2d_dz = np.zeros(stop+1)
coh02_fg_d2d_dz = np.zeros(stop+1)
coh12_fg_d2d_dz = np.zeros(stop+1)
coh01_td = np.zeros(stop+1)
coh02_td = np.zeros(stop+1)
coh12_td = np.zeros(stop+1)
coh01_td_dx = np.zeros(stop+1)
coh02_td_dx = np.zeros(stop+1)
coh12_td_dx = np.zeros(stop+1)
coh01_td_dy = np.zeros(stop+1)
coh02_td_dy = np.zeros(stop+1)
coh12_td_dy = np.zeros(stop+1)
coh01_td_dz = np.zeros(stop+1)
coh02_td_dz = np.zeros(stop+1)
coh12_td_dz = np.zeros(stop+1)
coh01_c2v = np.zeros(stop+1)
coh02_c2v = np.zeros(stop+1)
coh12_c2v = np.zeros(stop+1)
coh01_c2v_dx = np.zeros(stop+1)
coh02_c2v_dx = np.zeros(stop+1)
coh12_c2v_dx = np.zeros(stop+1)
coh01_c2v_dy = np.zeros(stop+1)
coh02_c2v_dy = np.zeros(stop+1)
coh12_c2v_dy = np.zeros(stop+1)
coh01_c2v_dz = np.zeros(stop+1)
coh02_c2v_dz = np.zeros(stop+1)
coh12_c2v_dz = np.zeros(stop+1)
coh01_d2d = np.zeros(stop+1)
coh02_d2d = np.zeros(stop+1)
coh12_d2d = np.zeros(stop+1)
coh01_d2d_dx = np.zeros(stop+1)
coh02_d2d_dx = np.zeros(stop+1)
coh12_d2d_dx = np.zeros(stop+1)
coh01_d2d_dy = np.zeros(stop+1)
coh02_d2d_dy = np.zeros(stop+1)
coh12_d2d_dy = np.zeros(stop+1)
coh01_d2d_dz = np.zeros(stop+1)
coh02_d2d_dz = np.zeros(stop+1)
coh12_d2d_dz = np.zeros(stop+1)


pop01 = np.zeros(stop+1)
pop02 = np.zeros(stop+1)
pop12 = np.zeros(stop+1)
pop01_fg_c2v = np.zeros(stop+1)
pop02_fg_c2v = np.zeros(stop+1)
pop12_fg_c2v = np.zeros(stop+1)
pop01_fg_d2d = np.zeros(stop+1)
pop02_fg_d2d = np.zeros(stop+1)
pop12_fg_d2d = np.zeros(stop+1)
pop01_td = np.zeros(stop+1)
pop02_td = np.zeros(stop+1)
pop12_td = np.zeros(stop+1)
pop01_c2v = np.zeros(stop+1)
pop02_c2v = np.zeros(stop+1)
pop12_c2v = np.zeros(stop+1)
pop01_d2d = np.zeros(stop+1)
pop02_d2d = np.zeros(stop+1)
pop12_d2d = np.zeros(stop+1)
pop01_cs = np.zeros(stop+1)
pop02_cs = np.zeros(stop+1)
pop12_cs = np.zeros(stop+1)
pop01_c2v9 = np.zeros(stop+1)
pop02_c2v9 = np.zeros(stop+1)
pop12_c2v9 = np.zeros(stop+1)
pop01_d2d9 = np.zeros(stop+1)
pop02_d2d9 = np.zeros(stop+1)
pop12_d2d9 = np.zeros(stop+1)
pop01_cs9 = np.zeros(stop+1)
pop02_cs9 = np.zeros(stop+1)
pop12_cs9 = np.zeros(stop+1)


pop_nac01 = np.zeros(stop+1)
pop_nac02 = np.zeros(stop+1)
pop_nac12 = np.zeros(stop+1)
pop_nac01_fg_c2v = np.zeros(stop+1)
pop_nac02_fg_c2v = np.zeros(stop+1)
pop_nac12_fg_c2v = np.zeros(stop+1)
pop_nac01_fg_d2d = np.zeros(stop+1)
pop_nac02_fg_d2d = np.zeros(stop+1)
pop_nac12_fg_d2d = np.zeros(stop+1)
pop_nac01_td = np.zeros(stop+1)
pop_nac02_td = np.zeros(stop+1)
pop_nac12_td = np.zeros(stop+1)
pop_nac01_c2v = np.zeros(stop+1)
pop_nac02_c2v = np.zeros(stop+1)
pop_nac12_c2v = np.zeros(stop+1)
pop_nac01_d2d = np.zeros(stop+1)
pop_nac02_d2d = np.zeros(stop+1)
pop_nac12_d2d = np.zeros(stop+1)


tdm01x = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm21xf.txt')
tdm01y = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm21yf.txt')
tdm01z = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm21zf.txt')
tdm02x = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm31xf.txt')
tdm02y = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm31yf.txt')
tdm02z = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm31zf.txt')
tdm12x = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm32xf.txt')
tdm12y = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm32yf.txt')
tdm12z = np.loadtxt('/home/users/cayo/atto-dynamics-CH4/tdm32zf.txt')

sv1x = tdm01x.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv1y = tdm01y.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv1z = tdm01z.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv2x = tdm02x.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv2y = tdm02y.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv2z = tdm02z.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv3x = tdm12x.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv3y = tdm12y.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
sv3z = tdm12z.reshape(nq1*nq2, order = 'F') #Vectorizing the subgrig
#
#for t in tqdm(inputs):
#    dsname = "real-coh-{:06}".format(t)
#    fname = '{}.h5'.format(dsname)
#    file = h5.File(fname, 'r')
#    aux2 = np.asarray(file[dsname])
#    z = aux2.reshape(3, int(nq2), int(nq1))
#    
#    s = 1
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#    
#    coh01_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
#    coh01_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
#    coh01_td [t-1] = cg[73-1,115-1] #Making a subgrid
#    coh01_c2v_dx[t-1] = cg[95-1,115-1] * tdm01x[95-1,115-1] #Making a subgrid
#    coh01_d2d_dx[t-1] = cg[73-1,99-1]  * tdm01x[73-1,99-1] #Making a subgrid
#    coh01_td_dx [t-1] = cg[73-1,115-1] * tdm01x[73-1,115-1] #Making a subgrid
#    coh01_c2v_dy[t-1] = cg[95-1,115-1] * tdm01y[95-1,115-1] #Making a subgrid
#    coh01_d2d_dy[t-1] = cg[73-1,99-1]  * tdm01y[73-1,99-1] #Making a subgrid
#    coh01_td_dy [t-1] = cg[73-1,115-1] * tdm01y[73-1,115-1] #Making a subgrid
#    coh01_c2v_dz[t-1] = cg[95-1,115-1] * tdm01z[95-1,115-1] #Making a subgrid
#    coh01_d2d_dz[t-1] = cg[73-1,99-1]  * tdm01z[73-1,99-1] #Making a subgrid
#    coh01_td_dz [t-1] = cg[73-1,115-1] * tdm01z[73-1,115-1] #Making a subgrid
#    for i in np.arange(nq1*nq2):
#        coh01[t-1] = coh01[t-1] + cg0[i]
#        coh01_dx[t-1] = coh01_dx[t-1] + cg0[i] * sv1x[i]
#        coh01_dy[t-1] = coh01_dy[t-1] + cg0[i] * sv1y[i]
#        coh01_dz[t-1] = coh01_dz[t-1] + cg0[i] * sv1z[i]
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dx = np.zeros(n1c2v*n2c2v)
#    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dy = np.zeros(n1c2v*n2c2v)
#    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dz = np.zeros(n1c2v*n2c2v)
#    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1c2v*n2c2v):
#        coh01_fg_c2v[t-1]    = coh01_fg_c2v[t-1] + sv[i]
#        coh01_fg_c2v_dx[t-1] = coh01_fg_c2v_dx[t-1] + sv[i] * sv_dx[i]
#        coh01_fg_c2v_dy[t-1] = coh01_fg_c2v_dy[t-1] + sv[i] * sv_dy[i]
#        coh01_fg_c2v_dz[t-1] = coh01_fg_c2v_dz[t-1] + sv[i] * sv_dz[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dx = np.zeros(n1d2d*n2d2d)
#    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dy = np.zeros(n1d2d*n2d2d)
#    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dz = np.zeros(n1d2d*n2d2d)
#    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1d2d*n2d2d):
#        coh01_fg_d2d[t-1]    = coh01_fg_d2d[t-1] + sv[i]
#        coh01_fg_d2d_dx[t-1] = coh01_fg_d2d_dx[t-1] + sv[i] * sv_dx[i]
#        coh01_fg_d2d_dy[t-1] = coh01_fg_d2d_dy[t-1] + sv[i] * sv_dy[i]
#        coh01_fg_d2d_dz[t-1] = coh01_fg_d2d_dz[t-1] + sv[i] * sv_dz[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    
#    s = 2
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#    
#    coh02_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
#    coh02_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
#    coh02_td [t-1] = cg[73-1,115-1] #Making a subgrid
#    coh02_c2v_dx[t-1] = cg[95-1,115-1] * tdm02x[95-1,115-1] #Making a subgrid
#    coh02_d2d_dx[t-1] = cg[73-1,99-1]  * tdm02x[73-1,99-1] #Making a subgrid
#    coh02_td_dx [t-1] = cg[73-1,115-1] * tdm02x[73-1,115-1] #Making a subgrid
#    coh02_c2v_dy[t-1] = cg[95-1,115-1] * tdm02y[95-1,115-1] #Making a subgrid
#    coh02_d2d_dy[t-1] = cg[73-1,99-1]  * tdm02y[73-1,99-1] #Making a subgrid
#    coh02_td_dy [t-1] = cg[73-1,115-1] * tdm02y[73-1,115-1] #Making a subgrid
#    coh02_c2v_dz[t-1] = cg[95-1,115-1] * tdm02z[95-1,115-1] #Making a subgrid
#    coh02_d2d_dz[t-1] = cg[73-1,99-1]  * tdm02z[73-1,99-1] #Making a subgrid
#    coh02_td_dz [t-1] = cg[73-1,115-1] * tdm02z[73-1,115-1] #Making a subgrid
#    for i in np.arange(nq1*nq2):
#        coh02[t-1]    = coh02[t-1] + cg0[i]
#        coh02_dx[t-1] = coh02_dx[t-1] + cg0[i] * sv2x[i]
#        coh02_dy[t-1] = coh02_dy[t-1] + cg0[i] * sv2y[i]
#        coh02_dz[t-1] = coh02_dz[t-1] + cg0[i] * sv2z[i]
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm02x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dx = np.zeros(n1c2v*n2c2v)
#    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    suby = tdm02y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dy = np.zeros(n1c2v*n2c2v)
#    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    subz = tdm02z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dz = np.zeros(n1c2v*n2c2v)
#    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1c2v*n2c2v):
#        coh02_fg_c2v[t-1]    = coh02_fg_c2v[t-1] + sv[i]
#        coh02_fg_c2v_dx[t-1] = coh02_fg_c2v_dx[t-1] + sv[i] * sv_dx[i]
#        coh02_fg_c2v_dy[t-1] = coh02_fg_c2v_dy[t-1] + sv[i] * sv_dy[i]
#        coh02_fg_c2v_dz[t-1] = coh02_fg_c2v_dz[t-1] + sv[i] * sv_dz[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm02x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dx = np.zeros(n1d2d*n2d2d)
#    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    suby = tdm02y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dy = np.zeros(n1d2d*n2d2d)
#    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    subz = tdm02z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dz = np.zeros(n1d2d*n2d2d)
#    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1d2d*n2d2d):
#        coh02_fg_d2d[t-1]    = coh02_fg_d2d[t-1] + sv[i]
#        coh02_fg_d2d_dx[t-1] = coh02_fg_d2d_dx[t-1] + sv[i] * sv_dx[i]
#        coh02_fg_d2d_dy[t-1] = coh02_fg_d2d_dy[t-1] + sv[i] * sv_dy[i]
#        coh02_fg_d2d_dz[t-1] = coh02_fg_d2d_dz[t-1] + sv[i] * sv_dz[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    
#    s = 3
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#    
#    coh12_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
#    coh12_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
#    coh12_td [t-1] = cg[73-1,115-1] #Making a subgrid
#    coh12_c2v_dx[t-1] = cg[95-1,115-1] * tdm12x[95-1,115-1] #Making a subgrid
#    coh12_d2d_dx[t-1] = cg[73-1,99-1]  * tdm12x[73-1,99-1] #Making a subgrid
#    coh12_td_dx [t-1] = cg[73-1,115-1] * tdm12x[73-1,115-1] #Making a subgrid
#    coh12_c2v_dy[t-1] = cg[95-1,115-1] * tdm12y[95-1,115-1] #Making a subgrid
#    coh12_d2d_dy[t-1] = cg[73-1,99-1]  * tdm12y[73-1,99-1] #Making a subgrid
#    coh12_td_dy [t-1] = cg[73-1,115-1] * tdm12y[73-1,115-1] #Making a subgrid
#    coh12_c2v_dz[t-1] = cg[95-1,115-1] * tdm12z[95-1,115-1] #Making a subgrid
#    coh12_d2d_dz[t-1] = cg[73-1,99-1]  * tdm12z[73-1,99-1] #Making a subgrid
#    coh12_td_dz [t-1] = cg[73-1,115-1] * tdm12z[73-1,115-1] #Making a subgrid
#    for i in np.arange(nq1*nq2):
#        coh12[t-1] = coh12[t-1] + cg0[i]
#        coh12_dx[t-1] = coh12_dx[t-1] + cg0[i] * sv1x[i]
#        coh12_dy[t-1] = coh12_dy[t-1] + cg0[i] * sv1y[i]
#        coh12_dz[t-1] = coh12_dz[t-1] + cg0[i] * sv1z[i]
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm12x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dx = np.zeros(n1c2v*n2c2v)
#    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    suby = tdm12y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dy = np.zeros(n1c2v*n2c2v)
#    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    subz = tdm12z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dz = np.zeros(n1c2v*n2c2v)
#    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1c2v*n2c2v):
#        coh12_fg_c2v[t-1]    = coh12_fg_c2v[t-1] + sv[i]
#        coh12_fg_c2v_dx[t-1] = coh12_fg_c2v_dx[t-1] + sv[i] * sv_dx[i]
#        coh12_fg_c2v_dy[t-1] = coh12_fg_c2v_dy[t-1] + sv[i] * sv_dy[i]
#        coh12_fg_c2v_dz[t-1] = coh12_fg_c2v_dz[t-1] + sv[i] * sv_dz[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm12x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dx = np.zeros(n1d2d*n2d2d)
#    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    suby = tdm12y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dy = np.zeros(n1d2d*n2d2d)
#    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    subz = tdm12z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dz = np.zeros(n1d2d*n2d2d)
#    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1d2d*n2d2d):
#        coh12_fg_d2d[t-1]    = coh12_fg_d2d[t-1] + sv[i]
#        coh12_fg_d2d_dx[t-1] = coh12_fg_d2d_dx[t-1] + sv[i] * sv_dx[i]
#        coh12_fg_d2d_dy[t-1] = coh12_fg_d2d_dy[t-1] + sv[i] * sv_dy[i]
#        coh12_fg_d2d_dz[t-1] = coh12_fg_d2d_dz[t-1] + sv[i] * sv_dz[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    ######################################################################### 
#
#np.savetxt('coh-fg-01.dat',coh01)
#np.savetxt('coh-fg-02.dat',coh02)
#np.savetxt('coh-fg-12.dat',coh12)
#np.savetxt('coh-fg-dx-01.dat',coh01_dx)
#np.savetxt('coh-fg-dx-02.dat',coh02_dx)
#np.savetxt('coh-fg-dx-12.dat',coh12_dx)
#np.savetxt('coh-fg-dy-01.dat',coh01_dy)
#np.savetxt('coh-fg-dy-02.dat',coh02_dy)
#np.savetxt('coh-fg-dy-12.dat',coh12_dy)
#np.savetxt('coh-fg-dz-01.dat',coh01_dz)
#np.savetxt('coh-fg-dz-02.dat',coh02_dz)
#np.savetxt('coh-fg-dz-12.dat',coh12_dz)
#np.savetxt('coh-c2v-region-01.dat',coh01_fg_c2v)
#np.savetxt('coh-c2v-region-02.dat',coh02_fg_c2v)
#np.savetxt('coh-c2v-region-12.dat',coh12_fg_c2v)
#np.savetxt('coh-d2d-region-01.dat',coh01_fg_d2d)
#np.savetxt('coh-d2d-region-02.dat',coh02_fg_d2d)
#np.savetxt('coh-d2d-region-12.dat',coh12_fg_d2d)
#np.savetxt('coh-c2v-region-dx-01.dat',coh01_fg_c2v_dx)
#np.savetxt('coh-c2v-region-dx-02.dat',coh02_fg_c2v_dx)
#np.savetxt('coh-c2v-region-dx-12.dat',coh12_fg_c2v_dx)
#np.savetxt('coh-d2d-region-dx-01.dat',coh01_fg_d2d_dx)
#np.savetxt('coh-d2d-region-dx-02.dat',coh02_fg_d2d_dx)
#np.savetxt('coh-d2d-region-dx-12.dat',coh12_fg_d2d_dx)
#np.savetxt('coh-c2v-region-dy-01.dat',coh01_fg_c2v_dy)
#np.savetxt('coh-c2v-region-dy-02.dat',coh02_fg_c2v_dy)
#np.savetxt('coh-c2v-region-dy-12.dat',coh12_fg_c2v_dy)
#np.savetxt('coh-d2d-region-dy-01.dat',coh01_fg_d2d_dy)
#np.savetxt('coh-d2d-region-dy-02.dat',coh02_fg_d2d_dy)
#np.savetxt('coh-d2d-region-dy-12.dat',coh12_fg_d2d_dy)
#np.savetxt('coh-c2v-region-dz-01.dat',coh01_fg_c2v_dz)
#np.savetxt('coh-c2v-region-dz-02.dat',coh02_fg_c2v_dz)
#np.savetxt('coh-c2v-region-dz-12.dat',coh12_fg_c2v_dz)
#np.savetxt('coh-d2d-region-dz-01.dat',coh01_fg_d2d_dz)
#np.savetxt('coh-d2d-region-dz-02.dat',coh02_fg_d2d_dz)
#np.savetxt('coh-d2d-region-dz-12.dat',coh12_fg_d2d_dz)
#np.savetxt('coh-c2v-point-01.dat',coh01_c2v)
#np.savetxt('coh-c2v-point-02.dat',coh02_c2v)
#np.savetxt('coh-c2v-point-12.dat',coh12_c2v)
#np.savetxt('coh-d2d-point-01.dat',coh01_d2d)
#np.savetxt('coh-d2d-point-02.dat',coh02_d2d)
#np.savetxt('coh-d2d-point-12.dat',coh12_d2d)
#np.savetxt('coh-td-point-01.dat',coh01_td)
#np.savetxt('coh-td-point-02.dat',coh02_td)
#np.savetxt('coh-td-point-12.dat',coh12_td)
#np.savetxt('coh-c2v-point-01-dx.dat',coh01_c2v_dx)
#np.savetxt('coh-c2v-point-02-dx.dat',coh02_c2v_dx)
#np.savetxt('coh-c2v-point-12-dx.dat',coh12_c2v_dx)
#np.savetxt('coh-d2d-point-01-dx.dat',coh01_d2d_dx)
#np.savetxt('coh-d2d-point-02-dx.dat',coh02_d2d_dx)
#np.savetxt('coh-d2d-point-12-dx.dat',coh12_d2d_dx)
#np.savetxt('coh-td-point-01-dx.dat',coh01_td_dx)
#np.savetxt('coh-td-point-02-dx.dat',coh02_td_dx)
#np.savetxt('coh-td-point-12-dx.dat',coh12_td_dx)
#np.savetxt('coh-c2v-point-01-dy.dat',coh01_c2v_dy)
#np.savetxt('coh-c2v-point-02-dy.dat',coh02_c2v_dy)
#np.savetxt('coh-c2v-point-12-dy.dat',coh12_c2v_dy)
#np.savetxt('coh-d2d-point-01-dy.dat',coh01_d2d_dy)
#np.savetxt('coh-d2d-point-02-dy.dat',coh02_d2d_dy)
#np.savetxt('coh-d2d-point-12-dy.dat',coh12_d2d_dy)
#np.savetxt('coh-td-point-01-dy.dat',coh01_td_dy)
#np.savetxt('coh-td-point-02-dy.dat',coh02_td_dy)
#np.savetxt('coh-td-point-12-dy.dat',coh12_td_dy)
#np.savetxt('coh-c2v-point-01-dz.dat',coh01_c2v_dz)
#np.savetxt('coh-c2v-point-02-dz.dat',coh02_c2v_dz)
#np.savetxt('coh-c2v-point-12-dz.dat',coh12_c2v_dz)
#np.savetxt('coh-d2d-point-01-dz.dat',coh01_d2d_dz)
#np.savetxt('coh-d2d-point-02-dz.dat',coh02_d2d_dz)
#np.savetxt('coh-d2d-point-12-dz.dat',coh12_d2d_dz)
#np.savetxt('coh-td-point-01-dz.dat',coh01_td_dz)
#np.savetxt('coh-td-point-02-dz.dat',coh02_td_dz)
#np.savetxt('coh-td-point-12-dz.dat',coh12_td_dz)
#
#
#outputname = 'coherence_on_C2v_point.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_c2v, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_c2v, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_c2v, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_point.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_d2d, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_d2d, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_d2d, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_Td_point.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_td, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_td, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_td, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_C2v_point_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_c2v_dx, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_c2v_dx, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_c2v_dx, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_point_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_d2d_dx, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_d2d_dx, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_d2d_dx, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#outputname = 'coherence_on_Td_point_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_td_dx, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_td_dx, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_td_dx, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_C2v_point_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_c2v_dy, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_c2v_dy, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_c2v_dy, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_point_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_d2d_dy, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_d2d_dy, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_d2d_dy, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#outputname = 'coherence_on_Td_point_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_td_dy, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_td_dy, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_td_dy, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'coherence_on_C2v_point_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_c2v_dz, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_c2v_dz, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_c2v_dz, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_point_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_d2d_dz, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_d2d_dz, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_d2d_dz, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#outputname = 'coherence_on_Td_point_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_td_dz, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_td_dz, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_td_dz, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#
#outputname = 'coherence_full_grid.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_full_grid_dipx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_dx, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_dx, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel(r'Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_dx, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#outputname = 'coherence_full_grid_dipy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_dy, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_dy, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_dy, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#outputname = 'coherence_full_grid_dipz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_dz, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_dz, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_dz, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_C2v_region.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_c2v, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_c2v, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_c2v, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_region.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_d2d, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_d2d, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_d2d, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'coherence_on_C2v_region_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_c2v_dx, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_c2v_dx, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_c2v_dx, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_region_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_d2d_dx, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_d2d_dx, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_d2d_dx, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#
#outputname = 'coherence_on_C2v_region_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_c2v_dy, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_c2v_dy, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_c2v_dy, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_region_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_d2d_dy, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_d2d_dy, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_d2d_dy, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#
#
#outputname = 'coherence_on_C2v_region_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_c2v_dz, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_c2v_dz, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_c2v_dz, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'coherence_on_D2d_region_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, coh01_fg_d2d_dz, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, coh02_fg_d2d_dz, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, coh12_fg_d2d_dz, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#time_m = np.arange(-tf,tf+tstep,tstep)
#coh01_fg_c2v_dx_m             = np.zeros(2*nf+1)
#coh01_fg_c2v_dy_m             = np.zeros(2*nf+1)
#coh01_fg_c2v_dz_m             = np.zeros(2*nf+1)
#coh01_fg_d2d_dx_m             = np.zeros(2*nf+1)
#coh01_fg_d2d_dy_m             = np.zeros(2*nf+1)
#coh01_fg_d2d_dz_m             = np.zeros(2*nf+1)
#coh01_c2v_dx_m             = np.zeros(2*nf+1)
#coh01_c2v_dy_m             = np.zeros(2*nf+1)
#coh01_c2v_dz_m             = np.zeros(2*nf+1)
#coh01_d2d_dx_m             = np.zeros(2*nf+1)
#coh01_d2d_dy_m             = np.zeros(2*nf+1)
#coh01_d2d_dz_m             = np.zeros(2*nf+1)
#coh01_td_dx_m             = np.zeros(2*nf+1)
#coh01_td_dy_m             = np.zeros(2*nf+1)
#coh01_td_dz_m             = np.zeros(2*nf+1)
#coh01_dx_m             = np.zeros(2*nf+1)
#coh01_dy_m             = np.zeros(2*nf+1)
#coh01_dz_m             = np.zeros(2*nf+1)
#coh01_fg_c2v_dx_m[ 0:nf]      = np.flip(coh01_fg_c2v_dx[1:nf+1])
#coh01_fg_c2v_dx_m[nf :2*nf+1] =         coh01_fg_c2v_dx
#coh01_fg_c2v_dy_m[ 0:nf]      = np.flip(coh01_fg_c2v_dy[1:nf+1])
#coh01_fg_c2v_dy_m[nf:2*nf+1]  =         coh01_fg_c2v_dy
#coh01_fg_c2v_dz_m[ 0:nf]      = np.flip(coh01_fg_c2v_dz[1:nf+1])
#coh01_fg_c2v_dz_m[nf:2*nf+1]  =         coh01_fg_c2v_dz
#coh01_fg_d2d_dx_m[ 0:nf]      = np.flip(coh01_fg_d2d_dx[1:nf+1])
#coh01_fg_d2d_dx_m[nf :2*nf+1] =         coh01_fg_d2d_dx
#coh01_fg_d2d_dy_m[ 0:nf]      = np.flip(coh01_fg_d2d_dy[1:nf+1])
#coh01_fg_d2d_dy_m[nf:2*nf+1]  =         coh01_fg_d2d_dy
#coh01_fg_d2d_dz_m[ 0:nf]      = np.flip(coh01_fg_d2d_dz[1:nf+1])
#coh01_fg_d2d_dz_m[nf:2*nf+1]  =         coh01_fg_d2d_dz
#coh01_c2v_dx_m[ 0:nf]         = np.flip(coh01_c2v_dx[1:nf+1])
#coh01_c2v_dx_m[nf :2*nf+1]    =         coh01_c2v_dx
#coh01_c2v_dy_m[ 0:nf]         = np.flip(coh01_c2v_dy[1:nf+1])
#coh01_c2v_dy_m[nf:2*nf+1]     =         coh01_c2v_dy
#coh01_c2v_dz_m[ 0:nf]         = np.flip(coh01_c2v_dz[1:nf+1])
#coh01_c2v_dz_m[nf:2*nf+1]     =         coh01_c2v_dz
#coh01_d2d_dx_m[ 0:nf]         = np.flip(coh01_d2d_dx[1:nf+1])
#coh01_d2d_dx_m[nf :2*nf+1]    =         coh01_d2d_dx
#coh01_d2d_dy_m[ 0:nf]         = np.flip(coh01_d2d_dy[1:nf+1])
#coh01_d2d_dy_m[nf:2*nf+1]     =         coh01_d2d_dy
#coh01_d2d_dz_m[ 0:nf]         = np.flip(coh01_d2d_dz[1:nf+1])
#coh01_d2d_dz_m[nf:2*nf+1]     =         coh01_d2d_dz
#coh01_td_dx_m[ 0:nf]          = np.flip(coh01_td_dx[1:nf+1])
#coh01_td_dx_m[nf :2*nf+1]     =         coh01_td_dx
#coh01_td_dy_m[ 0:nf]          = np.flip(coh01_td_dy[1:nf+1])
#coh01_td_dy_m[nf:2*nf+1]      =         coh01_td_dy
#coh01_td_dz_m[ 0:nf]          = np.flip(coh01_td_dz[1:nf+1])
#coh01_td_dz_m[nf:2*nf+1]      =         coh01_td_dz
#coh01_dx_m[ 0:nf]             = np.flip(coh01_dx[1:nf+1])
#coh01_dx_m[nf :2*nf+1]        =         coh01_dx
#coh01_dy_m[ 0:nf]             = np.flip(coh01_dy[1:nf+1])
#coh01_dy_m[nf:2*nf+1]         =         coh01_dy
#coh01_dz_m[ 0:nf]             = np.flip(coh01_dz[1:nf+1])
#coh01_dz_m[nf:2*nf+1]         =         coh01_dz
#
#coh02_fg_c2v_dx_m             = np.zeros(2*nf+1)
#coh02_fg_c2v_dy_m             = np.zeros(2*nf+1)
#coh02_fg_c2v_dz_m             = np.zeros(2*nf+1)
#coh02_fg_d2d_dx_m             = np.zeros(2*nf+1)
#coh02_fg_d2d_dy_m             = np.zeros(2*nf+1)
#coh02_fg_d2d_dz_m             = np.zeros(2*nf+1)
#coh02_c2v_dx_m             = np.zeros(2*nf+1)
#coh02_c2v_dy_m             = np.zeros(2*nf+1)
#coh02_c2v_dz_m             = np.zeros(2*nf+1)
#coh02_d2d_dx_m             = np.zeros(2*nf+1)
#coh02_d2d_dy_m             = np.zeros(2*nf+1)
#coh02_d2d_dz_m             = np.zeros(2*nf+1)
#coh02_td_dx_m             = np.zeros(2*nf+1)
#coh02_td_dy_m             = np.zeros(2*nf+1)
#coh02_td_dz_m             = np.zeros(2*nf+1)
#coh02_dx_m             = np.zeros(2*nf+1)
#coh02_dy_m             = np.zeros(2*nf+1)
#coh02_dz_m             = np.zeros(2*nf+1)
#coh02_fg_c2v_dx_m[ 0:nf]      = np.flip(coh02_fg_c2v_dx[1:nf+1])
#coh02_fg_c2v_dx_m[nf :2*nf+1] =         coh02_fg_c2v_dx
#coh02_fg_c2v_dy_m[ 0:nf]      = np.flip(coh02_fg_c2v_dy[1:nf+1])
#coh02_fg_c2v_dy_m[nf:2*nf+1]  =         coh02_fg_c2v_dy
#coh02_fg_c2v_dz_m[ 0:nf]      = np.flip(coh02_fg_c2v_dz[1:nf+1])
#coh02_fg_c2v_dz_m[nf:2*nf+1]  =         coh02_fg_c2v_dz
#coh02_fg_d2d_dx_m[ 0:nf]      = np.flip(coh02_fg_d2d_dx[1:nf+1])
#coh02_fg_d2d_dx_m[nf :2*nf+1] =         coh02_fg_d2d_dx
#coh02_fg_d2d_dy_m[ 0:nf]      = np.flip(coh02_fg_d2d_dy[1:nf+1])
#coh02_fg_d2d_dy_m[nf:2*nf+1]  =         coh02_fg_d2d_dy
#coh02_fg_d2d_dz_m[ 0:nf]      = np.flip(coh02_fg_d2d_dz[1:nf+1])
#coh02_fg_d2d_dz_m[nf:2*nf+1]  =         coh02_fg_d2d_dz
#coh02_c2v_dx_m[ 0:nf]         = np.flip(coh02_c2v_dx[1:nf+1])
#coh02_c2v_dx_m[nf :2*nf+1]    =         coh02_c2v_dx
#coh02_c2v_dy_m[ 0:nf]         = np.flip(coh02_c2v_dy[1:nf+1])
#coh02_c2v_dy_m[nf:2*nf+1]     =         coh02_c2v_dy
#coh02_c2v_dz_m[ 0:nf]         = np.flip(coh02_c2v_dz[1:nf+1])
#coh02_c2v_dz_m[nf:2*nf+1]     =         coh02_c2v_dz
#coh02_d2d_dx_m[ 0:nf]         = np.flip(coh02_d2d_dx[1:nf+1])
#coh02_d2d_dx_m[nf :2*nf+1]    =         coh02_d2d_dx
#coh02_d2d_dy_m[ 0:nf]         = np.flip(coh02_d2d_dy[1:nf+1])
#coh02_d2d_dy_m[nf:2*nf+1]     =         coh02_d2d_dy
#coh02_d2d_dz_m[ 0:nf]         = np.flip(coh02_d2d_dz[1:nf+1])
#coh02_d2d_dz_m[nf:2*nf+1]     =         coh02_d2d_dz
#coh02_td_dx_m[ 0:nf]          = np.flip(coh02_td_dx[1:nf+1])
#coh02_td_dx_m[nf :2*nf+1]     =         coh02_td_dx
#coh02_td_dy_m[ 0:nf]          = np.flip(coh02_td_dy[1:nf+1])
#coh02_td_dy_m[nf:2*nf+1]      =         coh02_td_dy
#coh02_td_dz_m[ 0:nf]          = np.flip(coh02_td_dz[1:nf+1])
#coh02_td_dz_m[nf:2*nf+1]      =         coh02_td_dz
#coh02_dx_m[ 0:nf]             = np.flip(coh02_dx[1:nf+1])
#coh02_dx_m[nf :2*nf+1]        =         coh02_dx
#coh02_dy_m[ 0:nf]             = np.flip(coh02_dy[1:nf+1])
#coh02_dy_m[nf:2*nf+1]         =         coh02_dy
#coh02_dz_m[ 0:nf]             = np.flip(coh02_dz[1:nf+1])
#coh02_dz_m[nf:2*nf+1]         =         coh02_dz
#
#coh12_fg_c2v_dx_m             = np.zeros(2*nf+1)
#coh12_fg_c2v_dy_m             = np.zeros(2*nf+1)
#coh12_fg_c2v_dz_m             = np.zeros(2*nf+1)
#coh12_fg_d2d_dx_m             = np.zeros(2*nf+1)
#coh12_fg_d2d_dy_m             = np.zeros(2*nf+1)
#coh12_fg_d2d_dz_m             = np.zeros(2*nf+1)
#coh12_c2v_dx_m             = np.zeros(2*nf+1)
#coh12_c2v_dy_m             = np.zeros(2*nf+1)
#coh12_c2v_dz_m             = np.zeros(2*nf+1)
#coh12_d2d_dx_m             = np.zeros(2*nf+1)
#coh12_d2d_dy_m             = np.zeros(2*nf+1)
#coh12_d2d_dz_m             = np.zeros(2*nf+1)
#coh12_td_dx_m             = np.zeros(2*nf+1)
#coh12_td_dy_m             = np.zeros(2*nf+1)
#coh12_td_dz_m             = np.zeros(2*nf+1)
#coh12_dx_m             = np.zeros(2*nf+1)
#coh12_dy_m             = np.zeros(2*nf+1)
#coh12_dz_m             = np.zeros(2*nf+1)
#coh12_fg_c2v_dx_m[ 0:nf]      = np.flip(coh12_fg_c2v_dx[1:nf+1])
#coh12_fg_c2v_dx_m[nf :2*nf+1] =         coh12_fg_c2v_dx
#coh12_fg_c2v_dy_m[ 0:nf]      = np.flip(coh12_fg_c2v_dy[1:nf+1])
#coh12_fg_c2v_dy_m[nf:2*nf+1]  =         coh12_fg_c2v_dy
#coh12_fg_c2v_dz_m[ 0:nf]      = np.flip(coh12_fg_c2v_dz[1:nf+1])
#coh12_fg_c2v_dz_m[nf:2*nf+1]  =         coh12_fg_c2v_dz
#coh12_fg_d2d_dx_m[ 0:nf]      = np.flip(coh12_fg_d2d_dx[1:nf+1])
#coh12_fg_d2d_dx_m[nf :2*nf+1] =         coh12_fg_d2d_dx
#coh12_fg_d2d_dy_m[ 0:nf]      = np.flip(coh12_fg_d2d_dy[1:nf+1])
#coh12_fg_d2d_dy_m[nf:2*nf+1]  =         coh12_fg_d2d_dy
#coh12_fg_d2d_dz_m[ 0:nf]      = np.flip(coh12_fg_d2d_dz[1:nf+1])
#coh12_fg_d2d_dz_m[nf:2*nf+1]  =         coh12_fg_d2d_dz
#coh12_c2v_dx_m[ 0:nf]         = np.flip(coh12_c2v_dx[1:nf+1])
#coh12_c2v_dx_m[nf :2*nf+1]    =         coh12_c2v_dx
#coh12_c2v_dy_m[ 0:nf]         = np.flip(coh12_c2v_dy[1:nf+1])
#coh12_c2v_dy_m[nf:2*nf+1]     =         coh12_c2v_dy
#coh12_c2v_dz_m[ 0:nf]         = np.flip(coh12_c2v_dz[1:nf+1])
#coh12_c2v_dz_m[nf:2*nf+1]     =         coh12_c2v_dz
#coh12_d2d_dx_m[ 0:nf]         = np.flip(coh12_d2d_dx[1:nf+1])
#coh12_d2d_dx_m[nf :2*nf+1]    =         coh12_d2d_dx
#coh12_d2d_dy_m[ 0:nf]         = np.flip(coh12_d2d_dy[1:nf+1])
#coh12_d2d_dy_m[nf:2*nf+1]     =         coh12_d2d_dy
#coh12_d2d_dz_m[ 0:nf]         = np.flip(coh12_d2d_dz[1:nf+1])
#coh12_d2d_dz_m[nf:2*nf+1]     =         coh12_d2d_dz
#coh12_td_dx_m[ 0:nf]          = np.flip(coh12_td_dx[1:nf+1])
#coh12_td_dx_m[nf :2*nf+1]     =         coh12_td_dx
#coh12_td_dy_m[ 0:nf]          = np.flip(coh12_td_dy[1:nf+1])
#coh12_td_dy_m[nf:2*nf+1]      =         coh12_td_dy
#coh12_td_dz_m[ 0:nf]          = np.flip(coh12_td_dz[1:nf+1])
#coh12_td_dz_m[nf:2*nf+1]      =         coh12_td_dz
#coh12_dx_m[ 0:nf]             = np.flip(coh12_dx[1:nf+1])
#coh12_dx_m[nf :2*nf+1]        =         coh12_dx
#coh12_dy_m[ 0:nf]             = np.flip(coh12_dy[1:nf+1])
#coh12_dy_m[nf:2*nf+1]         =         coh12_dy
#coh12_dz_m[ 0:nf]             = np.flip(coh12_dz[1:nf+1])
#coh12_dz_m[nf:2*nf+1]         =         coh12_dz
#
#
#
#
#sigma = time[nf]/4
#x0 = 0
#const = 1
#
## Multiply by a gaussian to ensure smooth decay to zero at the borders
#
#coh01_fg_c2v_dx_m = coh01_fg_c2v_dx_m * gaussian(const,time_m,x0,sigma)
#coh01_fg_c2v_dy_m = coh01_fg_c2v_dy_m * gaussian(const,time_m,x0,sigma)
#coh01_fg_c2v_dz_m = coh01_fg_c2v_dz_m * gaussian(const,time_m,x0,sigma)
#coh01_fg_d2d_dx_m = coh01_fg_d2d_dx_m * gaussian(const,time_m,x0,sigma)
#coh01_fg_d2d_dy_m = coh01_fg_d2d_dy_m * gaussian(const,time_m,x0,sigma)
#coh01_fg_d2d_dz_m = coh01_fg_d2d_dz_m * gaussian(const,time_m,x0,sigma)
#coh01_c2v_dx_m    = coh01_c2v_dx_m    * gaussian(const,time_m,x0,sigma)
#coh01_c2v_dy_m    = coh01_c2v_dy_m    * gaussian(const,time_m,x0,sigma)
#coh01_c2v_dz_m    = coh01_c2v_dz_m    * gaussian(const,time_m,x0,sigma)
#coh01_d2d_dx_m    = coh01_d2d_dx_m    * gaussian(const,time_m,x0,sigma)
#coh01_d2d_dy_m    = coh01_d2d_dy_m    * gaussian(const,time_m,x0,sigma)
#coh01_d2d_dz_m    = coh01_d2d_dz_m    * gaussian(const,time_m,x0,sigma)
#coh01_td_dx_m     = coh01_td_dx_m     * gaussian(const,time_m,x0,sigma)
#coh01_td_dy_m     = coh01_td_dy_m     * gaussian(const,time_m,x0,sigma)
#coh01_td_dz_m     = coh01_td_dz_m     * gaussian(const,time_m,x0,sigma)
#coh01_dx_m        = coh01_dx_m        * gaussian(const,time_m,x0,sigma)
#coh01_dy_m        = coh01_dy_m        * gaussian(const,time_m,x0,sigma)
#coh01_dz_m        = coh01_dz_m        * gaussian(const,time_m,x0,sigma)
#coh02_fg_c2v_dx_m = coh02_fg_c2v_dx_m * gaussian(const,time_m,x0,sigma)
#coh02_fg_c2v_dy_m = coh02_fg_c2v_dy_m * gaussian(const,time_m,x0,sigma)
#coh02_fg_c2v_dz_m = coh02_fg_c2v_dz_m * gaussian(const,time_m,x0,sigma)
#coh02_fg_d2d_dx_m = coh02_fg_d2d_dx_m * gaussian(const,time_m,x0,sigma)
#coh02_fg_d2d_dy_m = coh02_fg_d2d_dy_m * gaussian(const,time_m,x0,sigma)
#coh02_fg_d2d_dz_m = coh02_fg_d2d_dz_m * gaussian(const,time_m,x0,sigma)
#coh02_c2v_dx_m    = coh02_c2v_dx_m    * gaussian(const,time_m,x0,sigma)
#coh02_c2v_dy_m    = coh02_c2v_dy_m    * gaussian(const,time_m,x0,sigma)
#coh02_c2v_dz_m    = coh02_c2v_dz_m    * gaussian(const,time_m,x0,sigma)
#coh02_d2d_dx_m    = coh02_d2d_dx_m    * gaussian(const,time_m,x0,sigma)
#coh02_d2d_dy_m    = coh02_d2d_dy_m    * gaussian(const,time_m,x0,sigma)
#coh02_d2d_dz_m    = coh02_d2d_dz_m    * gaussian(const,time_m,x0,sigma)
#coh02_td_dx_m     = coh02_td_dx_m     * gaussian(const,time_m,x0,sigma)
#coh02_td_dy_m     = coh02_td_dy_m     * gaussian(const,time_m,x0,sigma)
#coh02_td_dz_m     = coh02_td_dz_m     * gaussian(const,time_m,x0,sigma)
#coh02_dx_m        = coh02_dx_m        * gaussian(const,time_m,x0,sigma)
#coh02_dy_m        = coh02_dy_m        * gaussian(const,time_m,x0,sigma)
#coh02_dz_m        = coh02_dz_m        * gaussian(const,time_m,x0,sigma)
#coh12_fg_c2v_dx_m = coh12_fg_c2v_dx_m * gaussian(const,time_m,x0,sigma)
#coh12_fg_c2v_dy_m = coh12_fg_c2v_dy_m * gaussian(const,time_m,x0,sigma)
#coh12_fg_c2v_dz_m = coh12_fg_c2v_dz_m * gaussian(const,time_m,x0,sigma)
#coh12_fg_d2d_dx_m = coh12_fg_d2d_dx_m * gaussian(const,time_m,x0,sigma)
#coh12_fg_d2d_dy_m = coh12_fg_d2d_dy_m * gaussian(const,time_m,x0,sigma)
#coh12_fg_d2d_dz_m = coh12_fg_d2d_dz_m * gaussian(const,time_m,x0,sigma)
#coh12_c2v_dx_m    = coh12_c2v_dx_m    * gaussian(const,time_m,x0,sigma)
#coh12_c2v_dy_m    = coh12_c2v_dy_m    * gaussian(const,time_m,x0,sigma)
#coh12_c2v_dz_m    = coh12_c2v_dz_m    * gaussian(const,time_m,x0,sigma)
#coh12_d2d_dx_m    = coh12_d2d_dx_m    * gaussian(const,time_m,x0,sigma)
#coh12_d2d_dy_m    = coh12_d2d_dy_m    * gaussian(const,time_m,x0,sigma)
#coh12_d2d_dz_m    = coh12_d2d_dz_m    * gaussian(const,time_m,x0,sigma)
#coh12_td_dx_m     = coh12_td_dx_m     * gaussian(const,time_m,x0,sigma)
#coh12_td_dy_m     = coh12_td_dy_m     * gaussian(const,time_m,x0,sigma)
#coh12_td_dz_m     = coh12_td_dz_m     * gaussian(const,time_m,x0,sigma)
#coh12_dx_m        = coh12_dx_m        * gaussian(const,time_m,x0,sigma)
#coh12_dy_m        = coh12_dy_m        * gaussian(const,time_m,x0,sigma)
#coh12_dz_m        = coh12_dz_m        * gaussian(const,time_m,x0,sigma)
#
#factor = 20
#
#fft_coh_dx_01        = np.absolute(np.fft.fft(       coh01_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dx_02        = np.absolute(np.fft.fft(       coh02_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dx_12        = np.absolute(np.fft.fft(       coh12_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dy_01        = np.absolute(np.fft.fft(       coh01_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dy_02        = np.absolute(np.fft.fft(       coh02_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dy_12        = np.absolute(np.fft.fft(       coh12_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dz_01        = np.absolute(np.fft.fft(       coh01_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dz_02        = np.absolute(np.fft.fft(       coh02_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_dz_12        = np.absolute(np.fft.fft(       coh12_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dx_01 = np.absolute(np.fft.fft(coh01_fg_c2v_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dx_02 = np.absolute(np.fft.fft(coh02_fg_c2v_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dx_12 = np.absolute(np.fft.fft(coh12_fg_c2v_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dy_01 = np.absolute(np.fft.fft(coh01_fg_c2v_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dy_02 = np.absolute(np.fft.fft(coh02_fg_c2v_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dy_12 = np.absolute(np.fft.fft(coh12_fg_c2v_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dz_01 = np.absolute(np.fft.fft(coh01_fg_c2v_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dz_02 = np.absolute(np.fft.fft(coh02_fg_c2v_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_c2v_dz_12 = np.absolute(np.fft.fft(coh12_fg_c2v_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dx_01 = np.absolute(np.fft.fft(coh01_fg_d2d_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dx_02 = np.absolute(np.fft.fft(coh02_fg_d2d_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dx_12 = np.absolute(np.fft.fft(coh12_fg_d2d_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dy_01 = np.absolute(np.fft.fft(coh01_fg_d2d_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dy_02 = np.absolute(np.fft.fft(coh02_fg_d2d_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dy_12 = np.absolute(np.fft.fft(coh12_fg_d2d_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dz_01 = np.absolute(np.fft.fft(coh01_fg_d2d_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dz_02 = np.absolute(np.fft.fft(coh02_fg_d2d_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_fg_d2d_dz_12 = np.absolute(np.fft.fft(coh12_fg_d2d_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dx_01    = np.absolute(np.fft.fft(   coh01_c2v_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dx_02    = np.absolute(np.fft.fft(   coh02_c2v_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dx_12    = np.absolute(np.fft.fft(   coh12_c2v_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dy_01    = np.absolute(np.fft.fft(   coh01_c2v_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dy_02    = np.absolute(np.fft.fft(   coh02_c2v_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dy_12    = np.absolute(np.fft.fft(   coh12_c2v_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dz_01    = np.absolute(np.fft.fft(   coh01_c2v_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dz_02    = np.absolute(np.fft.fft(   coh02_c2v_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_c2v_dz_12    = np.absolute(np.fft.fft(   coh12_c2v_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dx_01    = np.absolute(np.fft.fft(   coh01_d2d_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dx_02    = np.absolute(np.fft.fft(   coh02_d2d_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dx_12    = np.absolute(np.fft.fft(   coh12_d2d_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dy_01    = np.absolute(np.fft.fft(   coh01_d2d_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dy_02    = np.absolute(np.fft.fft(   coh02_d2d_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dy_12    = np.absolute(np.fft.fft(   coh12_d2d_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dz_01    = np.absolute(np.fft.fft(   coh01_d2d_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dz_02    = np.absolute(np.fft.fft(   coh02_d2d_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_d2d_dz_12    = np.absolute(np.fft.fft(   coh12_d2d_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dx_01     = np.absolute(np.fft.fft(    coh01_td_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dx_02     = np.absolute(np.fft.fft(    coh02_td_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dx_12     = np.absolute(np.fft.fft(    coh12_td_dx_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dy_01     = np.absolute(np.fft.fft(    coh01_td_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dy_02     = np.absolute(np.fft.fft(    coh02_td_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dy_12     = np.absolute(np.fft.fft(    coh12_td_dy_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dz_01     = np.absolute(np.fft.fft(    coh01_td_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dz_02     = np.absolute(np.fft.fft(    coh02_td_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#fft_coh_td_dz_12     = np.absolute(np.fft.fft(    coh12_td_dz_m, 2*nf*factor+1)*tstep)[0:nf*factor+1]
#
#dnu = 1 / ( 2*tf )
#fre = ( (np.arange(-factor*nf,factor*nf+1)) * dnu /( 2*nf*factor/(2*nf+1)) )[nf*factor:2*nf*factor+1]
#freq_ev = fre*27.211396*2*np.pi
#np.savetxt('fft_freq_ev.dat',freq_ev)
#
#
#np.savetxt('fft-coh-fg-dx-01.dat',fft_coh_dx_01)
#np.savetxt('fft-coh-fg-dx-02.dat',fft_coh_dx_02)
#np.savetxt('fft-coh-fg-dx-12.dat',fft_coh_dx_12)
#np.savetxt('fft-coh-fg-dy-01.dat',fft_coh_dy_01)
#np.savetxt('fft-coh-fg-dy-02.dat',fft_coh_dy_02)
#np.savetxt('fft-coh-fg-dy-12.dat',fft_coh_dy_12)
#np.savetxt('fft-coh-fg-dz-01.dat',fft_coh_dz_01)
#np.savetxt('fft-coh-fg-dz-02.dat',fft_coh_dz_02)
#np.savetxt('fft-coh-fg-dz-12.dat',fft_coh_dz_12)
#np.savetxt('fft-coh-c2v-region-dx-01.dat',fft_coh_fg_c2v_dx_01)
#np.savetxt('fft-coh-c2v-region-dx-02.dat',fft_coh_fg_c2v_dx_02)
#np.savetxt('fft-coh-c2v-region-dx-12.dat',fft_coh_fg_c2v_dx_12)
#np.savetxt('fft-coh-d2d-region-dx-01.dat',fft_coh_fg_d2d_dx_01)
#np.savetxt('fft-coh-d2d-region-dx-02.dat',fft_coh_fg_d2d_dx_02)
#np.savetxt('fft-coh-d2d-region-dx-12.dat',fft_coh_fg_d2d_dx_12)
#np.savetxt('fft-coh-c2v-region-dy-01.dat',fft_coh_fg_c2v_dy_01)
#np.savetxt('fft-coh-c2v-region-dy-02.dat',fft_coh_fg_c2v_dy_02)
#np.savetxt('fft-coh-c2v-region-dy-12.dat',fft_coh_fg_c2v_dy_12)
#np.savetxt('fft-coh-d2d-region-dy-01.dat',fft_coh_fg_d2d_dy_01)
#np.savetxt('fft-coh-d2d-region-dy-02.dat',fft_coh_fg_d2d_dy_02)
#np.savetxt('fft-coh-d2d-region-dy-12.dat',fft_coh_fg_d2d_dy_12)
#np.savetxt('fft-coh-c2v-region-dz-01.dat',fft_coh_fg_c2v_dz_01)
#np.savetxt('fft-coh-c2v-region-dz-02.dat',fft_coh_fg_c2v_dz_02)
#np.savetxt('fft-coh-c2v-region-dz-12.dat',fft_coh_fg_c2v_dz_12)
#np.savetxt('fft-coh-d2d-region-dz-01.dat',fft_coh_fg_d2d_dz_01)
#np.savetxt('fft-coh-d2d-region-dz-02.dat',fft_coh_fg_d2d_dz_02)
#np.savetxt('fft-coh-d2d-region-dz-12.dat',fft_coh_fg_d2d_dz_12)
#np.savetxt('fft-coh-c2v-point-01-dx.dat',fft_coh_c2v_dx_01)
#np.savetxt('fft-coh-c2v-point-02-dx.dat',fft_coh_c2v_dx_02)
#np.savetxt('fft-coh-c2v-point-12-dx.dat',fft_coh_c2v_dx_12)
#np.savetxt('fft-coh-d2d-point-01-dx.dat',fft_coh_d2d_dx_01)
#np.savetxt('fft-coh-d2d-point-02-dx.dat',fft_coh_d2d_dx_02)
#np.savetxt('fft-coh-d2d-point-12-dx.dat',fft_coh_d2d_dx_12)
#np.savetxt('fft-coh-td-point-01-dx.dat',fft_coh_td_dx_01)
#np.savetxt('fft-coh-td-point-02-dx.dat',fft_coh_td_dx_02)
#np.savetxt('fft-coh-td-point-12-dx.dat',fft_coh_td_dx_12)
#np.savetxt('fft-coh-c2v-point-01-dy.dat',fft_coh_c2v_dy_01)
#np.savetxt('fft-coh-c2v-point-02-dy.dat',fft_coh_c2v_dy_02)
#np.savetxt('fft-coh-c2v-point-12-dy.dat',fft_coh_c2v_dy_12)
#np.savetxt('fft-coh-d2d-point-01-dy.dat',fft_coh_d2d_dy_01)
#np.savetxt('fft-coh-d2d-point-02-dy.dat',fft_coh_d2d_dy_02)
#np.savetxt('fft-coh-d2d-point-12-dy.dat',fft_coh_d2d_dy_12)
#np.savetxt('fft-coh-td-point-01-dy.dat',fft_coh_td_dy_01)
#np.savetxt('fft-coh-td-point-02-dy.dat',fft_coh_td_dy_02)
#np.savetxt('fft-coh-td-point-12-dy.dat',fft_coh_td_dy_12)
#np.savetxt('fft-coh-c2v-point-01-dz.dat',fft_coh_c2v_dz_01)
#np.savetxt('fft-coh-c2v-point-02-dz.dat',fft_coh_c2v_dz_02)
#np.savetxt('fft-coh-c2v-point-12-dz.dat',fft_coh_c2v_dz_12)
#np.savetxt('fft-coh-d2d-point-01-dz.dat',fft_coh_d2d_dz_01)
#np.savetxt('fft-coh-d2d-point-02-dz.dat',fft_coh_d2d_dz_02)
#np.savetxt('fft-coh-d2d-point-12-dz.dat',fft_coh_d2d_dz_12)
#np.savetxt('fft-coh-td-point-01-dz.dat',fft_coh_td_dz_01)
#np.savetxt('fft-coh-td-point-02-dz.dat',fft_coh_td_dz_02)
#np.savetxt('fft-coh-td-point-12-dz.dat',fft_coh_td_dz_12)
#
#
#
#
#outputname = 'fft_coherence_on_C2v_point_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_c2v_dx_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_c2v_dx_02, 'b', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_c2v_dx_12, 'r', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_D2d_point_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_d2d_dx_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_d2d_dx_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_d2d_dx_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_Td_point_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_td_dx_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_td_dx_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_td_dx_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'fft_coherence_on_C2v_point_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_c2v_dy_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_c2v_dy_02, 'b', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_c2v_dy_12, 'r', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_D2d_point_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_d2d_dy_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_d2d_dy_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_d2d_dx_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_Td_point_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_td_dy_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_td_dy_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_td_dy_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'fft_coherence_on_C2v_point_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_c2v_dz_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_c2v_dz_02, 'b', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_c2v_dz_12, 'r', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_D2d_point_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_d2d_dz_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_d2d_dz_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_d2d_dz_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_Td_point_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_td_dz_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_td_dz_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_td_dz_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'fft_coherence_on_C2v_region_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_fg_c2v_dx_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_fg_c2v_dx_02, 'b', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_fg_c2v_dx_12, 'r', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_D2d_region_dx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_fg_d2d_dx_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_fg_d2d_dx_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_fg_d2d_dx_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'fft_coherence_on_C2v_region_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_fg_c2v_dy_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_fg_c2v_dy_02, 'b', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_fg_c2v_dy_12, 'r', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_D2d_region_dy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_fg_d2d_dy_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_fg_d2d_dy_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_fg_d2d_dx_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#outputname = 'fft_coherence_on_C2v_region_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_fg_c2v_dz_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_fg_c2v_dz_02, 'b', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_fg_c2v_dz_12, 'r', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_on_D2d_region_dz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_fg_d2d_dz_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_fg_d2d_dz_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_fg_d2d_dz_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#
#
#
#outputname = 'fft_coherence_full_grid_dipx.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_dx_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_dx_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_x$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_dx_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_full_grid_dipy.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_dy_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_dy_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_y$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_dy_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'fft_coherence_full_grid_dipz.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
#ax[2].plot(freq_ev, fft_coh_dz_01, 'k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Frequency', fontsize=15)
#ax[2].legend(['St0-St1'])
#ax[2].set_xlim([0,22.5])
#ax[2].locator_params(axis='x', nbins=14)
#
#ax[1].plot(freq_ev, fft_coh_dz_02, 'b', lw=1.5)
#ax[1].set_ylabel('FFT Coherence $\cdot \mu_z$', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St0-St2'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,22.5])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(freq_ev, fft_coh_dz_12, 'r', lw=1.5)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St1-St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,22.5])
#ax[0].axes.get_xaxis().set_visible(False)
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')


t = 1
dsname = "time-pop-{:06}".format(t)
fname = '{}.h5'.format(dsname)
file = h5.File(fname, 'r')
aux2 = np.asarray(file[dsname])
norma = np.linalg.norm(aux2, ord=1)
#norma = 1
for t in tqdm(inputs):
    dsname = "time-pop-{:06}".format(t)
    fname = '{}.h5'.format(dsname)
    file = h5.File(fname, 'r')
    aux2 = np.asarray(file[dsname])
    aux2 = aux2 / np.linalg.norm(aux2, ord=1)
    z = aux2.reshape(3, int(nq2), int(nq1))
    
    s = 1
    cg0 = aux2[(s-1)*ng:s*ng]
    cg = cg0.reshape(nq2,nq1)
    cg = np.transpose(cg)
    
    pop01_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
    pop01_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
    pop01_td [t-1] = cg[73-1,115-1] #Making a subgrid
    pop01_cs [t-1] = cg[60-1,133-1] #Making a subgrid
    for i in np.arange(nq1*nq2):
        pop01[t-1] = pop01[t-1] + cg0[i]
    #########################################################################
    #c2v region
    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv = np.zeros(n1c2v*n2c2v)
    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dx = np.zeros(n1c2v*n2c2v)
    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dy = np.zeros(n1c2v*n2c2v)
    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dz = np.zeros(n1c2v*n2c2v)
    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig

    for i in np.arange(n1c2v*n2c2v):
        pop01_fg_c2v[t-1]    = pop01_fg_c2v[t-1] + sv[i]
    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
    #########################################################################
    #D2d region
    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv = np.zeros(n1d2d*n2d2d)
    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dx = np.zeros(n1d2d*n2d2d)
    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dy = np.zeros(n1d2d*n2d2d)
    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dz = np.zeros(n1d2d*n2d2d)
    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig

    for i in np.arange(n1d2d*n2d2d):
        pop01_fg_d2d[t-1]    = pop01_fg_d2d[t-1] + sv[i]
    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
    #########################################################################
    #c2v region 9 points
    sub = cg[q1ic2v9-1:q1fc2v9,q2ic2v9-1:q2fc2v9] #Making a subgrid
    sv = np.zeros(n1c2v9*n2c2v9)
    sv = sub.reshape(n1c2v9*n2c2v9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1c2v9*n2c2v9):
        pop01_c2v9[t-1]    = pop01_c2v9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    #d2d region 9 points
    sub = cg[q1id2d9-1:q1fd2d9,q2id2d9-1:q2fd2d9] #Making a subgrid
    sv = np.zeros(n1d2d9*n2d2d9)
    sv = sub.reshape(n1d2d9*n2d2d9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1d2d9*n2d2d9):
        pop01_d2d9[t-1]    = pop01_d2d9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    #cs region 9 points
    sub = cg[q1ics9-1:q1fcs9,q2ics9-1:q2fcs9] #Making a subgrid
    sv = np.zeros(n1cs9*n2cs9)
    sv = sub.reshape(n1cs9*n2cs9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1cs9*n2cs9):
        pop01_cs9[t-1]    = pop01_cs9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    
    s = 2
    cg0 = aux2[(s-1)*ng:s*ng]
    cg = cg0.reshape(nq2,nq1)
    cg = np.transpose(cg)
    
    pop02_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
    pop02_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
    pop02_td [t-1] = cg[73-1,115-1] #Making a subgrid
    pop02_cs [t-1] = cg[60-1,133-1] #Making a subgrid
    for i in np.arange(nq1*nq2):
        pop02[t-1]    = pop02[t-1] + cg0[i]
    #########################################################################
    #c2v region
    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv = np.zeros(n1c2v*n2c2v)
    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dx = np.zeros(n1c2v*n2c2v)
    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dy = np.zeros(n1c2v*n2c2v)
    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dz = np.zeros(n1c2v*n2c2v)
    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig

    for i in np.arange(n1c2v*n2c2v):
        pop02_fg_c2v[t-1]    = pop02_fg_c2v[t-1] + sv[i]
    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
    #########################################################################
    #D2d region
    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv = np.zeros(n1d2d*n2d2d)
    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dx = np.zeros(n1d2d*n2d2d)
    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dy = np.zeros(n1d2d*n2d2d)
    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dz = np.zeros(n1d2d*n2d2d)
    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig

    for i in np.arange(n1d2d*n2d2d):
        pop02_fg_d2d[t-1]    = pop02_fg_d2d[t-1] + sv[i]
    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
    #########################################################################
    #c2v region 9 points
    sub = cg[q1ic2v9-1:q1fc2v9,q2ic2v9-1:q2fc2v9] #Making a subgrid
    sv = np.zeros(n1c2v9*n2c2v9)
    sv = sub.reshape(n1c2v9*n2c2v9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1c2v9*n2c2v9):
        pop02_c2v9[t-1]    = pop02_c2v9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    #d2d region 9 points
    sub = cg[q1id2d9-1:q1fd2d9,q2id2d9-1:q2fd2d9] #Making a subgrid
    sv = np.zeros(n1d2d9*n2d2d9)
    sv = sub.reshape(n1d2d9*n2d2d9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1d2d9*n2d2d9):
        pop02_d2d9[t-1]    = pop02_d2d9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    #cs region 9 points
    sub = cg[q1ics9-1:q1fcs9,q2ics9-1:q2fcs9] #Making a subgrid
    sv = np.zeros(n1cs9*n2cs9)
    sv = sub.reshape(n1cs9*n2cs9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1cs9*n2cs9):
        pop02_cs9[t-1]    = pop02_cs9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    
    s = 3
    cg0 = aux2[(s-1)*ng:s*ng]
    cg = cg0.reshape(nq2,nq1)
    cg = np.transpose(cg)
    
    pop12_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
    pop12_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
    pop12_td [t-1] = cg[73-1,115-1] #Making a subgrid
    pop12_cs [t-1] = cg[60-1,133-1] #Making a subgrid
    for i in np.arange(nq1*nq2):
        pop12[t-1]    = pop12[t-1] + cg0[i]
    #########################################################################
    #c2v region
    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv = np.zeros(n1c2v*n2c2v)
    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dx = np.zeros(n1c2v*n2c2v)
    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dy = np.zeros(n1c2v*n2c2v)
    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
    sv_dz = np.zeros(n1c2v*n2c2v)
    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig

    for i in np.arange(n1c2v*n2c2v):
        pop12_fg_c2v[t-1]    = pop12_fg_c2v[t-1] + sv[i]
    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
    #########################################################################
    #D2d region
    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv = np.zeros(n1d2d*n2d2d)
    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dx = np.zeros(n1d2d*n2d2d)
    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dy = np.zeros(n1d2d*n2d2d)
    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
    sv_dz = np.zeros(n1d2d*n2d2d)
    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig

    for i in np.arange(n1d2d*n2d2d):
        pop12_fg_d2d[t-1]    = pop12_fg_d2d[t-1] + sv[i]
    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
    ######################################################################### 
    #c2v region 9 points
    sub = cg[q1ic2v9-1:q1fc2v9,q2ic2v9-1:q2fc2v9] #Making a subgrid
    sv = np.zeros(n1c2v9*n2c2v9)
    sv = sub.reshape(n1c2v9*n2c2v9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1c2v9*n2c2v9):
        pop12_c2v9[t-1]    = pop12_c2v9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    #d2d region 9 points
    sub = cg[q1id2d9-1:q1fd2d9,q2id2d9-1:q2fd2d9] #Making a subgrid
    sv = np.zeros(n1d2d9*n2d2d9)
    sv = sub.reshape(n1d2d9*n2d2d9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1d2d9*n2d2d9):
        pop12_d2d9[t-1]    = pop12_d2d9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################
    #cs region 9 points
    sub = cg[q1ics9-1:q1fcs9,q2ics9-1:q2fcs9] #Making a subgrid
    sv = np.zeros(n1cs9*n2cs9)
    sv = sub.reshape(n1cs9*n2cs9, order = 'F') #Vectorizing the subgrig
    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
    
    for i in np.arange(n1cs9*n2cs9):
        pop12_cs9[t-1]    = pop12_cs9[t-1] + sv[i]
    del(sv,sub,sv1)
    #########################################################################

np.savetxt('population-fg-0.dat',pop01)
np.savetxt('population-fg-1.dat',pop02)
np.savetxt('population-fg-2.dat',pop12)
#np.savetxt('population-c2v-region-0.dat',pop01_fg_c2v)
#np.savetxt('population-c2v-region-1.dat',pop02_fg_c2v)
#np.savetxt('population-c2v-region-2.dat',pop12_fg_c2v)
#np.savetxt('population-d2d-region-0.dat',pop01_fg_d2d)
#np.savetxt('population-d2d-region-1.dat',pop02_fg_d2d)
#np.savetxt('population-d2d-region-2.dat',pop12_fg_d2d)
np.savetxt('population-c2v-point-0.dat',pop01_c2v)
np.savetxt('population-c2v-point-1.dat',pop02_c2v)
np.savetxt('population-c2v-point-2.dat',pop12_c2v)
np.savetxt('population-d2d-point-0.dat',pop01_d2d)
np.savetxt('population-d2d-point-1.dat',pop02_d2d)
np.savetxt('population-d2d-point-2.dat',pop12_d2d)
np.savetxt('population-cs-point-0.dat',pop01_cs)
np.savetxt('population-cs-point-1.dat',pop02_cs)
np.savetxt('population-cs-point-2.dat',pop12_cs)
np.savetxt('population-c2v-9points-0.dat',pop01_c2v9)
np.savetxt('population-c2v-9points-1.dat',pop02_c2v9)
np.savetxt('population-c2v-9points-2.dat',pop12_c2v9)
np.savetxt('population-d2d-9points-0.dat',pop01_d2d9)
np.savetxt('population-d2d-9points-1.dat',pop02_d2d9)
np.savetxt('population-d2d-9points-2.dat',pop12_d2d9)
np.savetxt('population-cs-9points-0.dat',pop01_cs9)
np.savetxt('population-cs-9points-1.dat',pop02_cs9)
np.savetxt('population-cs-9points-2.dat',pop12_cs9)


outputname = 'population_on_C2v_point.png'
fig, ax = plt.subplots(3, figsize=(6,4))
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
ax[2].plot(time*t_au, pop01_c2v, 'k', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[2].set_ylabel('', fontsize=15)
ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
ax[2].legend(['St0'])
ax[2].set_xlim([0,time[-1]*t_au])
ax[2].locator_params(axis='x', nbins=14)
# ax[2].axes.get_xaxis().set_visible(False)

ax[1].plot(time*t_au, pop02_c2v, 'b', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[1].set_xlim(0,time[-1]*t_au)
ax[1].set_ylabel('population', fontsize=15)
ax[1].set_xlabel('', fontsize=15)
ax[1].legend(['St1'])
ax[1].locator_params(axis='x', nbins=12)
ax[1].set_xlim([0,time[-1]*t_au])
ax[1].axes.get_xaxis().set_visible(False)

ax[0].plot(time*t_au, pop12_c2v, 'r', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[0].set_xlim(0,time[-1]*t_au)
ax[0].set_ylabel('', fontsize=15)
ax[0].set_xlabel('', fontsize=15)
ax[0].legend(['St2'])
ax[0].locator_params(axis='x', nbins=12)
ax[0].set_xlim([0,time[-1]*t_au])
ax[0].axes.get_xaxis().set_visible(False)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')


outputname = 'population_on_D2d_point.png'
fig, ax = plt.subplots(3, figsize=(6,4))
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
ax[2].plot(time*t_au, pop01_d2d, 'k', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[2].set_ylabel('', fontsize=15)
ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
ax[2].legend(['St0'])
ax[2].set_xlim([0,time[-1]*t_au])
ax[2].locator_params(axis='x', nbins=14)
# ax[2].axes.get_xaxis().set_visible(False)

ax[1].plot(time*t_au, pop02_d2d, 'b', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[1].set_xlim(0,time[-1]*t_au)
ax[1].set_ylabel('population', fontsize=15)
ax[1].set_xlabel('', fontsize=15)
ax[1].legend(['St1'])
ax[1].locator_params(axis='x', nbins=12)
ax[1].set_xlim([0,time[-1]*t_au])
ax[1].axes.get_xaxis().set_visible(False)

ax[0].plot(time*t_au, pop12_d2d, 'r', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[0].set_xlim(0,time[-1]*t_au)
ax[0].set_ylabel('', fontsize=15)
ax[0].set_xlabel('', fontsize=15)
ax[0].legend(['St2'])
ax[0].locator_params(axis='x', nbins=12)
ax[0].set_xlim([0,time[-1]*t_au])
ax[0].axes.get_xaxis().set_visible(False)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')


outputname = 'population_on_Td_point.png'
fig, ax = plt.subplots(3, figsize=(6,4))
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
ax[2].plot(time*t_au, pop01_td, 'k', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[2].set_ylabel('', fontsize=15)
ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
ax[2].legend(['St0'])
ax[2].set_xlim([0,time[-1]*t_au])
ax[2].locator_params(axis='x', nbins=14)
# ax[2].axes.get_xaxis().set_visible(False)

ax[1].plot(time*t_au, pop02_td, 'b', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[1].set_xlim(0,time[-1]*t_au)
ax[1].set_ylabel('population', fontsize=15)
ax[1].set_xlabel('', fontsize=15)
ax[1].legend(['St1'])
ax[1].locator_params(axis='x', nbins=12)
ax[1].set_xlim([0,time[-1]*t_au])
ax[1].axes.get_xaxis().set_visible(False)

ax[0].plot(time*t_au, pop12_td, 'r', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[0].set_xlim(0,time[-1]*t_au)
ax[0].set_ylabel('', fontsize=15)
ax[0].set_xlabel('', fontsize=15)
ax[0].legend(['St2'])
ax[0].locator_params(axis='x', nbins=12)
ax[0].set_xlim([0,time[-1]*t_au])
ax[0].axes.get_xaxis().set_visible(False)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')


outputname = 'population_full_grid.png'
fig, ax = plt.subplots(3, figsize=(6,4))
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
ax[2].plot(time*t_au, pop01, 'k', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[2].set_ylabel('', fontsize=15)
ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
ax[2].legend(['St0'])
ax[2].set_xlim([0,time[-1]*t_au])
ax[2].locator_params(axis='x', nbins=14)
# ax[2].axes.get_xaxis().set_visible(False)

ax[1].plot(time*t_au, pop02, 'b', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[1].set_xlim(0,time[-1]*t_au)
ax[1].set_ylabel('population', fontsize=15)
ax[1].set_xlabel('', fontsize=15)
ax[1].legend(['St1'])
ax[1].locator_params(axis='x', nbins=12)
ax[1].set_xlim([0,time[-1]*t_au])
ax[1].axes.get_xaxis().set_visible(False)

ax[0].plot(time*t_au, pop12, 'r', lw=1.5)
#ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
ax[0].set_xlim(0,time[-1]*t_au)
ax[0].set_ylabel('', fontsize=15)
ax[0].set_xlabel('', fontsize=15)
ax[0].legend(['St2'])
ax[0].locator_params(axis='x', nbins=12)
ax[0].set_xlim([0,time[-1]*t_au])
ax[0].axes.get_xaxis().set_visible(False)

fig.tight_layout()
fig.savefig(outputname)
plt.close('all')


#outputname = 'population_on_C2v_region.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop01_fg_c2v, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop02_fg_c2v, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop12_fg_c2v, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'population_on_D2d_region.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop01_fg_d2d, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop02_fg_d2d, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop12_fg_d2d, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
plt.close('all')




#for t in tqdm(inputs):
#    dsname = "diff-pop-{:06}".format(t)
#    fname = '{}.h5'.format(dsname)
#    file = h5.File(fname, 'r')
#    aux2 = np.asarray(file[dsname])
#    aux2 = aux2 / norma
#    z = aux2.reshape(3, int(nq2), int(nq1))
#    
#    s = 1
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#    
#    pop_nac01_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
#    pop_nac01_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
#    pop_nac01_td [t-1] = cg[73-1,115-1] #Making a subgrid
#    for i in np.arange(nq1*nq2):
#        pop_nac01[t-1] = pop_nac01[t-1] + cg0[i]
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dx = np.zeros(n1c2v*n2c2v)
#    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dy = np.zeros(n1c2v*n2c2v)
#    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dz = np.zeros(n1c2v*n2c2v)
#    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1c2v*n2c2v):
#        pop_nac01_fg_c2v[t-1]    = pop_nac01_fg_c2v[t-1] + sv[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dx = np.zeros(n1d2d*n2d2d)
#    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dy = np.zeros(n1d2d*n2d2d)
#    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dz = np.zeros(n1d2d*n2d2d)
#    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1d2d*n2d2d):
#        pop_nac01_fg_d2d[t-1]    = pop_nac01_fg_d2d[t-1] + sv[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    
#    s = 2
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#    
#    pop_nac02_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
#    pop_nac02_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
#    pop_nac02_td [t-1] = cg[73-1,115-1] #Making a subgrid
#    for i in np.arange(nq1*nq2):
#        pop_nac02[t-1]    = pop_nac02[t-1] + cg0[i]
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dx = np.zeros(n1c2v*n2c2v)
#    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dy = np.zeros(n1c2v*n2c2v)
#    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dz = np.zeros(n1c2v*n2c2v)
#    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1c2v*n2c2v):
#        pop_nac02_fg_c2v[t-1]    = pop_nac02_fg_c2v[t-1] + sv[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dx = np.zeros(n1d2d*n2d2d)
#    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dy = np.zeros(n1d2d*n2d2d)
#    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dz = np.zeros(n1d2d*n2d2d)
#    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1d2d*n2d2d):
#        pop_nac02_fg_d2d[t-1]    = pop_nac02_fg_d2d[t-1] + sv[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    
#    s = 3
#    cg0 = aux2[(s-1)*ng:s*ng]
#    cg = cg0.reshape(nq2,nq1)
#    cg = np.transpose(cg)
#    
#    pop_nac12_c2v[t-1] = cg[95-1,115-1] #Making a subgrid
#    pop_nac12_d2d[t-1] = cg[73-1,99-1] #Making a subgrid
#    pop_nac12_td [t-1] = cg[73-1,115-1] #Making a subgrid
#    for i in np.arange(nq1*nq2):
#        pop_nac12[t-1] = pop_nac12[t-1] + cg0[i]
#    #########################################################################
#    #c2v region
#    sub = cg[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv = np.zeros(n1c2v*n2c2v)
#    sv = sub.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dx = np.zeros(n1c2v*n2c2v)
#    sv_dx = subx.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dy = np.zeros(n1c2v*n2c2v)
#    sv_dy = suby.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1ic2v-1:q1fc2v,q2ic2v-1:q2fc2v] #Making a subgrid
#    sv_dz = np.zeros(n1c2v*n2c2v)
#    sv_dz = subz.reshape(n1c2v*n2c2v, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1c2v*n2c2v):
#        pop_nac12_fg_c2v[t-1]    = pop_nac12_fg_c2v[t-1] + sv[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    #########################################################################
#    #D2d region
#    sub = cg[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv = np.zeros(n1d2d*n2d2d)
#    sv = sub.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    sv1 = sv / (np.linalg.norm(sv, ord=1)) #Normalizing the subgrid vector
#    
#    subx = tdm01x[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dx = np.zeros(n1d2d*n2d2d)
#    sv_dx = subx.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    suby = tdm01y[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dy = np.zeros(n1d2d*n2d2d)
#    sv_dy = suby.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#    subz = tdm01z[q1id2d-1:q1fd2d,q2id2d-1:q2fd2d] #Making a subgrid
#    sv_dz = np.zeros(n1d2d*n2d2d)
#    sv_dz = subz.reshape(n1d2d*n2d2d, order = 'F') #Vectorizing the subgrig
#
#    for i in np.arange(n1d2d*n2d2d):
#        pop_nac12_fg_d2d[t-1]    = pop_nac12_fg_d2d[t-1] + sv[i]
#    del(sv,sub,sv1,sv_dx,sv_dy,sv_dz,subx,suby,subz)
#    ######################################################################### 
#
#
#np.savetxt('population-nac-fg-0.dat',pop_nac01)
#np.savetxt('population-nac-fg-1.dat',pop_nac02)
#np.savetxt('population-nac-fg-2.dat',pop_nac12)
#np.savetxt('population-nac-c2v-region-0.dat',pop_nac01_fg_c2v)
#np.savetxt('population-nac-c2v-region-1.dat',pop_nac02_fg_c2v)
#np.savetxt('population-nac-c2v-region-2.dat',pop_nac12_fg_c2v)
#np.savetxt('population-nac-d2d-region-0.dat',pop_nac01_fg_d2d)
#np.savetxt('population-nac-d2d-region-1.dat',pop_nac02_fg_d2d)
#np.savetxt('population-nac-d2d-region-2.dat',pop_nac12_fg_d2d)
#np.savetxt('population-nac-c2v-point-0.dat',pop_nac01_c2v)
#np.savetxt('population-nac-c2v-point-1.dat',pop_nac02_c2v)
#np.savetxt('population-nac-c2v-point-2.dat',pop_nac12_c2v)
#np.savetxt('population-nac-d2d-point-0.dat',pop_nac01_d2d)
#np.savetxt('population-nac-d2d-point-1.dat',pop_nac02_d2d)
#np.savetxt('population-nac-d2d-point-2.dat',pop_nac12_d2d)
#
#
#
#outputname = 'population_NAC_on_C2v_point.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop_nac01_c2v, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop_nac02_c2v, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop_nac12_c2v, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'population_NAC_on_D2d_point.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop_nac01_d2d, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop_nac02_d2d, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop_nac12_d2d, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'population_NAC_on_Td_point.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop_nac01_td, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop_nac02_td, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop_nac12_td, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'population_NAC_full_grid.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop_nac01, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop_nac02, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop_nac12, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'population_NAC_on_C2v_region.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop_nac01_fg_c2v, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop_nac02_fg_c2v, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop_nac12_fg_c2v, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
#
#outputname = 'population_NAC_on_D2d_region.png'
#fig, ax = plt.subplots(3, figsize=(6,4))
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh02+0.005, 'b', time*t_au, coh12+0.01, 'r', lw=1.5)
#ax[2].plot(time*t_au, pop_nac01_fg_d2d, 'k', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[2].set_ylabel('', fontsize=15)
#ax[2].set_xlabel('Time in femtoseconds', fontsize=15)
#ax[2].legend(['St0'])
#ax[2].set_xlim([0,time[-1]*t_au])
#ax[2].locator_params(axis='x', nbins=14)
## ax[2].axes.get_xaxis().set_visible(False)
#
#ax[1].plot(time*t_au, pop_nac02_fg_d2d, 'b', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[1].set_xlim(0,time[-1]*t_au)
#ax[1].set_ylabel('population', fontsize=15)
#ax[1].set_xlabel('', fontsize=15)
#ax[1].legend(['St1'])
#ax[1].locator_params(axis='x', nbins=12)
#ax[1].set_xlim([0,time[-1]*t_au])
#ax[1].axes.get_xaxis().set_visible(False)
#
#ax[0].plot(time*t_au, pop_nac12_fg_d2d, 'r', lw=1.5)
##ax.plot(time*t_au, coh01, 'k', time*t_au, coh011, '--k', lw=1.5)
#ax[0].set_xlim(0,time[-1]*t_au)
#ax[0].set_ylabel('', fontsize=15)
#ax[0].set_xlabel('', fontsize=15)
#ax[0].legend(['St2'])
#ax[0].locator_params(axis='x', nbins=12)
#ax[0].set_xlim([0,time[-1]*t_au])
#ax[0].axes.get_xaxis().set_visible(False)
#
#fig.tight_layout()
#fig.savefig(outputname)
#plt.close('all')
#
