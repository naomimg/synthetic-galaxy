#!/usr/bin/python
# Populate a galaxy with clouds and turn on,off extra non-circular velocities
# Based on simulate_pop.py written for McClure-Griffiths et al (2012)
# N McClure-Griffiths 10 June 2016


import os, commands, re
import pylab as pylab
import math as math
from astropy.io import fits
from astropy import wcs
import pylab
import scipy as scipy
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
  
def sigmar(R):
    R0 = 500.0 # in pc
    # Constant sigmaR in units of kpc^-2
    sigmaR = 500.0
    # Or sigmaR that is exponential
    sigmaR = 500.0*exp(-R/R0)

def rthetaz_pop():
# Creates a 1.e3 x 1.e3 x 1.e3 array of random numbers
    randa=scipy.random.random_sample(size=(10000,3))

# IAU Galactic constants
    R0 = 8500.0
    V0 = 220.0

# R is the radius, phi is the azimuth 
# Populate the cylinder
    phi = -math.pi + 2*math.pi*randa[:,0]
    R = 25000.*randa[:,1]
    z = -100. + 200*randa[:,2]

    # Transform the the (R, phi, z) into cartesian.  To be consistent with l,b standards:
    # coordinates (X,Y,Z), centred on the Galactic Centre. +Y should point towards the Sun
    # +X should point in the direction of rotation, and
    # +Z should point above the plane (positive b, North Galactic Pole).
    # Defined this way, phi is zero at the Sun-Centre line and increases in the direction of Galactic rotation.
    x = R * np.sin(phi) 
    y = R * np.cos(phi) 
    z = z
    
    print "Min x: %.2f   Max x: %.2f" % (min(x),max(x))
    print "Min y: %.2f   Max y: %.2f" % (min(y),max(y))
    print "Min z: %.2f   Max z: %.2f" % (min(z),max(z))

    # Make a spiral of the model used for velocity field
    X1 = np.arange(-10000, 10000, 100)
    Y1 = np.arange(-10000, 10000, 100)
    X1, Y1 = np.meshgrid(X1, Y1)
    R1 =np.sqrt(X1**2 + Y1**2)
    Z1 = 1/np.tan(math.pi/2) * R1

    # To Transform to the Galactic positions Y is along the Sun-GC line increasing away from the GC
    yprime=R0-1.0*y  # Shift to the position of the Sun
    d = np.sqrt(yprime**2 + x**2 + z**2)
       
    lat = np.degrees(np.arcsin(z/d))
    lon = np.degrees(np.arctan2(x,yprime))
    for i in range(len(R)):
        if R[i] >= R0 * abs(np.sin(np.radians(lon[i]))) -10. and R[i] <= R0 * abs(np.sin(np.radians(lon[i]))) +10.0:
            print " %.2f  %.2f %.2f %.2f %.2f" % (R[i], lon[i], d[i], x[i],y[i])
            print lon[i],np.degrees(phi[i])


    # Get the LSR velocity
    # NOTE phi is in radians whereas lon, and lat are in degrees
    vel = vlsr(R,z,d,lon,lat,phi)   # Standard Galactic rotation
    #vel = vnoncirc(R,z,d,lon,lat,phi)
    #vel = vel2(R,z,d,lon,lat,phi)    # Burton & Liszt (1978) tilted disk

    # Plot the distribution of points on a 3d plot
    fig = pylab.figure(2)
    pylab.clf()
    ax = Axes3D(fig)
    ax.scatter(x,y,z,marker='o',s=40, c=vel)
    surf = ax.plot_surface(X1, Y1, Z1,linewidth=0, antialiased=False)

    ax.set_xlabel('X (pc)')
    ax.set_ylabel('Y (pc)')
    ax.set_zlabel('Z (pc)')

    # Plot the x, y disk
    pylab.figure(1)
    pylab.clf()
    CS=pylab.scatter(x,y,s=20,c=vel)
    pylab.clim(-250,250)
    pylab.set_cmap('jet')
    CB=pylab.colorbar(CS)
    pylab.xlim(-25000,25000)
    pylab.ylim(-25000,25000)
    #pylab.contour(x,y,xi)
    pylab.xlabel("X (pc)")
    pylab.minorticks_on()
    pylab.ylabel("Y (pc)")
    pylab.savefig("sim_xyv.pdf")

    # Plot the l,b,v clouds
    pylab.figure(2)
    pylab.clf()
    CS=pylab.scatter(lon,lat,s=40,c=vel)
    pylab.clim(-250,250)
    pylab.set_cmap('jet')
    CB=pylab.colorbar(CS)
    pylab.xlim(180,-180)
    pylab.ylim(-5,5)
    pylab.xlabel("GALACTIC Longitude (deg)")
    pylab.ylabel("GALACTIC Latitude (deg)")
    pylab.savefig("sim_lb.pdf")

    return lon,lat,vel

# Rotation plus expansion from Liszt & Burton (1992)
def theta_vel(R):
    print np.shape(R)
    v=[]
    for i in range(len(R)):
        rkpc = R[i]/1000.
        a = 184 * rkpc**0.1
        if R[i]<400.:
            b = 190.0/400. * R[i]
        elif R[i]>=400. and R[i]< 900.0:
            b = 190.0
        else:
            b = 190.0 - 190.0/1400. * R[i]
        v.append(a + b)
        #        print a, b, v[i]
    vel=np.array(v)
    return vel

def vlsr(R,z,d,lon,lat,phi):
    V0 = 220.0 # km/s
    R0= 8500.0 # kpc
    vw = 0. # A wind-velocity with respect to the Galactic Centre
    vr = 0.0 # Default radial velocity component
    sinb=np.sin(np.radians(lat))
    cosb=np.cos(np.radians(lat))
    sinl=np.sin(np.radians(lon))
    cosl=np.cos(np.radians(lon))
    print min(phi), max(phi), min(np.radians(lon))
    coslphi=np.cos(np.radians(lon)+phi)
    print min(coslphi), max(coslphi)
    vz = vw * z/R
    vz=np.zeros(len(z))
#    vtheta = 221.6 - 0.44 * R/1000. # Fich, Blitz & Stark
    vtheta = V0

    vl=[]
    for i in range(len(R)):
        coeff = R0 * sinl[i] * (vtheta/R[i] - V0/R0)
        vsun = V0 * sinl[i]*cosb[i]
# Try experimenting with some radial motions
        if abs(sinl[i]) >0.45 and abs(sinl[i])<0.8:        
             vr = 0.
             vtheta = V0+10.
#        else:
#             vr = 0.0
#             vtheta = V0
        #vr = vw * (1 - (z[i]/R[i])**2)**0.5  # A wind model type vr
        #v = -1* vr * cosb[i]*coslphi[i]
#        if R[i] >= R0 * abs(sinl[i]) -2. and R[i] <= R0 * abs(sinl[i]) +2.:
#            print R[i], lon[i], d[i]
#            print lon[i],np.degrees(phi[i]),coslphi[i]
        v = (R0 * sinl[i] * (vtheta/R[i] - V0/R0) - vr * coslphi[i]) * cosb[i] + vz[i] * sinb[i]

        #print "%.2f %.2f %.2f %.2f %.2f %.2f %.2f" % (lon[i], lat[i], d[i], coeff, vr*coslphi[i], vz[i]* sinb[i], v)
        vl.append(v)
    vlsr=np.array(vl)
    return vl

def vnoncirc(R,z,d,lon,lat,phi):
# Attempt to put in Burton (1971) formalism for spiral potential
# Here phi is equivalent to Burton's theta and eta is Burton's phi, which relates to the spiral pattern
    V0 = 220.0 # km/s
    R0=8500.0 # kpc
    vw = 0. # A z-velocity with respect to the Galactic Centre
    sinb=np.sin(np.radians(lat))
    cosb=np.cos(np.radians(lat))
    sinl=np.sin(np.radians(lon))
    cosl=np.cos(np.radians(lon))
    coslphi=np.cos(np.radians(lon)+phi)
    cos90lphi=np.cos(np.radians(90.0) - (np.radians(lon)*phi))

#   Define some parameters related to spiral streaming
    t1 = -1e-5          # tan i = t1*rg +t0 pitch angle definition
    t0 = 0.35
    rnaught = 2000.          # Starting position for the spiral function
    pitch = 10.0             # Pitch angle in Roberts & Haussman
    itan = np.tan(np.radians(pitch))
    j=5 

    ar = 15.0   # magntitude of the radial stream
    at = 15.0   # magnitude of the azimuthal streaming
    vz = vw * z/R
    vz=np.zeros(len(z))
#    vtheta = 221.6 - 0.44 * R/1000. # Fich, Blitz & Stark
    vtheta = V0
 
    vl=[]
    for i in range(len(R)):
        coeff = R0 * sinl[i] * (vtheta/R[i] - V0/R0)
        vsun = V0 * sinl[i]*cosb[i]
        eta =  2 * math.log ((1 + (R[i]/rnaught)**j)) / (j*itan)  # Roberts & Haussman 
        #eta = -2. / t0 * math.log((R[i]/(t1*R[i] + t0)) * ((t1*rnaught + t0)/rnaught))
        ur = ar * np.cos (2*phi[i] - eta) 
        ut = at * np.sin (2*phi[i] - eta)
        #print ut, ur 

        #v = -1* vr * cosb[i]*coslphi[i]
        #v = (R0 * sinl[i] * (vtheta/R[i] - V0/R0) - vr * coslphi[i]) * cosb[i] + vz[i] * sinb[i]
        v = (R0 * sinl[i] * (vtheta/R[i] - V0/R0) + ut * cos90lphi[i] + ur * coslphi[i]) * cosb[i] + vz[i] * sinb[i]
        #print "%.2f %.2f %.2f %.2f %.2f %.2f %.2f" % (lon[i], lat[i], d[i], coeff, vr*coslphi[i], vz[i]* sinb[i], v)
        vl.append(v)

    vlsr=np.array(vl)
    return vl    

def new_vLSR(l,v_lsr):  # Use Reid et al (2014) values
    import numpy as np
    Uo_IAU = 10.27            # km/s precessed to J2000
    Vo_IAU = 15.32
    Wo_IAU =  7.74
# Modern Uo, Vo, Wo values from Reid et al (2014)
    Uo = 10.00                # km/s
    Vo =  5.25
    Wo =  7.17

# Assume b=0, so
    cos_b = 1.0
    sin_b = 0.0
    v_newlsr=np.zeros(len(v_lsr))
    v_out = np.zeros(len(v_lsr))
    for i in range(len(v_lsr)):
        cos_l = np.cos(np.radians(l[i]))
        sin_l = np.sin(np.radians(l[i]))
        v_helio = v_lsr[i] - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b -  Wo_IAU*sin_b
        v_newlsr[i] = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b +  Wo*sin_b

    return v_newlsr


# Use Burton & Liszt (1978) to define an inclined disk with its own motions
def vel2(R,z,d,lon,lat,phi):
    hd = 0.100 # kpc
    rd0 = 1.5 #kpc
    alpha = 13.5 #deg
    iang = 70 #deg
    dkpc = d/1000.
    rkpc = R /1000.
    cosalpha = np.cos(np.radians(alpha))
    sinalpha = np.sin(np.radians(alpha))
    cosi =  np.cos(np.radians(iang))
    sini =  np.sin(np.radians(iang))
    sinb=np.sin(np.radians(lat))
    cosb=np.cos(np.radians(lat))
    sinl=np.sin(np.radians(lon))
    cosl=np.cos(np.radians(lon))
    coslphi=np.cos(np.radians(lon+phi))
    r0=8.5 # kpc
    theta0=230.
    vw=150.0
    vw=0.0
    vz = vw * z/R
    zd = r0 * cosi + dkpc * (sini*(sinb*cosalpha + cosb*sinl*sinalpha) - cosi*cosb*cosl)
    rd = np.sqrt(dkpc**2 + r0**2 - 2*r0*dkpc*cosl*cosb - zd**2)

    v=[]
    for i in range(len(R)):
        if rd[i] < rd0 and abs(zd[i])<(3*hd):  # The B&L78 double condition
        #if rd[i] < 30*rd0 and abs(zd[i])<(30*hd):     # Adjusted the conditions to include all clouds
            pid = 170 * (1 - math.exp(-1*rd[i]/0.07))
            if rd[i] <=0.85:
                thetad = 180*(1-math.exp(-1*rd[i]/0.2))
            else:
                thetad = 180 * (1-math.exp(-1*(1.7-rd[i])/0.2))
                print "In the tilted disk", lon[i],lat[i],R[i]
            tmp=pid * (rd[i]*rd[i] - r0*(r0 - dkpc[i]*cosb[i]*cosl[i] - zd[i]*cosi)) / (rd[i] * dkpc[i]) - thetad * r0 * sini*(sinb[i]*sinalpha - cosb[i]*sinl[i]*cosalpha)/rd[i] - theta0 *sinl[i]*cosb[i]
        else:
            tmp = r0 * (theta0/rkpc[i] - theta0/r0) * sinl[i]*cosb[i]   # Standard rotation
                                                                        # Or add in a outward velocity component, vw
            vr = vw * (1 - (z[i]/R[i])**2)**0.5
            tmp = (r0 * sinl[i] * (theta0/rkpc[i] - theta0/r0) - vr * coslphi[i]) * cosb[i] + vz[i] * sinb[i]

            #        if lat[i]<0:
            #tmp += 10.*coslphi[i]*cosb[i] + 3000. * sinb[i]
            #      else:
            #tmp += 10.*coslphi[i]*cosb[i]  + 3000. * sinb[i]
        v.append(tmp)           
        #        print lon[i], phi[i],coslphi[i]
    vel=np.array(v)    
    return vel

# Function to analyse the simulated population for selection effects.  Here we pass in the variable vw, which is the wind velocity.
def analyse():
    j=0
    # Get the model population
    #lon,lat,v=xyz_population()
    lon,lat,vel=rthetaz_pop()

    # Correct the LSR velocities in terminal velocity curve for Reid et al values
    v=new_vLSR(lon,vel)


# Let's add some peculiar motions
    Vs = 10.0
    Us = 10.0
    sin_l = np.sin(np.radians(lon))
    cos_l = np.cos(np.radians(lon))
    cos_b = np.cos(np.radians(lat))

    v2 = v + (Vs*sin_l + Us*cos_l)*cos_b 

    for i in range(len(lon)):
        if sin_l[i]<=0.45 or sin_l[i]>0.8:
            Vs = -10.0
            Us = -10.0
            v2[i] = v[i] + (Vs*sin_l[i] + Us*cos_l[i])*cos_b[i] 

    # Find out some coordinate info about the cube
    vmin = -300.
    vmax = +300.

# Let's try to find the terminal velocity
# To do this we first interpolate the points on to a regular x-y grid, then work along a given longitude
# Combine the lists into one zipped array
#    zipped = zip(lon,lat,vel)
# Sort
#    newzip=sorted(zipped,key=lambda x: x[0])

# Bin the data
    lmid,maxv,minv=bin_data(lon,v,270)
    vt=np.zeros(len(lmid))

    for i in range(len(lmid)):        
        if lmid[i]<=0.0 and lmid[i]>=-90.0:
             vt[i] = minv[i]
        elif lmid[i]>=0.0 and lmid[i]<=90.0:
             vt[i] = maxv[i]
        else:
            vt[i]=0.0



# Plot the l-v diagram
 
    pylab.figure(3)
    pylab.clf()
    
    pylab.plot(lon,vel,linestyle='None',marker='+')
    pylab.plot(lmid,vt,marker='o',linestyle='None')

    pylab.xlim(90,-90)
    pylab.ylim(-200,200)
    pylab.xlabel("GALACTIC Longitude (deg)")
    pylab.ylabel("GALACTIC Velocity (km/s)")
    pylab.title("Simulated L-V diagram")
    pylab.savefig("sim_lv.pdf")

    pylab.figure(4)
    pylab.clf()
    sinl=np.sin(lmid*3.14159/180.)
    pylab.plot(abs(sinl[vt>0]),abs(vt[vt>0]),'xb')
    pylab.plot(abs(sinl[vt<0]),abs(vt[vt<0]),'xr')
    pylab.xlim(0.3,1.0)
    pylab.ylim(0.,140.)
    pylab.xlabel("sin(lon)")
    pylab.ylabel("V_t (km/s)")
    return lmid,vt



def getT(i,l,b):
   ''' Returns the distance from the origin to a point on the 
   plane through (1,0,0) inclined angle i degrees to the z-axis in the direction
   (l,b) (azimuth,elevation).
   Written by D. McConnell (on my request!)
   '''
   ir = np.radians(i)
   lr = np.radians(l)
   br = np.radians(b)
   A = np.cos(ir)**2
   C = np.cos(ir)*np.sin(ir)
   v1 = np.cos(lr)*np.cos(br)
   v2 = np.sin(lr)*np.cos(br)
   v3 = np.sin(br)
   t = A/(A*v1 + C*v3)
   x = v1*t
   y = v2*t
   z = v3*t
   R = math.sqrt(x*x+y*y+z*z)
   return t

# Function to bin data along the x-axis given a list of x and y values
def bin_data(x, y, nbin):
# Sort the data on the basis of x
    a=zip(x,y)
    a.sort()
    #  Now unzip it to get new x and y values
    x,y=zip(*a)
# Bin the velocities, returning the number of values and the lower edge
    nhist, edges = np.histogram(x,nbin)
    
# Find the middle value
    lower = np.resize(edges, len(edges)-1)    
    xmid = lower + 0.5*np.diff(edges)
# Go through the values and average
    maxy=[]
    miny=[]
    #stdey=[]
    #rmsy=[]
    j=0
    for val in nhist:
        k=j+val
        if val==0:
            maxy.append(0.0)
            miny.append(0.0)
            #rmsy.append(0.0)
# Return the mean and the 1.65 * standard error (std dev/sqrt(n)), which indicates 90% confidence interval
        else:
            #            print j, k, val, a[j:k]
            miny.append(np.min(y[j:k]))
            maxy.append(np.max(y[j:k]))
            #stdey.append(1.65*np.std(y[j:k])/np.sqrt(len(y[j:k])))
        j=k

    return xmid,maxy,miny

