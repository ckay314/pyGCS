import numpy as np
from numpy.linalg import inv 
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import cm
import glob
import sys


pi = np.pi
dtor = pi/180.

def rotx(vec, ang):
# Rotate a 3D vector by ang (input in degrees) about the x-axis
    ang *= dtor
    yout = np.cos(ang) * vec[1] - np.sin(ang) * vec[2]
    zout = np.sin(ang) * vec[1] + np.cos(ang) * vec[2]
    return [vec[0], yout, zout]

def roty(vec, ang):
# Rotate a 3D vector by ang (input in degrees) about the y-axis
    ang *= dtor
    xout = np.cos(ang) * vec[0] + np.sin(ang) * vec[2]
    zout =-np.sin(ang) * vec[0] + np.cos(ang) * vec[2]
    return [xout, vec[1], zout]

def rotz(vec, ang):
# Rotate a 3D vector by ang (input in degrees) about the y-axis
	ang *= dtor
	xout = np.cos(ang) * vec[0] - np.sin(ang) * vec[1]
	yout = np.sin(ang) * vec[0] + np.cos(ang) * vec[1]
	return [xout, yout, vec[2]]

def SPH2CART(sph_in):
    r = sph_in[0]
    colat = (90. - sph_in[1]) * dtor
    lon = sph_in[2] * dtor
    x = r * np.sin(colat) * np.cos(lon)
    y = r * np.sin(colat) * np.sin(lon)
    z = r * np.cos(colat)
    return [x, y, z]
    
    


def cmecloud(ang, hin, nleg, ncirc, k, ncross, hIsLeadingEdge=True):
    # This generates a horizontal GCS shape with nose along x axis and axis in xy plane
    h = hin
    # convert from distance of nose to Thernisien h (length of leg)
    if hIsLeadingEdge: h=hin*(1.-k)*np.cos(ang)/(1.+np.sin(ang))
    
    # Compute the shell points from axis and radius
    axisPTS, crossrads, betas = shellSkeleton(ang, h, nleg, ncirc, k)
    nbp = axisPTS.shape[0]
    theta = np.linspace(0, 360*(1-1./ncross) , ncross, endpoint=True)*dtor
    
    # Put things into massive arrays to avoid real for loops
    thetaMEGA = np.array([theta]*nbp).reshape([1,-1])
    crMEGA = np.array([[crossrads[i]]*ncross for i in range(nbp)]).reshape([-1])
    betaMEGA = np.array([[betas[i]]*ncross for i in range(nbp)]).reshape([-1])
    axisMEGA = np.array([[axisPTS[i]]*ncross for i in range(nbp)]).reshape([-1,3])
    
    # Calc the cross section vect in xyz
    radVec = crMEGA*(np.array([np.sin(thetaMEGA)*np.sin(betaMEGA), np.sin(thetaMEGA)*np.cos(betaMEGA), np.cos(thetaMEGA)]))     # Add to the axis to get the full shell  
    shell = np.transpose(radVec).reshape([ncross*nbp,3])+axisMEGA
        
    return np.array(shell)
        

def shellSkeleton(alpha, h, nleg, ncirc, k):
    # Determine the xyz position of axis, cross section radius, and beta angle
    gamma = np.arcsin(k)
    
    # Calculate the leg axis
    hrange = np.linspace(1, h, nleg)
    leftLeg = np.zeros([nleg,3])
    leftLeg[:,1] = -np.sin(alpha)*hrange
    leftLeg[:,0] = np.cos(alpha)*hrange
    rightLeg = np.zeros([nleg,3])
    rightLeg[:,1] = np.sin(alpha)*hrange
    rightLeg[:,0] = np.cos(alpha)*hrange
    rLeg = np.tan(gamma) * np.sqrt(rightLeg[:,1]**2 + rightLeg[:,0]**2)
    
    legBeta = np.ones(nleg) * -alpha
    
    rightCirc = np.zeros([ncirc, 3])
    leftCirc = np.zeros([ncirc, 3])
    
    # Calculate the circle axis
    beta = np.linspace(-alpha, pi/2 , ncirc, endpoint=True)
    b = h/np.cos(alpha) # b thernisien
    rho = h*np.tan(alpha) 
    
    X0 = (rho+b*k**2*np.sin(beta))/(1-k**2)
    rc = np.sqrt((b**2*k**2-rho**2)/(1-k**2)+X0**2)
        
    rightCirc[:,1] = X0*np.cos(beta) 
    rightCirc[:,0] = b+X0*np.sin(beta)
    leftCirc[:,1] = -rightCirc[:,1]
    leftCirc[:,0] = rightCirc[:,0]
     
    # Group into a list
    # radius of cross section
    crossrads = np.zeros(2*(nleg+ncirc)-3)
    crossrads[:nleg] = rLeg[:nleg]
    crossrads[-nleg:] = rLeg[:nleg][::-1]
    crossrads[nleg:nleg+ncirc-1] = rc[1:]
    crossrads[nleg+ncirc-1:-nleg] = rc[1:-1][::-1]
    
    # beta angle
    betas = np.zeros(2*(nleg+ncirc)-3)
    betas[:nleg] = legBeta[:nleg]
    betas[-nleg:] = pi-legBeta[:nleg][::-1]
    betas[nleg:nleg+ncirc-1] = beta[1:]
    betas[nleg+ncirc-1:-nleg] = pi-beta[1:-1][::-1]
    
    # xyz of axis
    axisPTS = np.zeros([2*(nleg+ncirc)-3,3]) 
    axisPTS[:nleg,:] = rightLeg[:nleg,:]
    axisPTS[-nleg:,:] = leftLeg[:nleg][::-1,:]
    axisPTS[nleg:nleg+ncirc-1,:] = rightCirc[1:,:]
    axisPTS[nleg+ncirc-1:-nleg,:] = leftCirc[1:-1][::-1,:]

    return axisPTS, crossrads, betas
    
    
def getGCS(CMElon, CMElat, CMEtilt, height, k, ang, nleg=5, ncirc=20, ncross=30):
    cloud = cmecloud(ang*dtor, height, nleg, ncirc, k, ncross) 
    # in order (of parens) rotx to tilt, roty by -lat, rotz by lon
    #cloud = np.transpose(rotz(roty(rotx(np.transpose(cloud), CMEtilt),-CMElat),CMElon))
    cXYZ = np.transpose(cloud) 
    cXYZ = rotz(roty(rotx(cXYZ, CMEtilt), -CMElat), CMElon)       
    # all the projection stuff is gone, return in Stony Cartesian
        
    return np.transpose(cXYZ)    
    
def getShell(CMElon, CMElat, CMEtilt, height, k, ang, nleg=0, ncross=30, ncirc=15):
    # Define  in model coords
    phis   = np.linspace(0, 2*np.pi, ncross)
    thetas = np.linspace(-np.pi/2 , np.pi/2, ncirc)
    wid1   = height * np.tan(ang*dtor) / (1 + np.tan(ang*dtor))
    wid2   = k * wid1
    xyz    = [] 
    for i in phis:
        for j in thetas:
            xyz.append(np.array([height-wid1+wid1*np.cos(j), wid1*np.sin(j)*np.sin(i), wid2*np.sin(j)*np.cos(i)]))
    # Rotate into stonyhurst heliographic
    stonyXYZ = rotz(roty(rotx(np.transpose(xyz), CMEtilt), -CMElat), CMElon) 
    return np.transpose(stonyXYZ)
    
    