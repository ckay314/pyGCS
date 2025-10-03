# |---------------------------------------------------------------------|
# |------ Things that can be edited to customize the GUI interface -----|
# |---------------------------------------------------------------------|

# Dictionary for the colors of the wireframe objects by type
# Defaults are an attempt to be color blind friends on black/white backgrounds
# GCS: standard green, Torus: teal, Spheres/Ellipses: blue, slab: orange
colorDict = {'GCS':'#9AE630', 'Torus': '#38D5BE', 'Sphere':'#2B7FFF', 'Half Sphere':'#2B7FFF', 'Ellipsoid':'#2B7FFF', 'Half Ellipsoid':'#2B7FFF', 'Slab':'#FF8904'}

# Extra colors to use if plotting more than one of the same type of WF
# order is maroon, purple, yellow
bonusColors = ['#FF2056', '#9810FA', '#FFD230']

# Number of points in the grid
# GCS:[nLeg,nCirc,nAxis], Torus:[nAxis,nCirc], Sphere/Ell[nTheta,nPhi], Slab [nx, ny, nz]
# Attempting to keep theta to lon angle and phi to lat angle throughout WOMBAT
# Everything is defined wrt to CK's Theoryland (which is Cartesian with 'nose' along x-axis
# and largest non-radial width in the z/vertical direction for antisymmetric shapes)
gridDict = {'GCS':[5,15,30], 'Torus':[25,15], 'Sphere':[50,25], 'Half Sphere':[25,25], 'Ellipsoid':[50,25], 'Half Ellipsoid':[25,25], 'Slab':[20,20,20]}




# |---------------------------------------------------------------------|
# |--------- Things that can be edited to customize the defaults -------|
# |------------- (Only change these if you understand them!) -----------|
# |---------------------------------------------------------------------|
# The default ranges for each of the parameter options
# These will be more COR appropriate than HI appropriate but will
# switch if the code see HI observations
rngDict = {'Height (Rs)':[1,25], 'Lon (deg)':[-180,180], 'Lat (deg)':[-90,90], 'Tilt (deg)':[-90,90], 'AW (deg)':[0,90], 'kappa':[0,1], 'AW_FO (deg)':[0,90], 'AW_EO (deg)':[0,90], 'deltaAx':[0,2], 'deltaCS':[0,2], 'ecc1':[0,1], 'ecc2':[0,1], 'Roll (deg)':[-90,90], 'Yaw (deg)':[-90,90], 'Pitch (deg)':[-90,90], 'Lx (Rs)':[0,25], 'Ly (Rs)':[0,25], 'Lz (Rs)':[0,25]}

# The HI values (only includes ones that change)
rngDictHI = {'Height (Rs)':[1,215], 'Lx (Rs)':[0,215], 'Ly (Rs)':[0,215], 'Lz (Rs)':[0,215]}

# The default values for each parameter (again COR appropriate)
defDict = {'Height (Rs)':10, 'Lon (deg)':0, 'Lat (deg)':0, 'Tilt (deg)':0, 'AW (deg)':30, 'AW_FO (deg)':40, 'AW_EO (deg)':15, 'kappa':0.3, 'deltaAx':1, 'deltaCS':1, 'ecc1':0.8, 'ecc2':0.5, 'Roll (deg)':0, 'Yaw (deg)':0, 'Pitch (deg)':0, 'Lx (Rs)':10, 'Ly (Rs)':4, 'Lz (Rs)':10}

# The HI values
defDictHI = {'Height (Rs)':50, 'Lx (Rs)':50, 'Ly (Rs)':10, 'Lz (Rs)':50}



# |---------------------------------------------------------------------|
# |---------------------------------------------------------------------|
# |------------- Things that probably shouldn't be touched -------------|
# |---------------------------------------------------------------------|
# |---------------------------------------------------------------------|

# imports
import numpy as np
import sys

from pyGCS import getGCS

# Dictionary for the number of points per WF type
npDict = {'GCS': 6, 'Torus': 8, 'Sphere':4, 'Half Sphere':4, 'Ellipsoid':7, 'Half Ellipsoid':7, 'Slab':9}

# Useful globals
global dtor, radeg, pi
dtor  = np.pi / 180.
radeg = 180. / np.pi
pi    = np.pi 

#|-------------------------------|
#|---- Geometry helper funcs ----|
#|-------------------------------|
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




# |---------------------------------------------------------------------|
# |--------- Parameters for each wireframe and their order -------------|
# |---------------------------------------------------------------------|
# GCS parameters - height, lon, lat, tilt, AW, kappa,  (6)
# Half torus - height, lon, lat,  tilt, AW, deltaAx, deltaCS (7) [add rot about nose?]
# Sphere/Half sphere - height, lon, lat,  AW (4)
# Ellipsoid/Half ellipsoid - height, lon, lat,  tilt, AW, epp (6)
# Slab - height, lon, lat, roll (~tilt), yaw (~lon), pitch (~lat),  Lx, Ly, Lz,  (9)

class wireframe():
    "Container to hold all the parameters and the points for a wireframe object"
    def __init__(self, WFtype):
        self.WFtype   = WFtype
        if type(WFtype) != type(None):
            if WFtype not in npDict:
                sys.exit('Unrecognized wireframe type. Exiting now... bye')
            self.nParams  = npDict[WFtype]
            self.WFcolor  = colorDict[WFtype]
            self.labels   = np.empty(self.nParams)
            self.ranges   = np.array([np.empty(2) for i in range(self.nParams)])
            self.params   = np.empty(self.nParams)
            self.gPoints  = None # set up the number of points in the grid
            self.points   = None # the WF grid points in theoryland coords
        
        if WFtype != None:
            # Set up labels
            self.setLabels()
            # Set up default ranges/parameters/grid points 
            self.setDefs()
            # Set up the number of grid points based on WF type
            self.gPoints = gridDict[WFtype]
        
            # Calc points - get the WF structure in theoryland
            self.getPoints()
        
    # |-----------------------------------------------------------------|
    # |--- Set up the text that will go above the sliders/text boxes ---|
    # |-----------------------------------------------------------------|
    def setLabels(self):
        # Set up the text that will go above the sliders/text boxes
        WFtype = self.WFtype
        if WFtype == 'GCS':
            self.labels = np.array(['Height (Rs)', 'Lon (deg)', 'Lat (deg)', 'Tilt (deg)', 'AW (deg)', 'kappa'])
        elif WFtype == 'Torus':
            self.labels = np.array(['Height (Rs)', 'Lon (deg)', 'Lat (deg)', 'Tilt (deg)', 'AW_FO (deg)', 'AW_EO (deg)', 'deltaAx', 'deltaCS'])
        elif WFtype in ['Sphere', 'Half Sphere']:
            self.labels = np.array(['Height (Rs)', 'Lon (deg)', 'Lat (deg)', 'AW (deg)'])
        elif WFtype in ['Ellipsoid', 'Half Ellipsoid']:
            self.labels = np.array(['Height (Rs)', 'Lon (deg)', 'Lat (deg)', 'Tilt (deg)', 'AW (deg)', 'ecc1', 'ecc2'])
        elif WFtype == 'Slab':
            self.labels = np.array(['Height (Rs)', 'Lon (deg)', 'Lat (deg)', 'Roll (deg)', 'Yaw (deg)', 'Pitch (deg)', 'Lx (Rs)', 'Ly (Rs)', 'Lz (Rs)'])
            
    
    # |-----------------------------------------------------------------|
    # |----------- Set up the min/max values for the sliders -----------|
    # |-------------------  and their default values -------------------|
    # |-----------------------------------------------------------------|
    def setDefs(self):
        # Start with COR defaults, will adjust if we include HI obs
        for i in range(self.nParams): 
            label = self.labels[i]
            self.ranges[i] = rngDict[label]
            self.params[i] = defDict[label]
        

    # |-----------------------------------------------------------------|
    # |--------- Calculate the shape in Stonyhurst Cartestian ----------|
    # |-----------------------------------------------------------------|
    def getPoints(self):
        WFtype = self.WFtype
        ps = self.params
        gps = self.gPoints
        
        # |---------------------------|
        # |------ GCS Wireframe ------|
        # |---------------------------|
        if WFtype == 'GCS': 
           
           # Refine param names to reuse existing code
           lon, lat, tilt, hin, k, alpha = ps[1], ps[2], ps[3], ps[0], ps[5], ps[4]*dtor
           nleg, ncirc, ncross = gps[0], gps[1], gps[2]
           
           # Assume we have leading edge height, which is not Thernisien h
           h=hin*(1.-k)*np.cos(alpha)/(1.+np.sin(alpha))
           
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
           
           # Compute the shell points from axis and radius
           nbp = axisPTS.shape[0]
           theta = np.linspace(0, 360*(1-1./ncross) , ncross, endpoint=True)*dtor
    
           # Put things into massive arrays to avoid real for loops
           thetaMEGA = np.array([theta]*nbp).reshape([1,-1])
           crMEGA = np.array([[crossrads[i]]*ncross for i in range(nbp)]).reshape([-1])
           betaMEGA = np.array([[betas[i]]*ncross for i in range(nbp)]).reshape([-1])
           axisMEGA = np.array([[axisPTS[i]]*ncross for i in range(nbp)]).reshape([-1,3])
    
           # Calc the cross section vect in xyz
           radVec = crMEGA*(np.array([np.sin(thetaMEGA)*np.sin(betaMEGA), np.sin(thetaMEGA)*np.cos(betaMEGA), np.cos(thetaMEGA)]))     # Add to the axis to get the full shell  
           shell = np.array(np.transpose(radVec).reshape([ncross*nbp,3])+axisMEGA)
           
           # Convert from theoryland to StonyCart
           cXYZ = np.transpose(shell) 
           self.points = np.transpose(rotz(roty(rotx(cXYZ, tilt), -lat), lon))  
           
        # |---------------------------|
        # |---------  Torus ----------|
        # |---------------------------|
        elif WFtype in ['Torus']:
            # Just a normal elliptical torus, not the OSPREI modified axis
            # First need to convert AWs and deltas in to lenth units
            AW, AWp = ps[4]*dtor, ps[5]*dtor
            deltaAx, deltaCS = ps[6], ps[7]
            rp = np.tan(AWp) / (1 + deltaCS * np.tan(AWp)) * ps[0]
            rr = deltaCS * rp
            Lp = (np.tan(AW) * (ps[0] - rr) - rr) / (1 + deltaAx * np.tan(AW)) 
            Lr = deltaAx * Lp
            
            x0 = ps[0] - rr - Lr
            
            # Set up parametric angles
            # theta is along axis, phi along CS
            thetas = np.linspace(-pi/2, pi/2, gps[1], endpoint=True)
            phis   = np.linspace(0, 2*pi, gps[0], endpoint=True)
            
            # Make the massive arrays
            thetaMEGA = np.array([thetas]*gps[0]).reshape([-1])
            phiMEGA = np.array([[phis[i]]*gps[1] for i in range(gps[0])]).reshape([-1])
            
            
            # Calcuate the surface
            xs = x0 + np.cos(thetaMEGA)*(Lr + rr * np.cos(phiMEGA))
            ys = rp * np.sin(phiMEGA)
            zs = np.sin(thetaMEGA)*(Lp + rr * np.cos(phiMEGA))
            
            xyz = np.array([xs, ys, zs])
            
            # Convert from theoryland to StonyCart
            self.points = np.transpose(rotz(roty(rotx(xyz, ps[3]), -ps[2]), ps[1])) 
            

        # |---------------------------|
        # |------ Spherioid F/H ------|
        # |---------------------------|
        elif WFtype in ['Sphere', 'Half Sphere']:
            if WFtype == 'Sphere':
                thetas = np.linspace(-pi, pi, gps[0], endpoint=True)
            else:
                thetas = np.linspace(-pi/2, pi/2, gps[0], endpoint=True)
            phis   = np.linspace(0, pi, gps[1], endpoint=True)
            
            # H is height of front, not center so gotta do some math
            h, aw = ps[0], ps[3]*dtor
            # get the radius of the sphere
            r = h * np.tan(aw) / (1 + np.tan(aw))
            
            # Make the massive arrays
            thetaMEGA = np.array([thetas]*gps[1]).reshape([-1])
            phiMEGA = np.array([[phis[i]]*gps[0] for i in range(gps[1])]).reshape([-1])
            # Make a sphere in cartesian
            sphere = r*np.array([np.sin(phiMEGA)*np.cos(thetaMEGA), np.sin(phiMEGA)*np.sin(thetaMEGA), np.cos(phiMEGA)])
            
            # Adjust it along the x axis
            sphere[0] += h - r
            
            # Convert from theoryland to StonyCart
            self.points = np.transpose(rotz(roty(sphere, -ps[2]), ps[1]))  
            
        # |---------------------------|
        # |------ Ellipsoid F/H ------|
        # |---------------------------|
        elif WFtype in ['Ellipsoid', 'Half Ellipsoid']: 
            if WFtype == 'Ellipsoid':
                thetas = np.linspace(-pi, pi, gps[0], endpoint=True)
            else:
                thetas = np.linspace(-pi/2, pi/2, gps[0], endpoint=True)
            phis   = np.linspace(0, pi, gps[1], endpoint=True)
            
            # H is height of front, not center so gotta do some math
            h, aw, ep, er = ps[0], ps[4]*dtor, ps[5], ps[6]
            
            # get the radius of the sphere
            r = h * np.tan(aw) / (1 + er*np.tan(aw))
            rr = er * r
            rp = ep * r
            
            # Make the massive arrays
            thetaMEGA = np.array([thetas]*gps[1]).reshape([-1])
            phiMEGA = np.array([[phis[i]]*gps[0] for i in range(gps[1])]).reshape([-1])
            # Make a sphere in cartesian
            ell = np.array([rr*np.sin(phiMEGA)*np.cos(thetaMEGA) + h-rr, rp*np.sin(phiMEGA)*np.sin(thetaMEGA), r*np.cos(phiMEGA)])
            
            # Convert from theoryland to StonyCart
            self.points = np.transpose(rotz(roty(rotx(ell, ps[3]), -ps[2]), ps[1])) 
            
        
        # |---------------------------|
        # |---------- Slab -----------|
        # |---------------------------|
        elif WFtype == 'Slab':
            h, lon, lat  = ps[0], ps[1], ps[2]
            roll, yaw, pitch = ps[3], ps[4], ps[5]
            Lx, Ly, Lz = ps[6], ps[7], ps[8]
            nx, ny, nz = gps[0], gps[1], gps[2]
            
            # Define the ranges
            xrs = np.linspace(-Lx/2,Lx/2, nx, endpoint=True)
            yrs = np.linspace(-Ly/2,Ly/2, ny, endpoint=True)
            zrs = np.linspace(-Lz/2,Lz/2, nz, endpoint=True)
            
            # Make each of the six faces and appends to long lists
            totPts = 2*(nx*ny + ny*nz + nz*nx)
            xs, ys, zs = [], [], []
            
            # front
            xs.append(np.array([xrs[-1]]*ny*nz).reshape(-1))
            ys.append(np.array([[yrs]*nz]).reshape(-1))
            zs.append(np.array([[[zrs[i]]*ny for i in range(nz)]]).reshape(-1))
            # back
            xs.append(np.array([xrs[0]]*ny*nz).reshape(-1))
            ys.append(np.array([[yrs]*nz]).reshape(-1))
            zs.append(np.array([[[zrs[i]]*ny for i in range(nz)]]).reshape(-1))
            # left
            ys.append(np.array([yrs[0]]*nx*nz).reshape(-1))
            xs.append(np.array([[xrs]*nz]).reshape(-1))
            zs.append(np.array([[[zrs[i]]*nx for i in range(nz)]]).reshape(-1))
            # right
            ys.append(np.array([yrs[-1]]*nx*nz).reshape(-1))
            xs.append(np.array([[xrs]*nz]).reshape(-1))
            zs.append(np.array([[[zrs[i]]*nx for i in range(nz)]]).reshape(-1))
            # bottom
            zs.append(np.array([zrs[0]]*nx*ny).reshape(-1))
            xs.append(np.array([[xrs]*ny]).reshape(-1))
            ys.append(np.array([[[yrs[i]]*nx for i in range(ny)]]).reshape(-1))
            # top
            zs.append(np.array([zrs[-1]]*nx*ny).reshape(-1))
            xs.append(np.array([[xrs]*ny]).reshape(-1))
            ys.append(np.array([[[yrs[i]]*nx for i in range(ny)]]).reshape(-1))
            
            xs = np.array(xs).reshape(-1)
            ys = np.array(ys).reshape(-1)
            zs = np.array(zs).reshape(-1)
            
            
            xyz = np.array([xs, ys, zs])
            
            # Rot by roll about x
            xyz = rotx(xyz,roll)
            # Rot by yaw about z
            xyz = rotz(xyz,yaw)
            # Rot by pitch about y
            xyz = roty(xyz,pitch)
            
            # add in the distance
            xyz[0] += h
            xyz = np.array(xyz)
            
            # Move in lat/lon
            self.points = np.transpose(rotz(roty(xyz, -lat), lon))  
            

            
           

        
'''wfGCS = wireframe('Torus')
#print (wfGCS.labels)
#print (wfGCS.ranges)
#print (wfGCS.params)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

pts = np.transpose(wfGCS.points)
ax.scatter(pts[0], pts[1], pts[2])
ax.set_xlim(-15,15)
ax.set_ylim(-15,15)
ax.set_zlim(-15,15)
plt.show()'''

