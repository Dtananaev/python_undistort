import DistHandler 
import ErrorProp

import numpy as np
import math
import matplotlib.pyplot as plt

def getCalib1():

    # Setting our parameters ... example copy
    fx=1460.2096
    fy=fx*1.009552
    cx=959.96967
    cy=636.4443
    k0=-0.3050208
    k1=0.13015677
    k2=0.0042051626
    k3=-0.000007903274
    k4=-0.0211596
    cols=1912
    rows=1273

    return [cols,rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] 


def getUndistMapsNumeric(cols,rows,step, fx, fy, cx, cy, k0, k1, k2, k3, k4):

    # Instance of undist model
    dh = DistHandler.DistHandler(fx, fy, k0, k1, k2, k3, k4)

    # Boostrap new array
    colsred=cols/step+1
    rowsred=rows/step+1

    mapx=np.zeros((rowsred,colsred))
    mapy=np.zeros((rowsred,colsred))

    # Make the subsampled image
    idy=0
    for y in range(0,rows,step):
        idx=0
        for x in range(0,cols,step):

            # pix -> z=1
            x_dirty_z1=(x-cx)/fx
            y_dirty_z1=(y-cy)/fy

            # dist -> undist, functionality for getting initial value from prev. optimization 
            [x_clean_z1, y_clean_z1] = dh.undist(x_dirty_z1, y_dirty_z1)

            # z=1 -> pix
            mapx[idy,idx]=x_clean_z1*fx+cx
            mapy[idy,idx]=y_clean_z1*fy+cy

            idx+=1
        idy+=1

    return [mapx, mapy]


def getUndistMapsPoly(cols,rows,step, fx, fy, cx, cy, k0, k1, k2, k3, k4):

    # Instance of undist model
    dh = DistHandler.DistHandler(fx, fy, k0, k1, k2, k3, k4)

    # Boostrap new array
    colsred=cols/step+1
    rowsred=rows/step+1

    mapx=np.zeros((rowsred,colsred))
    mapy=np.zeros((rowsred,colsred))

    # Make the subsampled image
    idy=0
    for y in range(0,rows,step):
        idx=0
        for x in range(0,cols,step):

            # pix -> z=1
            x_dirty_z1=(x-cx)/fx
            y_dirty_z1=(y-cy)/fy

            # dist -> undist, functionality for getting initial value from prev. optimization 
            [x_clean_z1, y_clean_z1] = dh.undistInvPoly(x_dirty_z1, y_dirty_z1)

            # z=1 -> pix
            mapx[idy,idx]=x_clean_z1*fx+cx
            mapy[idy,idx]=y_clean_z1*fy+cy

            idx+=1
        idy+=1

    return [mapx, mapy]

# Ordering as follows
# 00 01
# 10 11
def interpolateMap(sparsemap, step, cols, rows):
    
    # Boostrap new array
    fullmap=np.zeros((rows,cols))

    # Helpers
    colsred=np.size(sparsemap,1)
    rowsred=np.size(sparsemap,0)

    for y in range(0,rows-step):
        
        lowidy=int(math.floor(y/step))
        
        for x in range(0,cols-step):
        
            lowidx=int(math.floor(x/step))
        
            f00=sparsemap[lowidy,lowidx]   # y=0 x=0
            f01=sparsemap[lowidy,lowidx+1]   # y=0 x=1  
            f10=sparsemap[lowidy+1,lowidx]    # y=1 x=0
            f11=sparsemap[lowidy+1,lowidx+1]  # y=1 x=1  
            
            xf=float(x)/float(step)-float(lowidx)
            yf=float(y)/float(step)-float(lowidy)
            
            fullmap[y,x] = f00*(1.0-yf)*(1.0-xf) + f10*yf*(1.0-xf) + f01*(1.0-yf)*xf + f11*xf*yf  

    return fullmap


def evluateInvMappingPoly():
    
    # Load Calib
    [cols, rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] = getCalib1()

    # Map holding undistorted coordinates 
    step = 20
    [mapx, mapy] = getUndistMapsNumeric(cols, rows, step, fx, fy, cx, cy, k0, k1, k2, k3, k4)
    [gt_mapx, gt_mapy] = getUndistMapsNumeric(cols, rows, 1, fx, fy, cx, cy, k0, k1, k2, k3, k4)

    print(np.size(mapx,1), cols)
    print(np.size(mapx,0), rows)
   
    fullmapx_inter=interpolateMap(mapx, step, cols, rows)
    fullmapy_inter=interpolateMap(mapy, step, cols, rows)

    print(np.size(fullmapx_inter,1), cols)
    print(np.size(fullmapx_inter,0), rows)

    print(np.size(gt_mapx,1), cols)
    print(np.size(gt_mapx,0), rows)

    #plt.imshow(fullmapy_inter[:-step, :-step] - gt_mapy[:np.size(fullmapy_inter,0)-step,:np.size(fullmapy_inter,1)-step])
    plt.imshow(fullmapx_inter[:-step, :-step] - gt_mapx[:np.size(fullmapx_inter,0)-step,:np.size(fullmapx_inter,1)-step])
    plt.show()
    #print(fullmapy-inter)

def evaluateAnalytic():

    # Load Calib
    [cols, rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] = getCalib1()

    # Get the clean maps numeric
    [gt_mapx, gt_mapy] = getUndistMapsNumeric(cols, rows, 1, fx, fy, cx, cy, k0, k1, k2, k3, k4)
    
    # Get the clean maps analytic
    [a_mapx, a_mapy] = getUndistMapsPoly(cols,rows, 1, fx, fy, cx, cy, k0, k1, k2, k3, k4)

    plt.imshow(a_mapx - gt_mapx)
    plt.colorbar()
    plt.title('Displacement x')
    plt.contour(a_mapx - gt_mapx, levels=[-0.1,0.1])
    plt.show()

def evluateBipolarInterpolation():
    
    # Load Calib
    [cols, rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] = getCalib1()

    # Map holding undistorted coordinates 
    step =20 
    [mapx, mapy] = getUndistMapsNumeric(cols, rows, step, fx, fy, cx, cy, k0, k1, k2, k3, k4)
    [gt_mapx, gt_mapy] = getUndistMapsNumeric(cols, rows, 1, fx, fy, cx, cy, k0, k1, k2, k3, k4)

    print(np.size(mapx,1), cols)
    print(np.size(mapx,0), rows)
   
    fullmapx_inter=interpolateMap(mapx, step, cols, rows)
    fullmapy_inter=interpolateMap(mapy, step, cols, rows)

    print(np.size(fullmapx_inter,1), cols)
    print(np.size(fullmapx_inter,0), rows)

    print(np.size(gt_mapx,1), cols)
    print(np.size(gt_mapx,0), rows)

    plt.imshow(fullmapy_inter[:-step, :-step] - gt_mapy[:np.size(fullmapy_inter,0)-step,:np.size(fullmapy_inter,1)-step])
    #plt.imshow(fullmapx_inter[:-step, :-step] - gt_mapx[:np.size(fullmapx_inter,0)-step,:np.size(fullmapx_inter,1)-step])
    plt.colorbar()
    plt.title('Displacement y')
    plt.show()
    #print(fullmapy-inter)

def main_():

    # Ours
    #k0=-0.3050208
   # k1=0.13015677
   # k4=-0.0211596

    # Theirs
    k0=0.0001532
    k1=-0.00000009656
    k4=0.00000000007245

    b0=-k0
    b1=3.0*k0*k0-k1
    b4=8.0*k0*k1-12.0*k0*k0*k0-k4

    print(b0,b1,b4)

    # coefficents seem ok? Yes 



def checkFisheyeProjection():

    [cols,rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] = getCalib1()

    ep=ErrorProp.ErrorProp(fx,fy,cx,cy,k0,k1,k2,k3,k4)

    # Define test point #
    X=-1.0
    Y=1.0
    Z=1.0        
    
    norm=math.sqrt(X*X+Y*Y+Z*Z);
    X/=norm
    Y/=norm
    Z/=norm

    [x,y]=ep.XYZ2xy(X,Y,Z)

    print ( X,Y,Z, " -> ", x, y )

    [Xnew,Ynew,Znew, theta,phi]=ep.xy2XYZ(x,y)

    print ( x,y, " -> ", Xnew, Ynew, Znew )
    print ( " diff ", Xnew-X, Ynew-Y, Znew-Z )
    print ( x,y, " -> ",  theta/math.pi*180.0, phi/math.pi*180.0 )

   # [x_z1,y_z1,theta,phi]=ep.thetaProjection(X,Y,Z)
    
   # print ( X,Y,Z, " -> ", x_z1, y_z1, theta, phi )

   # [ X_new, Y_new, Z_new ]= ep.invThetaProjection(x_z1, y_z1) 

   # print ('\n', x_z1, y_z1, " -> ", X_new, Y_new, Z_new)


def testErrorProp():

    # Loading Data
    [cols,rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] = getCalib1()

    # Setting up propagation module
    ep=ErrorProp.ErrorProp(fx,fy,cx,cy,k0,k1,k2,k3,k4)

    # Get the error
    #x_dirty=1912.0/2.0
    #y_dirty=1273.0/2.0
    
    x_dirty=cx+1000.0
    y_dirty=cy+1000.0#20.10
    
    sigmax=1.0
    sigmay=1.0

    [sigmatheta, sigmaphi]=ep.getError(x_dirty,y_dirty,sigmax,sigmay)
    print('sigmatheta=', sigmatheta/math.pi*180.0)
    print('sigmaphi=', sigmaphi/math.pi*180.0)


def main():
   
    testErrorProp()
    #checkFisheyeProjection()
    #evaluateAnalytic()
    #evluateBipolarInterpolation()

def main_m():


    [cols,rows, fx, fy, cx, cy, k0, k1, k2, k3, k4] = getCalib1()
    
    # Set up the Undistortion object
    dh = DistHandler.DistHandler(fx, fy, k0, k1, k2, k3, k4)
    
    # Define clean pix
    x_clean=100.0;
    y_clean=100.0;

    # Transform to z=1
    x_clean_z1=(x_clean-cx)/fx
    y_clean_z1=(y_clean-cy)/fy
    
    # Undist print result
    [x_dirty_z1, y_dirty_z1] = dh.dist(x_clean_z1, y_clean_z1)

    # Transform into image
    x_dirty=x_dirty_z1*fx+cx
    y_dirty=y_dirty_z1*fy+cy

    # Print result
    print( x_clean, " -> ", x_dirty)
    print( y_clean, " -> ", y_dirty)

    # Distort the thing
    [newx_clean_z1, newy_clean_z1] = dh.undist(x_dirty_z1,y_dirty_z1)
    [newx_clean_z1_dir, newy_clean_z1_dir] = dh.undistInvPoly(x_dirty_z1,y_dirty_z1)

    print( x_dirty_z1, " -> ", newx_clean_z1)
    print( y_dirty_z1, " -> ", newy_clean_z1)
    
    print( x_dirty_z1, " -> ", newx_clean_z1_dir)
    print( y_dirty_z1, " -> ", newy_clean_z1_dir)

    # Transform into image
    newx_clean=newx_clean_z1*fx+cx
    newy_clean=newy_clean_z1*fy+cy
    
    newx_clean_dir=newx_clean_z1_dir*fx+cx
    newy_clean_dir=newy_clean_z1_dir*fy+cy

    print("Iterative solution:", newx_clean)
    print("Iterative solution:", newy_clean)

    print("Direct solution:", newx_clean_dir)
    print("Direct solution:", newy_clean_dir)

    # TODO: Verify the transformation with c++ code, still there is some strange stuff regarding k_i #
    # Done. works. Double check with Can coeff C++ model implementation

    # Actually what we are interested in: for a n arificial 

    # Also make model for polynomimals test accuracy throughout image

    # Now: can we save memory by sub-sampled grid and bilinear interpolation?


if __name__ == "__main__":
    main()
