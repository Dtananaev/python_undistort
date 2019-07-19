import numpy as np 
import DistHandler 
from PIL import Image

def getCalibrationParameters():
    fx = 324.68185
    fy = 324.68185
    cx = 643.0355
    cy = 376.54645
    k0 = 0.045771178
    k1 = 0.0
    k2 = 0.0
    k3 = 0.0 
    k4 = 0.0
    return fx, fy, cx, cy, k0, k1, k2, k3, k4


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
            [x_clean_z1, y_clean_z1] = dh.undist(x_dirty_z1, y_dirty_z1)
            # z=1 -> pix
            mapx[idy,idx]=x_clean_z1*fx+cx
            mapy[idy,idx]=y_clean_z1*fy+cy

            idx+=1
        idy+=1

    return [mapx, mapy]


def main():
    imagename = "image.png"
    fx, fy, cx, cy, k0, k1, k2, k3, k4 = getCalibrationParameters()
    img = Image.open(imagename)
    img = np.array(img)
    cols = img.shape[1]
    rows =  img.shape[0]
    step = 1
    maptable = getUndistMapsPoly(cols, rows, step, fx, fy, cx, cy, k0, k1, k2, k3, k4 )
    rgbArray = np.zeros((rows,cols,3), 'uint8')
    maptable = np.array(maptable)
    print(rgbArray.shape)


    for y in range(0,rows):
        for x in range(0,cols):
            idx = int(maptable[0,y,x])
            idy = int(maptable[1,y,x])
            #print(idx)
            #print(idy)
            rgbArray[y,x,0] = img[idy, idx, 0]
            rgbArray[y,x,1] = img[idy, idx, 1]
            rgbArray[y,x,2] = img[idy, idx, 2]
    img = Image.fromarray(rgbArray)
    img.save('myimg.png')

main()
