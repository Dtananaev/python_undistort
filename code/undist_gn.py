
class DistHandler(object):

    
    def __init__(self,fx,fy,k0,k1,k2,k3,k4):
        
        # Distortion params
        self.k0=k0
        self.k1=k1
        self.k2=k2
        self.k3=k3
    
        # Helpers break criteria
        self.fx2_scale=1.0/fx/fx;
        self.fy2_scale=1.0/fy/fy;
        self.maxit=20 
        self.maxidelta=0.01*0.01 #[pix*pix]

    # This assumes normalized coordinates (z=1) 
    def undist( x_dirty, y_dirty, k0, k1, k2, k3, k4):

        # Helpers
        x2=x_dirty*x_dirty
        y2=y_dirty*y_dirty
        r2=x2+y2
        r4=r2*r2
        r6=r4*r2

        # Crunch 
        x_clean=x_dirty*(1.0 + k0*r2 + k1*r4 + k4*r6) + 2.0*k2*x_dirty*y_dirty + k3*(r2+2.0*x2)
        y_clean=y_dirty*(1.0 + k0*r2 + k1*r4 + k4*r6) + 2.0*k3*x_dirty*y_dirty + k2*(r2+2.0*y2) 

        return [x_clean, y_clean]

    # This assumes normalized coordinates (z=1) 
    def dist( x_clean, y_clean, k0, k1, k2, k3, k4):

        # Initial values 
        x_dirty=x_clean
        y_dirty=y_clean

        # Compute derivatives at working point 
        numit=0
        delta=maxdelta+1.0
        while (numit<self.maxit) and (delta>self.maxdelta):

            [x_clean_bar,y_clean_bar]=undist(x_dirty, y_dirty, k0, k1, k2, k3, k4)
        
            deltax=(x_clean_bar-x_clean)
            deltay=(y_clean_bar-y_clean)
        
            xdirty += deltax
            ydirty += deltay
        
            # Compute the delta in pixel coordinates && increment
            delta= deltax*deltax*self.fx2_scale+deltay*deltay*fy2_scale
            numit+=1

        return [x_dirty, y_dirty]
    


