
class DistHandler(object):

    
    def __init__(self,fx,fy,k0,k1,k2,k3,k4):
        
        # Distortion params
        self.k0=k0
        self.k1=k1
        self.k2=k2
        self.k3=k3
        self.k4=k4
    
        # Inverse params for analytical inverse
        self.b0=-k0
        self.b1=3.0*k0*k0-k1
        self.b2=8.0*k0*k1-12.0*k0*k0*k0-k4
        self.b3=55.0*k0*k0*k0*k0+10.0*k0*k4-55.0*k0*k0*k1+5.0*k1*k1
        self.b4=-273.0*k0*k0*k0*k0*k0+364.0*k0*k0*k0*k1-78.0*k0*k1*k1-78.0*k0*k0*k4+12.0*k1*k4
        self.b5=1428.0*k0*k0*k0*k0*k0*k0-2380.0*k0*k0*k0*k0*k1+840.0*k0*k0*k1*k1-35.0*k1*k1*k1+560.0*k0*k0*k0*k4-210.0*k0*k1*k4+7.0*k4*k4
        self.b6=-7752.0*pow(k0,7)+15504.0*pow(k0,5)*k1-7752.0*k0*k0*k0*k1*k1+816.0*k0*k1*k1*k1-3876.0*pow(k0,4)*k4+2448*k0*k0*k1*k4-136.0*k1*k1*k4-136*k0*k4*k4 
        self.b7=43263.0*pow(k0,8)-100947.0*pow(k0,6)*k1+65835.0*pow(k0,4)*k1*k1-11970.0*k0*k0*k1*k1*k1+285.0*pow(k1,4)+26334.0*pow(k0,5)*k4-23940.0*k0*k0*k0*k1*k4+3420.0*k0*k1*k1*k4+1710.0*k0*k0*k4*k4-171.0*k1*k4*k4

        # Helpers break criteria
        self.fx2_scale=fx*fx;
        self.fy2_scale=fy*fy;
        self.maxit=40 
        self.maxdelta=0.0001 # [pix*pix]


    # This assumes normalized coordinates (z=1) 
    def dist(self, x_clean, y_clean):

        # Helpers
        x2=x_clean*x_clean
        y2=y_clean*y_clean
        r2=x2+y2
        r4=r2*r2
        r6=r4*r2

        # Crunch 
        x_dirty=x_clean*(1.0 + self.k0*r2 + self.k1*r4 + self.k4*r6) + 2.0*self.k2*x_clean*y_clean + self.k3*(r2+2.0*x2)
        y_dirty=y_clean*(1.0 + self.k0*r2 + self.k1*r4 + self.k4*r6) + 2.0*self.k3*x_clean*y_clean + self.k2*(r2+2.0*y2) 

        return [x_dirty, y_dirty]

    # This assumes normalized coordinates (z=1) 
    def undist(self, x_dirty, y_dirty):

        # Initial values 
        x_clean=x_dirty
        y_clean=y_dirty

        # Compute derivatives at working point 
        numit=0
        delta=self.maxdelta+1.0

        # Iteratively get dirty params 
        while (numit<self.maxit) and (delta>self.maxdelta):

            [x_dirty_bar,y_dirty_bar]=self.dist(x_clean, y_clean)
            
            deltax=(x_dirty_bar-x_dirty)
            deltay=(y_dirty_bar-y_dirty)
        
            x_clean -= deltax
            y_clean -= deltay
        
            # Compute the delta in pix coordinates && increment
            delta= deltax*deltax*self.fx2_scale + deltay*deltay*self.fy2_scale
            numit+=1

        return [x_clean, y_clean]
    
    # Analitical inverse model
    # This assumes normalized coordinates (z=1) 
    # If you ask me: kind of poor paper writing, tangential stuff not modeled, why does tan stuff even appear?
    def undistInvPoly(self, x_dirty, y_dirty):
        
        r2=x_dirty*x_dirty+y_dirty*y_dirty
        r4=r2*r2
        r6=r2*r4
        r8=r2*r6
        r10=r8*r2
        r12=r10*r2
        r14=r12*r2
        r16=r14*r2

        rad=1.0+self.b0*r2+self.b1*r4+self.b2*r6+self.b3*r8+self.b4*r10+r12*self.b5+self.b6*r14 +self.b7*r16
            
        x_clean=x_dirty*rad
        y_clean=y_dirty*rad

        return [x_clean, y_clean]





