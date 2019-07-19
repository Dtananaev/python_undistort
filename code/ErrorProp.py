import DistHandler
import math

class ErrorProp(object):

    def __init__(self,fx,fy,cx,cy,k0,k1,k2,k3,k4):

        # Intrinsics
        self.fx=fx
        self.fy=fy
        self.cx=cx
        self.cy=cy

        # Distortion Parameters
        self.k0=k0
        self.k1=k1
        self.k2=k2
        self.k3=k3
        self.k4=k4

        self.dh=DistHandler.DistHandler(fx,fy,k0,k1,k2,k3,k4)

    # Input is arbitary 3d ppoint or ray in cam COS
    def thetaProjection(self, X, Y, Z):

        theta=0.0
        norm=math.sqrt(pow(X,2)+pow(Y,2)+pow(Z,2))
        if norm<1e-9:
            theta=0.0
        else:
            theta=math.acos(Z/norm)
      
        # Compute phi ... just for debug purposes 
        #phi=math.tan(Y/X)

        # Now scale coordinates
        norm=math.sqrt(pow(X,2)+pow(Y,2))
        if(fabs(norm)>1e-9):
            x_z1=X/norm*theta
            y_z1=Y/norm*theta
        else:
            x_z1=0
            y_z1=0

        #return [x_z1, y_z1, theta, phi]
        return [x_z1, y_z1]


    # Input undistorted x,y  output object points in cam coords (unit vector)
    def invThetaProjection(self,x_z1,y_z1):

        # Theta
        theta=math.sqrt(pow(x_z1,2)+pow(y_z1,2))

        # Phi not needed in this representation for 
        phi=math.atan2(y_z1,x_z1)

        # Scale x and y on unit length
        X=x_z1/theta
        Y=y_z1/theta

        # Z Coordinate
        Z=1.0/math.tan(theta)

        # Normalize vector
        norm=math.sqrt(pow(X,2)+pow(Y,2)+pow(Z,2))
        X/=norm
        Y/=norm
        Z/=norm

        return [X, Y, Z, theta, phi]
        
        
    # K*x
    def projection(self,x_z1,y_z1):
        
        x=x_z1*self.fx + self.cx
        y=y_z1*self.fy + self.cy

        return [x,y]


    # Kinv*x
    def unprojection(self,x,y):
        
        x_z1=(x-self.cx)/self.fx
        y_z1=(y-self.cy)/self.fy

        return [x_z1,y_z1]

    def XYZ2xy(self,X,Y,Z):
        
        # Fisheye
        [x_z1_clean, y_z1_clean]=self.thetaProjection(X,Y,Z)

        # Distortion
        [x_z1_dirty, y_z1_dirty]=self.dh.dist(x_z1_clean, y_z1_clean)

        # Pinhole projection (more scaling)
        [x, y]=self.projection(x_z1_dirty, y_z1_dirty)

        return [x, y]
        
    def xy2XYZ(self, x, y):
        
        # Inverse pinhole (more scaling)
        [x_z1_dirty, y_z1_dirty]=self.unprojection(x, y)
        
        # Distortion
        [x_z1_clean, y_z1_clean]=self.dh.undist(x_z1_dirty, y_z1_dirty)

        # Fisheye
        [X, Y, Z, theta, phi]=self.invThetaProjection(x_z1_clean,y_z1_clean)

        return [X, Y, Z, theta, phi]

    # This takes distorted x,y on z=1 (normalized)
    # Outputs the corresponding angular angle of \phi and \theta 
    def getError(self, x_dirty, y_dirty, sigmax, sigmay):

        # Transform dirty coordinates to theta and phi ... I guess till here it should be debugged
        [X, Y, Z, theta, phi]=self.xy2XYZ( x_dirty, y_dirty)

        print('Objectspace entities')
        print('X Y Z', X,Y,Z)
        print('theta phi', theta/math.pi*180.0,phi/math.pi*180.0)

        # Helpers
        theta2=theta*theta
        theta3=theta2*theta
        theta4=theta3*theta
        theta5=theta4*theta
        theta6=theta5*theta
        theta7=theta6*theta
        cosphi=math.cos(phi)
        sinphi=math.sin(phi)
        cosphi2=cosphi*cosphi
        sinphi2=sinphi*sinphi

        # Compute partial derivatives ... could be further simplified  
        df_dt = self.fx * ( cosphi * ( 1.0 + 3.0*self.k0*theta2 + 5.0*self.k1*theta4 + 7.0*self.k4*theta6 ) \
                + 2.0*theta * ( 2.0*self.k2*cosphi*sinphi + self.k3*( 1.0 + 2.0*cosphi2  ) ) ) # :=a 
        df_dp = self.fx * ( -sinphi * ( theta + self.k0*theta3 + self.k1*theta5 + self.k4*theta7 ) \
                + theta2 * ( 2.0*self.k2*(cosphi2 - sinphi2 ) - self.k3 *4.0 *cosphi *sinphi ) ) # :=b

        dg_dt = self.fy * ( sinphi * ( 1.0 + 3.0*self.k0*theta2 + 5.0*self.k1*theta4 + 7.0*self.k4*theta6) \
                + 2.0*theta * ( self.k2*(1.0 + 2.0*sinphi2) + 2.0*self.k3*sinphi*cosphi ) ) # :=c
        dg_dp = self.fy * ( cosphi * (theta + self.k0 * theta3 + self.k1*theta5 + self.k4*theta7 ) \
                + theta2 *( 4.0 *self.k2*cosphi*sinphi + 2.0*self.k3*(cosphi2-sinphi2) ) ) # :=d

        print('df_dt', df_dt)
        print('df_dp', df_dp)
        print('dg_dt', dg_dt)
        print('dg_dp', dg_dp)
        
        # Buggy: We get square root of negative terms. How to proceed here? Is phi and theta captured correctly? 
            
        print(1.0 - (dg_dt*dg_dt) / (df_dt*df_dt) )
        print (dg_dp*dg_dp - dg_dt*dg_dt / df_dt/df_dt * df_dp*df_dp )

        print(1.0 - df_dp*df_dp)

        # Do the error propagation ... could be further simplified  
        #if ( math.fabs(1.0 - dg_dt*dg_dt/df_dt/df_dt)<1e-4):
            #dummy1 = math.sqrt( dg_dp*dg_dp - dg_dt*dg_dt / df_dt/df_dt * df_dp*df_dp )
            sigmaphi = math.sqrt(math.fabs(( sigmay*sigmay - dg_dt*dg_dt/df_dt/df_dt*sigmax*sigmax ) / ( dg_dp*dg_dp - dg_dt*dg_dt / df_dt/df_dt * df_dp*df_dp )))
            sigmatheta = math.sqrt ( math.fabs( (sigmax*sigmax - ( df_dp*df_dp * sigmaphi*sigmaphi ) ) / ( df_dt*df_dt ) ) )

        return [ sigmatheta, sigmaphi ]  
        
