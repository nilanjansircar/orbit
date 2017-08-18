#!/usr/bin/env python

# Library Imports
import  numpy as np
import matplotlib.pyplot as plt
import sys, math
##############


# Initial Conditions and parameters
delta_t= 1.0/500 #proper time increment
phi_h= 0.0 # For simplicity we will always set phi at tau = 1/2 to zero
r_h = float(raw_input("Please enter initial radial position (in units of GM):  " )) # radius at tau = 1/2
otype=raw_input("Do you want to find circular orbit? Y/N:  ")
	
if otype == 'Y' or otype=='y':
	while r_h<=3:
		r_h = float(raw_input("For cylcic orbit r>3 GM. Please enter the radius again:  " ))
    
	vr_h= 0.0 # radial velocity at tau = 1/2
	l= r_h/math.sqrt((r_h-3))# angular momentum evaluated assuming r=constant
	vphi_h=l/(r_h ** 2) # Calculating Initial angular velocity
	total_t=2* math.pi*(r_h**2)/l # total time is set to time period in circular orbit
	print("Time period for the orbit %f GM" % total_t)
else:
	vr_h=float(raw_input("Please enter the initial radial velocity:  ")) # radial velocity at tau = 1/2
	if r_h>3.0:
		vphi_h_f=float(raw_input("Please enter initial angular velocity as a fraction of that of the circular orbit at the same radius:  "))
		vphi_h=vphi_h_f/(r_h * math.sqrt(r_h-3)) # Calculate the initial angular velocity from the input in previous line.
	else:
		vphi_h=float(raw_input("Please enter initial angular velocity (in units of (GM)^(-1)):  "))
	l=vphi_h *(r_h ** 2) # angular momentum
	total_t = float(raw_input("Please input total proper time (in units of GM):  ")) 
	


NI = int(total_t/delta_t) # number of itterations
 

################


rdata = [r_h-vr_h * delta_t/2,r_h+vr_h * delta_t/2] # Initiating radial position data set at t0 and t1=t0 + delta_t
phidata= [phi_h - vphi_h * delta_t/2, phi_h + vphi_h * delta_t/2] # Initiating angular position data set at t0 and t1=t0 + delta_t

if r_h<=0 or rdata[0]<=0 or rdata[1]<=0: sys.exit("You are at the singularity!!! Try some other initial condition")

################## Solving the difference equation to get the solution
i=1 # counter
while i < NI:
	rnew = 2 * rdata[i] - rdata[i-1] + (delta_t**2) * (-1/(rdata[i]**2)+(l**2)/(rdata[i]**3)-3.0 *(l**2)/(rdata[i]**4)) #Calculate r(t_(i+1)) from r(t_i) and r(t_(i-1))
	if rnew<0: # Check if r is close to the singularity, break the loop if it has reached. Also r cannot be less than zero.
		print("You have reached the sugilarity at proper time =%f GM" %((i+1)*delta_t))
		break 
	rdata.append(rnew)
	phinew=phidata[i] + 4 * l* delta_t/((rdata[i]+rdata[i+1])**2) # Calculate phi(t_(i+1)) from phi(t_i) , r(t_i) , rnew=r(t_(i+1))
	phidata.append(phinew)
	i+=1
######################################################

#Ceate numpy arrays
rdataN=np.array(rdata)
phidataN=np.array(phidata)
##########################################
xdata=rdataN * np.cos(phidataN) # Covert the trajectory from polar to Cartesian.
ydata=rdataN * np.sin(phidataN) # Covert the trajectory from polar to Cartesian.
# Ploting in cartesian coordinates
plt.plot(xdata,ydata,linestyle='dotted',color='black')
plt.plot(xdata[0],ydata[0],'r^',markersize=10)
plt.plot(0,0,'ro',markersize=20)
plt.ylabel("y (GM)")
plt.xlabel("x (GM)")
plt.title('Orbit in Schwarzschild Space time')
plt.show()
