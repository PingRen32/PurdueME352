# This is a script for ME 352 Fall 2017 Project II
# Set the origin at I_12 (input point for torque)
# All calculations are kept in this single file for easier presentation for Lab report
import numpy as np
from math import sin, cos, sqrt
from numpy import arcsin, arccos
import matplotlib.pyplot as plt

# Input Given Values (mm)
R_11,R_3,R_2 = 110,200,100
R_15,R_4 = 150,250
# Input Initial Angles
# Use t for theta
t_11,t_c = 270,150
t_15,t_51 = 90,180
# Given Point P location
R_AP,R_BP = 200, 150
# Given link 2 angular velocity
omega_2 = 25 # CCW for positive k direction
# Define resolution for solution
resolution = 0.01 * 3.1415 / 180

# Shift between degree angles to radius
def degreetoradius(degree):
    radius = degree * 3.1415 / 180
    return radius
def radiustodegree(radius):
    degree = radius * 180 /3.1415
    return degree

# Turn all degree input to radius for numpy
t_11 = degreetoradius(t_11)
t_c = degreetoradius(t_c)
t_15 = degreetoradius(t_15)
t_51 = degreetoradius(t_51)

# Calculate Initial Values @ t_2 = 0
R_51 = sqrt((R_4 * R_4)-(R_15*R_15))- R_2
t_4 = 3.1415 - arcsin(R_15 / R_4)
R_c = 105
t_3 = 3.1415926+0.28934558
R_AB = R_4 # share same value (just different name for easier use)
# Calculate P point relative location
t_ABP = arccos((R_AB**2+R_BP**2-R_AP**2)/(2*R_BP*R_AB))
t_BAP = arccos((R_AB**2+R_AP**2-R_BP**2)/(2*R_AB*R_AP))

# Initialize storage
storage_R_c,storage_R_51,storage_t_3,storage_t_4 = np.zeros(360),np.zeros(360),np.zeros(360),np.zeros(360)
firstorder_s_R_51, firstorder_s_R_c, firstorder_s_t_3, firstorder_s_t_4 = np.zeros(360),np.zeros(360),np.zeros(360),np.zeros(360)
secondorder_s_R_51,secondorder_s_R_c,secondorder_s_t_3,secondorder_s_t_4 = np.zeros(360),np.zeros(360),np.zeros(360),np.zeros(360)
storage_X_p,storage_Y_p = np.zeros(360),np.zeros(360)
storage_X_p_1,storage_Y_p_1 = np.zeros(360),np.zeros(360)
storage_X_p_2,storage_Y_p_2 = np.zeros(360),np.zeros(360)
storage_V_p,storage_A_p = np.zeros(360),np.zeros(360)
storage_U_t_i,storage_U_n_i,storage_P_p = np.zeros(360),np.zeros(360),np.zeros(360)
storage_U_t_j,storage_U_n_j = np.zeros(360),np.zeros(360)

# Solve for each single theta_2
for t_2 in range(0,360,1):
    i = t_2
    # Calculate for Position
    t_2 = degreetoradius(t_2)
    delta_R_c, delta_t_3, delta_R_51, delta_t_4 = 1,1,1,1
    while (abs(delta_t_3) >= resolution) & (abs(delta_R_c) >= resolution):
        a1 = [[cos(t_c),R_3*sin(t_3)],[sin(t_c),-R_3*cos(t_3)]]
        b1 = [[R_2*cos(t_2)+R_3*cos(t_3)-R_c*cos(t_c)],[R_2*sin(t_2)+R_3*sin(t_3)-R_c*sin(t_c)-R_11]]
        res1 = np.linalg.solve(a1,b1)
        delta_R_c = res1.item(0)
        delta_t_3 = res1.item(1)
        t_3 += delta_t_3
        R_c += delta_R_c

    while (abs(delta_t_4) >= resolution) & (abs(delta_R_51) >= resolution):
        a2 = [[cos(t_51),R_4*sin(t_4)],[sin(t_51),-R_4*cos(t_4)]]
        b2 = [[R_2*cos(t_2)+R_4*cos(t_4)-R_51*cos(t_51)],[R_2*sin(t_2)+R_4*sin(t_4)-R_51*sin(t_51)-R_15]]  
        res2 = np.linalg.solve(a2,b2)
        delta_R_51 = res2.item(0)
        delta_t_4 = res2.item(1)
        t_4 += delta_t_4
        R_51 += delta_R_51
    
    storage_R_51[i] = R_51
    storage_R_c[i] = R_c
    storage_t_3[i] = radiustodegree(t_3)
    storage_t_4[i] = radiustodegree(t_4)
    
    # Calculation for first and second order kinematic coefficients
    # First set of calculation
    a3 = [[cos(t_c),R_3*sin(t_3)],[sin(t_c),-R_3*cos(t_3)]]
    b3 = [[-R_2*sin(t_2)],[R_2*cos(t_2)]]
    res3 = np.linalg.solve(a3,b3)
    firstorder_R_c = res3.item(0)
    firstorder_t_3 = res3.item(1)    
    
    a4 = [[cos(t_c),R_3*sin(t_3)],[sin(t_c),-R_3*cos(t_3)]]
    b4 = [[-R_2*cos(t_2)-R_3*cos(t_3)*firstorder_t_3*firstorder_t_3],[-R_2*sin(t_2)-R_3*sin(t_3)*firstorder_t_3*firstorder_t_3]]
    res4 = np.linalg.solve(a4,b4)
    secondorder_R_c = res4.item(0)
    secondorder_t_3 = res4.item(1) 

    # Second set of calculation 
    a5 = [[cos(t_51),R_4*sin(t_4)],[sin(t_51),-R_4*cos(t_4)]]
    b5 = [[-R_2*sin(t_2)],[R_2*cos(t_2)]]
    res5 = np.linalg.solve(a5,b5)
    firstorder_R_51 = res5.item(0)
    firstorder_t_4 = res5.item(1)    
    
    a6 = [[cos(t_51),R_4*sin(t_4)],[sin(t_51),-R_4*cos(t_4)]]
    b6 = [[-R_2*cos(t_2)-R_4*cos(t_4)*firstorder_t_4*firstorder_t_4],[-R_2*sin(t_2)-R_4*sin(t_4)*firstorder_t_4*firstorder_t_4]]
    res6 = np.linalg.solve(a6,b6)
    secondorder_R_51 = res6.item(0)
    secondorder_t_4 = res6.item(1)   
      
    # Storage the first order to string for plotting    
    firstorder_s_R_51[i] = firstorder_R_51
    firstorder_s_R_c[i] = firstorder_R_c
    firstorder_s_t_3[i] = firstorder_t_3
    firstorder_s_t_4[i] = firstorder_t_4  
    # Storage the second order to string for plotting
    secondorder_s_R_51[i] = secondorder_R_51
    secondorder_s_R_c[i] = secondorder_R_c
    secondorder_s_t_3[i] = secondorder_t_3
    secondorder_s_t_4[i] = secondorder_t_4
    
    # After getting all coefficients, proceed to point P
    t_AP = t_4 - t_BAP
    storage_X_p[i] = R_2*cos(t_2)+R_AP*cos(t_AP)
    storage_Y_p[i] = R_2*sin(t_2)+R_AP*sin(t_AP)
    storage_X_p_1[i] = -R_2*sin(t_2)-R_AP*sin(t_AP)*firstorder_t_4
    storage_Y_p_1[i] = R_2*cos(t_2)+R_AP*cos(t_AP)*firstorder_t_4
    storage_X_p_2[i] = -R_2*cos(t_2)-R_AP*cos(t_AP)*(firstorder_t_4**2)-R_AP*sin(t_AP)*secondorder_t_4
    storage_Y_p_2[i] = -R_2*sin(t_2)-R_AP*sin(t_AP)*(firstorder_t_4**2)+R_AP*cos(t_AP)*secondorder_t_4
    
    # Storage the velocity of point P as a scaler instead of vector
    storage_V_p[i] = sqrt(storage_X_p_1[i]**2+storage_Y_p_1[i]**2) / 1000 * omega_2
    # Storage the acceleration of point P as scaler
    # In this case, there is no angular acceleration only angular
    storage_A_p[i] = sqrt(storage_X_p_2[i]**2+storage_Y_p_2[i]**2) / 1000 / 1000 *(omega_2**2)
    R_F_1 = sqrt(storage_X_p_1[i]**2+storage_Y_p_1[i]**2)
    storage_U_n_i[i] = storage_Y_p_1[i] / R_F_1
    storage_U_n_j[i] = storage_X_p_1[i] / R_F_1
    storage_U_t_i[i] = storage_X_p_1[i] / R_F_1
    storage_U_t_j[i] = storage_Y_p_1[i] / R_F_1
    # Radius of curvature of the path of P
    storage_P_p[i] = (R_F_1**3)/(storage_X_p_1[i]*storage_Y_p_2[i]-storage_Y_p_1[i]*storage_X_p_2[i])/1000
    

DAT = np.column_stack((firstorder_s_R_51,firstorder_s_R_c,firstorder_s_t_3,firstorder_s_t_4 \
                      ,secondorder_s_R_51,secondorder_s_R_c,secondorder_s_t_3,secondorder_s_t_4 \
                      ,storage_X_p,storage_X_p_1,storage_X_p_2 \
                      ,storage_Y_p,storage_Y_p_1,storage_Y_p_2 \
                      ,storage_V_p,storage_A_p \
                      ,storage_U_n_i,storage_U_n_j \
                      ,storage_U_t_i,storage_U_t_j))
np.savetxt('results.csv', DAT, delimiter=",", fmt="%s") 

c_distance = max(storage_R_c)-min(storage_R_c)
print('The distance movement of C is (mm):',c_distance)

plt.figure(1)
# Plot the Position, First Order and Second Order
plt.subplot(4,3,1)
plt.plot(storage_R_51,'g')
plt.plot(storage_R_c,'r')
plt.title('Position of R_51 green and R_c red')
plt.grid(True)

plt.subplot(4,3,2)
plt.plot(storage_t_3,'g')
plt.plot(storage_t_4,'r')
plt.title('Angles of theta 3 green and 4 red')
plt.grid(True)

plt.subplot(4,3,3)
plt.ylim(-2,4)
plt.plot(firstorder_s_t_3,'r')
plt.plot(secondorder_s_t_3,'b')
plt.title('First and second order t_3 red blue')
plt.grid(True)

plt.subplot(4,3,4)
plt.ylim(-2,4)
plt.plot(firstorder_s_t_4,'g')
plt.plot(secondorder_s_t_4,'y')
plt.title('First and second order t_4 green yellow')
plt.grid(True)

plt.subplot(4,3,5)
plt.ylim(-300,500)
plt.plot(firstorder_s_R_51,'r')
plt.plot(secondorder_s_R_51,'b')
plt.title('First and second order R_51 red blue')
plt.grid(True)

plt.subplot(4,3,6)
plt.ylim(-300,500)
plt.plot(firstorder_s_R_c,'g')
plt.plot(secondorder_s_R_c,'y')
plt.title('First and second order R_c green yellow')
plt.grid(True)

# Plot the P point position, velocity and acceleration
plt.subplot(4,3,7)
plt.plot(storage_X_p,storage_Y_p)
plt.title('Positions of point P (Moving Trace)')
plt.grid(True)

plt.subplot(4,3,8)
plt.ylim(-200,300)
plt.plot(storage_X_p,'r')
plt.plot(storage_X_p_1,'g')
plt.plot(storage_X_p_2,'b')
plt.title('Position,first,second order KE of P in X rgb')
plt.grid(True)

plt.subplot(4,3,9)
plt.ylim(-200,300)
plt.plot(storage_Y_p,'r')
plt.plot(storage_Y_p_1,'g')
plt.plot(storage_Y_p_2,'b')
plt.title('Position,first,second order KE of P in Y rgb')
plt.grid(True)

plt.subplot(4,3,10)
plt.plot(storage_V_p)
plt.title('The velocity (scaler) of point P (m/s)')
plt.grid(True)

plt.subplot(4,3,11)
plt.ylim(-10,50)
plt.plot(storage_A_p)
plt.title('The acceleration (scaler) of point P (m^2/s)')
plt.grid(True)

plt.subplot(4,3,12)

plt.plot(storage_P_p)
plt.title('The radius of curvature of the path (m)')
plt.grid(True)

plt.show()