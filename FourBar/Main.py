# This is written for ME 352 Fall 2017 HW02 2.33 part d
# Could be used as a universal solver for single simple 4 bar linkage
import numpy as np
from math import sin, cos
# Initial data
R_1, R_2, R_3, R_4 = 0,0,0,0
theta_1, theta_2, theta_3, theta_4 = 0,0,0,0
resolution = 0.01 * 3.1415 / 180
delta_theta_3, delta_theta_4 = 1, 1
# Take input for initial state could done by code but easier if already done
R_1 = float(input('Input R1:')) # In this problem, R1 = 14
R_2 = float(input('Input R2:')) # In this problem, R2 = 7
R_3 = float(input('Input R3:')) # In this problem, R3 = 10
R_4 = float(input('Input R4:')) # In this problem, R4 = 8
theta_1 = 3.14 * float(input('Input theta_1(degree):')) / 180 # In this problem, theta1 = 0
theta_2 = 3.14 * float(input('Input theta_2(degree):')) / 180 # In this problem, theta2 = 60
theta_3 = 3.14 * float(input('Input estimated theta_3(degree):')) / 180 # In this problem, theta3 is close to 11
theta_4 = 3.14 * float(input('Input estimated theta_4(degree):')) / 180 # In this problem, theta4 is close to 95
# Use turn to VLE then into first-order Taylor-series
# Then into matrixs
while (abs(delta_theta_3) >= resolution) & (abs(delta_theta_4) >= resolution):
    a = [[R_3 * sin(theta_3),-R_4 * sin(theta_4)],[-R_3 * cos(theta_3),R_4 * cos(theta_4)]]
    b = [[R_2 * cos(theta_2) + R_3 * cos(theta_3) - R_4 * cos(theta_4) - R_1],[R_2 * sin(theta_2) + R_3 * sin(theta_3) - R_4 * sin(theta_4)]]  
    # Get solution for delta_theta_3 and delta_theta_4
    x = np.linalg.solve(a,b)
    delta_theta_3 = x.item(0)
    delta_theta_4 = x.item(1)
    theta_3 += delta_theta_3
    theta_4 += delta_theta_4
    print('\n Delta 3,4 in rad', delta_theta_3, delta_theta_4)
# Turn radius back to degrees
theta_3 = theta_3 * 180 / 3.14
theta_4 = theta_4 * 180 / 3.14
print('\n Final Result', theta_3, theta_4)



