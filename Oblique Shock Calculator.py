#NOTE FOR INSTRUCTORS:
#Please ensure that you have the scipy module installed on your computer using pip
#My code will not work otherwise as it uses fsolve

from math import *
from scipy.optimize import fsolve


#Inputs:
#1. Gamma, angles in degrees
#2. Incoming Mach Number
#3. Turning Angle in degrees
#4. Flag designating weak or strong solutions
#     False = weak
#     True = strong

#Outputs:
#1. Wave Angle in degrees
#2. Normal Mach number of incoming Mach number
#3. Normal Mach number of outgoing Mach number
#4. Outgoing Mach Number
#5. ratio of outgoing density to incoming density, ro2/ro1
#6. ratio of outgoing pressure to incoming pressure, P2/P1
#7. ratio of outgoing temperature to incoming temperature, T2/T1
#8. ratio of outgoing stagnation pressure to incoming stagnation pressure, Po2/Po1
def CalculateShockProperties(gamma, incoming_mach_num, turning_angle, flag=False):
    #Treats as normal shock if turning angle is 0.
    if turning_angle==0:
        #if turning angle is 0, beta = mach angle = arcsin(1/M1)
        beta=asin(1/incoming_mach_num)*180/pi

        #Calculates outgoing mach number
        outgoing_mach_num=sqrt( (1+((gamma-1)/2)*incoming_mach_num**2)/(gamma*incoming_mach_num**2-(gamma-1)/2) )

        #Calculates ratio of outgoing density to incoming density
        ro2_to_ro1= ( (gamma+1) * incoming_mach_num**2)/(2+ (gamma-1)*incoming_mach_num**2)

        #Calculates ratio of outgoing pressure to incoming pressure
        p2_to_p1=1 + (2*gamma*(incoming_mach_num**2 -1))/(gamma+1)

        #Calculates ratio of outgoing temperature to incoming temperature
        t2_to_t1=(1+((gamma-1)/2)*incoming_mach_num**2)/(1+((gamma-1)/2)*outgoing_mach_num**2)

        #Calculates ratio of outgoing stagnation pressure to incoming stagnation pressure
        po2_to_po1=pow(1/t2_to_t1,(gamma/(gamma-1)))*p2_to_p1

        return beta,outgoing_mach_num,ro2_to_ro1,p2_to_p1,t2_to_t1,po2_to_po1

    else:
        #Calculates beta using the fsolve function. Weak solution when initial guess is 5, strong when initial guess is 100
        if not flag:
            theta_beta_mach=lambda beta:tan(turning_angle*pi/180) - 2 * ( 1 / tan(beta[0]*pi/180) ) * ( (( incoming_mach_num)**2 * ( sin(beta[0]*pi/180) )**2 -1) / ( (incoming_mach_num)**2 * (gamma + cos(2*beta[0]*pi/180)) + 2))
            beta=fsolve(theta_beta_mach, [5])[0]
        else:
            theta_beta_mach=lambda beta:tan(turning_angle*pi/180) - 2 * ( 1 / tan(beta[0]*pi/180) ) * ( (( incoming_mach_num)**2 * ( sin(beta[0]*pi/180) )**2 -1) / ( (incoming_mach_num)**2 * (gamma + cos(2*beta[0]*pi/180)) + 2))
            beta=fsolve(theta_beta_mach, [100])[0]

    #Calculates Mn1, Mn2, outgoing Mach number
    Mn1=incoming_mach_num*sin(beta*pi/180)
    
    Mn2=sqrt( (1+((gamma-1)/2)*Mn1**2)/(gamma*Mn1**2-(gamma-1)/2))
    
    outgoing_mach_num=Mn2/sin((beta-turning_angle)*pi/180)

    #Calculates ratio of outgoing density to incoming density
    ro2_to_ro1= ( (gamma+1) * Mn1**2)/(2+ (gamma-1)*Mn1**2)

    #Calculates ratio of outgoing pressure to incoming pressure
    p2_to_p1=1 + (2*gamma*(Mn1**2 -1))/(gamma+1)

    #Calculates ratio of outgoing temperature to incoming temperature
    t2_to_t1=(1+((gamma-1)/2)*incoming_mach_num**2)/(1+((gamma-1)/2)*outgoing_mach_num**2)

    #Calculates ratio of outgoing stagnation pressure to incoming stagnation pressure
    po2_to_po1=pow(1/t2_to_t1,(gamma/(gamma-1)))*p2_to_p1
    
    
    return beta,outgoing_mach_num,ro2_to_ro1,p2_to_p1,t2_to_t1,po2_to_po1


#This function solves an entire column of the HW3 Problem 2 Table:
#INPUTS:
#1. Total Turning Angle in degrees
#2. Turning Angle 1 in degrees
#3. Ma - Incoming Mach Number
#4. gamma

def SolveHW3Num2(total_turning_angle, turning_angle_1, Ma, gamma=1.4):
    theta_2=total_turning_angle-turning_angle_1

    beta,M1,ro2_to_ro1,p1_to_pa,t2_to_t1,po1_to_poa=CalculateShockProperties(gamma, Ma,turning_angle_1)

    print(f'M1 = {M1:.4f}')
    print(f'Po1/Poa (Stagnation Pressure Ratio) = {po1_to_poa:.4f}')
    print(f'P1/Pa (Pressure Ratio) = {p1_to_pa:.4f}')
    
    beta,M2,ro2_to_ro1,p2_to_p1,t2_to_t1,po2_to_po1=CalculateShockProperties(gamma, M1,theta_2)

    print(f'M2 = {M2:.4f}')
    print(f'Po2/Po1 (Stagnation Pressure Ratio) = {po2_to_po1:.4f}')
    print(f'P2/P1 (Pressure Ratio) = {p2_to_p1:.4f}')

    beta,M3,ro2_to_ro1,p3_to_p2,t2_to_t1,po3_to_po2=CalculateShockProperties(gamma, M2,0)

    print(f'M3 = {M3:.4f}')
    print(f'Po3/Po2 (Stagnation Pressure Ratio) = {po3_to_po2:.4f}')
    print(f'P3/P2 (Pressure Ratio) = {p3_to_p2:.4f}')

    po3_to_poa=po3_to_po2*po2_to_po1*po1_to_poa
    p3_to_pa=p3_to_p2*p2_to_p1*p1_to_pa

    print(f'Po3/Poa = {po3_to_poa:.4f}')
    print(f'P3/Pa = {p3_to_pa:.4f}')
    
    

def main():
    #This flag determines if you want to solve HW3 problem 2.
    #   True = Output CalculateShockProperies values AND Output HW3 Problem 2 answers
    #   False = Only output CalculateShockProperties values
    flag=False

    #TO TEST: 
    #1. Set Incoming Mach Number, gamma, and turning angle
    #2. Run the code
    #It should output the correct values in the form seen below
    incoming_mach_num=1.4009
    gamma=1.4
    turning_angle=0
    beta,outgoing_mach_num,ro2_to_ro1,p2_to_p1,t2_to_t1,po2_to_po1=CalculateShockProperties(gamma, incoming_mach_num,turning_angle)

    #Print Statements
    print(f'Beta (Wave Angle) = {beta:.4f}')
    print(f'Outgoing Mach Number = {outgoing_mach_num:.4f}')
    print(f'ro2/ro1 (Density Ratio) = {ro2_to_ro1:.4f}')
    print(f'P2/P1 (Pressure Ratio) = {p2_to_p1:.4f}')
    print(f'T2/T1 (Temperature Ratio) = {t2_to_t1:.4f}')
    print(f'Po2/Po1 (Stagnation Pressure Ratio) = {po2_to_po1:.4f}')
    
    #############################

    #start Solving problem 2
    if flag:
        total_turning_angles=[30,30,30,25,20]
        turning_angle_1s=[15,14,10,12,9]
        incoming_mach_num=2.6

        for i in range(len(total_turning_angles)):
            print()
            print('####################################')
            print()
            SolveHW3Num2(total_turning_angles[i],turning_angle_1s[i],incoming_mach_num)

if __name__ == '__main__':
    main()