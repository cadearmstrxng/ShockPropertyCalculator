from ObliqueShockCalculator import CalculateShockProperties

def main():
    gamma=float(input("Enter Specific Heat Ratio: "))
    incom_mach_num = float(input("Enter Incoming Mach Number: "))
    theta=float(input("Enter Turning Angle: "))
    flag=input("Do you want weak solution? (yes/no)")
    if flag=="no":
        flag = True
    else:
        flag = False

    
    beta,outgoing_mach_num,ro2_to_ro1,p2_to_p1,t2_to_t1,po2_to_po1=CalculateShockProperties(gamma, incom_mach_num,theta,flag)

    #Print Statements
    print(f'Beta (Wave Angle) = {beta:.4f}')
    print(f'Outgoing Mach Number = {outgoing_mach_num:.4f}')
    print(f'ro2/ro1 (Density Ratio) = {ro2_to_ro1:.4f}')
    print(f'P2/P1 (Pressure Ratio) = {p2_to_p1:.4f}')
    print(f'T2/T1 (Temperature Ratio) = {t2_to_t1:.4f}')
    print(f'Po2/Po1 (Stagnation Pressure Ratio) = {po2_to_po1:.4f}')


if __name__ == '__main__':
    main()
