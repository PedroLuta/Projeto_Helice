import numpy as nup
import math

#def inverse(a,c,n): 
#
#    L = nup.zeros( (n, n) )
#    U = nup.zeros( (n, n) )
#    b = nup.zeros( (n) )
#    d = nup.zeros( (n) )
#    x = nup.zeros( (n) )
#
#    for k in range(0,n-1):                         
#        for i in range(k+1,n):  
#            coeff = a[i][k]/a[k][k]
#            L[i][k] = coeff
#            for j in range(k+1,n):                          
#                a[i][j] = a[i][j] - (coeff*a[k][j])
#
#    for i in range(0,n):                                    
#        L[i][i] = 1.0
#    
#    for j in range(0,n):                                    
#        for i in range(0,n):                                
#            U[i][j] = a[i][j]
#    
#    for k in range(0,n):                                    
#        b[k] = 1.0
#        d[0] = b[0]
#        for i in range(1,n):                                
#            d[i] = b[i]
#            for j in range(0,i):                            
#                d[i] = d[i] - (L[i][j]*d[j])
#        x[n - 1] = d[n - 1]/U[n - 1][n - 1]
#        for i in range(n-2,-1,-1):                          
#            x[i] = d[i]
#            for j in range (n-1,i,-1):                      
#                x[i] = x[i] - (U[i][j]*x[j])
#            x[i] = x[i]/U[i][i]
#        for i in range(0,n):                                
#            c[i][k] = x[i]
#        b[k] = 0.0
#
#def airfoilgen(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo, coeffs = False):
#    try:
#        coeff_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\codigos_0605\Coeff.txt' 
#        coord_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\codigos_0605\airfoil.txt' 
#
#        a_up = nup.zeros( (6, 12) )
#        a_lo = nup.zeros( (6, 12) )
#        a_inv_up = nup.zeros( (6, 6) )
#        a_inv_lo = nup.zeros( (6, 6) )
#        cup = nup.zeros( (6) )
#        clo = nup.zeros( (6) )
#        coe_up = nup.zeros( (6) )
#        coe_lo = nup.zeros( (6) )
#
#        pi = nup.pi
#        np = 50
#
#        x = nup.zeros((np))
#        z_up = nup.zeros((np))
#        z_lo = nup.zeros((np))
#
#        for j in range(0,6):    
#            a_up[3][j]=1.0
#            a_up[0][j]=xup**(j + 0.5) 
#            a_up[4][j]=j + 0.5 
#            a_up[1][j]=(j + 0.5)*xup**(j - 0.5) 
#            a_up[2][j]=(j - 0.5)*(j + 0.5)*xup**(j - 1.5) 
#            a_up[5][j]=0.0
#
#            a_lo[3][j]=1.0
#            a_lo[0][j]=xlo**(j + 0.5) 
#            a_lo[4][j]=j + 0.5 
#            a_lo[1][j]=(j + 0.5)*xlo**(j - 0.5) 
#            a_lo[2][j]=(j - 0.5)*(j + 0.5)*xlo**(j - 1.5) 
#            a_lo[5][j]=0.0
#
#        a_up[5][0]=1.0
#        a_lo[5][0]=1.0
#
#        inverse(a_up,a_inv_up,6)
#        inverse(a_lo,a_inv_lo,6)
#
#        cup[3]=zte+0.5*dzte 
#        cup[0]=zup 
#        cup[4]=math.tan((ate-0.5*bte)*pi/180.0) 
#        cup[1]=0.0 
#        cup[2]=zxxup 
#        cup[5]=(2.0*rup)**0.5 
#
#        clo[3]=zte-0.5*dzte 
#        clo[0]=zlo 
#        clo[4]=math.tan((ate+0.5*bte)*pi/180.0) 
#        clo[1]=0.0 
#        clo[2]=zxxlo 
#        clo[5]=-(2.0*rlo)**0.5 
#
#        if coeffs:
#            with open(coeff_file, 'w') as coefficient:
#                coe_temp_up = 0.0
#                coe_temp_lo = 0.0
#                for i in range(0,6):                                    
#                    for j in range(0,6):                                
#                        coe_temp_up = coe_temp_up + a_inv_up[i][j]*cup[j]
#                        coe_temp_lo = coe_temp_lo + a_inv_lo[i][j]*clo[j]
#                    coe_up[i] = coe_temp_up
#                    coe_lo[i] = coe_temp_lo
#                    coefficient.write(str(coe_up[i]) + '\n')                   
#                    coefficient.write(str(coe_lo[i]) + '\n\n')                   
#                    coe_temp_up = 0.0
#                    coe_temp_lo = 0.0
#        else:
#            coe_temp_up = 0.0
#            coe_temp_lo = 0.0
#            for i in range(0,6):                                    
#                for j in range(0,6):                                
#                    coe_temp_up = coe_temp_up + a_inv_up[i][j]*cup[j]
#                    coe_temp_lo = coe_temp_lo + a_inv_lo[i][j]*clo[j]
#                coe_up[i] = coe_temp_up
#                coe_lo[i] = coe_temp_lo  
#                print(coe_up[i]) 
#                print(coe_lo[i])               
#                coe_temp_up = 0.0
#                coe_temp_lo = 0.0
#
#
#        with open(coord_file, 'w') as coord:
#            z_temp_up = 0.0
#            z_temp_lo = 0.0
#            for i in range(0,np):                                   
#                theta = (180.0/(np - 1))*(np - i - 1)
#                x[i] = 0.5 - (0.5*math.cos(theta*pi/180.0))
#                for j in range(0,6):                                
#                    z_temp_up = z_temp_up + coe_up[j]*x[i]**(j + 0.5)
#                    z_temp_lo = z_temp_lo + coe_lo[j]*x[i]**(j + 0.5)
#                z_up[i] = z_temp_up
#                z_lo[i] = z_temp_lo
#                z_temp_up = 0.0
#                z_temp_lo = 0.0
#                coord.write(str(x[i]) + ' ' + str(z_up[i]) + '\n')              
#
#            for i in range(np-2,-1,-1):
#                coord.write(str(x[i]) + ' ' + str(z_lo[i]) + '\n')
#    except:
#        return False
#    return True

def get_curve_coeffs(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo):
    Cup = nup.zeros((6, 6))
    Cup[0] = [1, 1, 1, 1, 1, 1]
    Cup[1] = [xup**0.5, xup**1.5, xup**2.5, xup**3.5, xup**4.5, xup**5.5]
    Cup[2] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    Cup[3] = [0.5*(xup**(-0.5)), 1.5*(xup**0.5), 2.5*(xup**1.5), 3.5*(xup**2.5), 4.5*(xup**3.5), 5.5*(xup**4.5)]
    Cup[4] = [-0.25*(xup**(-1.5)), 0.75*(xup**(-0.5)), 3.75*(xup**0.5), 8.75*(xup**1.5), 15.75*(xup**2.5), 24.75*(xup**3.5)]
    Cup[5] = [1, 0, 0, 0, 0, 0]

    bup = nup.array(([zte + (dzte/2)], [zup], [nup.tan(nup.radians(ate - (bte/2)))], [0], [zxxup], [(2*rup)**0.5]))

    Clo = nup.zeros((6, 6))
    Clo[0] = [1, 1, 1, 1, 1, 1]
    Clo[1] = [xlo**0.5, xlo**1.5, xlo**2.5, xlo**3.5, xlo**4.5, xlo**5.5]
    Clo[2] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    Clo[3] = [0.5*(xlo**(-0.5)), 1.5*(xlo**0.5), 2.5*(xlo**1.5), 3.5*(xlo**2.5), 4.5*(xlo**3.5), 5.5*(xlo**4.5)]
    Clo[4] = [-0.25*(xlo**(-1.5)), 0.75*(xlo**(-0.5)), 3.75*(xlo**0.5), 8.75*(xlo**1.5), 15.75*(xlo**2.5), 24.75*(xlo**3.5)]
    Clo[5] = [1, 0, 0, 0, 0, 0]

    blo = nup.array(([zte - (dzte/2)], [zlo], [nup.tan(nup.radians(ate + (bte/2)))], [0], [zxxlo], [-(2*rlo)**0.5]))

    aup = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(Cup), bup)))
    alo = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(Clo), blo)))

    return aup, alo

def airfoilgen1(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo):
    coord_file = r'C:\Users\PEDRO\Desktop\IC DE HÉLICE\codigos_0605\airfoil1.txt'
    np = 50

    aup, alo = get_curve_coeffs(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo)

    with open(coord_file, 'w') as inp:
        for x in nup.linspace(1, 0, np + 1):
            inp.write(f'{x} {aup[0]*(x**0.5) + aup[1]*(x**1.5) + aup[2]*(x**2.5) + aup[3]*(x**3.5) + aup[4]*(x**4.5) + aup[5]*(x**5.5)}\n')
        for x in nup.linspace(0, 1, np + 1):
            inp.write(f'{x} {alo[0]*(x**0.5) + alo[1]*(x**1.5) + alo[2]*(x**2.5) + alo[3]*(x**3.5) + alo[4]*(x**4.5) + alo[5]*(x**5.5)}\n')

def airfoilgen_params(PARSEC_params, coeffs = False):
    ate = PARSEC_params[0]
    bte = PARSEC_params[1]
    zte = PARSEC_params[2]
    dzte = PARSEC_params[3]
    rup = PARSEC_params[4]
    xup = PARSEC_params[5]
    zup = PARSEC_params[6]
    zxxup = PARSEC_params[7]
    rlo = PARSEC_params[8]
    xlo = PARSEC_params[9]
    zlo = PARSEC_params[10]
    zxxlo = PARSEC_params[11]
    return airfoilgen1(ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo, coeffs)

#run = airfoilgen(-6.4434644618472925 , 9.4527710459963057, 0.001, 0.003, 8.10992437339230765E-003 , 0.43083556083047875 , 6.29208143868678071E-002, -0.75, 8.49891911271616340E-003,  0.34404601144765617, -5.89010630429923354E-002, 0.75, coeffs = False)
run = airfoilgen1(-6.4434644618472925 , 9.4527710459963057, 0.001, 0.003, 8.10992437339230765E-003 , 0.43083556083047875 , 6.29208143868678071E-002, -0.75, 8.49891911271616340E-003,  0.34404601144765617, -5.89010630429923354E-002, 0.75)

#-6.4434644618472925 , 9.4527710459963057, 2.33147293329238892E-004, -1.89922749996185303E-004, 8.10992437339230765E-003 , 0.43083556083047875 , 6.29208143868678071E-002, -0.42542379733615032, 8.49891911271616340E-003,  0.34404601144765617, -5.89010630429923354E-002, 0.70572080361993572
#ate, bte, zte, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo
