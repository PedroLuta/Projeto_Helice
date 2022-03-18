from importing import *



def compatibility_function(point1, point2, value1, value2, der1, der2):
    M = nup.zeros((4, 4))
    M[0] = [point2**3, point2**2, point2, 1]
    M[1] = [point1**3, point1**2, point1, 1]
    M[2] = [3*(point2**2), 2*point2, 1, 0]
    M[3] = [3*(point1**2), 2*point1, 1, 0]
    
    rhs = nup.array(([value2], [value1], [der2], [der1]))

    coeffs = nup.squeeze(nup.transpose(nup.matmul(nup.linalg.inv(M), rhs)))

    return coeffs

#print(compatibility_function(0.1, 0.3, 0, 30, 0, 0))