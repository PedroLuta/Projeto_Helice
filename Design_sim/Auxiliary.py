from importing import *

def area_under_curve(x, curve):
    area = 0
    for i in range(1, len(x)):
        area += (x[i] - x[i - 1])*(curve[i] + curve[i - 1])/2
    return area

def linear_interpolate(x0, x1, y0, y1, x):
    if x1 - x0 == 0:
        return (y0 + y1)/2
    return y0 + ((x - x0)*(y1 - y0)/(x1 - x0))

def truncate_float(float, decimal_places):
    float *= 10**decimal_places
    float = int(float)
    float /= 10**decimal_places
    return float

def calculate_most_eff_alpha(a_list, cl_list, cd_list):
    a, cl, cd = 0, 0, 1
    clcd_remember = 0
    for i in range(len(a_list)):
        if cd_list[i] == 0:
            continue
        clcd_try = cl_list[i]/cd_list[i]
        if clcd_try > clcd_remember:
            a, cl, cd = a_list[i], cl_list[i], cd_list[i]
            clcd_remember = clcd_try
    return a, cl, cd