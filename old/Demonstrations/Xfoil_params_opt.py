from xfoil_interface import get_curve_fileMethod
import matplotlib.pyplot as plt
from PARSEC_functions import *

alpha, cl, cd = get_curve_fileMethod(50000, afile = "airfoils\\clarky.txt")#, N_panels = 342)
alpha2, cl2, cd2 = get_curve_fileMethod(50000, afile = "airfoils\\clarky.txt", N_panels = 90, P_bunch = 1, Te_Le_ratio = 0.15, Refined_Le_ratio = 0.2, \
            top_refined_init = 1, top_refined_end = 1, \
                bot_refined_init = 1, bot_refined_end = 1)#, N_panels = 342)
alpha3, cl3, cd3 = get_curve_fileMethod(50000, afile = "airfoils\\clarky.txt", N_panels = 342, P_bunch = 1, Te_Le_ratio = 0.15, Refined_Le_ratio = 0.2, \
            top_refined_init = 1, top_refined_end = 1, \
                bot_refined_init = 1, bot_refined_end = 1)#, N_panels = 342)
print(alpha)
print(cl)
print(cd)
plt.plot(alpha, cl, label = "Optimized")
plt.plot(alpha2, cl2, label = "same panels, not optimized")
plt.plot(alpha3, cl3, label = "correct")
plt.legend()
plt.show()
plt.plot(alpha, cd, label = "Optimized")
plt.plot(alpha2, cd2, label = "same panels, not optimized")
plt.plot(alpha3, cd3, label = "correct")
plt.legend()
plt.show()
plt.plot(cl, cd, label = "Optimized")
plt.plot(cl2, cd2, label = "same panels, not optimized")
plt.plot(cl3, cd3, label = "correct")
plt.legend()
plt.show()