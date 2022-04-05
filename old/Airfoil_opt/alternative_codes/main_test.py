from PARSEC_functions import *
from xfoil_interface import *

Cl = 0
Cd = 0.2
coord_file = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\airfoils\airfoil.txt'
i = 0
while (Cl < 0.8) or (Cd > 0.05):
    while True:
        params = rand_params_gen()
        params[2] = 0
        run_check = write_from_params_vec(params, coord_file = coord_file, np = 50) 
        aup, alo = coeffs_from_params_vec(params)
        if check_valid(aup, alo, params):
            break
    Cl, Cd, _ = get_properties_fixed_alpha(200000, 0)
    i += 1
    if i%100 == 0:
        print(f"Currently on airfoil {i}")
#print(Cl)
print(f'Cd = {Cd}')
#print(f"It took {i} tries!")
print(f'Cl = {Cl}')
plot_file(coord_file)