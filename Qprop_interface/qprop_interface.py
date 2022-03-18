import os
import numpy as nup
import matplotlib.pyplot as plt

def run_spec(propfile, motorfile, runfile, outfile):
    os.system(f'cmd /c qprop {propfile} {motorfile} {runfile} > {outfile}\n')

def read_output(outfile, skip_lines = 17):
    simulation_vector = []
    read_counter = 0
    with open(outfile, 'r') as o:
        while True:
            if read_counter < skip_lines:
                content = o.readline()
                read_counter += 1
            else:
                content = o.readline()
                if content == '':
                    break
                if content == '\n':
                    continue
                content = content.replace('\n', '').split()
                for i in content:
                    f = float(i)
                    content[content.index(i)] = f
                simulation_vector.append(content)
    return simulation_vector

def build_constant_power_vector(info_vec, max_power = 700, power_index = 15, v_index = 0):
    V_vector = []
    for line in info_vec:
        v = line[v_index]
        if v not in V_vector:
            V_vector.append(v)
        else:
            break
    
    main_vec = []
    for v in V_vector:
        iter_vec = []
        for i in range(len(info_vec)):
            if info_vec[i][v_index] == v:
                iter_vec.append(info_vec[i])
        compare_vector = [abs(line[power_index] - max_power) for line in iter_vec]
        minimum_diff = min(compare_vector)
        minimum_diff_index = compare_vector.index(minimum_diff)
        main_vec.append(iter_vec[minimum_diff_index])
    return main_vec

name = "helice"
max_radius = 302.26
Blades = 4
Cl0 = 0.55
Cla = 7.73
Clmin = 0
Clmax = 1.5
Cd0 = 0.02
Cd2u = 0.04
Cd2l = 0.055
ClCd0 = 0.8
Re_ref = 100000
Re_exp = -0.5
R_fac = 0.001
C_fac = 0.001
B_fac = 1
Radd = 0
Cadd = 0
Badd = 0

r_vec = [45.339, 60.452, 75.565, 90.678, 105.791, 120.904, 136.017, 151.13, 166.243, 181.356, 196.469, 211.582, 226.695, 241.808, 256.921, 272.034, 287.147, 302.26]
chord_vec = [32.566, 43.586, 54.159, 60.574, 62.065, 59.571, 57.048, 54.201, 51.175, 48.157, 45.247, 42.464, 39.738, 36.920, 33.772, 29.974, 25.122, 2.720]
beta_vec = [29.830, 40.206, 40.206, 35.160, 31.121, 27.846, 25.154, 22.910, 21.017, 19.401, 18.009, 16.797, 15.735, 14.796, 13.961, 13.213, 12.540, 11.931]

propfile = "propfile_test.txt"
motorfile = "motorfile.txt"
runfile = "runfile.txt"
outfile = "outputs\out.txt"

with open(propfile, 'w') as inp:
    inp.write(f'{name}\n')
    inp.write(f'{Blades} {max_radius}\n')
    inp.write(f'{Cl0} {Cla}\n')
    inp.write(f'{Clmin} {Clmax}\n')
    inp.write(f'{Cd0} {Cd2u} {Cd2l} {ClCd0}\n')
    inp.write(f"{Re_ref} {Re_exp}\n")
    inp.write(f'{R_fac} {C_fac} {B_fac}\n')
    inp.write(f'{Radd} {Cadd} {Badd}\n')
    for i in range(len(r_vec)):
        inp.write(f'{r_vec[i]} {chord_vec[i]} {beta_vec[i]}\n')
run_spec(propfile, motorfile, runfile, outfile)
simulation_vector = read_output(outfile)
cte_power_vec = build_constant_power_vector(simulation_vector)
cte_power_vec = nup.array(cte_power_vec).T.tolist()
prop_results = {
        "Velocities": cte_power_vec[0],
        "rpm": cte_power_vec[1],
        "DBeta": cte_power_vec[2],
        "Thrust": cte_power_vec[3],
        "Torque": cte_power_vec[4],
        "Power_at_shaft": cte_power_vec[5],
        "Volts": cte_power_vec[6],
        "Amps": cte_power_vec[7],
        "effmot": cte_power_vec[8],
        "effprop": cte_power_vec[9],
        "adv": cte_power_vec[10],
        "CT": cte_power_vec[11],
        "CP": cte_power_vec[12],
        "DV": cte_power_vec[13],
        "eff": cte_power_vec[14],
        "Power": cte_power_vec[15],
        "Pprop": cte_power_vec[16],
        "cl_avg": cte_power_vec[17],
        "cd_avg": cte_power_vec[18],
    }

with open("results.txt", 'w') as inp:
    for key in prop_results:
        inp.write(f'{key}:\n')
        inp.write(f'{prop_results[key]}\n\n')
