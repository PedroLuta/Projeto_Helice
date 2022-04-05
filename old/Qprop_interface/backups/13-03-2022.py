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

def exclude_informations(info_vec, exclude_index_vec):
    while len(exclude_index_vec) > 0:
        biggest_exclude_index = 0
        for i in range(len(exclude_index_vec)):
            if exclude_index_vec[i] > exclude_index_vec[biggest_exclude_index]:
                biggest_exclude_index = i
        for line in info_vec:
            line.pop(exclude_index_vec[biggest_exclude_index])
        exclude_index_vec.pop(biggest_exclude_index)
    return info_vec

def build_constant_power_vector(info_vec, max_power = 700):
    V_vector = []
    for line in info_vec:
        v = line[0]
        if v not in V_vector:
            V_vector.append(v)
    
    main_vec = []
    for v in V_vector:
        iter_vec = []
        for i in range(len(info_vec)):
            if info_vec[i][0] == v:
                iter_vec.append(info_vec[i])
        #print(iter_vec)
        compare_vector = [abs(line[2] - max_power) for line in iter_vec]
        minimum_diff = min(compare_vector)
        minimum_diff_index = compare_vector.index(minimum_diff)
        main_vec.append(iter_vec[minimum_diff_index])
    return main_vec



Blades = 4
Cl0 = 0.55
Cla = 7.73
Clmin = 0
Clmax = 1.5
Cd0 = 0.02
Cd2u = 0.04
Cd2l = 0.055
ClCd0 = 0.8

propfile_test = "propfile_test.txt"
with open(propfile_test, 'w') as inp:
    inp.write('helice\n')
    inp.write(f'{Blades} {302.26}\n')
    inp.write(f'{Cl0} {Cla}\n')
    inp.write(f'{Clmin} {Clmax}\n')
    inp.write(f'{Cd0} {Cd2u} {Cd2l} {ClCd0}\n')
    inp.write("100000 -0.5\n0.001 0.001 1.0\n0. 0. 0.\n45.339 32.56690738 29.830716\n60.452 43.5860406 40.20649549\n75.565 54.15905935 40.20649549\n90.678 60.5742611 35.16022374 \n105.791 62.06594684 31.12177578\n120.904	 59.57144197	27.84693006\n136.017	 57.04897529	25.15422488\n151.13	 54.20132447	22.9103745\n166.243	 51.17559149	21.01722617\n181.356	 48.15736994	19.4018361\n196.469	 45.24793082	18.00934099\n211.582	 42.46422254	16.79793758\n226.695	 39.73887089	15.73535684\n241.808	 36.92017908	14.79637473\n256.921	 33.77212774	13.96103831\n272.034	 29.97437487	13.21338756\n287.147	 25.1222559	12.54052339\n302.26	 2.72034	11.93191978")

propfile = "propfile_simple.txt"
motorfile = "motorfile.txt"
runfile = "runfile.txt"
outfile = "outputs\out.txt"

run_spec(propfile, motorfile, runfile, outfile)
simulation_vector = read_output(outfile)
important_info_vector = exclude_informations(simulation_vector, [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18]) #Apenas V, Tração e Potência

cte_power_vec = build_constant_power_vector(important_info_vector) 
cte_power_vec = exclude_informations(cte_power_vec, [2]) #Apenas V e Tração
cte_power_vec = nup.array(cte_power_vec).T.tolist()

# Power_limit = 700
# v_vec = []
# rpm_vec = []
# t_vec = []
# q_vec = []
# em_vec = []
# ep_vec = []
# p_vec = []

# ii = 0
# for i in simulation_vector:
#     V, rpm, Thrust, Torque, effmot, effprop, Power_measured = i[0], i[1], i[3], i[4], i[8], i[9], i[15]
#     v_vec.append(V)
#     rpm_vec.append(rpm)
#     t_vec.append(Thrust)
#     q_vec.append(Torque)
#     em_vec.append(effmot)
#     ep_vec.append(effprop)
#     p_vec.append(Power_measured)
#     #print(Power_measured)

# first = True
# placeholder_vec = []
# ideal_combinations = []
# last_V = 999
# number_of_velocities = 36
# vmin = 0
# vmax = 35
# first_swap = True

# i = 0
# while i < len(p_vec): #for Power in p_vec:
#     #index = p_vec.index(Power)
#     current_V = v_vec[i]
    
#     if first:
#         placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])
#         previous_V = current_V
#         first = False
#         continue

#     if (current_V < previous_V) and first_swap:
#         ideal_power = 0
#         last_V = previous_V
#         ii = 0
#         while ii < len(placeholder_vec): 
#             current_power = placeholder_vec[ii][4]
#             if (current_power < Power_limit) and (abs(current_power - Power_limit) < abs(ideal_power - Power_limit)) and (placeholder_vec[ii][2] > 0):
#                 ideal_power = current_power
#                 index2 = placeholder_vec.index(placeholder_vec[ii])
#             ii += 1
#         if (placeholder_vec[index2][4] < Power_limit):
#             #print(placeholder_vec[index2])
#             ideal_combinations.append(placeholder_vec[index2])
#         placeholder_vec.clear()
#         placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])
#         first_swap = False

#     elif current_V == last_V:
#         ideal_power = 0
#         ii = 0
#         while ii < len(placeholder_vec):
#             current_power = placeholder_vec[ii][4]
#             if (abs(current_power) < Power_limit) and (abs(current_power - Power_limit) < abs(ideal_power - Power_limit)) and (placeholder_vec[ii][2] > 0):
#                 ideal_power = current_power
#                 index2 = placeholder_vec.index(placeholder_vec[ii])
#             ii += 1
#         if (placeholder_vec[index2][4] < Power_limit):
#             #print(placeholder_vec[index2])
#             ideal_combinations.append(placeholder_vec[index2])
#         placeholder_vec.clear()
#         placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])


#     placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])
#     previous_V = current_V
    
#     i += 1

# v_dist = nup.linspace(vmin, vmax, number_of_velocities)
# #print(v_dist)
# #print(v_dist)
# iii = 0
# distribution = nup.zeros((number_of_velocities - 6, 5))
# #print(distribution)
# #first = True
# #print(ideal_combinations)
# previous_V = vmin - 1
# velocity = vmin
# index_v = 0
# while iii < len(ideal_combinations):
#     if (ideal_combinations[iii][0] != v_dist[index_v]) and (ideal_combinations[iii][0] != previous_V):
#         #first = False
#         iiii = 0
#         while iiii < len(ideal_combinations[iii]): #for i in ideal_combinations[iii]: 
#             distribution[index_v][iiii] = ideal_combinations[iii-1][iiii]
#             iiii += 1
#         index_v += 1
#         previous_V = ideal_combinations[iii][0]
#     iii += 1
# #print(distribution)


# v_vector = []
# rpm_vector = []
# t_vector = []
# q_vector = []
# iiiii = 0
# while iiiii < len(distribution):
#     v_vector.append(distribution[iiiii][0])
#     rpm_vector.append(distribution[iiiii][1])
#     t_vector.append(distribution[iiiii][2])
#     q_vector.append(distribution[iiiii][3])
#     iiiii += 1

# plt.plot(v_vector, t_vector, label = 'Thrust')
# plt.legend()
# plt.show()
# plt.plot(v_vector, q_vector, label = 'Torque')
# plt.legend()
# plt.show()
# plt.plot(v_vector, rpm_vector, label = 'rpm')
# plt.legend()
# plt.show()

# #print(v_vector)
# #print(rpm_vector)
# #print(t_vector)
# #print(q_vector)