import os
import numpy as nup
import matplotlib.pyplot as plt

def run_default(outfile):
    #outfile = 'out.txt'
    os.system(f'cmd /c qprop propfile.txt motorfile.txt runfile.txt > {outfile}\n')

def run_spec(propfile, motorfile, runfile, outfile):
    os.system(f'cmd /c qprop {propfile} {motorfile} {runfile} > {outfile}\n')

def read_output(outfile):
    simulation_vector = []
    read_counter = 0
    with open(outfile, 'r') as o:
        while True:
            if read_counter < 17:
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
    #os.remove(outfile)


#run_default(r"outputs\out.txt")
run_spec("propfile.txt", "motorfile.txt", "runfile.txt", r"outputs\out.txt")
simulation_vector = read_output(r"outputs\out.txt")
#print(simulation_vector)

Power_limit = 700
v_vec = []
rpm_vec = []
t_vec = []
q_vec = []
em_vec = []
ep_vec = []
p_vec = []

ii = 0
for i in simulation_vector:
    V, rpm, Thrust, Torque, effmot, effprop, Power_measured = i[0], i[1], i[3], i[4], i[8], i[9], i[15]
    v_vec.append(V)
    rpm_vec.append(rpm)
    t_vec.append(Thrust)
    q_vec.append(Torque)
    em_vec.append(effmot)
    ep_vec.append(effprop)
    p_vec.append(Power_measured)
    #print(Power_measured)

first = True
placeholder_vec = []
ideal_combinations = []
last_V = 999
number_of_velocities = 36
vmin = 0
vmax = 35
first_swap = True

i = 0
while i < len(p_vec): #for Power in p_vec:
    #index = p_vec.index(Power)
    current_V = v_vec[i]
    
    if first:
        placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])
        previous_V = current_V
        first = False
        continue

    if (current_V < previous_V) and first_swap:
        ideal_power = 0
        last_V = previous_V
        ii = 0
        while ii < len(placeholder_vec): 
            current_power = placeholder_vec[ii][4]
            if (current_power < Power_limit) and (abs(current_power - Power_limit) < abs(ideal_power - Power_limit)) and (placeholder_vec[ii][2] > 0):
                ideal_power = current_power
                index2 = placeholder_vec.index(placeholder_vec[ii])
            ii += 1
        if (placeholder_vec[index2][4] < Power_limit):
            #print(placeholder_vec[index2])
            ideal_combinations.append(placeholder_vec[index2])
        placeholder_vec.clear()
        placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])
        first_swap = False

    elif current_V == last_V:
        ideal_power = 0
        ii = 0
        while ii < len(placeholder_vec):
            current_power = placeholder_vec[ii][4]
            if (abs(current_power) < Power_limit) and (abs(current_power - Power_limit) < abs(ideal_power - Power_limit)) and (placeholder_vec[ii][2] > 0):
                ideal_power = current_power
                index2 = placeholder_vec.index(placeholder_vec[ii])
            ii += 1
        if (placeholder_vec[index2][4] < Power_limit):
            #print(placeholder_vec[index2])
            ideal_combinations.append(placeholder_vec[index2])
        placeholder_vec.clear()
        placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])


    placeholder_vec.append([v_vec[i], rpm_vec[i], t_vec[i], q_vec[i], p_vec[i]])
    previous_V = current_V
    
    i += 1

v_dist = nup.linspace(vmin, vmax, number_of_velocities)
#print(v_dist)
#print(v_dist)
iii = 0
distribution = nup.zeros((number_of_velocities - 6, 5))
#print(distribution)
#first = True
#print(ideal_combinations)
previous_V = vmin - 1
velocity = vmin
index_v = 0
while iii < len(ideal_combinations):
    if (ideal_combinations[iii][0] != v_dist[index_v]) and (ideal_combinations[iii][0] != previous_V):
        #first = False
        iiii = 0
        while iiii < len(ideal_combinations[iii]): #for i in ideal_combinations[iii]: 
            distribution[index_v][iiii] = ideal_combinations[iii-1][iiii]
            iiii += 1
        index_v += 1
        previous_V = ideal_combinations[iii][0]
    iii += 1
#print(distribution)


v_vector = []
rpm_vector = []
t_vector = []
q_vector = []
iiiii = 0
while iiiii < len(distribution):
    v_vector.append(distribution[iiiii][0])
    rpm_vector.append(distribution[iiiii][1])
    t_vector.append(distribution[iiiii][2])
    q_vector.append(distribution[iiiii][3])
    iiiii += 1

plt.plot(v_vector, t_vector, label = 'Thrust')
plt.legend()
plt.show()
plt.plot(v_vector, q_vector, label = 'Torque')
plt.legend()
plt.show()
plt.plot(v_vector, rpm_vector, label = 'rpm')
plt.legend()
plt.show()

#print(v_vector)
#print(rpm_vector)
#print(t_vector)
#print(q_vector)