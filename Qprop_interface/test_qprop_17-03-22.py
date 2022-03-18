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

# def exclude_informations(info_vec, exclude_index_vec):
#     while len(exclude_index_vec) > 0:
#         biggest_exclude_index = 0
#         for i in range(len(exclude_index_vec)):
#             if exclude_index_vec[i] > exclude_index_vec[biggest_exclude_index]:
#                 biggest_exclude_index = i
#         for line in info_vec:
#             line.pop(exclude_index_vec[biggest_exclude_index])
#         exclude_index_vec.pop(biggest_exclude_index)
#     return info_vec

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


Cl0min = 0
Cl0max = 0.8
Clamin = 4
Clamax = 10
ClCd0min = 0
ClCd0max = 0.8

valores = 3

Cl0_vec = nup.linspace(Cl0min, Cl0max, valores)
Cla_vec = nup.linspace(Clamin, Clamax, valores)
ClCd0_vec = nup.linspace(ClCd0min, ClCd0max, valores)

# Blades = 4
Cl0 = 0.55
Cla = 7.73
# Clmin = 0
# Clmax = 1.5
Cd0 = 0.02
Cd2u = 0.04
Cd2l = 0.055
ClCd0 = 0.8

propfile = "propfile_test.txt"
motorfile = "motorfile.txt"
runfile = "runfile.txt"
outfile = "outputs\out.txt"

biggest_area = 0
for Cl0 in Cl0_vec:
    #for Cla in Cla_vec:
    for ClCd0 in ClCd0_vec:    
        if Cl0 > ClCd0:
            continue
        with open(propfile, 'w') as inp:
            inp.write('helice\n')
            inp.write(f'{4} {302.26}\n')
            inp.write(f'{Cl0} {Cla}\n')
            inp.write(f'{0} {1.5}\n')
            inp.write(f'{Cd0} {Cd2u} {Cd2l} {ClCd0}\n')
            inp.write("100000 -0.5\n0.001 0.001 1.0\n0. 0. 0.\n45.339 32.56690738 29.830716\n60.452 43.5860406 40.20649549\n75.565 54.15905935 40.20649549\n90.678 60.5742611 35.16022374 \n105.791 62.06594684 31.12177578\n120.904	 59.57144197	27.84693006\n136.017	 57.04897529	25.15422488\n151.13	 54.20132447	22.9103745\n166.243	 51.17559149	21.01722617\n181.356	 48.15736994	19.4018361\n196.469	 45.24793082	18.00934099\n211.582	 42.46422254	16.79793758\n226.695	 39.73887089	15.73535684\n241.808	 36.92017908	14.79637473\n256.921	 33.77212774	13.96103831\n272.034	 29.97437487	13.21338756\n287.147	 25.1222559	12.54052339\n302.26	 2.72034	11.93191978")
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
        plt.plot(prop_results["Velocities"], prop_results["Thrust"], label = f"Cl0 = {round(Cl0, 2)}, ClCd0 = {round(ClCd0, 2)}")
        # area = nup.trapz(prop_results["Thrust"], prop_results["Velocities"])
        # if area > biggest_area:
        #     biggest_area = area
        #     best_combination = [Cl0, Cla, ClCd0]
        #     print(best_combination)
plt.legend()
plt.show()

Cl0, Cla, ClCd0 = best_combination[0], best_combination[1], best_combination[2]

with open(propfile, 'w') as inp:
    inp.write('helice\n')
    inp.write(f'{4} {302.26}\n')
    inp.write(f'{Cl0} {Cla}\n')
    inp.write(f'{0} {1.5}\n')
    inp.write(f'{Cd0} {Cd2u} {Cd2l} {ClCd0}\n')
    inp.write("100000 -0.5\n0.001 0.001 1.0\n0. 0. 0.\n45.339 32.56690738 29.830716\n60.452 43.5860406 40.20649549\n75.565 54.15905935 40.20649549\n90.678 60.5742611 35.16022374 \n105.791 62.06594684 31.12177578\n120.904	 59.57144197	27.84693006\n136.017	 57.04897529	25.15422488\n151.13	 54.20132447	22.9103745\n166.243	 51.17559149	21.01722617\n181.356	 48.15736994	19.4018361\n196.469	 45.24793082	18.00934099\n211.582	 42.46422254	16.79793758\n226.695	 39.73887089	15.73535684\n241.808	 36.92017908	14.79637473\n256.921	 33.77212774	13.96103831\n272.034	 29.97437487	13.21338756\n287.147	 25.1222559	12.54052339\n302.26	 2.72034	11.93191978")
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
