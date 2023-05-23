import numpy as np
import matplotlib.pyplot as plt

Kv = 450    #rpm/Volt
i0 = 1.4    #Amp
Rm = 0.5      #Ohms
Vb = 22.2   #Volt

omega_vec = np.linspace(500, 10000, 101)


for Rm in [0.05, 0.1, 0.25, 0.5, 1]:
    i_vec = [(Vb - (omega/Kv))/Rm for omega in omega_vec]
    Q_vec = [(i - i0)/Kv for i in i_vec]
    eff_mot = []
    omega_plot = []
    for j in range(len(omega_vec)):
        if (i_vec[j] < 0) or (i_vec[j] > 45):
            pass
        else:
            eff = Q_vec[j]*omega_vec[j]/(i_vec[j]*Vb)
            if eff > 0:
                eff_mot.append(eff)
                omega_plot.append(omega_vec[j])
    plt.plot(omega_plot, eff_mot, label = f'Rm = {Rm}')
    plt.title("Efficiency by variation of Rm")

# for Kv in [100, 150, 200, 250, 300]:
#     i_vec = [(Vb - (omega/Kv))/Rm for omega in omega_vec]
#     Q_vec = [(i - i0)/Kv for i in i_vec]
#     eff_mot = []
#     omega_plot = []
#     for j in range(len(omega_vec)):
#         if (i_vec[j] < 0) or (i_vec[j] > 45):
#             pass
#         else:
#             eff = Q_vec[j]*omega_vec[j]/(i_vec[j]*Vb)
#             if eff > 0:
#                 eff_mot.append(eff)
#                 omega_plot.append(omega_vec[j])
#     plt.plot(omega_plot, eff_mot, label = f'Kv = {Kv}')
#     plt.title("Efficiency by variation of Kv")

# for i0 in [0.1, 0.5, 1, 2, 5]:
#     i_vec = [(Vb - (omega/Kv))/Rm for omega in omega_vec]
#     Q_vec = [(i - i0)/Kv for i in i_vec]
#     eff_mot = []
#     omega_plot = []
#     for j in range(len(omega_vec)):
#         if (i_vec[j] < 0) or (i_vec[j] > 45):
#             pass
#         else:
#             eff = Q_vec[j]*omega_vec[j]/(i_vec[j]*Vb)
#             if eff > 0:
#                 eff_mot.append(eff)
#                 omega_plot.append(omega_vec[j])
#     plt.plot(omega_plot, eff_mot, label = f'i0 = {i0}')
#     plt.title("Efficiency by variation of i0")

# i_vec = [(Vb - (omega/Kv))/Rm for omega in omega_vec]
# Q_vec = [(i - i0)/Kv for i in i_vec]
# eff_mot = []
# omega_plot = []
# for j in range(len(omega_vec)):
#     if (i_vec[j] < 0) or (i_vec[j] > 45):
#         pass
#         print("here")
#     else:
#         eff = Q_vec[j]*omega_vec[j]/(i_vec[j]*Vb)
#         if eff > 0:
#             eff_mot.append(eff)
#             omega_plot.append(omega_vec[j])
# plt.plot(omega_plot, eff_mot)

plt.legend()
plt.show()