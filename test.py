import random
import numpy as np
import matplotlib.pyplot as plt

r_R = np.linspace(0, 1, 101, endpoint = True)
# ChordVec_m_R = [-1]
# i = 0
# while i < 50:
#     a = random.uniform(-1, 5)
#     b = random.uniform(-5, 5)
#     c = random.uniform(-5, 5)



#     ChordSquaredVec_m2_R2 = [(a*(r**2)) + (b*r) + c for r in r_R]
#     if any(chord < 0 for chord in ChordSquaredVec_m2_R2):
#         continue
#     ChordVec_m_R = [ChordSquared**0.5 for ChordSquared in ChordSquaredVec_m2_R2]
#     if any(chord < 0.05 for chord in ChordVec_m_R) or any(chord > 0.5 for chord in ChordVec_m_R):
#         continue
#     plt.plot(r_R, ChordVec_m_R)
#     i += 1

# print(a)
# print(b)
# print(c)
# plt.ylim([0,1])
# plt.show()
Ctip_adim = 0.1
Croot_adim = 0.2
Cmax_adim = 0.3
rmax_adim = 0.7

# b = (((Croot_adim**2)*((rmax_adim**2) - 1)) - ((Ctip_adim**2)*(rmax_adim**2)) + (Cmax_adim**2))/(rmax_adim*(1 - rmax_adim))
# a = (Ctip_adim**2) - (Croot_adim**2) - b
# c = Croot_adim**2

# plt.plot(r_R, [((a*(r**2)) + (b*r) + c)**0.5 for r in r_R])
# # plt.ylim([0,1])
# print([((a*(r**2)) + (b*r) + c)**0.5 for r in r_R])
i = 0
while i < 50:
    Ctip_adim = random.uniform(0, 0.3)
    Croot_adim = random.uniform(0, 0.3)
    Cmax_adim = random.uniform(0, 0.5)
    rmax_adim = random.uniform(0, 1)

    a = -(2*Cmax_adim**2*rmax_adim - Cmax_adim**2 + Croot_adim**2*rmax_adim**2 - 2*Croot_adim**2*rmax_adim + Croot_adim**2 - Ctip_adim**2*rmax_adim**2)/(rmax_adim**2*(rmax_adim - 1)**2)
    b = -(- 3*Cmax_adim**2*rmax_adim**2 + Cmax_adim**2 - 2*Croot_adim**2*rmax_adim**3 + 3*Croot_adim**2*rmax_adim**2 - Croot_adim**2 + 2*Ctip_adim**2*rmax_adim**3)/(rmax_adim**2*(rmax_adim - 1)**2)
    c = -(3*Cmax_adim**2*rmax_adim - 2*Cmax_adim**2 + Croot_adim**2*rmax_adim**3 - 3*Croot_adim**2*rmax_adim + 2*Croot_adim**2 - Ctip_adim**2*rmax_adim**3)/(rmax_adim*(rmax_adim - 1)**2)
    d = Croot_adim**2

    DerivativeZero = c/(2*(d**0.5))
    if DerivativeZero < 0:
        continue
    DerivativeOne = (3*a + 2*b + c)/(2*((a + b + c + d)**0.5))
    if DerivativeOne > 0:
        continue

    ChordVec_m_R = [((a*(r**3)) + (b*(r**2)) + (c*r) + d)**0.5 for r in r_R]
    plt.plot(r_R, ChordVec_m_R)
    i += 1

plt.show()