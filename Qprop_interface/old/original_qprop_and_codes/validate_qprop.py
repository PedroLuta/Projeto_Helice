from importing import *
import simulate
import pandas as pd
import Auxiliary

airfoil = 'airfoils\\NACA3408.txt'

Diameter_inches = 25
Blades = 4
V = -12
rpm = 2972

#25 - 4 pás
r_vector = [0.0635, 0.0889, 0.1143, 0.1397, 0.1651, 0.1905, 0.2159, 0.2413, 0.2667, 0.2921, 0.3175]
chord_vector = [0.05071150166862077, 0.05996753007121807, 0.05768077702596391, 0.05309911111279459, 0.04806618715033626, 0.04319872897387827, 0.03854172371102357, 0.0338110311345318, 0.028384754248340135, 0.020880684085393363, 0]
#Passo Subida (rpm 2972)
Betas = [42.44888894494408, 30.848360995261963, 24.853497285567453, 21.104093484163368, 18.54829284926484, 16.695575820390268, 16.305115895281034, 15.197140576615206, 14.304711284941176, 13.571271764540654, 0.0]
#Passo Cruzeiro (rpm 2514)
#Betas = [56.4330418227288, 44.069283032288254, 35.98489615919879, 30.651892949738684, 26.90518496935885, 24.138218716458734, 22.015349201040877, 20.335142236221486, 18.97358372365689, 17.84910553692154, 0.0]



# #25 polegadas - 2 pás
# r_vector = [0.0635, 0.0889, 0.1143, 0.1397, 0.1651, 0.1905, 0.2159, 0.2413, 0.2667, 0.2921, 0.3175]
# chord_vector = [0.08160809858217317, 0.07423668712269352, 0.067, 0.06, 0.05339666808125843, 0.04604578660344783, 0.03950339075434698, 0.03325316476044051, 0.02668162100820269, 0.01674304220145237, 0]
# #Passo Subida (rpm 4286)
# Betas = [29.63962979020654, 22.596023769994304, 18.587210475126916, 16.066195084289763, 14.3363999744356, 13.075764141121162, 12.115938259093284, 11.360511377492532, 10.750483995557822, 11.249241739697437, 0.0]
# #Passo Cruzeiro (rpm 3720)
# #Betas = [40.44944262534906, 31.091650863508967, 25.421451096371094, 21.750332105850088, 19.191016850850033, 17.30783042271711, 15.865076039425599, 14.724712131248488, 14.804884583865626, 14.040466329251814, 0.0]



# #15 polegadas - 4 pás
# r_vector = [0.0381,  0.05334, 0.06858, 0.08382, 0.09906, 0.1143, 0.12954, 0.14478, 0.16002, 0.17526, 0.1905]
# chord_vector = [0.06659883387308367, 0.06863299926342242, 0.06393414017715361, 0.06467690732222628, 0.05767946309057695, 0.051200641527185675, 0.04519809922627309, 0.03925268831976199, 0.032617247706923265, 0.02373794829137861, 0]
# #Passo Subida (rpm 6436)
# #Betas = [36.351841448433746, 26.503867344756923, 21.070894405175547, 17.833167336508097, 15.691751178004537, 14.170019712032808, 13.031786143258593, 12.146344263340685, 12.444947409666149, 11.861770633582532, 0.0]
# #Passo Cruzeiro (rpm 5574)
# Betas = [49.427640737217004, 36.736812404531065, 29.126668649203484, 24.427237989674524, 21.261657858729194, 18.988297831086896, 17.276978077127964, 16.95933327149443, 15.880847937561002, 14.9979270691836, 0.0]


    
# #15 polegadas - 2 pás
# r_vector = [0.0381,  0.05334, 0.06858, 0.08382, 0.09906, 0.1143,  0.12954, 0.14478, 0.16002, 0.17526, 0.1905]
# chord_vector = [0.06050352544917628, 0.055206206755384166, 0.04487421657582313, 0.03704451874511294, 0.03094411215586708, 0.02598765131199184, 0.02178152418310504, 0.017896321000794975, 0.014038324212474444, 0.009673340166681949, 0]
# #Passo Subida (rpm 9192):
# #Betas = [25.434677738037525, 19.23709181002249, 15.83379545403429, 13.751526383372685, 12.346584404246064, 11.333672703781438, 10.568138415158703, 9.968500916780577, 9.486080819409555, 9.089707102809784, 0.0] 
# #Passo Cruzeiro (rpm 8031):
# Betas = [34.968835426342736, 26.29964016644934, 21.33537111487524, 18.247406744072748, 16.145962326722596, 14.623419628507577, 13.469201570990142, 12.563291618109913, 11.833544487970457, 12.236603711125246, 0.0]




R = Diameter_inches*0.0254/2
dT_vector, dQ_vector, r_vector, Re_vector, WA_vector, Cl_vector, Cd_vector = simulate.qprop_fixed_pitch(V, rpm*2*pi/60, Blades, R, r_vector, Betas, chord_vector, airfoil = airfoil)

# array = [[x/R for x in r_vector], r_vector, Cl_vector, Cd_vector, Re_vector, WA_vector]
# df = pd.DataFrame(array).T
# df.to_excel(excel_writer = "excel.xlsx")