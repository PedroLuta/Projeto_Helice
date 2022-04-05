#default imports
import numpy as nup
import math
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import messagebox 

#functions import
from xfoil_interface import *
from induction import *
from PARSEC_functions import *


# Interface Config
root= tk.Tk()
root.title('Otimizador de Incidências')
canvas1 = tk.Canvas(root, width = 450, height = 250, bg = 'lightgray')
canvas1.pack()
# Fim de Interface Config

def chord_distribution(x, factor = 1):
    if x > 0.4:
        return factor*(-1.8598247908*(x**4) + 4.7243249915*(x**3) - 4.4099344990*(x**2) + 1.5997758465*x + 0.0076143297)
    else:
        return factor*(-775.7226922959*(x**6) + 1170.2815273553*(x**5) - 660.9381028488*(x**4) + 161.9663843811*(x**3) - 13.9599756954*(x**2) + 0.1950071591*x + 0.1005241407)

#main starts here
def main():
    pi = nup.pi
    euler = nup.e
    #operating conditions
    rho = simpledialog.askfloat("Input", "Densidade do ar [kg/(m^3)]", parent=root, minvalue=0.0, maxvalue=3.0) #[kg/(m^3)]
    kvisc = simpledialog.askfloat("Input", "Viscosidade cinemática [kg/(m*s)]", parent=root, minvalue=0.0, maxvalue=1) #[kg/(m*s)]
    rpm = simpledialog.askfloat("Input", "rpm [1/(s*60)]", parent=root, minvalue=0.0, maxvalue=10000.0)
    vi = simpledialog.askfloat("Input", "Velocidade da aeronave [m/s]", parent=root, minvalue=0.0, maxvalue=100.0) #[m/s]
    radps = (rpm*2*pi)/60 #[1/s]

    #propeller characterization
    Blades = simpledialog.askinteger("Input", "Quantas pás?", parent=root, minvalue=0, maxvalue=10)
    sections = simpledialog.askinteger("Input", "Quantas seções utilizar para o cálculo?", parent=root, minvalue=0, maxvalue=100) #more = more precision
    p = simpledialog.askfloat("Input", "Posição radial (em porcentagem 0-1) do primeiro perfil", parent=root, minvalue=0.0, maxvalue=1.0) #radius (percentage) at which there's no more geometric constraints regarding the hub
    R = simpledialog.askfloat("Input", "Raio da hélice [m]", parent=root, minvalue=0.0, maxvalue=1.0)
    a1 = simpledialog.askfloat("Input", "alpha inicial", parent=root, minvalue=0.0, maxvalue=10.0)
    a2 = simpledialog.askfloat("Input", "alpha final", parent=root, minvalue=a1, maxvalue=25.0)
    astep = simpledialog.askfloat("Input", "passo de alpha", parent=root, minvalue=0.0, maxvalue=a2 - a1)

    #initialization
    Thrust = 0
    Torque = 0
    Beta_vector = []
    radius_vector = []
    Thrust_vector = []
    Torque_vector = []
    first = True

    airfoil_file = 'airfoil_clarky.txt'

    for rr in nup.linspace(p*R, R, num = sections):
        dT = 0
        dQ = 0
        chord = R*chord_distribution(rr/R)
        Vr = radps*rr
        V = ((Vr**2)+(vi**2))**0.5
        Re = ((V*chord)/kvisc) 
        alpha, Cl, Cd, _ = communicate_range_flexible(Re, a1, a2, astep, afile = airfoil_file, timeout = 25)
        try:
            dT, dQ, phi, _, _ = jitted_induction(radps, rr, Cl, Cd, Blades, rho, R, chord, vi)
        except:
            dT, dQ, phi, _, _ = 0, 0, 0, 0, 0 
        if first:
            dQ_remember = dQ
            dT_remember = dT
            r_remember = rr
            first = False
        elif (dT == 0 and dQ == 0) or math.isnan(dT) or math.isnan(dQ):
            pass
        else:
            Beta_vector.append(alpha + math.degrees(phi))
            radius_vector.append(rr)
            Torque_vector.append(dQ)
            Thrust_vector.append(dT)
            Thrust += ((dT + dT_remember)/2)*(rr - r_remember)
            Torque += ((dQ + dQ_remember)/2)*(rr - r_remember)
            dQ_remember = dQ
            dT_remember = dT
            r_remember = rr

    print(Thrust)
    print(Torque)
    print(Beta_vector)
    print(radius_vector)




browseButton = tk.Button(text='Run', command=main, bg='green', fg='white', font=('helvetica', 12, 'bold'))
canvas1.create_window(300, 150, window=browseButton)
root.mainloop() #Chama a Interface e abre a janela com o botão