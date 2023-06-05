import math
import numpy as np
import matplotlib.pyplot as plt
import os
from tkinter import simpledialog
pi = np.pi
euler = np.e

def simple_pitch_inches2(r_vector, Pitch_inches): #only pitch in inches
    Pitch = Pitch_inches*0.0254
    return simple_pitch(r_vector, Pitch)

def simple_pitch(r_vector, Pitch): #All values in meters
    Beta_dist = []
    for rr in r_vector:
        Beta_dist.append(math.degrees(math.atan(Pitch/(2*pi*rr))))
    return Beta_dist

def apply_chord(x, y, chord):
    x_copy = x[:]
    y_copy = y[:]
    for i in range(len(x)):
        x_copy[i] *= chord
        y_copy[i] *= chord
    return x_copy, y_copy

def apply_offset_x(x, offset):
    x_copy = x[:]
    for i in range(len(x)):
        x_copy[i] += offset
    return x_copy

def apply_offset_y(y, offset):
    for i in range(len(y)):
        y[i] += offset
    return y

def apply_rotation_deg(x, y, angle):
    return apply_rotation_rad(x, y, math.radians(angle))

def apply_rotation_rad(x, y, angle):
    rot_mat = [[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]]
    x_rot = np.array([])
    y_rot = np.array([])
    for i in range(len(x)):
        case = np.array([[x[i]], [y[i]]])
        rotated = np.matmul(rot_mat, case)
        x_rot = np.append(x_rot, rotated[0])
        y_rot = np.append(y_rot, rotated[1])
    return x_rot, y_rot

def apply_thickness_on_camber(yl, yu, camber, thickness):
    thick_u = []
    thick_l = []
    for i in range(len(camber)):
        thick_u.append(yu[i] - camber[i])
        thick_l.append(camber[i] - yl[i])
    thick_l_after = [y*thickness for y in thick_l]
    thick_u_after = [y*thickness for y in thick_u]
    yu_after = []
    yl_after = []
    for i in range(len(camber)):
        yl_after.append(camber[i] - thick_l_after[i])
        yu_after.append(camber[i] + thick_u_after[i])
    return yl_after, yu_after

def apply_thickness_on_camber_easy(x, y, thickness):
    xl, yl, xu, yu = split_upper_lower(x, y)
    xu, camber_line = calculate_camber(xl, yl, xu, yu)
    yl, yu = apply_thickness_on_camber(yl, yu, camber_line, thickness)
    xu_copy = xu.copy()
    xu_reverse = xu_copy[::-1]
    yu_reverse = yu[::-1]
    xu = np.delete(xu, 0)
    yl = np.delete(yl, 0)
    x = np.append(xu_reverse, xu)
    y = np.append(yu_reverse, yl)
    return x, y


#auxiliary functions
def read_foil(airfoil_file, header_lines):
    x, y = np.loadtxt(airfoil_file, skiprows = header_lines, unpack = True)
    return x, y

def count_header_lines(file):
    with open(file, 'r') as inp:
        i = 0
        while True:
            try:
                content = inp.readline()
                cont = content.split()
                _ = float(cont[0])
                break
            except:
                i += 1
    return i

def split_upper_lower(xx, yy):
    x = np.array(xx.copy())
    y = np.array(yy.copy())
    i = 0
    diff_remember = math.inf
    while i < len(x):
        if x[i] > diff_remember:
            index0 = i - 1
            if len(x)%2 == 0:
                x_lower = x[i-1:]
                y_lower = y[i-1:]
                x_upper = np.flip(x[:i-1])
                y_upper = np.flip(y[:i-1])
                break
            else:
                x_lower = x[i-1:]
                y_lower = y[i-1:]
                x_upper = np.flip(x[:i])
                y_upper = np.flip(y[:i])
                break
        diff_remember = x[i]
        i += 1 
    delta = y[index0]
    y_upper = y_upper - delta
    y_lower = y_lower - delta
    delta = x[index0]
    x_upper = x_upper - delta
    x_lower = x_lower - delta
    return x_lower.tolist(), y_lower.tolist(), x_upper.tolist(), y_upper.tolist()

def calculate_camber(xl, yl, xu, yu):
    len_l = len(xl)
    len_u = len(xu)
    camber_line = []
    if len_l < len_u:
        for i in range(len_l):
            camber_line.append((yu[i] + yl[i])/2)
        return xl, camber_line
    else:
        for i in range(len_u):
            camber_line.append((yu[i] + yl[i])/2)
        return xu, camber_line


#====================================COMECO DO CODIGO===============================
Diametro_pol = 16   #int(input("Diametro [pol]: "))
Diametro_mm = Diametro_pol*25.4
Raio_mm = Diametro_mm/2
Passo_pol = 9       #int(input("Passo [pol]: "))
Passo_mm = Passo_pol*25.4

NumberOfStations = 21

RVector_adim = np.linspace(0, 1, NumberOfStations, endpoint = True).tolist()
# print(RVector_adim)

#====================================DEFINICAO DA HELICE===============================
#----Hélice T-Motor
OriginalAirfoil = [[1.00000, 0.99096, 0.97206, 0.94900, 0.92481, 0.90076, 0.87681, 0.85272, 0.82847, 0.80422, 0.78008, 0.75605, 0.73201, 0.70781, 0.68342, 0.65893, 0.63455, 0.61032, 0.58616, 0.56187, 0.53743, 0.51301, 0.48883, 0.46489, 0.44096, 0.41692, 0.39280, 0.36878, 0.34494, 0.32127, 0.29764, 0.27397, 0.25030, 0.22674, 0.20341, 0.18037, 0.15778, 0.13578, 0.11441, 0.09367, 0.07366, 0.05484, 0.03855, 0.02605, 0.01737, 0.01156, 0.00756, 0.00479, 0.00295, 0.00188, 0.00153, 0.00192, 0.00299, 0.00473, 0.00725, 0.01079, 0.01565, 0.02198, 0.03109, 0.04543, 0.06552, 0.08839, 0.11204, 0.13580, 0.16036, 0.18534, 0.21016, 0.23512, 0.26026, 0.28539, 0.31047, 0.33526, 0.35998, 0.38478, 0.40974, 0.43491, 0.46016, 0.48511, 0.50969, 0.53421, 0.55901, 0.58411, 0.60918, 0.63391, 0.65830, 0.68263, 0.70725, 0.73225, 0.75743, 0.78244, 0.80723, 0.83190, 0.85658, 0.88126, 0.90579, 0.92988, 0.95279, 0.97373, 0.99115, 1.00000], \
                   [0.00453, 0.00658, 0.01101, 0.01656, 0.02243, 0.02813, 0.03355, 0.03877, 0.04383, 0.04871, 0.05336, 0.05774, 0.06181, 0.06564, 0.06926, 0.07271, 0.07595, 0.07892, 0.08157, 0.08393, 0.08605, 0.08795, 0.08955, 0.09077, 0.09159, 0.09202, 0.09211, 0.09186, 0.09122, 0.09015, 0.08861, 0.08663, 0.08421, 0.08139, 0.07813, 0.07441, 0.07021, 0.06546, 0.06006, 0.05395, 0.04711, 0.03970, 0.03233, 0.02577, 0.02054, 0.01628, 0.01249, 0.00897, 0.00568, 0.00261, -0.0002, -0.0030, -0.0055, -0.0076, -0.0090, -0.0097, -0.0098, -0.0100, -0.0097, -0.0090, -0.0077, -0.0064, -0.0053, -0.0037, -0.0016, 0.00047, 0.00249, 0.00436, 0.00619, 0.00810, 0.01000, 0.01177, 0.01334, 0.01473, 0.01598, 0.01714, 0.01829, 0.01936, 0.02023, 0.02082, 0.02119, 0.02145, 0.02165, 0.02172, 0.02154, 0.02103, 0.02025, 0.01928, 0.01826, 0.01717, 0.01594, 0.01453, 0.01291, 0.01111, 0.00913, 0.00697, 0.00454, 0.00170, -0.0013, -0.0031]]
ChordVector_adim = []   #c/R
Beta_deg = []
Offset_adim = []        #1/c
FatorEspessura = []
for i in range(len(RVector_adim)):

    #Chord
    if i <= 8:
        ChordVector_adim.append(-775.7226922959*RVector_adim[i]**6 + 1170.2815273553*RVector_adim[i]**5 - 660.9381028488*RVector_adim[i]**4 + 161.9663843811*RVector_adim[i]**3 - 13.9599756954*RVector_adim[i]**2 + 0.1950071591*RVector_adim[i] + 0.1005241407)
    elif i <= 19:
        ChordVector_adim.append(-1.8598247908*RVector_adim[i]**4 + 4.7243249915*RVector_adim[i]**3 - 4.4099344990*RVector_adim[i]**2 + 1.5997758465*RVector_adim[i] + 0.0076143297)
    else:
        ChordVector_adim.append(0.009)

    #Beta
    if i == 0:
        Beta_deg.append(0)
    elif i == 1:
        Beta_deg.append(-0.0063*Passo_pol**2 + 0.475*Passo_pol + 0.2)
    elif i == 2:
        Beta_deg.append(-0.0208*Passo_pol**2 + 1.3333*Passo_pol - 4E-14)
    elif i == 3:
        Beta_deg.append(-0.0376*Passo_pol**2 + 2.4821*Passo_pol)
    elif i == 4:
        Beta_deg.append(math.degrees(math.atan(Passo_mm/(pi*2*Raio_mm*RVector_adim[5]))) - 7)
    else:
        Beta_deg.append(math.degrees(math.atan(Passo_mm/(pi*2*Raio_mm*RVector_adim[i]))))

    #Offset
    if i <= 6:
        Offset_adim.append(27252.63302612300000000000*RVector_adim[i]**6 - 24869.73093736170000000000*RVector_adim[i]**5 + 8475.72818899154000000000*RVector_adim[i]**4 - 1330.68601571768000000000*RVector_adim[i]**3 + 98.21577155683180000000*RVector_adim[i]**2 - 2.82604525527859000000*RVector_adim[i])
    elif i <= 11:
        Offset_adim.append(-101.51524423435300000000*RVector_adim[i]**5 + 179.22220443189100000000*RVector_adim[i]**4 - 115.59856522083200000000*RVector_adim[i]**3 + 30.86231640912590000000*RVector_adim[i]**2 - 2.37040664989036000000*RVector_adim[i])
    elif i <= 17:
        Offset_adim.append(-40.70481866598120000000*RVector_adim[i]**6 + 175.21347379684400000000*RVector_adim[i]**5 - 293.00339770317000000000*RVector_adim[i]**4 + 238.77772212028500000000*RVector_adim[i]**3 - 95.69989842176430000000*RVector_adim[i]**2 + 15.33286256343120000000*RVector_adim[i])
    else:
        Offset_adim.append(5183.29989910027000000000*RVector_adim[i]**3 - 13989.28422993000000000000*RVector_adim[i]**2 + 12572.05843030770000000000*RVector_adim[i] - 3762.21499505409000000000)

    #Fator Espessura
    C = ChordVector_adim[i]*Raio_mm
    if i == 0:
        FatorEspessura.append(1.85)
    elif i == 1:
        FatorEspessura.append(1.85)
    elif i == 2:
        FatorEspessura.append(1.85)
    elif i == 3:
        FatorEspessura.append(1.85)
    elif i == 4:
        FatorEspessura.append(1.25)
    else:
        FatorEspessura.append(1.0)

#====================================COMECO DOS MACROS===============================
outfile = simpledialog.askstring("Filename", "Coloque o nome do arquivo de saída: ")
sketchfile = outfile + "_(esboço).swb"
outfile = outfile + "_(curvas).swb"     #nomeia os arquivos de saída
try:
    os.remove(outfile)
    os.remove(sketchfile)   #tenta deletar os arquivos se já existirem
except:
    pass                    #caso o processo de deletar os arquivos der algum erro (qualquer erro), ele ignora, não tenta de novo e continua o código

#====================================ESCREVER MACROS===============================
with open(sketchfile, 'w') as s: #abre o objeto como um arquivo a ser escrito ('w')
    with open(outfile, 'w') as o:
        s.write("Dim swApp As Object\n\nDim Part As Object\nDim skPoint As Object\nDim skSegment As Object\nDim boolstatus As Boolean\nDim longstatus As Long, longwarnings As Long\nDim myRefPlane As Object\n\nSub main()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n")   #burocracia do vba para o sketch file
        o.write("Dim swApp As Object\n\nDim Part As Object\nDim boolstatus As Boolean\nDim longstatus As Long, longwarnings As Long\n\nSub main()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n") #Burocracia do vba para o outfile
        #Nas duas linhas acima, apenas escreve em cada arquivo o cabeçalho do código em VBA
        
        curves_counter = 0
        lista = []
        for i in range(len(RVector_adim)):
            XCoords = OriginalAirfoil[0].copy()
            YCoords = OriginalAirfoil[1].copy()

            XCoords, YCoords = apply_thickness_on_camber_easy(XCoords, YCoords, FatorEspessura[i])
            XCoords = apply_offset_x(XCoords, -0.5)
            XCoords = apply_offset_x(XCoords, Offset_adim[i])
            XCoords, YCoords = apply_chord(XCoords, YCoords, ChordVector_adim[i]*Raio_mm)
            XCoords, YCoords = apply_rotation_deg(XCoords, YCoords, Beta_deg[i])
            FinalAirfoil = [[RVector_adim[i]*Raio_mm]*len(XCoords), XCoords, YCoords]

            o.write("Part.InsertCurveFileBegin\n")      #Como cada arquivo é uma curva, começa a inserir a curva
            for j in range(len(FinalAirfoil[0])):
                o.write(f"boolstatus = Part.InsertCurveFilePoint({FinalAirfoil[0][j]}, {FinalAirfoil[1][j]}, {FinalAirfoil[2][j]})\n")  #Escreva o ponto a ser inserido
            o.write("boolstatus = Part.InsertCurveFileEnd()\n") #Termine o processo de inserção de curva

            curves_counter += 1                         #Conte o número de curvas inseridas 

            if (curves_counter%10 == 0):                    #Divida os processos no VBA pq ele reclama quando o arquivo é muito grande
                o.write(f"\nCall main{int((curves_counter/10) + 1)}()\nEnd Sub" + '\n'*5 + f"Sub main{int((curves_counter/10) + 1)}()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n")

            a,b = FinalAirfoil[1][0], FinalAirfoil[2][0]
            c,d = FinalAirfoil[1][-1], FinalAirfoil[2][-1]
            lista.extend([-b, a, 0, -d, c, 0])

        counter = 1
        x1 = 0
        while counter <= curves_counter: #Para cada curva, crie um plano
            s.write(f"\nboolstatus = Part.Extension.selectbyid2(\"Curva{counter}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nSet myRefPlane = Part.FeatureManager.InsertRefPlane(4, 0, 0, 0, 0, 0)") #escreve a criação de plano com base na curva que acabou de ser criada  
            counter += 1    
        counter = 1
        x1 = 0
        while counter <= curves_counter: #Para cada curva, crie e feche o esboço
            if (counter%10 == 0):
                s.write(f"\nCall main{int((counter/10) + 1)}()\nEnd Sub" + '\n'*5 + f"Sub main{int((counter/10) + 1)}()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n")
            s.write(f"\n\n\n\nboolstatus = Part.Extension.SelectByID2(\"Plano{counter}\", \"PLANE\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchManager.InsertSketch True\nboolstatus = Part.Extension.selectbyid2(\"Curva{counter}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.SketchManager.SketchUseEdge3(False, False)\nSet skSegment = Part.SketchManager.CreateLine({lista[x1]}, {lista[x1 + 1]}, {lista[x1 + 2]}, {lista[x1 + 3]}, {lista[x1 + 4]}, {lista[x1 + 5]})\nSet skSegment = Part.SketchManager.CreateCenterLine(-0.03, -0.05, 0, -0.04, -0.06, 0)\nSet skPoint = Part.SketchManager.CreatePoint(0, 0, 0)\n\nboolstatus = Part.Extension.SelectByID2(\"Point1\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point2\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point3\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point4\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point5\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point6\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)")
            s.write("\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Line1\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\n\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgFIXED\"\nboolstatus = Part.Extension.SelectByID2(\"Line1\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgFIXED\"\nboolstatus = Part.Extension.SelectByID2(\"Line1\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)")
            s.write("\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgPARALLEL\"\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgTANGENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point6\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point6\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point5\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point5\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nPart.SketchManager.InsertSketch True")
            #s.write("\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.EditDelete\nPart.SketchManager.InsertSketch True")
            counter += 1
            x1 = x1 + 6

        s.write("\nEnd Sub")
        o.write("\nEnd Sub") #Escreva para os macros em VBA terminarem
# for i in range(len(RVector_adim)):
#     XCoords = OriginalAirfoil[0].copy()
#     YCoords = OriginalAirfoil[1].copy()

#     XCoords, YCoords = apply_thickness_on_camber_easy(XCoords, YCoords, FatorEspessura[i])
#     XCoords = apply_offset_x(XCoords, -0.5)
#     XCoords = apply_offset_x(XCoords, Offset_adim[i])
#     XCoords, YCoords = apply_chord(XCoords, YCoords, ChordVector_adim[i]*Raio_mm)
#     XCoords, YCoords = apply_rotation_deg(XCoords, YCoords, Beta_deg[i])
#     FinalAirfoil = [[RVector_adim[i]*Raio_mm]*len(XCoords), XCoords, YCoords]
#     with open(f'{i + 1}.txt', 'w') as OutFile:
#         for j in range(len(XCoords)):
#             OutFile.write(f'{FinalAirfoil[0][j]},{FinalAirfoil[1][j]},{FinalAirfoil[2][j]}\n')
    # ax.plot3D(FinalAirfoil[0], FinalAirfoil[1], FinalAirfoil[2])
# plt.show()
