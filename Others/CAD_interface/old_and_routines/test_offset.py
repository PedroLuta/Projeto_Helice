import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
import numpy as nup
import os

# Interface Config
root= tk.Tk()
root.title('TXT To VBA Formatter')
canvas1 = tk.Canvas(root, width = 300, height = 150, bg = 'lightgray')
canvas1.pack()
# Fim de Interface Config

def apply_offset(x, y, offset, invert = True):
    x_delta = nup.max(x) - nup.min(x)
    y_delta = nup.max(y) - nup.min(y)
    V = nup.array([[0, 0]])
    i = 1
    while i < len(x) - 1:
        tanvec_x = (x[i + 1] - x[i - 1])/2
        tanvec_y = ((y[i + 1] - y[i - 1])/2)*(x_delta/y_delta)
        if invert:
            tanvec_x = -tanvec_x
            tanvec_y = -tanvec_y
        normvec_x = -tanvec_y
        normvec_y = tanvec_x
        placeholder = ((normvec_x**2) + (normvec_y**2))**0.5
        normvec_x /= placeholder
        normvec_y /= placeholder
        V = nup.append(V, [[normvec_x, normvec_y]], axis = 0)
        i += 1
    x_offset = nup.array([])
    y_offset = nup.array([])
    i = 1
    while i < len(x) - 1:
        x_offset = nup.append(x_offset, x[i] + V[i][0]*offset)
        y_offset = nup.append(y_offset, y[i] + V[i][1]*offset/(x_delta/y_delta))
        i += 1
    return x_offset, y_offset

def Routine():
    outfile = simpledialog.askstring("Filename", "Coloque o nome do arquivo de saída: ")
    offset = simpledialog.askfloat("Offset", "Offset (em metros): ")
    if outfile == "":
        return
    sketchfile = outfile + " (esboço).swb"
    outfile = outfile + " (curvas).swb"
    try:
        os.remove(outfile)
        os.remove(sketchfile)   #tenta deletar os arquivos se já existirem
    except:
        pass                    #caso o processo de deletar os arquivos der algum erro (qualquer erro), ele ignora, não tenta de novo e continua o código
    with open(sketchfile, 'w') as s:
        with open(outfile, 'w') as o:
            s.write("Dim swApp As Object\n\nDim Part As Object\nDim skPoint As Object\nDim skSegment As Object\nDim boolstatus As Boolean\nDim longstatus As Long, longwarnings As Long\nDim myRefPlane As Object\n\nSub main()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n")   #burocracia do vba para o sketch file
            o.write("Dim swApp As Object\n\nDim Part As Object\nDim boolstatus As Boolean\nDim longstatus As Long, longwarnings As Long\n\nSub main()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n") #Burocracia do vba para o outfile
            curves_counter = 0
            airfoil_files = filedialog.askopenfilenames()
            for foil in airfoil_files: 
                xxx = nup.array([])
                yyy = nup.array([])
                zzz = nup.array([]) 
                o.write("Part.InsertCurveFileBegin\n")
                with open(foil, 'r') as i:
                    while True:
                        content = i.readline()
                        if content == "":
                            o.write("boolstatus = Part.InsertCurveFileEnd()\n")
                            o.write("Part.InsertCurveFileBegin\n") 
                            xxx, yyy = apply_offset(xxx, yyy, offset)
                            i = 0
                            while i < len(xxx):
                                o.write(f"boolstatus = Part.InsertCurveFilePoint({zzz[i]}, {xxx[i]}, {yyy[i]})\n")
                                i += 1
                            o.write("boolstatus = Part.InsertCurveFileEnd()\n")
                            break
                        content.replace("\n", "")
                        x = content.split(",")
                        x[0] = float(x[0])/1000            #divide o valor por 1000 porque vba considera todas as medidas como metro
                        x[1] = float(x[1])/1000
                        x[2] = float(x[2])/1000
                        o.write(f"boolstatus = Part.InsertCurveFilePoint({x[0]}, {x[1]}, {x[2]})\n")
                        zzz = nup.append(zzz, x[0])
                        xxx = nup.append(xxx, x[1])
                        yyy = nup.append(yyy, x[2])
                    
                    curves_counter += 2
                if (curves_counter%10 == 0):
                    o.write(f"\nCall main{int((curves_counter/10) + 1)}()\nEnd Sub" + '\n'*5 + f"Sub main{int((curves_counter/10) + 1)}()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n")
            counter = 1
            x1 = 0
            while counter <= curves_counter:
                s.write(f"\nboolstatus = Part.Extension.selectbyid2(\"Curva{counter}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nSet myRefPlane = Part.FeatureManager.InsertRefPlane(4, 0, 0, 0, 0, 0)") #escreve a criação de plano com base na curva que acabou de ser criada  
                counter += 2    
            counter = 1
            x1 = 0
            while counter <= curves_counter + 2:
                s.write(f"\n\n\n\nboolstatus = Part.Extension.SelectByID2(\"Plano{int(counter/2)}\", \"PLANE\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchManager.InsertSketch True\nboolstatus = Part.Extension.selectbyid2(\"Curva{int(counter - 2)}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.SketchManager.SketchUseEdge3(False, False)\nboolstatus = Part.Extension.selectbyid2(\"Curva{int(counter - 1)}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.SketchManager.SketchUseEdge3(False, False)\nPart.SketchManager.InsertSketch True")
                counter += 2
                x1 = x1 + 6
            s.write("\nEnd Sub")
            o.write("\nEnd Sub")






# Interface Config e funcionalidade do botão
browseButton = tk.Button(text='Import TXT File', command=Routine, bg='green', fg='white', font=('helvetica', 12, 'bold'))
canvas1.create_window(150, 75, window=browseButton)
# Fim de Interface Config e funcionalidade do botão, aqui a única função do botão é chamar a função Routine()

root.mainloop() #Chama a Interface e abre a janela com o botão