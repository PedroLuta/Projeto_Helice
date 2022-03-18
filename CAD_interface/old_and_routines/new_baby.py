import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
import os

# Interface Config
root= tk.Tk()
root.title('TXT To VBA Formatter')
canvas1 = tk.Canvas(root, width = 300, height = 150, bg = 'lightgray')
canvas1.pack()
# Fim de Interface Config

def Routine():
    outfile = simpledialog.askstring("Filename", "Coloque o nome do arquivo de saída: ")
    if outfile == "":
        return
    outfile = outfile.replace(".txt", "")
    sketchfile = outfile + " (sketch).swb"
    outfile = outfile + " (vba).swb"
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
            lista = []
            while True:
                airfoil_file = filedialog.askopenfilename()
                if airfoil_file == "":
                    o.write("\nEnd Sub")
                    break
                else:
                    o.write("Part.InsertCurveFileBegin\n")
                    with open(airfoil_file, 'r') as i:
                        first = True
                        a,b = 0,0
                        c,d = 1,1
                        while True:
                            content = i.readline()
                            if content == "":
                                lista.extend([-b, a, 0, -d, c, 0])
                                break
                            content.replace("\n", "")
                            x = content.split(",")
                            x[0] = float(x[0])/1000            #divide o valor por 1000 porque vba considera todas as medidas como metro
                            x[1] = float(x[1])/1000
                            x[2] = float(x[2])/1000
                            if first:
                                a = x[1]
                                b = x[2]
                                first = False
                            c = x[1]
                            d = x[2]
                            o.write(f"boolstatus = Part.InsertCurveFilePoint({x[0]}, {x[1]}, {x[2]})\n")
                        o.write("boolstatus = Part.InsertCurveFileEnd()\n")
                        curves_counter += 1
                if curves_counter%10 == 0:
                    o.write(f"\nCall main{int((curves_counter/10) + 1)}()\nEnd Sub" + '\n'*5 + f"Sub main{int((curves_counter/10) + 1)}()\n\nSet swApp = Application.SldWorks\n\nSet Part = swApp.ActiveDoc\n\n\n")
            counter = 1
            x1 = 0
            while counter <= curves_counter:
                s.write(f"\nboolstatus = Part.Extension.selectbyid2(\"Curva{counter}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nSet myRefPlane = Part.FeatureManager.InsertRefPlane(4, 0, 0, 0, 0, 0)") #escreve a criação de plano com base na curva que acabou de ser criada  
                counter += 1    
            counter = 1
            x1 = 0
            while counter <= curves_counter:
                s.write(f"\n\n\n\nboolstatus = Part.Extension.SelectByID2(\"Plano{counter}\", \"PLANE\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchManager.InsertSketch True\nboolstatus = Part.Extension.selectbyid2(\"Curva{counter}\", \"REFERENCECURVES\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.SketchManager.SketchUseEdge3(False, False)\nSet skSegment = Part.SketchManager.CreateLine({lista[x1]}, {lista[x1 + 1]}, {lista[x1 + 2]}, {lista[x1 + 3]}, {lista[x1 + 4]}, {lista[x1 + 5]})\nSet skSegment = Part.SketchManager.CreateCenterLine(-0.03, -0.05, 0, -0.04, -0.06, 0)\nSet skPoint = Part.SketchManager.CreatePoint(0, 0, 0)\n\nboolstatus = Part.Extension.SelectByID2(\"Point1\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point2\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point3\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point4\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point5\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Point6\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)")
                s.write("\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Line1\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchConstraintsDelAll\n\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgFIXED\"\nboolstatus = Part.Extension.SelectByID2(\"Line1\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgFIXED\"\nboolstatus = Part.Extension.SelectByID2(\"Line1\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)")
                s.write("\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgPARALLEL\"\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgTANGENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point6\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point6\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point5\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Line2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nboolstatus = Part.Extension.SelectByID2(\"Point5\", \"SKETCHPOINT\", 0, 0, 0, False, 0, Nothing, 0)\nboolstatus = Part.Extension.SelectByID2(\"Spline2\", \"SKETCHSEGMENT\", 0, 0, 0, True, 0, Nothing, 0)\nPart.SketchAddConstraints \"sgCOINCIDENT\"\nPart.SketchManager.InsertSketch True")
                counter += 1
                x1 = x1 + 6
            s.write("\nEnd Sub")






# Interface Config e funcionalidade do botão
browseButton = tk.Button(text='Import TXT File', command=Routine, bg='green', fg='white', font=('helvetica', 12, 'bold'))
canvas1.create_window(150, 75, window=browseButton)
# Fim de Interface Config e funcionalidade do botão, aqui a única função do botão é chamar a função Routine()

root.mainloop() #Chama a Interface e abre a janela com o botão