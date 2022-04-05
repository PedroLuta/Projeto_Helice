from mosa.mosalib.mosa import Anneal
import subprocess as sp
from PARSEC_functions import *

population = {"ate":[], "bte":[], "dzte":[], "rup":[], "xup":[], "zup":[], "zxxup":[], "rlo":[], "xlo":[], "zlo":[], "zxxlo":[]} # "zte":[],

def func(solution):
    Rey, alpha, itr = 200000, 0, 250
    airfoil_file = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\airfoils\airfoil.txt'
    
    ate = solution["ate"][0]
    bte = solution["bte"][0]
    #zte = solution["zte"][0]
    dzte = solution["dzte"][0]
    rup = solution["rup"][0]
    xup = solution["xup"][0]
    zup = solution["zup"][0]
    zxxup = solution["zxxup"][0]
    rlo = solution["rlo"][0]
    xlo = solution["xlo"][0]
    zlo = solution["zlo"][0]
    zxxlo = solution["zxxlo"][0]
    params = [ate, bte, 0, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo]
    aup, alo = coeffs_from_params(ate, bte, 0, dzte, rup, xup, zup, zxxup, rlo, xlo, zlo, zxxlo)
    if not check_valid(aup, alo, params, checkWavy = True):
        return [0, 1000]
    write_check = write_from_params_vec(params, airfoil_file)
    if not write_check:
        return [0, 1000]

    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load airfoils\\airfoil.txt\n\noper\niter {itr}\nvisc {Rey}\na {alpha}\n\n\n\nquit \n', timeout = 2.5) 
    except:
        ps.kill() 
        return [0, 1000]
    
    list = output.split()
    i = 0
    id_remember = 0
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]

        if (list[i].lower() == "rms:") and (id_remember < itr):
            id_remember = int(list[i - 1])
        elif id_remember >= itr:
            ps.kill()
            return [0, 1000]
        
        if (i == len(list) - 1) and (id_remember < itr):
            ps.kill()
            return [-float(cl), float(cd)]
        i += 1
    ps.kill()
    return 0, 1000

xstep10 = { "ate":4.5,      "bte":4.5,  "dzte":0.005,   "rup":0.01,     "xup":0.05,     "zup":0.02,     "zxxup":0.5,    "rlo":0.01,    "xlo":0.05,     "zlo":0.01,      "zxxlo":1} #"zte":0.01,     
# xstep5 = {  "ate":2.25,     "bte":2.25, "dzte":0.0025,  "rup":0.005,    "xup":0.025,    "zup":0.01,   "zxxup":0.25,   "rlo":0.005,   "xlo":0.025,    "zlo":0.005,     "zxxlo":0.5}#"zte":0.005,    
# xstep1 = {  "ate":0.45,     "bte":0.45, "dzte":0.0005,  "rup":0.001,    "xup":0.005,    "zup":0.002,    "zxxup":0.05,   "rlo":0.001,   "xlo":0.005,    "zlo":0.001,     "zxxlo":0.1}#"zte":0.001,    

#"zte":[-0.05, 0.05], 
#"zte":0.5,                 

opt = Anneal()
opt.setarchiveconfig(archivesize=100,maxarchivereject=1000000)
opt.settemp(initemp=750,niter=1000,ntemp=20)
opt.setxconfig( xbounds         ={"ate":[-30, 15],  "bte":[0, 45],  "dzte":[0, 0.05],   "rup":[0, 0.1], "xup":[0, 0.5], "zup":[0.0001, 0.2],   "zxxup":[-5, 0],    "rlo":[0, 0.1], "xlo":[0, 0.5],     "zlo":[-0.1, -0.0001], "zxxlo":[0, 10]}, \
                selweight       ={"ate":1.0,        "bte":1.0,      "dzte":0.1,         "rup":1.0,      "xup":1.0,      "zup":1.0,              "zxxup":1.0,        "rlo":1.0,      "xlo":1.0,          "zlo":1.0,              "zxxlo":1.0}, \
                xstep           =xstep10)
                
                #xnel            ={"ate":1,          "bte":1,        "zte":1,                "dzte":1,           "rup":1,        "xup":1,        "zup":1,                "zxxup":1,          "rlo":1,        "xlo":1,        "zlo":1,                "zxxlo":1}, \
                #exchangeprob    ={"ate":1.0,        "bte":1.0,      "zte":1.0,              "dzte":1.0,         "rup":1.0,      "xup":1.0,      "zup":1.0,              "zxxup":1.0,        "rlo":1.0,      "xlo":1.0,      "zlo":1.0,              "zxxlo":1.0},\
                #xdistinct={"Component":True},
opt.setpopulation(population)

opt.evolve(func)

opt.printx()
opt.plotfront()