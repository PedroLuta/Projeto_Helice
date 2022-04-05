import numpy as nup
#from scipy import linalg
import math
import os
import subprocess as sp
#import time
from numba import njit, jit



def write_to_xfoil(x, y):
    y.stdin.write(x)

def reset_xfoil(y):
    y.stdin.write('\n\n\n')

#returns most efficient alpha of a fixed airfoil between the range (files)
def get_properties(Rey, a1, a2, astep):
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil_output.txt'
    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, universal_newlines = True, startupinfo = startupinfo, stdout = stout) # 

    write_to_xfoil('load airfoils\\airfoil.txt\n\n', ps)
    #write_to_xfoil(r'load C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\airfoils\airfoil.txt\n\n', ps)
    #write_to_xfoil('load airfoil_sunnysky.dat\n\n', ps)
    #write_to_xfoil('load airfoil_clarky.txt\n\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {int(round(Rey))}\n', ps)
    write_to_xfoil('seqp\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'aseq {a1} {a2} {astep}\n', ps)
    reset_xfoil(ps)
    write_to_xfoil('quit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = 5)
    except:
        #print("XFOIL timed out")
        ps.kill()
        return 0, 0, 0

    LDm = 0
    read_counter = 0
    try:
        with open(xfoil_outfile, 'r') as o:
            while True:
                if read_counter < 12:
                    content = o.readline()
                    read_counter += 1
                else:
                    content = o.readline()
                    if content == '':
                        break
                    cont = content.split()
                    alpha = float(cont[0])
                    Cl = float(cont[1])
                    Cd = float(cont[2])
                    if Cl/Cd > LDm:
                        LDm = Cl/Cd
                        Clm = Cl
                        Cdm = Cd
                        alpham = alpha
        try:
            os.remove(xfoil_outfile)
        except:
            pass
        try:
            ps.kill()
            return Clm, Cdm, alpham
        except UnboundLocalError:
            ps.kill()
            return 0, 0, 0

    except FileNotFoundError:
        ps.kill()
        return 0, 0, 0

#returns properties for a fixed airfoil at a specified alpha (file)
def get_properties_fixed_alpha(Rey, alpha):
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil_output.txt'
    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, startupinfo=startupinfo, universal_newlines = True, stdout = stout)

    #write_to_xfoil('load sunnysky_100_pontos.dat\n\n', ps)
    #write_to_xfoil('load clarky.txt\n\n', ps)
    write_to_xfoil('load airfoils\\airfoil.txt\n\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {Rey}\n', ps)
    write_to_xfoil('seqp\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'a {alpha}\n', ps)
    reset_xfoil(ps)
    write_to_xfoil('quit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = 10)
    except:
        #print("Timeout")
        ps.kill()
        return 0, 0, 0

    read_counter = 0
    try:
        with open(xfoil_outfile, 'r') as o:
            while True:
                if read_counter < 12:
                    content = o.readline()
                    read_counter += 1
                else:
                    content = o.readline()
                    if content == '':
                        break
                    cont = content.split()
                    alpha = float(cont[0])
                    Cl = float(cont[1])
                    Cd = float(cont[2])
        try:
            os.remove(xfoil_outfile)
        except:
            pass
        try:
            ps.kill()
            return Cl, Cd, alpha
        except UnboundLocalError:
            #print("Unbound Local")
            ps.kill()
            return 0, 0, 0

    except FileNotFoundError:
        #print("File Not Found")
        ps.kill()
        return 0, 0, 0

#returns cl, cd, cm of a fixed airfoil at a specified alpha (communicate)
def communicate_fixed(Rey, alpha):
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load airfoils\\airfoil_clarky.txt\n\noper\niter {itr}\nvisc {Rey}\na {alpha}\n\n\n\nquit \n', timeout = 25) 
    except: 
        ps.kill()
        return 0, 1000, 1000
    
    list = output.split()
    i = 0
    id_remember = 0
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
        if list[i].lower() == "cm":
            cm = list[i + 2]

        if (list[i].lower() == "rms:") and (id_remember < itr):
            id_remember = int(list[i - 1])
        elif id_remember >= itr:
            ps.kill()
            return 0, 1000, 1000
        
        if (i == len(list) - 1) and (id_remember < itr):
            ps.kill()
            return cl, cd, cm
        i += 1
    ps.kill()
    return 0, 1000, 1000

#returns cl, cd, cm of a specified airfoil at a specified alpha (communicate)
def communicate_fixed_flexible(Rey, alpha, afile = 'airfoils\\airfoil.txt'):
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\noper\niter {itr}\nvisc {Rey}\na {alpha}\n\n\n\nquit \n', timeout = 25) 
    except: 
        ps.kill()
        return 0, 1000, 1000
    
    list = output.split()
    i = 0
    id_remember = 0
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
        if list[i].lower() == "cm":
            cm = list[i + 2]

        if (list[i].lower() == "rms:") and (id_remember < itr):
            id_remember = int(list[i - 1])
        elif id_remember >= itr:
            ps.kill()
            return 0, 1000, 1000
        
        if (i == len(list) - 1) and (id_remember < itr):
            ps.kill()
            return cl, cd, cm
        i += 1
    ps.kill()
    return 0, 1000, 1000


def communicate_range(Rey, a1, a2, astep):
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load airfoils\\airfoil_clarky.txt\n\noper\niter {itr}\nvisc {Rey}\naseq {a1} {a2} {astep}\n\n\n\nquit \n', timeout = 25) 
    except: 
        ps.kill()
        return 0, 1000, 1000
    
    list = output.split()
    i = 0
    id_remember = 0
    cl_list = []
    cd_list = []
    cm_list = []
    a_list = []
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
        if list[i].lower() == "cm":
            cm = list[i + 2]
        if list[i].lower() == "a":
            a = list[i + 2]
        
        #skip 
        if id_remember < itr:
            skip = False

        if list[i].lower() == "rms:":
            if (id_remember >= itr) and (not skip):
                a_list.append(float(a))
                cl_list.append(0.0)
                cd_list.append(1000.0)
                cm_list.append(1000.0)
                skip = True
            elif int(list[i - 1]) < id_remember:
                a_list.append(float(a))
                cl_list.append(float(cl))
                cd_list.append(float(cd))
                cm_list.append(float(cm))
            id_remember = int(list[i - 1])
        
        if (i == len(list) - 1) and (id_remember < itr):
            a_list.append(float(a))
            cl_list.append(float(cl))
            cd_list.append(float(cd))
            cm_list.append(float(cm))
        
        i += 1

    i = 0
    clcd_remember = 0
    while i < len(a_list):
        a_try = a_list[i]
        cl_try = cl_list[i]
        cd_try = cd_list[i]
        cm_try = cm_list[i]
        if cl_try/cd_try > clcd_remember:
            a = a_try
            cl = cl_try
            cd = cd_try
            cm = cm_try
            clcd_remember = cl_try/cd_try
        i += 1

    ps.kill()
    return a, cl, cd, cm


def communicate_range_flexible(Rey, a1, a2, astep, afile = 'airfoils\\airfoil.txt', timeout = 2.5):
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\xfoil.exe'
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
    itr = 250

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\noper\niter {itr}\nvisc {Rey}\naseq {a1} {a2} {astep}\n\n\n\nquit \n', timeout = timeout) 
    except: 
        ps.kill()
        return 0, 0, 1000, 1000
    
    list = output.split()
    i = 0
    id_remember = 0
    cl_list = []
    cd_list = []
    cm_list = []
    a_list = []
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
        if list[i].lower() == "cm":
            cm = list[i + 2]
        if list[i].lower() == "a":
            a = list[i + 2]
        
        if id_remember < itr:
            skip = False

        if list[i].lower() == "rms:":
            if (id_remember >= itr) and (not skip):
                try:
                    a_list.append(float(a))
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                    cm_list.append(1000.0)
                except:
                    a_list.append(0.0)
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                    cm_list.append(1000.0)
                skip = True
            elif int(list[i - 1]) < id_remember:
                try:
                    a_list.append(float(a))
                    cl_list.append(float(cl))
                    cd_list.append(float(cd))
                    cm_list.append(float(cm))
                except:
                    a_list.append(0.0)
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                    cm_list.append(1000.0)
            id_remember = int(list[i - 1])
        
        if (i == len(list) - 1) and (id_remember < itr):
            try:
                a_list.append(float(a))
                cl_list.append(float(cl))
                cd_list.append(float(cd))
                cm_list.append(float(cm))
            except:
                a_list.append(0.0)
                cl_list.append(0.0)
                cd_list.append(1000.0)
                cm_list.append(1000.0)
        
        i += 1

    i = 0
    clcd_remember = 0
    while ((i < len(a_list)) and (i < len(cl_list)) and (i < len(cd_list)) and (i < len(cm_list))):
        a_try = a_list[i]
        cl_try = cl_list[i]
        cd_try = cd_list[i]
        cm_try = cm_list[i]
        try:
            clcd_try = cl_try/cd_try
        except:
            clcd_try = 0
        if clcd_try > clcd_remember:
            a = a_try
            cl = cl_try
            cd = cd_try
            cm = cm_try
            clcd_remember = cl_try/cd_try
        i += 1
    try:
        ps.kill()
        return a, cl, cd, cm
    except:
        ps.kill()
        return 0, 0, 1000, 1000