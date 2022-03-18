from importing import *



def write_to_xfoil(x, y):
    y.stdin.write(x)

def calculate_most_eff_alpha(a_list, cl_list, cd_list):
    if (len(a_list) != len(cl_list)) or (len(a_list) != len(cd_list)):
        return
    length = len(a_list)
    if length == 0:
        return 0, 0, 1
    clcd_remember = 0
    converged = False
    for i in range(length):
        clcd_try = cl_list[i]/cd_list[i]
        if clcd_try > clcd_remember:
            a, cl, cd = a_list[i], cl_list[i], cd_list[i]
            clcd_remember = clcd_try
            converged = True
    if not converged:
        return a_list[0], cl_list[0], cd_list[0]
    return a, cl, cd

def get_curve_com_default(Rey, a1, a2, astep, afile = 'airfoils\\airfoil.txt', timeout = 5, \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', itr = 250, M = 0):
    
    if a1 < 0 and a2 > 0:
        a_list1, cl_list1, cd_list1 = get_curve_com_default(Rey, 0, a1, -astep, afile = afile, timeout = timeout, xfoilpath = xfoilpath, itr = itr, M = M)
        a_list2, cl_list2, cd_list2 = get_curve_com_default(Rey, 0, a2, astep, afile = afile, timeout = timeout, xfoilpath = xfoilpath, itr = itr, M = M)
        if len(a_list1) > 0:
            a_list1.pop(0)
            cl_list1.pop(0)
            cd_list1.pop(0)
        a_list1.reverse()
        cl_list1.reverse()
        cd_list1.reverse()
        a_list1.extend(a_list2)
        cl_list1.extend(cl_list2)
        cd_list1.extend(cd_list2)
        return a_list1, cl_list1, cd_list1

    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\noper\niter {itr}\nvisc {Rey}\nM {M}\naseq {a1} {a2} {astep}\n\n\n\nquit \n', timeout = timeout) 
    except: 
        ps.kill()
        return [], [], []
    list = output.split()
    i = 0
    id_remember = 0
    cl_list = []
    cd_list = []
    a_list = []
    length = len(list)
    for i in range(length):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
        if list[i].lower() == "a":
            a = list[i + 2]
        
        if id_remember < itr:
            skip = False

        if list[i].lower() == "rms:":
            if (id_remember >= itr) and (not skip):
                skip = True
            elif (int(list[i - 1]) < id_remember) and (not skip):
                try:
                    a_list.append(float(a))
                    cl_list.append(float(cl))
                    cd_list.append(float(cd))
                except:
                    pass
                skip = True
            id_remember = int(list[i - 1])
        
        if (i == len(list) - 1) and (id_remember < itr):
            try:
                a_list.append(float(a))
                cl_list.append(float(cl))
                cd_list.append(float(cd))
            except:
                pass
    try:
        ps.kill()
        return a_list, cl_list, cd_list
    except:
        ps.kill()
        return [], [], []

def get_curve_file_default(Rey, a1, a2, astep, afile = 'airfoils\\airfoil.txt', \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', \
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_output.txt', itr = 250, timeout = 5, M = 0):
    
    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, universal_newlines = True, startupinfo = startupinfo, stdout = stout) # 

    write_to_xfoil(f'load {afile}\n\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {int(round(Rey))}\nM {M}\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'aseq {a1} {a2} {astep}\n', ps)
    write_to_xfoil(f'\n\nquit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = timeout)
    except:
        ps.kill()
        return [], [], []

    read_counter = 0
    Cl_dist = []
    Cd_dist = []
    a_dist = []
    try:
        with open(xfoil_outfile, 'r') as o:
            while True:
                if read_counter < 12:
                    o.readline()
                    read_counter += 1
                else:
                    content = o.readline()
                    if content == '':
                        break
                    cont = content.split()
                    a_dist.append(float(cont[0]))
                    Cl_dist.append(float(cont[1]))
                    Cd_dist.append(float(cont[2]))
        try:
            os.remove(xfoil_outfile)
        except:
            pass
        try:
            ps.kill()
            return a_dist, Cl_dist, Cd_dist
        except UnboundLocalError:
            ps.kill()
            return [], [], []

    except FileNotFoundError:
        ps.kill()
        return [], [], []










#returns properties for a fixed airfoil at a specified alpha (file)
def get_properties_fixed(Rey, alpha, afile = 'airfoils\\airfoil.txt', \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', \
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_output.txt', itr = 250, M = 0, \
        N_panels = 90, P_bunch = 2, Te_Le_ratio = 0.05, Refined_Le_ratio = 0.25, \
            top_refined_init = 0.01, top_refined_end = 0.8, \
                bot_refined_init = 1, bot_refined_end = 1):

    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, startupinfo=startupinfo, universal_newlines = True, stdout = stout)

    write_to_xfoil(f'load {afile}\n\nN {N_panels}\nP {P_bunch}\nT {Te_Le_ratio}\nXT {top_refined_init} {top_refined_end}\nXB {bot_refined_init} {bot_refined_end}\nR {Refined_Le_ratio}\n\n\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {Rey}\nM {M}\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'a {alpha}\n', ps)
    write_to_xfoil(f'\n\n\n quit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = 5)
    except:
        ps.kill()
        return 0, 0

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
            return Cl, Cd
        except UnboundLocalError:
            ps.kill()
            return 0, 0

    except FileNotFoundError:
        ps.kill()
        return 0, 0

#returns cl, cd, cm of a specified airfoil at a specified alpha (communicate)
def communicate_fixed(Rey, alpha, afile = 'airfoils\\airfoil.txt', \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', itr = 250, M = 0, \
        N_panels = 90, P_bunch = 2, Te_Le_ratio = 0.05, Refined_Le_ratio = 0.25, \
            top_refined_init = 0.01, top_refined_end = 0.8, \
                bot_refined_init = 1, bot_refined_end = 1):

    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\nN {N_panels}\nP {P_bunch}\nT {Te_Le_ratio}\nXT {top_refined_init} {top_refined_end}\nXB {bot_refined_init} {bot_refined_end}\nR {Refined_Le_ratio}\n\n\noper\niter {itr}\nvisc {Rey}\nM {M}\na {alpha}\n\n\n\nquit \n', timeout = 25) 
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

def get_curve_comMethod_ppar(Rey, a1, a2, astep, afile = 'airfoils\\airfoil.txt', timeout = 5, \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', itr = 250, M = 0, \
        N_panels = 90, P_bunch = 2, Te_Le_ratio = 0.05, Refined_Le_ratio = 0.25, \
            top_refined_init = 0.01, top_refined_end = 0.8, \
                bot_refined_init = 1, bot_refined_end = 1):
    
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\nppar\nN {N_panels}\nP {P_bunch}\nT {Te_Le_ratio}\nXT {top_refined_init} {top_refined_end}\nXB {bot_refined_init} {bot_refined_end}\nR {Refined_Le_ratio}\n\n\noper\niter {itr}\nvisc {Rey}\nM {M}\naseq {a1} {a2} {astep}\n\n\n\nquit \n', timeout = timeout) 
    except: 
        ps.kill()
        return [], [], []
    list = output.split()
    i = 0
    id_remember = 0
    cl_list = []
    cd_list = []
    a_list = []
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
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
                except:
                    a_list.append(0.0)
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                skip = True
            elif (int(list[i - 1]) < id_remember) and (not skip):
                try:
                    a_list.append(float(a))
                    cl_list.append(float(cl))
                    cd_list.append(float(cd))
                except:
                    a_list.append(0.0)
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                skip = True
            id_remember = int(list[i - 1])
        
        if (i == len(list) - 1) and (id_remember < itr):
            try:
                a_list.append(float(a))
                cl_list.append(float(cl))
                cd_list.append(float(cd))
            except:
                a_list.append(0.0)
                cl_list.append(0.0)
                cd_list.append(1000.0)
        i += 1
    try:
        ps.kill()
        return a_list, cl_list, cd_list
    except:
        ps.kill()
        return [], [], []

def get_curve_fileMethod_ppar(Rey, a1 = 0, a2 = 10, astep = 1, afile = 'airfoils\\airfoil.txt', \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', \
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_output.txt', itr = 250, timeout = 5, M = 0, \
        N_panels = 90, P_bunch = 2, Te_Le_ratio = 0.05, Refined_Le_ratio = 0.25, \
            top_refined_init = 0.01, top_refined_end = 0.8, \
                bot_refined_init = 1, bot_refined_end = 1):
    
    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, universal_newlines = True, startupinfo = startupinfo, stdout = stout) # 

    write_to_xfoil(f'load {afile}\n\n', ps)
    write_to_xfoil(f'ppar\n', ps)
    write_to_xfoil(f'N {N_panels}\n', ps)
    write_to_xfoil(f'P {P_bunch}\n', ps)
    write_to_xfoil(f'T {Te_Le_ratio}\n', ps)
    write_to_xfoil(f'XT {top_refined_init} {top_refined_end}\n', ps)
    write_to_xfoil(f'XB {bot_refined_init} {bot_refined_end}\n', ps)
    write_to_xfoil(f'R {Refined_Le_ratio}\n\n\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {int(round(Rey))}\nM {M}\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'aseq {a1} {a2} {astep}\n', ps)
    write_to_xfoil(f'\n\nquit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = timeout)
    except:
        ps.kill()
        return [], [], []

    read_counter = 0
    Cl_dist = []
    Cd_dist = []
    a_dist = []
    try:
        with open(xfoil_outfile, 'r') as o:
            while True:
                if read_counter < 12:
                    o.readline()
                    read_counter += 1
                else:
                    content = o.readline()
                    if content == '':
                        break
                    cont = content.split()
                    a_dist.append(float(cont[0]))
                    Cl_dist.append(float(cont[1]))
                    Cd_dist.append(float(cont[2]))
        try:
            os.remove(xfoil_outfile)
        except:
            pass
        try:
            ps.kill()
            return a_dist, Cl_dist, Cd_dist
        except UnboundLocalError:
            ps.kill()
            return [], [], []

    except FileNotFoundError:
        ps.kill()
        return [], [], []

def get_curve_comMethod_pcop(Rey, a1, a2, astep, afile = 'airfoils\\airfoil.txt', timeout = 5, \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', itr = 250, M = 0):
    
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, startupinfo = startupinfo, stdout = sp.PIPE) # 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\npcop\noper\niter {itr}\nvisc {Rey}\nM {M}\naseq {a1} {a2} {astep}\n\n\n\nquit \n', timeout = timeout) 
    except: 
        ps.kill()
        return [], [], []
    list = output.split()
    i = 0
    id_remember = 0
    cl_list = []
    cd_list = []
    a_list = []
    while i < len(list):
        if list[i].lower() == "cl":
            cl = list[i + 2]
        if list[i].lower() == "cd":
            cd = list[i + 2]
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
                except:
                    a_list.append(0.0)
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                skip = True
            elif (int(list[i - 1]) < id_remember) and (not skip):
                try:
                    a_list.append(float(a))
                    cl_list.append(float(cl))
                    cd_list.append(float(cd))
                except:
                    a_list.append(0.0)
                    cl_list.append(0.0)
                    cd_list.append(1000.0)
                skip = True
            id_remember = int(list[i - 1])
        
        if (i == len(list) - 1) and (id_remember < itr):
            try:
                a_list.append(float(a))
                cl_list.append(float(cl))
                cd_list.append(float(cd))
            except:
                a_list.append(0.0)
                cl_list.append(0.0)
                cd_list.append(1000.0)
        i += 1
    try:
        ps.kill()
        return a_list, cl_list, cd_list
    except:
        ps.kill()
        return [], [], []

def get_curve_fileMethod_pcop(Rey, a1 = 0, a2 = 10, astep = 1, afile = 'airfoils\\airfoil.txt', \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', \
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_output.txt', itr = 250, timeout = 5, M = 0):
    
    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, universal_newlines = True, startupinfo = startupinfo, stdout = stout) # 

    write_to_xfoil(f'load {afile}\n\n', ps)
    #write_to_xfoil(f'mdes\n', ps) ###################
    #write_to_xfoil(f'smoo\nsmoo\nsmoo\nsmoo\nsmoo\n\n', ps) ###################
    write_to_xfoil(f'pcop\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {int(round(Rey))}\nM {M}\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'aseq {a1} {a2} {astep}\n', ps)
    write_to_xfoil(f'\n\nquit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = timeout)
    except:
        ps.kill()
        return [], [], []

    read_counter = 0
    Cl_dist = []
    Cd_dist = []
    a_dist = []
    try:
        with open(xfoil_outfile, 'r') as o:
            while True:
                if read_counter < 12:
                    o.readline()
                    read_counter += 1
                else:
                    content = o.readline()
                    if content == '':
                        break
                    cont = content.split()
                    a_dist.append(float(cont[0]))
                    Cl_dist.append(float(cont[1]))
                    Cd_dist.append(float(cont[2]))
        try:
            os.remove(xfoil_outfile)
        except:
            pass
        try:
            ps.kill()
            return a_dist, Cl_dist, Cd_dist
        except UnboundLocalError:
            ps.kill()
            return [], [], []

    except FileNotFoundError:
        ps.kill()
        return [], [], []

def get_curve_fileMethod_N(Rey, a1 = 0, a2 = 10, astep = 1, afile = 'airfoils\\airfoil.txt', \
    xfoilpath = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil.exe', \
    xfoil_outfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_prop\xfoil_output.txt', itr = 250, timeout = 5, M = 0, N = 160):
    
    stout = 0
    startupinfo = sp.STARTUPINFO()
    startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    try:
        os.remove(xfoil_outfile)
    except:
        pass

    ps = sp.Popen(xfoilpath, stdin = sp.PIPE, universal_newlines = True, startupinfo = startupinfo, stdout = stout) # 

    write_to_xfoil(f'load {afile}\n\n', ps)
    write_to_xfoil(f'ppar\nN {N}\n\n\n', ps)
    write_to_xfoil('oper\n', ps)
    write_to_xfoil(f'iter {itr}\n', ps)
    write_to_xfoil(f'visc {int(round(Rey))}\nM {M}\n', ps)
    write_to_xfoil('pacc\n', ps)
    write_to_xfoil('xfoil_output.txt\n\n', ps)
    write_to_xfoil(f'aseq {a1} {a2} {astep}\n', ps)
    write_to_xfoil(f'\n\nquit\n', ps)

    ps.stdin.close()
    try:
        ps.wait(timeout = timeout)
    except:
        ps.kill()
        return [], [], []

    read_counter = 0
    Cl_dist = []
    Cd_dist = []
    a_dist = []
    try:
        with open(xfoil_outfile, 'r') as o:
            while True:
                if read_counter < 12:
                    o.readline()
                    read_counter += 1
                else:
                    content = o.readline()
                    if content == '':
                        break
                    cont = content.split()
                    a_dist.append(float(cont[0]))
                    Cl_dist.append(float(cont[1]))
                    Cd_dist.append(float(cont[2]))
        try:
            os.remove(xfoil_outfile)
        except:
            pass
        try:
            ps.kill()
            return a_dist, Cl_dist, Cd_dist
        except UnboundLocalError:
            ps.kill()
            return [], [], []

    except FileNotFoundError:
        ps.kill()
        return [], [], []

