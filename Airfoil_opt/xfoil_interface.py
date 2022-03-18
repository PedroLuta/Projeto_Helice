from importing import *



def write_to_xfoil(x, y):
    y.stdin.write(x)

def calculate_most_eff_alpha(a_list, cl_list, cd_list):
    # if (len(a_list) != len(cl_list)) or (len(a_list) != len(cd_list)):
    #     return
    if len(a_list) == 0:
        return 0, 0, 1

    clcd_remember = 0
    # converged = False
    for i in range(len(a_list)):
        if cd_list[i] == 0:
            continue
        clcd_try = cl_list[i]/cd_list[i]
        if clcd_try > clcd_remember:
            a, cl, cd = a_list[i], cl_list[i], cd_list[i]
            clcd_remember = clcd_try
            # converged = True
    # if not converged:
    #     return a_list[0], cl_list[0], cd_list[0]
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

    #startupinfo = sp.STARTUPINFO()
    #startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW

    ps = sp.Popen(xfoilpath, universal_newlines = True, stdin = sp.PIPE, stdout = sp.PIPE) #  startupinfo = startupinfo, 

    try:
        output, _ = ps.communicate(input = f'load {afile}\n\nplop\ng f\n\noper\niter {itr}\nvisc {Rey}\nM {M}\naseq {a1} {a2} {astep}\n\n\n\nquit \n', timeout = timeout) 
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

