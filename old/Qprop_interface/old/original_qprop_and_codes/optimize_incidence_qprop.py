
def read_propfile(propfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\qprop\propfile.txt'):
    with open(propfile, 'r') as inp:
        i = 1
        r_vec = []
        c_vec = []
        B_vec = []
        while True:
            content = inp.readline()
            if content == "\n":
                continue
            if content == '':
                break
            split = content.replace("\n", '').split()
            if i == 1:
                header = content.replace("\n", '')
            elif i == 2:
                NBlades = int(split[0])
                R = float(split[1])
            elif i == 3:
                CL0 = float(split[0])
                CL_a = float(split[1])
            elif i == 4:
                CLmin = float(split[0])
                CLmax = float(split[1])
            elif i == 5:
                CD0 = float(split[0])
                CD2u = float(split[1])
                CD2l = float(split[2])
                CLCD0 = float(split[3])
            elif i == 6:
                REref = float(split[0])
                REexp = float(split[1])
            elif i == 7:
                Rfac = float(split[0])
                Cfac = float(split[1])
                Bfac = float(split[2])
            elif i == 8:
                Radd = float(split[0])
                Cadd = float(split[1])
                Badd = float(split[2])
            elif i == 9:
                pass
            else:
                r_vec.append(float(split[0]))
                c_vec.append(float(split[1]))
                B_vec.append(float(split[2]))
            i += 1
    return r_vec, c_vec, B_vec, [header, NBlades, R, CL0, CL_a, CLmin, CLmax, CD0, CD2u, CD2l, CLCD0, REref, REexp, Rfac, Cfac, Bfac, Radd, Cadd, Badd]

def read_motorfile(motorfile = r'C:\Users\PEDRO\Desktop\IC_DE_HÉLICE\git_control\git_prop\qprop\motorfile.txt', motor_type = 1):
    with open(motorfile, 'r') as inp:
        i = 1
        while True:
            content = inp.readline()
            if content == "\n":
                continue
            if content == '':
                break
            split = content.replace("\n", '').split()
            
            if motor_type == 1:
                if i == 1:
                    header = content.replace("\n", '')
                elif i == 2:
                    pass
                if i == 3:
                    Res = float(split[0])
                elif i == 4:
                    NLcurr = float(split[0])
                elif i == 5:
                    KV = float(split[0])
            elif motor_type == 2:
                pass
            else:
                pass
            i += 1

    if motor_type == 1:
        return header, Res, NLcurr, KV
    elif motor_type == 2:
        return 0
    else:
        return 0

r_vec, c_vec, B_vec, vector = read_propfile()
header, Res, NLcurr, KV = read_motorfile()

print(header)
print(Res)
print(NLcurr)
print(KV)
