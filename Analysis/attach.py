L1 = 0
L2 = 0
R1 = 0
R2 = 0

with open("att.dat", "r") as f:
    flines = f.readlines()

KtL = int(flines[0].split()[2])
KtR = 0

for i in range(len(flines)):
    if (int(flines[i].split()[2]) != KtL):
        KtR = int(flines[i].split()[2])

with open("MTnum.dat", "w") as h:
    h.write("0.000000" + "\t" + str(L1) + "\t" + str(L2) + "\t" + str(R1) + "\t" + str(R2) + "\t" + str(L1 + L2 + R1 + R2) + "\n")
    for i in range(len(flines)):
        if (int(flines[i].split()[2]) == KtL and flines[i].split()[3] == "ATTACHMENT" and int(flines[i].split()[1]) < 8250):
            L1 += 1
        if (int(flines[i].split()[2]) == KtL and flines[i].split()[3] == "ATTACHMENT" and int(flines[i].split()[1]) >= 8250):
            R1 += 1
        if (int(flines[i].split()[2]) == KtR and flines[i].split()[3] == "ATTACHMENT" and int(flines[i].split()[1]) < 8250):
            L2 += 1
        if (int(flines[i].split()[2]) == KtR and flines[i].split()[3] == "ATTACHMENT" and int(flines[i].split()[1]) >= 8250):
            R2 += 1

        if (int(flines[i].split()[2]) == KtL and flines[i].split()[3] == "DETACHMENT" and int(flines[i].split()[1]) < 8250):
            L1 -= 1
        if (int(flines[i].split()[2]) == KtL and flines[i].split()[3] == "DETACHMENT" and int(flines[i].split()[1]) >= 8250):
            R1 -= 1
        if (int(flines[i].split()[2]) == KtR and flines[i].split()[3] == "DETACHMENT" and int(flines[i].split()[1]) < 8250):
            L2 -= 1
        if (int(flines[i].split()[2]) == KtR and flines[i].split()[3] == "DETACHMENT" and int(flines[i].split()[1]) >= 8250):
            R2 -= 1

        h.write(flines[i].split()[0] + "\t" + str(L1) + "\t" + str(L2) + "\t" + str(R1) + "\t" + str(R2) + "\t" + str(L1 + L2 + R1 + R2) + "\n")

print(L2)

with open("MTnum.dat", "r") as h:
    hlines = h.readlines()

delta = 4.0

with open("STATUS.dat", "w") as g:
# 0 - no att; 1 - mero; 2 - mero-synt; 3 - synt; 4 - mero-amph; 5 - amph
    for i in range(len(hlines)):
        L = int(hlines[i].split()[1]) + int(hlines[i].split()[2])
        R = int(hlines[i].split()[3]) + int(hlines[i].split()[4])
        A1 = int(hlines[i].split()[1]) + int(hlines[i].split()[4])
        A2 = int(hlines[i].split()[2]) + int(hlines[i].split()[3])

        if (L == 0 and R == 0):
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "0\n")
        elif ((L == 0 and int(hlines[i].split()[3]) != 0 and int(hlines[i].split()[4]) != 0) or (R == 0 and int(hlines[i].split()[1]) != 0 and int(hlines[i].split()[2]) != 0)):
            print(R, L1, L2)
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "3\n")
        elif (L == 0 or R == 0):
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "0\n")
        elif (A1 == 0 or A2 == 0):
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "5\n")
        elif (A1/A2 < 1/delta or A1/A2 > delta):
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "4\n")
        elif (L/R < 1/delta or L/R > delta):
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "2\n")
        else:
            g.write(str(float(hlines[i].split()[0])/60.0) + "\t" + "1\n")
    

with open("STATUS.dat", "r") as g:
    glines = g.readlines()

with open("Amph.dat", "w") as k:
    for i in range(len(glines)):
        if (int(glines[i].split()[1]) != 5):
            k.write(glines[i].split()[0] + "\t" + "0\n")
        else:
            k.write(glines[i].split()[0] + "\t" + "1\n")
    
