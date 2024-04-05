with open("new.log", "r") as f:
    flines = f.readlines()
with open("force1.dat", "w") as h:
    for i in range(len(flines)):
        if (flines[i] == "KT#1 FORCES:	N_total	N_pull	N_push	F_pull	F_push\n"):
            h.write(flines[i + 1])
with open("force0.dat", "w") as g:
    for i in range(len(flines)):
        if (flines[i] == "KT#0 FORCES:	N_total	N_pull	N_push	F_pull	F_push\n"):
            g.write(flines[i + 1])
