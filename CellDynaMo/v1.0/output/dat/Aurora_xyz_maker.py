import numpy as np

#for i in range(1000):
#    s = np.random.uniform(0,1)
#    print(s) 

with open("aurora.dat", "r") as f:
    flines = f.readlines()
num = 0



with open("aurora.dat", "r") as f:
    for i in range(len(flines)):
        if (flines[i] == '*******************START*******************\n'):
            num += 1
            print(num)

sys = [0]*num
start = [0]*num

start[0] = 0

frame_id = 0
line_id = 0

bad_words = ['.']

sumA = 0

for i in range(len(flines)):
    line_id += 1
    if (flines[i] != '*******************START*******************\n' and flines[i] != '********************END********************\n' and flines[i] != 'x	y	z	cell_type	Aurora_B	Aurora_A\n'):
        if not any(bad_word in flines[i] for bad_word in bad_words):
            sumA += int(flines[i].split()[4])    
    
    if (flines[i] == '********************END********************\n'):
        print(sumA)
        sumA = 0
        sys[frame_id] = line_id
        if (frame_id < num - 1):
            start[frame_id + 1] = i
        frame_id += 1
        line_id = 0

print(sys)

NB = 3505*3
#NA = 1025680*2
NA = 0

pNB = NB
pNA = NA
Svsize = 250

XBOX = 2*8000 + Svsize
YBOX = 2*5000 + Svsize
ZBOX = 2*5000 + Svsize

mdN = NB + NA

x_coord_B = [0.0]*NB
y_coord_B = [0.0]*NB
z_coord_B = [0.0]*NB

x_coord_A = [0.0]*NA
y_coord_A = [0.0]*NA
z_coord_A = [0.0]*NA

new_i_B = 0
new_i_A = 0

#sys = 531445

with open("aurora.xyz", "w") as h:
    for k in range(num):
        for i in range(sys[k]):
            if (flines[start[k] + i] != '*******************START*******************\n' and flines[start[k] + i] != '********************END********************\n' and flines[start[k] + i] != 'x	y	z	cell_type	Aurora_B	Aurora_A\n'):
                if not any(bad_word in flines[start[k] + i] for bad_word in bad_words):
                    if (int(flines[start[k] + i].split()[4]) != 0):
                        for j in range(int(flines[start[k] + i].split()[4])):
                            print (new_i_B)
                            x_coord_B[new_i_B] = Svsize*np.random.uniform(0,1) + Svsize*(float(flines[start[k] + i].split()[0]) - 1.0) - XBOX/2.0
                            y_coord_B[new_i_B] = Svsize*np.random.uniform(0,1) + Svsize*(float(flines[start[k] + i].split()[1]) - 1.0) - YBOX/2.0
                            z_coord_B[new_i_B] = Svsize*np.random.uniform(0,1) + Svsize*(float(flines[start[k] + i].split()[2]) - 1.0) - ZBOX/2.0  
                            new_i_B += 1
                    if (int(flines[start[k] + i].split()[5]) != 0):
                        for j in range(int(flines[start[k] + i].split()[5])):
                            x_coord_A[new_i_A] = Svsize*np.random.uniform(0,1) + Svsize*(float(flines[start[k] + i].split()[0]) - 1.0) - XBOX/2.0
                            y_coord_A[new_i_A] = Svsize*np.random.uniform(0,1) + Svsize*(float(flines[start[k] + i].split()[1]) - 1.0) - YBOX/2.0
                            z_coord_A[new_i_A] = Svsize*np.random.uniform(0,1) + Svsize*(float(flines[start[k] + i].split()[2]) - 1.0) - ZBOX/2.0  
                            new_i_A += 1

            if (flines[start[k] + i] == '********************END********************\n'):
                new_i_B = 0
                new_i_A = 0

                h.write(str(int(mdN)) + "\n" + "aurora visualization\n")
                for j in range(NB):
                    h.write("AB\t" + str(float(x_coord_B[j])) + "\t" + str(float(y_coord_B[j])) + "\t" + str(float(z_coord_B[j])) + "\n")
                for j in range(NA):
                    h.write("AA\t" + str(float(x_coord_A[j])) + "\t" + str(float(y_coord_A[j])) + "\t" + str(float(z_coord_A[j])) + "\n")
                print("frame " + str(int(k)) + " was writen")
       
