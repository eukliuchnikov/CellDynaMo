import math 
import numpy

with open("KtPol.xyz", "r") as f:
    flines = f.readlines()

with open("ang_dist.dat", "w") as h:
    for i in range(len(flines)):
        if (flines[i].split()[0] == 'KTLC'):
            dx = float(flines[i].split()[1])-float(flines[i + 1].split()[1])
            dy = float(flines[i].split()[2])-float(flines[i + 1].split()[2])
            dz = float(flines[i].split()[3])-float(flines[i + 1].split()[3])
            dist = math.sqrt(dx*dx + dy*dy + dz*dz) + 2*362.5

            x0 = (float(flines[i].split()[1])+float(flines[i + 1].split()[1]))/2
            y0 = (float(flines[i].split()[2])+float(flines[i + 1].split()[2]))/2
            z0 = (float(flines[i].split()[3])+float(flines[i + 1].split()[3]))/2
            
#            dx0 = x0-float(flines[i + 2].split()[1])
#            dy0 = y0-float(flines[i + 2].split()[2])
#            dz0 = z0-float(flines[i + 2].split()[3])
            dist_equat = x0
            if (dist_equat < 0):
                dist_equat = -dist_equat

            dist_spind = math.sqrt(y0*y0 + z0*z0)

            dxK = float(flines[i].split()[1])-float(flines[i + 1].split()[1])
            dyK = float(flines[i].split()[2])-float(flines[i + 1].split()[2])
            dzK = float(flines[i].split()[3])-float(flines[i + 1].split()[3])

            dxP = float(flines[i + 2].split()[1])-float(flines[i + 3].split()[1])
            dyP = float(flines[i + 2].split()[2])-float(flines[i + 3].split()[2])
            dzP = float(flines[i + 2].split()[3])-float(flines[i + 3].split()[3])

            cos = (dxK*dxP + dyK*dyP + dzK*dzP)/(math.sqrt(dxK*dxK + dyK*dyK + dzK*dzK)*math.sqrt(dxP*dxP + dyP*dyP + dzP*dzP))
            angle = math.acos(cos)
            order = math.cos(2*angle)
            angle = angle*180/3.14
            if (angle > 90):
                angle = 180.0 - angle

            h.write(str(dist/1000) + "\t" + str(dist_equat/1000) + "\t" + str(dist_spind/1000) + "\t" + str(angle) + "\t" + str(order) + "\t" + str(dist) + "\n")
