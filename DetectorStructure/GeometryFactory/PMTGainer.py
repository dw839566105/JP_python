import random


f1 = open("PMT_Gain.xml","w")
f1.truncate()

f1.write("<extension>\n\n")

#f2 = open("QE_1828.dat","r")

#qe = []

#with open("QE_1828.dat") as f2:
#    f2.readline()
#    i = 0
#    for line in f2:
#        data = line.split()
#        qe.append([float(data[0]), float(data[1])])
#        i = i+1


for i in range(4274):
#for i in range(4590):
#for i in range(30):
    f1.write('<PMT name="PMT_'+str(i+1)+'">\n')
    gainfluc = random.gauss(1,0.1/3)
    f1.write('\t<property name="Gain" value="'+str("%.2f"%(155*gainfluc))+'"/>\n')
    f1.write('\t<property name="DarkNoiseRate" value="2000"/>\n')
    #f1.write('\t<property name="Model" model="HZC_XP1805"/>\n')
    f1.write('\t<property name="Model" model="Hamamatsu_R3600_02"/>\n')
    fluc = random.gauss(1,0.05/3)
    f1.write('\t<property name="QECorrection" value="'+str("%.4f"%fluc)+'"/>\n')
    f1.write('</PMT>\n\n')


f1.write("\n</extension>\n");

f1.close()
