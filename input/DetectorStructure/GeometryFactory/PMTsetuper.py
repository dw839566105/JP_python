import random


f1 = open("PMT_Setup.xml","w")
f1.truncate()

f1.write("<extension>\n\n")

f2 = open("QE_1828.dat","r")

qe = []

with open("QE_1828.dat") as f2:
    f2.readline()
    i = 0
    for line in f2:
        data = line.split()
        qe.append([float(data[0]), float(data[1])])
        i = i+1


for i in range(30):
    f1.write('<PMT name="PMT_'+str(i+1)+'">\n')
    f1.write('\t<property name="Gain" value="1700"/>\n')
    f1.write('\t<property name="CollectionEff" value="0.8"/>\n')
    f1.write('\t<property name="DarkNoiseRate" value="10000"/>\n')
    f1.write('\t<property name="QE" QEList="\n')
    fluc = random.gauss(1,0.05/3)
    for item in qe:
        item[1] = item[1]*fluc
        f1.write('\t\t'+str(item[0])+'\t'+str(item[1])+'\n')
    f1.write('\t"/>\n')
    f1.write('</PMT>\n\n')


f1.write("\n</extension>\n");

f1.close()
