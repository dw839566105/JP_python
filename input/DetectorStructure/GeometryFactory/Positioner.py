f1 = open("PMT_Position.xml","w");
f1.truncate()

fread = open("PMT_pos.dat", "r")
slines = fread.readlines()
del slines[0]
for a in slines:
   data = a.split()

   PMT_No = data[0]
   if PMT_No == '#': 
       f1.write('\r\n<!--')
       for b in data:
           f1.write(b+' ')
       f1.write('-->\r\n')
   else :
       x = data[1]
       y = data[2]
       z = data[3]
       Rx = data[4]
       Ry = data[5]
       Rz = data[6]
   
       f1.write('<physvol name="PMT_'+PMT_No+'">\r\n')
       f1.write('\t<file name="_20inPMT.gdml"/>\r\n')
       f1.write('\t<position x="'+x+'" y="'+y+'" z="'+z+'" unit="m"/>\r\n')
       f1.write('\t<rotation x="'+Rx+'" y="'+Ry+'" z="'+Rz+'" unit="deg"/>\r\n')
       f1.write('</physvol>\r\n\r\n')

fread.close()



f1.close()

