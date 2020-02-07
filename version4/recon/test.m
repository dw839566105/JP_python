E_total = [];

for i = 1:1:100
    A = h5read(['Recon',num2str(i),'.h5'],'/Recon');
    E_total = [E_total; A.E_sph];
   
end
 hist(E_total,100)