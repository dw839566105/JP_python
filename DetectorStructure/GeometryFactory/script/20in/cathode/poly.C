void poly(){
	ofstream ne("Newpoly.txt");
	const int num=152;
	double r[num],z[num],dz;
	for(int i=0;i<num;i=i+5){
		z[i]=i;
		dz=z[i];
		if(i<33)r[i]=sqrt(124*124-(124-dz)*(124-dz));
		else if(i<129)r[i]=41+sqrt(59*59-(73.2-dz)*(73.2-dz));
		else r[i]=77-sqrt(28*28-(152-dz)*(152-dz));
		ne<<"<rzpoint r=\""<<2.5*r[i]<<"\" z=\""<<375-2.5*z[i]<<"\"/>"<<endl;
	}
}

