void poly(){
	ofstream ne("Newpoly.txt");
	const int num=152;
	double r[num],z[num],dz;
	for(int i=0;i<num;i=i+5){
		z[i]=i-2;
		dz=z[i];
		if(i<33)r[i]=sqrt(126*126-(124-dz)*(124-dz));
		else if(i<129)r[i]=43+sqrt(60*60-(73.2-dz)*(73.2-dz));
		else r[i]=79-sqrt(28*28-(153-dz)*(153-dz));
		ne<<"<rzpoint r=\""<<2.5*r[i]<<"\" z=\""<<375-2.5*z[i]<<"\"/>"<<endl;
		cout<<z[i]<<" "<<r[i]<<endl;
	}
}

