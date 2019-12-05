void poly(){
	ofstream ne("Newpoly.txt");
	const int num=200;
	double z[num],r[num];
	for(int i=0;i<num;i=i+10){
		z[i]=i;
		if(i<30)r[i]=127;
		else if(i<117)r[i]=sqrt(86*86-(z[i]-30)*(z[i]-30))+41;
		else if(i<127)r[i]=(r[i-10]+41*2)/3;
		else if(i<157)r[i]=41;
		else r[i]=sqrt(41*41-(z[i]-157)*(z[i]-157));
		ne<<"<rzpoint r=\""<<r[i]<<"\" z=\""<<-z[i]-50<<"\"/>"<<endl;
	}
}
