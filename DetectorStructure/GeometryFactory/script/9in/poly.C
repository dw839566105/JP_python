void poly(){
	ifstream po("PMTpoly.txt");
	ofstream ne("Newpoly.txt");
	const int num=44;
	double r[num],z[num];
	for(int i=0;i<num;i++){
		po>>r[i]>>z[i];
		ne<<"<rzpoint r=\""<<2*r[i]<<"\" z=\""<<2*z[i]<<"\"/>"<<endl;
	}
}

