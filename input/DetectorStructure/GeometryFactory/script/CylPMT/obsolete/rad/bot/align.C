void align(){
	const double pi=3.14159;
	ofstream ou("Align.txt");
	double radius=830, pr=26, covf=2.0;
	int nring=radius/2/pr/covf;
	double r,x,y;
	for(int i=0;i<nring+1;i++){
		r=i*pr*covf*2;
		cout<<r<<endl;
		for(int j=0;j<int(2*pi*(r+1)/2/pr/covf)+1e-8;j++){
			ou<<"<physvol>"<<endl<<"<volumeref ref=\"PMTLog\"/>"<<endl;
			x=r*cos(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			y=r*sin(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			if(r==0){
				x=0;
				y=0;
			}
			ou<<"<position x=\""<<x<<"\" y=\""<<y<<"\" z=\""<<-1680./2<<"\" unit=\"cm\"/>"<<endl;
			ou<<"</physvol>"<<endl<<endl;
		}
	}
}

