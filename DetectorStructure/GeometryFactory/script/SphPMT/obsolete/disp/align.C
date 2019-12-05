void align(){
	const double pi=3.14159;
	ofstream ou("Align.txt");
	double radius=830, pr=26, covf=1.5, covradius=810;
	int nring=pi*radius/2/pr/covf;
	int PMTid=1;
	double r,x,y,z;
	for(int i=0;i<nring+1;i++){
		r=radius*sin(i*pi/nring);
		z=radius*cos(i*pi/nring);
		for(int j=0;j<int(2*pi*(r+1)/2/pr/covf);j++){
			PMTid++;
			x=r*cos(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			y=r*sin(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			if(y<0)continue;
			ou<<"<physvol>"<<endl<<"<volumeref ref=\"PMTLog\"/>"<<endl;
			ou<<"<position x=\""<<x<<"\" y=\""<<y<<"\" z=\""<<z<<"\" unit=\"cm\"/>"<<endl;
			if(z>0)ou<<"<rotation x=\""<<-pi+asin(y/sqrt(z*z+y*y))<<"\" y=\""<<asin(x/radius)<<"\" unit=\"rad\"/>"<<endl;
			else ou<<"<rotation x=\""<<-asin(y/sqrt(z*z+y*y))<<"\" y=\""<<asin(x/radius)<<"\" unit=\"rad\"/>"<<endl;
			ou<<"</physvol>"<<endl<<endl;
		}
	}
	cout<<PMTid<<endl;
	cout<<PMTid*pr*pr/(4*covradius*covradius)<<endl;
}

