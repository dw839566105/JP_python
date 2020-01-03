vid align(){
	const double pi=3.14159;
	ofstream ou("A.txt");
	double radius=650, pr=26, covf=2;
	int nring=pi*radius/2/pr/covf;
	int PMTid=1;
	double r,x,y,z;
	for(int i=0;i<nring+1;i++){
		r=radius*sin(i*pi/nring);
		z=radius*cos(i*pi/nring);
		cout<<r<<endl;
		for(int j=0;j<int(2*pi*(r+1)/2/pr/covf);j++){
			PMTid++;
			ou<<"<physvol name=\""<<Form("PMT_%03i",PMTid)<<"\">"<<endl;
			ou<<"<file name=\"_20inPMT.gdml\"/>"<<endl;
			cout<<r<<" "<<int(2*pi*(r+1)/2/pr/covf)<<endl;
			x=r*cos(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			y=r*sin(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			ou<<"<position x=\""<<x<<"\" y=\""<<y<<"\" z=\""<<z<<"\" unit=\"cm\"/>"<<endl;
			if(z>0)ou<<"<rotation x=\""<<-pi+asin(y/sqrt(z*z+y*y))<<"\" y=\""<<asin(x/radius)<<"\" unit=\"rad\"/>"<<endl;
			else ou<<"<rotation x=\""<<-asin(y/sqrt(z*z+y*y))<<"\" y=\""<<asin(x/radius)<<"\" unit=\"rad\"/>"<<endl;
			ou<<"</physvol>"<<endl<<endl;
		}
	}
}

