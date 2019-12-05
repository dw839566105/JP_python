void align(){
	const double pi=3.14159;
	ofstream ou("Align.txt");
	double radius=830, pr=26, covf=1.05, covradius=810;
	int nring=pi*radius/2/pr/covf;
	int PMTid=0;
	double r,x,y,z;
	ou<<"No	x	y	z	rx	ry	rz"<<endl;
	for(int i=0;i<nring+1;i++){
		r=radius*sin(i*pi/nring);
		z=radius*cos(i*pi/nring);
		for(int j=0;j<int(2*pi*(r+1)/2/pr/covf);j++){
			PMTid++;
			x=r*cos(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			y=r*sin(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);
			ou<<PMTid<<" "<<x/1e2<<" "<<y/1e2<<" "<<z/1e2<<" ";
			if(z>0)ou<<(-pi+asin(y/sqrt(z*z+y*y)))/pi*180<<" "<<(asin(x/radius))/pi*180<<" 0"<<endl;
			else ou<<(-asin(y/sqrt(z*z+y*y)))/pi*180<<" "<<(asin(x/radius))/pi*180<<" 0"<<endl;
		}
	}
	cout<<PMTid<<endl;
	cout<<PMTid*pr*pr/(4*covradius*covradius)<<endl;
}

