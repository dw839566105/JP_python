const double pi=3.14159265;
ofstream ou("Align.txt");
double radius=830, pr=26, covf=1.5,covradius=810;
int nrow=pi*radius/pr/covf, nli=radius/pr/covf;
int nPMT=0;

void align(){
	endcap(1);
	endcap(-1);
	barrel();
	cout<<"nPMT="<<nPMT<<endl;
	cout<<"Coverage="<<nPMT*pr*pr/(covradius*covradius*2+2*covradius*2*covradius)<<endl;
}

void endcap(int flag){
	if(flag>0)ou<<"<!--top endcap-->"<<endl;
	else ou<<"<!--bottom endcap-->"<<endl;
	double x,y;
	double rotx=(flag+1)/2*180;
	double z=radius*flag;
	double phi;
	for(int i=0;i<=int(radius/pr/covf);i++){
		for(int j=0;j<=int(radius/pr/covf);j++){
			x=(double(i)/int(radius/pr/covf)*2-1)*radius;
			y=(double(j)/int(radius/pr/covf)*2-1)*radius;
			if((x*x+y*y)>(radius-pr)*(radius-pr))continue;
			if(y<0)continue;
			ou<<"<physvol>"<<endl;
			ou<<"<volumeref ref=\"PMTLog\"/>"<<endl;
			ou<<"<position x=\""<<x<<"\" y=\""<<y<<"\" z=\""<<z<<"\" unit=\"cm\"/>"<<endl;
			ou<<"<rotation x=\""<<rotx<<"\" unit=\"deg\"/>"<<endl;
			ou<<"</physvol>"<<endl<<endl;
			nPMT++;
		}
	}
	cout<<endl;
}

void barrel(){
	ou<<"<!--barrel-->"<<endl;
	double x,y,z,phi;
	double rotx, roty;
	for(int i=0;i<=int((radius-pr*2.5)/pr/covf);i++){
		z=(double(i)/int((radius-pr*2.5)/pr/covf)*2-1)*(radius-pr*2.5);
		for(int j=0;j<int(pi*radius/pr/covf);j++){
			phi=(double(j)/int(pi*radius/pr/covf)*2-1)*pi;
			if(phi<0)continue;
			x=radius*cos(phi);
			y=radius*sin(phi);
			rotx=pi/2;
			roty=pi/2+phi;
			ou<<"<physvol>"<<endl;
			ou<<"<volumeref ref=\"PMTLog\"/>"<<endl;
			ou<<"<position x=\""<<x<<"\" y=\""<<y<<"\" z=\""<<z<<"\" unit=\"cm\"/>"<<endl;
			ou<<"<rotation x=\""<<rotx<<"\" y=\""<<roty<<"\" unit=\"rad\"/>"<<endl;
			ou<<"</physvol>"<<endl<<endl;
			nPMT++;
		}
	}
}

