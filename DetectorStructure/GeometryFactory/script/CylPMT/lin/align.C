const double pi=3.14159265;
ofstream ou("Align.txt");
double radius=830, pr=26, covf=1.05,covradius=810;
int nrow=pi*radius/pr/covf, nli=radius/pr/covf;
int nPMT=0;

void align(){
	ou<<"No	x	y	z	Rx	Ry	Rz"<<endl;
	barrel();
	endcap(1);
	endcap(-1);
	cout<<"nPMT="<<nPMT<<endl;
	cout<<"Coverage="<<nPMT*pr*pr/(covradius*covradius*2+2*covradius*2*covradius)<<endl;
}

void endcap(int flag){
//	if(flag>0)ou<<"#top_endcap#"<<endl;
//	else ou<<"#bottom_endcap#"<<endl;
	double x,y;
	double rotx=(flag+1)/2*180, roty=0, rotz=0;
	double z=radius*flag;
	for(int i=0;i<=int(radius/pr/covf);i++){
		for(int j=0;j<=int(radius/pr/covf);j++){
			x=(double(i)/int(radius/pr/covf)*2-1)*radius;
			y=(double(j)/int(radius/pr/covf)*2-1)*radius;
			if((x*x+y*y)>(radius-pr)*(radius-pr))continue;
			nPMT++;
			ou<<nPMT<<"	"<<x/1e2<<" "<<y/1e2<<" "<<z/1e2<<" "<<rotx<<" "<<roty<<" "<<rotz<<endl;
		}
	}
	cout<<"Endcap "<<flag<<" "<<nPMT<<endl;
}

void barrel(){
//	ou<<"#barrel#"<<endl;
	double x,y,z,phi;
	double rotx, roty,rotz=0;
	for(int i=0;i<=int((radius-pr*2.5)/pr/covf);i++){
		z=(double(i)/int((radius-pr*2.5)/pr/covf)*2-1)*(radius-pr*2.5);
		for(int j=0;j<int(pi*radius/pr/covf);j++){
			phi=(double(j)/int(pi*radius/pr/covf)*2-1)*pi;
			x=radius*cos(phi);
			y=radius*sin(phi);
			rotx=pi/2;
			roty=pi/2+phi;
			nPMT++;
			ou<<nPMT<<" "<<x/1e2<<" "<<y/1e2<<" "<<z/1e2<<" "<<rotx*180/pi<<" "<<roty*180/pi<<" "<<rotz*180/pi<<endl;
		}
	}
	cout<<"Barrel "<<nPMT<<endl;
}

