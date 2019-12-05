const double pi=3.14159265;
ofstream ou("PMT_Position.xml");
// 20 inch PMT largest radius unit:cm 
double pr =26 ;

// cyclinder radius unit :cm   
double radius= 800+55;

// cyclinder height unit :cm 
double height = 800+55 ;

// Sphere outer PMT cycle radius vs. Radius ratio(Acrylic Sphere tangent)  
double covf=1.05   ;

//covradius=810;

// nrow : cycle around cyclinder nli: numbers of lines put PMT  
int nrow=pi*radius/pr/covf, nli=radius/pr/covf;

int nPMT=0;


void Cylinder(){
    // cyclinder barrel  function 
     barrel();
    // cyclinder top endcap  function 
    endcap(1);
    // cyclinder bottom endtop function 
    endcap(-1);

    cout << nPMT << endl;
    //	cout<<"Coverage="<<nPMT*pr*pr/(covradius*covradius*2+2*covradius*2*covradius)<<endl;
}

void endcap(int flag)
{

    double x,y;
    double rotx=(flag+1)/2*180, roty=0, rotz=0;
    double z=height*flag;
    for(int i=0;i<=int(radius/pr/covf);i++){
        for(int j=0;j<=int(radius/pr/covf);j++){
            x=(double(i)/int(radius/pr/covf)*2-1)*radius;
            y=(double(j)/int(radius/pr/covf)*2-1)*radius;
            if((x*x+y*y)>(radius-pr)*(radius-pr))continue;
            nPMT++;
            ou<<"<physvol name = \"PMT_"<< nPMT << "\">"<< endl;
            ou<<"<file name = \"../DetectorStructure/PMTlib/Hamamatsu_R3600_02/Hamamatsu_R3600_02_Geo.gdml\"/>"<< endl;
            ou<< "<position x=\"" << x/1e2 << "\" y=\" "<< y/1e2 <<" \" z=\" "<< z/1e2<< "\" unit=\"m\"/>" << endl;    
            ou<< "<rotation x=\""<< rotx <<"\" y=\"" << roty <<"\" z=\""<< "0 \" unit=\"deg\"/>"<< endl;
            ou << "</physvol>" << endl ;

        }
    }
}
void barrel(){
    //	ou<<"#barrel#"<<endl;
    double x,y,z,phi;
    double rotx, roty,rotz=0;
    //for(int i=0;i<=int((radius-pr*2.5)/pr/covf);i++){
    for(int i=0;i<=int((height-pr*2.5)/pr/covf);i++){
        // z=(double(i)/int((radius-pr*2.5)/pr/covf)*2-1)*(radius-pr*2.5);
        z=(double(i)/int((height-pr*2.5)/pr/covf)*2-1)*(height-pr*2.5);
        for(int j=0;j<int(pi*radius/pr/covf);j++){
            phi=(double(j)/int(pi*radius/pr/covf)*2-1)*pi;
            x=radius*cos(phi);
            y=radius*sin(phi);
            rotx=pi/2;
            roty=pi/2+phi;
            nPMT++;
            ou<<"<physvol name = \"PMT_"<< nPMT << "\">"<< endl;
            ou<<"<file name = \"../DetectorStructure/PMTlib/Hamamatsu_R3600_02/Hamamatsu_R3600_02_Geo.gdml\"/>"<< endl;
            ou<< "<position x=\"" << x/1e2 << "\" y=\" "<< y/1e2 <<" \" z=\" "<< z/1e2<< "\" unit=\"m\"/>" << endl; 
            ou<< "<rotation x=\""<< rotx*180/pi <<"\" y=\"" << roty*180/pi <<"\" z=\""<< rotz*180/pi <<"\"  unit=\"deg\"/>"<< endl; 
            ou << "</physvol>" << endl ;
        }
    }
}
