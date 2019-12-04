void Sphere(){
    const double pi=3.14159;

    ofstream ou("PMT_Position.xml");

    // 20 inch PMT largest radius unit :cm  
    double pr =26 ;

    // Radius distacne between  20 inch PMT center  and sphere center unit: cm 
    double radius = 1500 ;

    // Sphere outer PMT cycle radius vs. Radius ratio(Acrylic Sphere tangent)  
    double covf=1.05   ;

    // PMT range rings around Acrylic Sphere 
    //int nring=pi*radius/2/pr/covf;
    int nring=pi*radius/2/pr/covf;

    // PMT ID 
    int PMTid=0;

    double r,x,y,z;
    double r_x,r_y,r_z;

     for(int i=0;i<nring+1;i++){

        // calculate every ring radius 
        r=radius*sin(i*pi/nring);

        // z distance from center 
        z=radius*cos(i*pi/nring);

        for(int j=0;j<int(2*pi*(r+1)/2/pr/covf);j++){

            PMTid++;

            x=r*cos(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);

            y=r*sin(2*pi/int(2*pi*(r+1)/2/pr/covf)*j);

            ou<<"<physvol name = \"PMT_"<< PMTid << "\">"<< endl; 

            ou<<"<file name = \"../DetectorStructure/PMTlib/Hamamatsu_R3600_02/Hamamatsu_R3600_02_Geo.gdml\"/>"<< endl; 

            ou<< "<position x=\"" << x/1e2 << "\" y=\" "<< y/1e2 <<" \" z=\" "<< z/1e2<< "\" unit=\"m\"/>" << endl; 

            if(z>0)
            {
                r_x = (-pi+asin(y/sqrt(z*z+y*y)))/pi*180 ; 
                r_y = (asin(x/radius))/pi*180 ; 
                r_z = 0 ;
                // ou<<(-pi+asin(y/sqrt(z*z+y*y)))/pi*180<<" "<<(asin(x/radius))/pi*180<<" 0"<<endl;
                ou<< "<rotation x=\""<< r_x <<"\" y=\"" << r_y <<"\" z=\""<< "0 \" unit=\"deg\"/>"<< endl; 
            }    
            else
            {
                r_x = (-asin(y/sqrt(z*z+y*y)))/pi*180 ; 
                r_y =   (asin(x/radius))/pi*180 ; 

                ou<< "<rotation x=\""<< r_x <<"\" y=\"" << r_y <<"\" z=\""<< "0 \" unit=\"deg\"/>"<< endl; 
            }
            ou << "</physvol>" << endl ;
        }
    }

    }

