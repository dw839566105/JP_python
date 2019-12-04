// Script to calculate the time constant from the measurement result

Double_t xmax = 100;

Double_t Waveform(Double_t *x,Double_t *par) 
{
	Double_t val;
	Double_t c = (2/TMath::Sqrt(2*TMath::Pi()))/((par[0]/3)+(par[1]/3));
	if(x[0]<0)
	  val = 0;
	if(x[0]>0 && x[0]<par[0])
	  val = c*TMath::Gaus(x[0], par[0], par[0]/3, kFALSE);
	else
	  val = c*TMath::Gaus(x[0], par[0], par[1]/3, kFALSE);

	return val;
}

Double_t Waveform1(Double_t *x,Double_t *par) 
{
	Double_t val;
	Double_t c = (par[0]+par[1])/(par[1]*par[1])*1.6E-19*100*6E6*1E12;
	if(x[0]<0)
	  val = 0;
	else
	  val = c*TMath::Exp(-x[0]/par[1])*(1-TMath::Exp(-x[0]/par[0]));

	return val;
}


void GetMsmFromConst(Double_t rConst, Double_t dConst, Double_t& rMsmt, Double_t& dMsmt)
{
	TF1* f1 = new TF1("twoexp","([0]+[1])/([1]*[1])*TMath::Exp(-x/[1])*(1-TMath::Exp(-x/[0]))", 0, xmax);
	f1->SetNpx(1000);
	f1->SetParameters(rConst, dConst);
	Double_t xOfmax = f1->GetMaximumX(0,xmax);
	Double_t ymax = f1->GetMaximum(0,xmax);
	Double_t tt1 = f1->GetX(0.1*ymax,0,xOfmax);
	Double_t tt2 = f1->GetX(0.9*ymax,0,xOfmax);
	Double_t tt3 = f1->GetX(0.9*ymax,xOfmax,xmax);
	Double_t tt4 = f1->GetX(0.1*ymax,xOfmax,xmax);
	rMsmt = tt2-tt1;
	dMsmt = tt4-tt3;


}

void Plot(Double_t rConst, Double_t dConst)
{
	/*TF1* f1 = new TF1("twoexp","([0]+[1])/([1]*[1])*TMath::Exp(-x/[1])*(1-TMath::Exp(-x/[0]))", 0, xmax);
	f1->SetNpx(1000);
	f1->SetParameters(r, d);
	Double_t xOfmax = f1->GetMaximumX(0,xmax);
	Double_t ymax = f1->GetMaximum(0,xmax);
	Double_t tt1 = f1->GetX(0.1*ymax,0,xOfmax);
	Double_t tt2 = f1->GetX(0.9*ymax,0,xOfmax);
	Double_t tt3 = f1->GetX(0.9*ymax,xOfmax,xmax);
	Double_t tt4 = f1->GetX(0.1*ymax,xOfmax,xmax);
	TH2F* h1 = new TH2F("wave", "Waveform", 40, 0, xmax, 40, 0, 1.1*ymax);
	h1->SetStats(kFALSE);
	h1->Draw();
	f1->Draw("same");

	TLine* a = new TLine(tt4,0,tt4,1.1*ymax);
	a->SetLineStyle(7);
	a->Draw();
	TLine* a1 = new TLine(tt1,0,tt1,1.1*ymax);
	a1->SetLineStyle(7);
	a1->Draw();
	TLine* a11 = new TLine(tt2,0,tt2,1.1*ymax);
	a11->SetLineStyle(7);
	a11->Draw();
	TLine* a12 = new TLine(tt3,0,tt3,1.1*ymax);
	a12->SetLineStyle(7);
	a12->Draw();
	TLine* a2 = new TLine(0,0.1*ymax,xmax,0.1*ymax);
	a2->SetLineStyle(7);
	a2->Draw();
	TLine* a3 = new TLine(0,0.9*ymax,xmax,0.9*ymax);
	a3->SetLineStyle(7);
	a3->Draw();

	cout<<"Rising time constant: "<<r<<", Rising time measurement: "<<tt2-tt1<<endl;
	cout<<"Decay time constant: "<<d<<", Decay time measurement: "<<tt4-tt3<<endl;*/
	Double_t rMsmt = 0;
	Double_t dMsmt = 0;
	GetMsmFromConst(rConst, dConst, rMsmt, dMsmt);
	//TF1* f1 = new TF1("twoexp","([0]+[1])/([1]*[1])*TMath::Exp(-x/[1])*(1-TMath::Exp(-x/[0]))", 0, xmax);
	TF1 *f1 = new TF1("fit",Waveform1,0,xmax,2);
	//TF1* f1 = new TF1("twoexp","TMath::Exp(-x/[0])-TMath::Exp(-x/[1])", 0, xmax);
	f1->SetNpx(1000);
	f1->SetParameters(rConst, dConst);
	Double_t ymax = f1->GetMaximum(0,xmax);
	//TH2F* h1 = new TH2F("wave", "Waveform", 40, 0, xmax, 40, 0, 1.1*ymax);
	//h1->SetStats(kFALSE);
	//h1->Draw();
	f1->Draw();
	Double_t tt1 = f1->GetX(0.5*ymax,0,rConst);
	Double_t tt2 = f1->GetX(0.5*ymax,rConst,xmax);
	cout<<"half: "<<tt2-tt1<<endl;


}

void SearchTimeConstant()
{
	Double_t rMsmt = 10;
	Double_t dMsmt = 17;

	//a,b = GetFromConst();
	//minimum (a-rMsmt)^2+(b-dMsmt)^2


	Plot(10, 15);


}
