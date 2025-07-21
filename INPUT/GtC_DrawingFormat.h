
int my_root_color[13]={1,    1,      2, kGreen+1,        4,   kOrange-3,  
	                   kPink+6,kAzure+7,  kTeal+4,kViolet-4,   kOrange-6,
	                   kGray+1,kRed-8 };
int my_root_color_10[10]={kGreen+1,        4,   kOrange-3,  kPink+6,kAzure+7,  kTeal+4,kViolet-4,   kOrange-6,
	                   kGray+1,kRed-8 };	         
int A_color_25[25]={kRed,kSpring,kAzure,kCyan,kTeal, kOrange,kViolet,1,kGray,kRed+2,   kPink+1,46,kOrange,kOrange-3,kCyan+3,    kSpring-6,38,41,kCyan+1,42,    kAzure+7,6,kViolet+6,9,38};
int Z_color_20[20]={0,0,0,0,1, kGray,2,kPink+1,46,kOrange-3,    5,3,kSpring-6,8,7,    kCyan+1,4,kAzure+7,6,kViolet+6};


void AxisFormat(TGraph* gr,const char* title= NULL,const char* x= NULL,const char* y= NULL,int color=1)
{
	// color 1=black 2=red 3=green 4=blue 5=yellow 6=pink 7=light blue 8=darkgreen 9=purple
	//color = (color%9==0)?9:color%9;
	gr->SetTitle(title);
	gr->GetXaxis()->SetTitle(x);
	gr->GetXaxis()->CenterTitle(1);
    gr->GetXaxis()->SetTitleOffset(0.8);
	gr->GetXaxis()->SetTitleSize(0.06);
	gr->GetXaxis()->SetLabelSize(0.05); //enlarge  axis label
	//gr->GetXaxis()->SetLabelFont(34);

	gr->GetYaxis()->SetTitle(y);
	gr->GetYaxis()->CenterTitle(1);
	gr->GetYaxis()->SetTitleOffset(0.8);
	gr->GetYaxis()->SetTitleSize(0.06);
	gr->GetYaxis()->SetLabelSize(0.04); //enlarge  axis label
	gr->SetLineColor(color);
	gr->SetMarkerColor(color);
	
	gr->SetLineWidth(2);
	gr->SetMarkerStyle(8);
	
	gr->SetMarkerSize(2);
}
void SetMyROOTColor(TGraph* gr,int color)
{
	gr->SetLineColor(my_root_color[color]);
	gr->SetMarkerColor(my_root_color[color]);
}

void AxisFormat(TGraphErrors *grerr,const char* title= NULL,const char* x= NULL,const char* y= NULL,int color=1)
{
	// color 1=black 2=red 3=green 4=blue 5=yellow 6=pink 7=light blue 8=darkgreen 9=purple
	grerr->SetTitle(title);
	grerr->GetXaxis()->SetTitle(x);
	grerr->GetXaxis()->CenterTitle(1);
    grerr->GetXaxis()->SetTitleOffset(0.8);
	grerr->GetXaxis()->SetTitleSize(0.05);
	grerr->GetXaxis()->SetLabelOffset(0.006); 
	grerr->GetXaxis()->SetLabelSize(0.04); //enlarge  axis label
	grerr->GetYaxis()->SetTitle(y);
	grerr->GetYaxis()->CenterTitle(1);
	grerr->GetYaxis()->SetTitleOffset(0.8);
	grerr->GetYaxis()->SetTitleSize(0.05);
	grerr->GetYaxis()->SetLabelOffset(0.006); 
	grerr->GetYaxis()->SetLabelSize(0.04); //enlarge  axis label
	grerr->SetLineColor(color);
	grerr->SetLineWidth(3);
	grerr->SetMarkerStyle(8);
	grerr->SetMarkerColor(color);
	grerr->SetMarkerSize(2);
}

void AxisFormat(TGraph2D* gr2d,const char* title= NULL,const char* x= NULL,const char* y= NULL,const char* z= NULL,int color=1)
{
	// color 1=black 2=red 3=green 4=blue 5=yellow 6=pink 7=light blue 8=darkgreen 9=purple
	//color = (color%9==0)?9:color%9;
	gr2d->SetTitle(title);
	gr2d->GetXaxis()->SetTitle(x);
	gr2d->GetXaxis()->CenterTitle(1);
    gr2d->GetXaxis()->SetTitleOffset(1.0);
	gr2d->GetXaxis()->SetTitleSize(0.03);
	gr2d->GetXaxis()->SetLabelSize(0.03); //enlarge  axis label
      
	gr2d->GetYaxis()->SetTitle(y);
	gr2d->GetYaxis()->CenterTitle(1);
	gr2d->GetYaxis()->SetTitleOffset(1.0);
	gr2d->GetYaxis()->SetTitleSize(0.03);
	gr2d->GetYaxis()->SetLabelSize(0.03); //enlarge  axis label
   
	gr2d->GetZaxis()->SetTitle(z);
	gr2d->GetZaxis()->CenterTitle(1);
	gr2d->GetZaxis()->SetTitleOffset(1.0);
	gr2d->GetZaxis()->SetTitleSize(0.03);
	gr2d->GetZaxis()->SetLabelSize(0.03); //enlarge  axis label
    
	gr2d->SetLineColor(my_root_color[color]);
	gr2d->SetLineWidth(2);
	gr2d->SetMarkerStyle(8);
	gr2d->SetMarkerColor(my_root_color[color]);
	gr2d->SetMarkerSize(3);
}

void AxisFormat(TH1* h,const char* title = NULL,const char* x=NULL,const char* y=NULL,int color=1)
{
    h->SetTitle(title);
	h->SetLineColor(color);
	h->SetLineWidth(2);
	h->GetXaxis()->SetTitle(x);
	h->GetXaxis()->CenterTitle(1);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.0);
	h->GetXaxis()->SetLabelSize(0.04);

	h->GetYaxis()->SetTitle(y);
	h->GetYaxis()->CenterTitle(1);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.0);
	h->GetYaxis()->SetLabelSize(0.04);
}
void AxisFormat(TMultiGraph* mg,const char* title= NULL,const char* x= NULL,const char* y= NULL)
{
	mg->SetTitle(title);
	mg->GetXaxis()->SetTitle(x);
	mg->GetXaxis()->CenterTitle(1);
	
    mg->GetXaxis()->SetTitleOffset(0.8);
	mg->GetXaxis()->SetTitleSize(0.06);
	mg->GetXaxis()->SetLabelSize(0.05); //enlarge  axis label
	mg->GetYaxis()->SetTitle(y);
	mg->GetYaxis()->CenterTitle(1);
	mg->GetYaxis()->SetTitleOffset(0.8);
	mg->GetYaxis()->SetTitleSize(0.06);
	mg->GetYaxis()->SetLabelSize(0.05); //enlarge  axis label
}

void Draw_one_TGraph(TGraph* gr,TString canvas_name)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	gr->Draw("ap");
}
void Draw_one_TGraph_yaxis_range(TGraph* gr,TString canvas_name, double y1, double y2)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	gr->Draw("ap");
	gr->GetYaxis()->SetRangeUser(y1,y2);
}
void Draw_one_TGraph_xaxis_range(TGraph* gr,TString canvas_name, double x1, double x2)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	gr->Draw("ap");
	gr->GetXaxis()->SetRangeUser(x1,x2);
}
void Draw_one_TGraph_xyaxis_range(TGraph* gr,TString canvas_name, double x1, double x2,double y1, double y2)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	gr->Draw("ap");
	gr->GetXaxis()->SetRangeUser(x1,x2);
	gr->GetYaxis()->SetRangeUser(y1,y2);

}
void Draw_one_TGraph_range(TGraph* gr,TString canvas_name, int opt,double x1, double x2,double y1, double y2)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	gr->Draw("ap");
	if(opt==0){}
	else if(opt==1){gr->GetXaxis()->SetRangeUser(x1,x2);}
	else if(opt==2){gr->GetYaxis()->SetRangeUser(y1,y2);}
	else if(opt==3){gr->GetXaxis()->SetRangeUser(x1,x2);gr->GetYaxis()->SetRangeUser(y1,y2);}	
	else {}

}
void Draw_one_histogram(TH1* h1,TString canvas_name)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	h1->Draw();
}
void Draw_one_histogram_range(TH1* h1,TString canvas_name,int opt,double x1, double x2,double y1, double y2)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,500);
	h1->Draw();
	if(opt==0){}
	else if(opt==1){h1->GetXaxis()->SetRangeUser(x1,x2);}
	else if(opt==2){h1->GetYaxis()->SetRangeUser(y1,y2);}
	else if(opt==3){h1->GetXaxis()->SetRangeUser(x1,x2);h1->GetYaxis()->SetRangeUser(y1,y2);}	
	else {}
}
void Draw_one_2Dhistogram(TH2* h2,TString canvas_name, int opt=0)
{
	TCanvas* c1 = new TCanvas(canvas_name,canvas_name,1000,800);
	h2->Draw("colz");
	//gPad->Update();
	if(opt==1)h2->GetXaxis()->SetNdivisions(505);
}