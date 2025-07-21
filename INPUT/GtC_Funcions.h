int ReadFileNames(string DIR_IN,string* FileNames,int opt1=1,int opt2=1)
{
    ofstream outfile;
    if(opt2)outfile.open(DIR_IN+"/ALL_FILENAMES_THIS_DIR.txt");
            
    string str;
            
    string Dir = DIR_IN; 
    int i=0;

    void* OpenFile = gSystem->OpenDirectory(Dir.c_str());  
    char *Direntry; 
    Direntry = ( (char *)gSystem->GetDirEntry(OpenFile) ) ;
    //get the names of files ,one by one,
    if(opt1)cout<<" ---------filenames in the directory:    "<<DIR_IN<<" :"<<endl;
    while(Direntry!=NULL)   
    {
        str=Direntry;
        if(str!="." &&str!=".."&&str!="ALL_FILENAMES_THIS_DIR.txt" &&str!="SHOW_ALL_FILENAMES.cpp" ) 
        {
            if(opt1)cout<<str<<endl;
            if(opt2)outfile<<str<<endl;
            FileNames[i]=str;
            i++;
        }
        Direntry = ( (char *)gSystem->GetDirEntry(OpenFile) ) ;
    }
    if(opt1)cout<<"---------- total:     "<<i<<endl;

    if(opt2)outfile.close();
    return i;   
}

// compare function used in qsort
int compare_double(const void * a, const void * b)
{
    if( ( (*(double*)a) - (*(double*)b) )<0 ) return -1 ;
    else if( ( (*(double*)a) - (*(double*)b) )>0 ) return 1 ;
    else return 0;
}
int compareIONSpecies (const void * a, const void * b)
{
    if( ( (*(IONSpecies*)a).AveT - (*(IONSpecies*)b).AveT)<0 ) return -1 ;
    else if( ( (*(IONSpecies*)a).AveT - (*(IONSpecies*)b).AveT)>0 ) return 1 ;
    else return 0;
}

int compareION (const void * a, const void * b)
{
    return ( (*(ION*)a).Species - (*(ION*)b).Species);
}
int compareION_C (const void * a, const void * b)
{  
    if( ( (*(ION*)a).C - (*(ION*)b).C)<0 ) return -1 ;
    else if( ( (*(ION*)a).C - (*(ION*)b).C)>0 ) return 1 ;
    else return 0;
}


Double_t average(Double_t *x, Int_t len)
{
    Double_t sum = 0;
    for (Int_t i = 0; i < len; i++) // 求和
        {sum += x[i];}
    return sum/len; // 得到平均值
}

Double_t variance(Double_t *x, Int_t len) /// variance of sample
{ 
    Double_t sum = 0;
    for (Int_t i = 0; i < len; i++) // 求和
        {sum += x[i];}
    Double_t aver = sum/len; // 得到平均值
    sum=0;
    for (Int_t i = 0; i < len; i++) // 求和 
        {sum += pow(x[i] - aver, 2);} 
    //return sqrt(sum/(len-1)); 
    return (sum/(len)); 
}
double GetVectorMean(vector<double>& v1)
{
    if(v1.size()>0)
    {
        double sum = std::accumulate(v1.begin(), v1.end(), 0.0);
        double mean = sum / v1.size();
        return mean;
    }
    else
    {
        return 0;
    }
    
}

double GetVectorStdDev(vector<double>& v1, double mean)
{
    if(v1.size()>1)
    {
        double sum=0;
        for(auto& i:v1){sum+= (i-mean)*(i-mean);}
        double StdDev = sqrt(sum / (v1.size()-1));
        return StdDev;
    }
    else
    {
        return 0;
    }
   
}
double Get_residue(TGraph* gr , TF1* f, TGraph* gr2)
{
    int n = gr->GetN();
    double x,y;
    double res=0;
    for(int i=0;i<n  ;i++)
    {
        gr->GetPoint(i,x,y);
        res+=pow(y-f->Eval(x),2);
        gr2->SetPoint(i,x,pow(y-f->Eval(x),2) );
    }
    return sqrt(res)/n;
}
double Get_residue(TGraphErrors* grerr , TF1* f, TGraphErrors* grerr2)
{
    int n = grerr->GetN();
    double x,y;
    double xerr,yerr;
    double res=0;
    for(int i=0;i<n  ;i++)
    {
        grerr->GetPoint(i,x,y);
        xerr=grerr->GetErrorX(i);
        yerr=grerr->GetErrorY(i);
        if(yerr==0){cout<<"yerr=0 at point"<<i<<" : ("<<x<<" , "<<y<<endl;}
        else {res+= pow(y-f->Eval(x),2)/(yerr*yerr);}
        grerr2->SetPoint(i,x,pow( y-f->Eval(x),2 ));
        grerr2->SetPointError(i,xerr,2*(y-f->Eval(x))*yerr);
    }
    return sqrt(res)/n;
}
double sigmaU_XY2(double U,double X,double SX,double Y,double SY)         //error for U=XY
{
    return U*U* ((SX*SX)/(X*X)+(SY*SY)/(Y*Y)); 
}

double ME_to_bare_nuc_mass(double ME, int A, int Z)
{//KeV
    double u  = 931494.102417;   //KeV
    double Me=510.9989;
    return A*u + ME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
}
void Print_fitfun_poln(TF1* f,int n)
{
    // 比如2阶多项式 pol2, n=2, 3 个系数
    double* p = new double[n] ;
   
    TString strtmp;
    f->GetParameters(p);
    for(int i=0;i<=n  ;i++)
    {
        cout<<" para: "<<i<<" = "<<p[i]<<" +- "<<f->GetParErrors()[i]<<endl;
    }
    delete []p;

}
TString Info_fitfun_pol3(TF1* f)
{
    double p[4];
    TString strtmp;
    f->GetParameters(p);
    return strtmp.Format(" %.5f + %.5f*x + %.5f*x^{2} + %.5f*x^{3}",p[0],p[1],p[2],p[3]);
}
TString Info_fitfun_pol2(TF1* f)
{
    double p[3];
    TString strtmp;
    f->GetParameters(p);
    return strtmp.Format(" %.5f + %.5f*x + %.5f*x^{2} ",p[0],p[1],p[2]);
}
TString Info_fitfun_pol1(TF1* f)
{
    double p[2];
    TString strtmp;
    f->GetParameters(p);
    return strtmp.Format(" %.5f + %.5f*x   ",p[0],p[1]);
}
TString Info_fitfun_pol0(TF1* f)
{
    double p[1];
    double perr[1];
    TString strtmp;
    f->GetParameters(p);
    perr[0] = f->GetParErrors()[0];
    return strtmp.Format(" %.5f +- %.5f   ",p[0],perr[0]);
}

void TGraph2D_to_outfile(TString outfile_name, TGraph2D* gr2d)
{
    ofstream outfile;
    outfile.open(outfile_name);
    int n = gr2d->GetN();
    double x=0,y=0,z=0;
    for(int i=0;i<n;i++)
    {
        gr2d->GetPoint(i,x,y,z);
        outfile<<x<<" "<<y<<" "<<z<<endl;
    }
    outfile.close();
    cout<<"-------TGraph2D_to_outfile : save "<<gr2d->GetTitle()<<" to "<<outfile_name<<endl;
}
void TGraph_to_outfile(TString outfile_name, TGraph* gr, int p=3)
{
    ofstream outfile;
    outfile.open(outfile_name);
    int n = gr->GetN();
    double x=0,y=0,z=0;
    for(int i=0;i<n;i++)
    {
        gr->GetPoint(i,x,y);
        if(p>3) outfile<<fixed<<setprecision(p)<<x<<" "<<y<<endl;
        else outfile<<x<<" "<<y<<endl;
    }
    outfile.close();
    cout<<"-------TGraph_to_outfile : save "<<gr->GetTitle()<<" to "<<outfile_name<<endl;
}
void TGraphErrors_to_outfile(TString outfile_name, TGraphErrors* grr)
{
    ofstream outfile;
    outfile.open(outfile_name);
    int n = grr->GetN();
    double x=0,y=0,yerr=0;
    for(int i=0;i<n;i++)
    {
        grr->GetPoint(i,x,y);
        yerr=grr->GetErrorY(i);
        outfile<<x<<" "<<y<<" "<<yerr<<endl;
    }
    outfile.close();
    cout<<"-------TGraphErrors_to_outfile() : save "<<grr->GetTitle()<<" to "<<outfile_name<<endl;
}
void TGraph_from_infile(TString infile_name, TGraph* gr)
{
    ifstream infile;
    infile.open(infile_name);
    int n=0 ;
    double x=0,y=0,z=0;
    while(infile>>x>>y)
    {
        gr->SetPoint(n++,x,y);
    }
    infile.close();
    cout<<"-------TGraph_from_infile "<<infile_name<<" read in lines: "<<n<<endl;
}
void TGraphErrors_from_infile(TString infile_name, TGraphErrors* grr)
{
    ifstream infile;
    infile.open(infile_name);
    int n=0 ;
    double x=0,y=0,yerr=0;
    while(infile>>x>>y>>yerr)
    {
        grr->SetPoint(n,x,y);
        grr->SetPointError(n,0,yerr);
        n++;
    }
    infile.close();
    cout<<"-------TGraphErrors_from_infile "<<infile_name<<" read in lines: "<<n<<endl;
}
void Hide_ErrorBar(TGraphErrors* grerr)
{
    grerr->SetLineColor(kWhite);
    grerr->SetLineWidth(0);
}
void Clear_ErrorBar(TGraphErrors* grerr)
{
    for(int i=0;i<grerr->GetN();i++)
    {
        grerr->SetPointError(i,0,0);
    }
}
void TGraph_to_TH1F(TGraph* gr,TH1F* h)
{
    int nn=gr->GetN();
    double xx,yy;
    for(int j=0;j<nn  ;j++)
    {
        gr->GetPoint(j,xx,yy);
        h->Fill(yy);
    }
}
    
// ================== two graph difference ========================
void Get_difference_of_two_TGraphErrors(TGraphErrors* grr1,TGraphErrors* grr2,TGraphErrors* grrd,double dd_limit)
{
    int N1 = grr1->GetN();
    int N2 = grr2->GetN();
    
    for(int i=0;i<N1  ;i++)
    {
        double tx1,ty1,tyerr1;
        
        grr1->GetPoint(i,tx1,ty1);
        tyerr1 =  grr1->GetErrorY(i);
        for(int j=0;j<N2  ;j++)
        {
            double tx2,ty2,tyerr2;
            grr2->GetPoint(j,tx2,ty2);
            tyerr2 =  grr2->GetErrorY(j);
            if( abs(tx1-tx2)<dd_limit )
            {
                grrd->SetPoint(grrd->GetN(),tx1,ty1-ty2);
                grrd->SetPointError(grrd->GetN()-1,0,sqrt(tyerr1*tyerr1+tyerr2*tyerr2) );

            }
        }
    }
}
//=================== grerr +- sigma to gr ==============================

void Grerr_sigma_to_gr(TGraphErrors* grerr_in, TGraph* gr,int opt)
{
    int n = grerr_in->GetN();
    double x,y,yerr;
    for(int i=0;i<n  ;i++)
    {
        grerr_in->GetPoint(i,x,y);
        yerr = grerr_in->GetErrorY(i);

        if(opt==1)gr->SetPoint(i,x,y+yerr);
        else if(opt==2)gr->SetPoint(i,x,y-yerr);
        
    }
    //Draw_one_TGraph(gr,strtmp.Format("c_gr_convert_check_%d",opt) );
    //gr->Print();
}

void FIND_TH1F_effective_L_R(TH1* h1, int nbins,int lb, int lb_n, double& L, double& R)
{   // lb: lower_bound. lb_n: continuous n bins that have counts more than lb.
    int L_i=1;
    int R_i=nbins;
    for(int i=1;i<nbins+1-lb_n-1  ;i++)
    {
        int count1=0;  // continuous n bins that have counts more than lb.
        int count2=0;  // monotonic bins
        for(int j=i;j<=i+lb_n  ;j++)
        {
            if(h1->GetBinContent(j) >=lb ){count1++;}
            if(h1->GetBinContent(j) <=h1->GetBinContent(j+1) ) {count2++;}
        }
        if(count1==lb_n&&count2==lb_n){ L_i = i;break;}
    }
     L=h1->GetBinCenter(L_i);
    if(L_i >= nbins/2){  L=h1->GetBinCenter(0); cout<<endl<<"FIND_TH1F_effective_L_R fail!! at L " <<endl;}
   

    for(int i=nbins-1;i>=1+lb_n+1  ;i--)
    {
        int count1=0;  // continuous n bins that have counts more than lb.
        int count2=0;  // monotonic bins
        for(int j=i;j>=i-lb_n+1  ;j--)
        {
            if(h1->GetBinContent(j) >=lb ){count1++;/*cout<<" count1 "<<j<<" "<<h1->GetBinContent(j)<<" "<<lb;*/}
            if(h1->GetBinContent(j) <=h1->GetBinContent(j-1) ) {count2++;/*cout<<" count2 "<<j<<" "<<h1->GetBinContent(j)<<" "<<h1->GetBinContent(j-1);*/}
            //cout<<endl<<j<<" "<<h1->GetBinContent(j)<<endl;
        }
        //cout<<count1<<" ,,, "<<count2<<endl;
        if(count1==lb_n&&count2==lb_n){ R_i = i;break;}
    }
    R=h1->GetBinCenter(R_i);
    if(R_i <= nbins/2){  R=h1->GetBinCenter(nbins-1); cout<<endl<<"FIND_TH1F_effective_L_R fail!! at R " <<endl;}
    

}


void Get_h1_percent_region(TH1* h1,double p, double& L, double& R)
{
   if(p<0||p>100){cout<<" error in Get_h1_percent_region !!! "<<endl<<" input p is out of 0~100"<<endl;}
   const Int_t nq = 10000;// 精度是 nq/100 %
   //double xq[100],yq[100];
   Double_t xq[nq];  // position where to compute the quantiles in [0,1]
   Double_t yq[nq];  // array to contain the quantiles
   for(int i=0;i<nq  ;i++)
   {
      xq[i] = float(i+1)/nq;
   }
   h1->GetQuantiles(nq,yq,xq);
   cout<<"debug!! in function: Get_h1_percent_region"<<endl;
   
   L=yq[int ( (100-p)/2/100*nq) ];
   R=yq[int(nq- (100-p)/2/100*nq ) ];

   cout<<"L[ "<<int ( (100-p)/2/100*nq)<<"]: "<<L<<endl;
   cout<<"R[ "<<int(nq- (100-p)/2/100*nq )<<"]: "<<R<<endl;
}


double ErrorFrom_UD_CENTER(double u, double d, double c )
{
    if( (u-c)/(c-d) >0 ) return abs(u-d)*0.5;
    else return  0.5* (abs(u-c)+abs(d-c));
}
double ErrorFrom_UD_CENTER_count(double u, double d, double c, int count1, int count2 )
{
    if( (u-c)/(c-d) >0 ) {count1++; return abs(u-d)*0.5; }
    else {count2++; return  0.5* (abs(u-c)+abs(d-c));}
}
void TGraph_shift(TGraph* gr1, TGraph* gr1_shift,double shift)
{
    double xtmp,ytmp=0;
    for(int i=0;i<gr1->GetN() ; i++)
    {
        gr1->GetPoint(i,xtmp,ytmp);
        gr1_shift->SetPoint(i,xtmp,ytmp+shift);
    }
}
//-------------------------------  vector ------------------------------------------

void GetErrorWeightedResult(vector<double> v1_m,  vector<double> v2_err, double& m, double& merr)
{
    if(v1_m.size()!=v2_err.size()) { cout<<" error in GetErrorWeightedResult , size of v1_m, v2_err not the same !!"<<endl;return; }
    int n = v1_m.size();
    vector<double> tmp_err_square;
    vector<double> tmp_err_square_m;
    for(auto &i:v2_err){ tmp_err_square.push_back(1.0/(i*i) ); }
    for(int i=0;i<n;i++){ tmp_err_square_m.push_back( tmp_err_square[i]*v1_m[i] ); }

    //for(auto &i:tmp_err_square){cout<<" tmp_err_square: "<<i;}   //debug
    //for(auto &i:tmp_err_square_m){cout<<" tmp_err_square_m: "<<i;}//debug

    double sum1 = std::accumulate(tmp_err_square_m.begin(), tmp_err_square_m.end(), 0.0);  // 最后一个必须是0.0 告诉accumulate 用double 而不是 int
    double sum2 = std::accumulate(tmp_err_square.begin(), tmp_err_square.end(), 0.0);
    //cout<<" sum1= "<<sum1<<" sum2= "<<sum2<<endl;//debug
    m = sum1/sum2;
    merr = sqrt(1/sum2);
}

double Get_IQR(vector<double> & v1)
{
    if(v1.size()==0){return 0;}
    sort(v1.begin(), v1.end());
    int n = v1.size();
    return v1[n*0.75] - v1[n*0.25];
}
double vector_StdDev(const std::vector<double>& vec) 
{
    if (vec.empty()) return 0;
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    double variance = 0.0;
    for (double x : vec) {variance += (x - mean) * (x - mean);}
    variance /= vec.size();
    return sqrt(variance);
}
double vector_average(const std::vector<double>& vec) 
{
    if (vec.empty()) return 0;
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    return mean;
}

void Get_subregion_StdDev_from_gr(TGraph* gr1, TH1F* h1,int h1_nbins, TGraphErrors* grr2)
{
    double x,y;
    vector<double>* vs = new vector<double>[h1_nbins];
    for(int i=0;i<gr1->GetN()  ;i++)
    {
        gr1->GetPoint(i,x,y);
        int ii = h1->FindBin(x)-1; // 0~N-1
        if(ii>=0&&ii<=h1_nbins-1)
        {
            vs[ii].push_back(y);
        }
    }
 
    //grr2 = new TGraphErrors();  // 在外面 new 再传进来
    for(int i=0;i<h1_nbins  ;i++)
    {
        if(vs[i].size()<1){continue;}
        double std = vector_StdDev(vs[i]);
        grr2->SetPoint(grr2->GetN(), h1->GetBinCenter(i) , std);
        grr2->SetPointError(grr2->GetN()-1, 0, std/(sqrt( 2*vs[i].size()-2 ) ) );
        //cout<<" debug in func: "<<i<<" std= "<<std<<endl;
    }

}

//-------------------------------------------------------------------------
void Do_2gaus(TH1F* h, double fmin, double fmax, double a1,double mu1, double sigma1,double a2,double mu2, double sigma2,
    double& fit_mu1,double& fit_mu1_err,double& fit_mu2,double& fit_mu2_err,
    double& fit_sigma1,double& fit_sigma1_err,double& fit_sigma2,double& fit_sigma2_err,
    bool to_show)
{
    TF1* fitfun_2gaus = new TF1("fit_2gaus isomer","gaus(0)+gaus(3)",fmin,fmax);
    fitfun_2gaus->SetParameters(a1,mu1, sigma1,a2,mu2,sigma2);
    //h->Fit(fitfun_2gaus , "Q");
    TFitResultPtr FitResult_h_mvq = h->Fit(fitfun_2gaus , "SQ");
    FitResult_h_mvq->Print("V");     // print full information of fit including covariance matrix 
    
    fit_mu1      = fitfun_2gaus->GetParameter(1);
    fit_sigma1   = fitfun_2gaus->GetParameter(2);
    fit_mu2    = fitfun_2gaus->GetParameter(4);
    fit_sigma2 = fitfun_2gaus->GetParameter(5);

    fit_mu1_err      = fitfun_2gaus->GetParErrors()[1];
    fit_sigma1_err   = fitfun_2gaus->GetParErrors()[2];
    fit_mu2_err    = fitfun_2gaus->GetParErrors()[4];
    fit_sigma2_err = fitfun_2gaus->GetParErrors()[5];

    TString h_title = h->GetTitle();
    cout<<" call Do_2gaus for h1: "<<h_title <<" ---- "<<endl;
    cout<<fixed<<setprecision(8)
        <<" mu1= "<<fit_mu1<<" +- "<<fit_mu1_err<<" sigma1= "<<fit_sigma1<<" +- "<<fit_sigma1_err<<endl
        <<" mu2= "<<fit_mu2<<" +- "<<fit_mu2_err<<" sigma2= "<<fit_sigma2<<" +- "<<fit_sigma2_err<<endl;  
    cout<<" chi2= "<<fitfun_2gaus->GetChisquare();

    if(to_show)
    {
        TCanvas *c_h_2gaus = new TCanvas("c_2gaus_"+h_title,"c_2gaus"+h_title, 1000,500);
        h->DrawClone();
    }
}


