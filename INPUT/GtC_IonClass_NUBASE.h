//--------------------------------IMP CSR IMS DToF GtC -------------------------------------------------------
//Institue of Modern Physics, Cooler Storage Ring group, Isochronous Mass Spectrometry, Double Time of Flight
// gammat-C method
// header IONCLASS
const Double_t V_c=0.2997924580;      //m/ns
const Double_t Vc_ms=299792458.0;      //m/s

const Double_t u  = 931494.102417;   //KeV
const Double_t mass_neutron = 939565.420 ;   // u+ 8071.318062(440)
const Double_t mass_proton = 938272.089 ;   // u+ 7288.971064(13) - 511
const Double_t mass_Hydrogen = 938783.0734810; //

#define ZN_MAX_Z 120
#define ZN_MAX_N 200
int ZN_n;
int ZN_AME[ZN_MAX_Z][ZN_MAX_N];  //对应于ISS序号--ground state -- NUBASE没有的是 -1



//10.3642686577806 = u[keV] / (V_c*V_c) *10^6
double v_to_gamma(double v)
{
	if(v<300000000.0&&v>0.3)return pow((1-(v*v/Vc_ms/Vc_ms)),-0.5);
	else if((v<0.3&&v>=0.0))return pow((1-(v*v/V_c/V_c)),-0.5);
	else {cout<<endl<<" error in v_to_gamma!!!! v unit is not m/s or m/ns "<<endl;return -1;}
}
//根据名字返回Z
int convert_name_to_z(TString name)
{
	if(name.CompareTo("n")==0 ) return 0;
	else if(name.CompareTo("H")==0 ) return 1;
	else if(name.CompareTo("He")==0 ) return 2;
	else if(name.CompareTo("Li")==0 ) return 3;
	else if(name.CompareTo("Be")==0 ) return 4;
	else if(name.CompareTo("B" )==0 ) return 5;
	else if(name.CompareTo("C" )==0 ) return 6;
	else if(name.CompareTo("N" )==0 ) return 7;
	else if(name.CompareTo("O" )==0 ) return 8;
	else if(name.CompareTo("F" )==0 ) return 9;
    
    else if(name.CompareTo("Ne")==0 ) return 10;else if(name.CompareTo("Na")==0 ) return 11;else if(name.CompareTo("Mg")==0 ) return 12;else if(name.CompareTo("Al")==0 ) return 13;else if(name.CompareTo("Si")==0 ) return 14;else if(name.CompareTo("P")==0 )  return 15;else if(name.CompareTo("S")==0 )  return 16;else if(name.CompareTo("Cl")==0 ) return 17;else if(name.CompareTo("Ar")==0 ) return 18;else if(name.CompareTo("K")==0 )  return 19;
    else if(name.CompareTo("Ca")==0 ) return 20;else if(name.CompareTo("Sc")==0 ) return 21;else if(name.CompareTo("Ti")==0 ) return 22;else if(name.CompareTo("V")==0 )  return 23;else if(name.CompareTo("Cr")==0 ) return 24;else if(name.CompareTo("Mn")==0 ) return 25;else if(name.CompareTo("Fe")==0 ) return 26;else if(name.CompareTo("Co")==0 ) return 27;else if(name.CompareTo("Ni")==0 ) return 28;else if(name.CompareTo("Cu")==0 ) return 29;
    else if(name.CompareTo("Zn")==0 ) return 30;else if(name.CompareTo("Ga")==0 ) return 31;else if(name.CompareTo("Ge")==0 ) return 32;else if(name.CompareTo("As")==0 ) return 33;else if(name.CompareTo("Se")==0 ) return 34;else if(name.CompareTo("Br")==0 ) return 35;else if(name.CompareTo("Kr")==0 ) return 36;else if(name.CompareTo("Rb")==0 ) return 37;else if(name.CompareTo("Sr")==0 ) return 38;else if(name.CompareTo("Y")==0 ) return 39;
    else if(name.CompareTo("Zr")==0 ) return 40;else if(name.CompareTo("Nb")==0 ) return 41;else if(name.CompareTo("Mo")==0 ) return 42;else if(name.CompareTo("Tc")==0 ) return 43;else if(name.CompareTo("Ru")==0 ) return 44;else if(name.CompareTo("Rh")==0 ) return 45;else if(name.CompareTo("Pd")==0 ) return 46;else if(name.CompareTo("Ag")==0 ) return 47;else if(name.CompareTo("Cd")==0 ) return 48;else if(name.CompareTo("In")==0 ) return 49;
    else if(name.CompareTo("Sn")==0 ) return 50;else if(name.CompareTo("Sb")==0 ) return 51;else if(name.CompareTo("Te")==0 ) return 52;else if(name.CompareTo("I")==0 ) return 53;else if(name.CompareTo("Xe")==0 ) return 54;else if(name.CompareTo("Cs")==0 ) return 55;else if(name.CompareTo("Ba")==0 ) return 56;else if(name.CompareTo("La")==0 ) return 57;else if(name.CompareTo("Ce")==0 ) return 58;else if(name.CompareTo("Pr")==0 ) return 59;
    else if(name.CompareTo("Nd")==0 ) return 60;else if(name.CompareTo("Pm")==0 ) return 61;else if(name.CompareTo("Sm")==0 ) return 62;else if(name.CompareTo("Eu")==0 ) return 63;else if(name.CompareTo("Gd")==0 ) return 64;else if(name.CompareTo("Tb")==0 ) return 65;else if(name.CompareTo("Dy")==0 ) return 66;else if(name.CompareTo("Ho")==0 ) return 67;else if(name.CompareTo("Er")==0 ) return 68;else if(name.CompareTo("Tm")==0 ) return 69;
    else if(name.CompareTo("Yb")==0 ) return 70;else if(name.CompareTo("Lu")==0 ) return 71;else if(name.CompareTo("Hf")==0 ) return 72;else if(name.CompareTo("Ta")==0 ) return 73;else if(name.CompareTo("W")==0 ) return 74;else if(name.CompareTo("Re")==0 ) return 75;else if(name.CompareTo("Os")==0 ) return 76;else if(name.CompareTo("Ir")==0 ) return 77;else if(name.CompareTo("Pt")==0 ) return 78;else if(name.CompareTo("Au")==0 ) return 79;
    else if(name.CompareTo("Hg")==0 ) return 80;else if(name.CompareTo("Tl")==0 ) return 81;else if(name.CompareTo("Pb")==0 ) return 82;else if(name.CompareTo("Bi")==0 ) return 83;else if(name.CompareTo("Po")==0 ) return 84;else if(name.CompareTo("At")==0 ) return 85;else if(name.CompareTo("Rn")==0 ) return 86;else if(name.CompareTo("Fr")==0 ) return 87;else if(name.CompareTo("Ra")==0 ) return 88;else if(name.CompareTo("Ac")==0 ) return 89;
    else if(name.CompareTo("Th")==0 ) return 90;else if(name.CompareTo("Pa")==0 ) return 91;else if(name.CompareTo("U")==0 ) return 92;else if(name.CompareTo("Np")==0 ) return 93;else if(name.CompareTo("Pu")==0 ) return 94;else if(name.CompareTo("Am")==0 ) return 95;else if(name.CompareTo("Cm")==0 ) return 96;else if(name.CompareTo("Bk")==0 ) return 97;else if(name.CompareTo("Cf")==0 ) return 98;else if(name.CompareTo("Es")==0 ) return 99;
    else if(name.CompareTo("Fm")==0 ) return 100;else if(name.CompareTo("Md")==0 ) return 101;else if(name.CompareTo("No")==0 ) return 102;else if(name.CompareTo("Lr")==0 ) return 103;else if(name.CompareTo("Rf")==0 ) return 104;else if(name.CompareTo("Db")==0 ) return 105;else if(name.CompareTo("Sg")==0 ) return 106;else if(name.CompareTo("Bh")==0 ) return 107;else if(name.CompareTo("Hs")==0 ) return 108;else if(name.CompareTo("Mt")==0 ) return 109;
    else if(name.CompareTo("Ds")==0 ) return 110;else if(name.CompareTo("Rg")==0 ) return 111;else if(name.CompareTo("Cn")==0 ) return 112;else if(name.CompareTo("Nh")==0 ) return 113;else if(name.CompareTo("Fl")==0 ) return 114;else if(name.CompareTo("Mc")==0 ) return 115;else if(name.CompareTo("Lv")==0 ) return 116;else if(name.CompareTo("Ts")==0 ) return 117;else if(name.CompareTo("Og")==0 ) return 118;
	

	else return -1;

}
//根据z返回名字
TString convert_z_to_name(int z)
{
	switch(z)
	{
		case 0:return"n" ;break;
		case 1:return"H" ;break;
		case 2:return"He" ;break;
		case 3:return"Li" ;break;
		case 4:return"Be" ;break;
		case 5:return"B" ;break;
		case 6:return"C" ;break;
		case 7:return"N" ;break;
		case 8:return"O" ;break;
		case 9:return"F" ;break;
		case 10:return"Ne" ;break;
		case 11:return"Na" ;break;
		case 12:return"Mg" ;break;
		case 13:return"Al" ;break;
		case 14:return"Si" ;break;
		case 15:return"P" ;break;
		case 16:return"S" ;break;
		case 17:return"Cl" ;break;
		case 18:return"Ar" ;break;
		case 19:return"K" ;break;
		case 20:return"Ca" ;break;
		case 21:return"Sc" ;break;
		case 22:return"Ti" ;break;
		case 23:return"V" ;break;
		case 24:return"Cr" ;break;
		case 25:return"Mn" ;break;
		case 26:return"Fe" ;break;
		case 27:return"Co" ;break;
		case 28:return"Ni" ;break;
		case 29:return"Cu" ;break;

		case 30:return"Zn" ;break;
		case 31:return"Ga" ;break;
		case 32:return"Ge" ;break;
		case 33:return"As" ;break;
		case 34:return"Se" ;break;
		case 35:return"Br" ;break;
		case 36:return"Kr" ;break;
		case 37:return"Rb" ;break;
		case 38:return"Sr" ;break;
		case 39:return"Y" ;break;
		                                    
		case 40:return"Zr" ;break;
		case 41:return"Nb" ;break;
		case 42:return"Mo" ;break;
		case 43:return"Tc" ;break;
		case 44:return"Ru" ;break;
		case 45:return"Rh" ;break;
		case 46:return"Pd" ;break;
		case 47:return"Ag" ;break;
		case 48:return"Cd" ;break;
		case 49:return"In" ;break;
		                              
		case 50:return"Sn" ;break;
		case 51:return"Sb" ;break;
		case 52:return"Te" ;break;
		case 53:return"I" ;break;
		case 54:return"Xe" ;break;
		case 55:return"Cs" ;break;
		case 56:return"Ba" ;break;
		case 57:return"La" ;break;
		case 58:return"Ce" ;break;
		case 59:return"Pr" ;break;
		                              
		case 60:return"Nd" ;break;
		case 61:return"Pm" ;break;
		case 62:return"Sm" ;break;
		case 63:return"Eu" ;break;
		case 64:return"Gd" ;break;
		case 65:return"Tb" ;break;
		case 66:return"Dy" ;break;
		case 67:return"Ho" ;break;
		case 68:return"Er" ;break;
		case 69:return"Tm" ;break;
		                            
		case 70:return"Yb" ;break;
		case 71:return"Lu" ;break;
		case 72:return"Hf" ;break;
		case 73:return"Ta" ;break;
		case 74:return"W" ;break;
		case 75:return"Re" ;break;
		case 76:return"Os" ;break;
		case 77:return"Ir" ;break;
		case 78:return"Pt" ;break;
		case 79:return"Au" ;break;
                                      
		case 80:return"Hg" ;break;
		case 81:return"Tl" ;break;
		case 82:return"Pb" ;break;
		case 83:return"Bi" ;break;
		case 84:return"Po" ;break;
		case 85:return"At" ;break;
		case 86:return"Rn" ;break;
		case 87:return"Fr" ;break;
		case 88:return"Ra" ;break;
		case 89:return"Ac" ;break;
		                                
		case 90:return"Th" ;break;
		case 91:return"Pa" ;break;
		case 92:return"U" ;break;
		case 93:return"Np" ;break;
		case 94:return"Pu" ;break;
		case 95:return"Am" ;break;
		case 96:return"Cm" ;break;
		case 97:return"Bk" ;break;
		case 98:return"Cf" ;break;
		case 99:return"Es" ;break;
                                                     
		case 100:return"Fm" ;break;
		case 101:return"Md" ;break;
		case 102:return"No" ;break;
		case 103:return"Lr" ;break;
		case 104:return"Rf" ;break;
		case 105:return"Db" ;break;
		case 106:return"Sg" ;break;
		case 107:return"Bh" ;break;
		case 108:return"Hs" ;break;
		case 109:return"Mt" ;break;
                                             
		case 110:return"Ds" ;break;
		case 111:return"Rg" ;break;
		case 112:return"Cn" ;break;
		case 113:return"Nh" ;break;
		case 114:return"Fl" ;break;
		case 115:return"Mc" ;break;
		case 116:return"Lv" ;break;
		case 117:return"Ts" ;break;
		case 118:return"Og" ;break;
		

		default:return"default";break;
	}
}
class ION
{
public:
	int ID;
	TString name;
	TString inject_filename;
	int ion_number;					//serial number in this injection
	int inject_number;              //appear in which injection
	int Species;  //from 0 to NSpecies
	int A;
	int Z;
	int C_region;

	double T,T_err,v,v_err;    // T :[ns] v [m/ns]
	double v_err2;
	// dA0--A5   dA1--A6
    double A1,A1err,A2,A2err,A3,A4,dA0,dA0err,dA1,dA1err;   //注意这里所有的Aerr 都是误差的平方！ [平方皮秒]
    long double A3err;
    double cov12,cov15,cov16,cov25,cov26,cov56;
	double C,C_err2;    //m
	double Bp_true;   //Tm
	double Bp;   //Tm
	double Bp_err2;   //Tm
	
	//======= 2017_58Ni =======
	double turn_middle;
	long double cov13,cov23,cov35,cov36;
	//======= ISF =======
	double cov_Brou_C;
	double cov_C_v;
	double dFr;//参考核误差影响因子
	double dFt;//目标核误差影响因子

	//double v_ns;    //m/ns
	//double v_err_ns;    //m/ns
	double gammat;
	double gammat_err;   //20230707
	//double gammat_v2;  // another way of calculating gt using Bp-C
	bool gt_abnormal;
	double time;
	double m_VE,m_VE_err;    // calculated final mass
	double mvq_v1;           // m/q v1

	bool Do_dA0_T_flag;

	
	ION():ID(-1),T(0),C(0),Bp(0),v(0),dFr(0), dFt(0) ,gammat(0),gammat_err(0),gt_abnormal(0),
	m_VE(0),m_VE_err(0),turn_middle(0)
	{};

	//Set
	
	void SetA(TString str)
	{
		A=10*( int(str[0])-48 ) + (int(str[1]) -48);
	}
	void SetZ(TString str)
	{
		while( (int(str[0])-48)<=9&&(int(str[0])-48)>=0 )
		{str.Remove(0,1);}
		Z=convert_name_to_z(str);
	}
	//---------- Nm =0 ------------ 
	double Calculate_only_return_T()
	{
		return (A1-dA1*0.5)*0.001; //[ns]
	}
	double Calculate_only_return_v(double L ,double ddT)
	{
		return L/(dA0*0.001+ddT); //[m/ns]
	}
	
	

	double Calculate_only_return_gamma()
	{
		return 1.0/sqrt(1-(v*v/V_c/V_c));
	}
	double Calculate_only_return_Bp(double vns, double m_ref)
	{
		return m_ref/Z*0.000001/(V_c*V_c)*sqrt( (1/(vns*vns)) - (1/(V_c*V_c)) );
	}
	void Calculate_gt(double ddT)
	{
		gammat = sqrt( (1.0/(1 - pow(v/V_c,2))) / (1 - A2*2/T*(dA0*0.001+ddT)/dA1 ) );
	}
	void Calculate_gt_from_a_with_err(double L,double ddT)
	{
		//a5=dA0 a6=dA1
		// all Aerr is [ps*ps]
		double vc_ps = V_c*0.001;  // speed of light [m/ps]
		double v_ps = v*0.001;    // time:[ps] length:[m]
		double ddt_ps = ddT*1000; 
		double ga = pow( 1-(v_ps/vc_ps)*(v_ps/vc_ps) , -0.5);
		double Da1,Da2,Da5,Da6=0;
		double S = A1*dA1 - 0.5*dA1*dA1;
		double B = 2*A2*L/S;
		double A = 1 - B/v;
		gammat = pow(1-2*A2*(dA0+ddt_ps)/(A1-dA1*0.5)/dA1 , -0.5 )*ga;

		Da1 = (-1.0)*ga*L/v*pow(A,-1.5)*A2*dA1*pow(S,-2);
		Da2 = ga*L/v*pow(A,-1.5)/S;
		Da5 = pow((1 - B/v_ps - v_ps*v_ps/vc_ps/vc_ps + B/vc_ps/vc_ps*v_ps),-1.5)
				*(B - 2*v_ps*v_ps*v_ps/vc_ps/vc_ps +B/vc_ps/vc_ps*v_ps*v_ps) *0.5 /L;
		Da6 = (-1.0)*ga*L/v*pow(A,-1.5)*A2*(A1-dA1)*pow(S,-2);
		gammat_err = Da1*Da1*A1err + Da2*Da2*A2err + Da5*Da5*dA0err + Da6*Da6*dA1err
				   + Da1*Da2*cov12 + Da1*Da5*cov15 + Da1*Da6*cov16 + Da2*Da5*cov25 + Da2*Da6*cov26 + Da5*Da6*cov56;
		gammat_err = sqrt(gammat_err);
	}
	void Calculate_gt_from_a_with_err_58Ni(double L,double ddT)
	{
		//a5=dA0 a6=dA1
		// all Aerr is [ps*ps]
		double vc_ps = V_c*0.001;  // speed of light [m/ps]
		double v_ps = v*0.001;    // time:[ps] length:[m]
		double ddt_ps = ddT*1000; 
		double ga = pow( 1-(v_ps/vc_ps)*(v_ps/vc_ps) , -0.5);
		double Da1,Da2,Da5,Da6=0;
		double S = A1*dA1 - 0.5*dA1*dA1;
		double B = 2*A2*L/S;
		double A = 1 - B/v;
		gammat = pow(1-(2*A2+6*A3*turn_middle)*(dA0+ddt_ps)/(A1-dA1*0.5)/dA1 , -0.5 )*ga;

		//未考虑A3err
		Da1 = (-1.0)*ga*L/v*pow(A,-1.5)*A2*dA1*pow(S,-2);
		Da2 = ga*L/v*pow(A,-1.5)/S;
		Da5 = pow((1 - B/v_ps - v_ps*v_ps/vc_ps/vc_ps + B/vc_ps/vc_ps*v_ps),-1.5)
				*(B - 2*v_ps*v_ps*v_ps/vc_ps/vc_ps +B/vc_ps/vc_ps*v_ps*v_ps) *0.5 /L;
		Da6 = (-1.0)*ga*L/v*pow(A,-1.5)*A2*(A1-dA1)*pow(S,-2);
		gammat_err = Da1*Da1*A1err + Da2*Da2*A2err + Da5*Da5*dA0err + Da6*Da6*dA1err
				   + Da1*Da2*cov12 + Da1*Da5*cov15 + Da1*Da6*cov16 + Da2*Da5*cov25 + Da2*Da6*cov26 + Da5*Da6*cov56;
		gammat_err = sqrt(gammat_err);
	}
	void Calculate_gt_from_a_with_err_58Ni_Nm(double L,double ddT, int i)
	{
		bool test = (i%2000 == 0);
		double vc_ps = V_c*0.001;  // speed of light [m/ps]
		double v_ps = v*0.001;    // [m/ps]
		double ddt_ps = ddT*1000;  //[ps]
		double ga = pow( 1-(v_ps/vc_ps)*(v_ps/vc_ps) , -0.5);  //gamma
		gammat = pow(1-(2*A2+6*A3*turn_middle)*(dA0+ddt_ps)/(A1-dA1*0.5)/dA1 , -0.5 )*ga;
		// gt = gamma* pow( 1 - (W*U/S) , -0.5 );
		double S = A1*dA1 - 0.5*dA1*dA1;   // a1 ,a6
		double W = 2*A2+6*A3*turn_middle;  //a2,a3
		double U = dA0+ddt_ps;			  //a5
		
		//if(test)cout<<"test gt :"<<gammat<<" || "<<ga* pow( 1-(W*U/S) , -0.5 )<<endl ;
		
		double tmp1 = -0.5*pow(gammat,3.0)/ga/ga;
		double Dgt_S = tmp1 * W * U / S / S ;
		double Dgt_W = tmp1 * (-1.0) * U / S;
		double Dgt_U = tmp1 * (-1.0) * W / S;

		double DS_A1 = dA1;
		double DS_A6 = A1-dA1;
		double DW_A2 = 2;
		double DW_A3 = 6*turn_middle;
		// chain rule of derivative 
		double Dgt_A1,Dgt_A2,Dgt_A3,Dgt_A5,Dgt_A6;
		Dgt_A1 = Dgt_S * DS_A1;
		Dgt_A2 = Dgt_W * DW_A2;
		Dgt_A3 = Dgt_W * DW_A3;
		Dgt_A5 = Dgt_U;
		Dgt_A6 = Dgt_S * DS_A6;

		gammat_err = Dgt_A1*Dgt_A1*A1err + Dgt_A2*Dgt_A2*A2err + Dgt_A3*Dgt_A3*A3err + Dgt_A5*Dgt_A5*dA0err + Dgt_A6*Dgt_A6*dA1err
				   + 2*Dgt_A1*Dgt_A2*cov12 + 2*Dgt_A1*Dgt_A3*cov13 + 2*Dgt_A1*Dgt_A5*cov15 + 2*Dgt_A1*Dgt_A6*cov16 
				   + 2*Dgt_A2*Dgt_A3*cov23 + 2*Dgt_A2*Dgt_A5*cov25 + 2*Dgt_A2*Dgt_A6*cov26 
				   + 2*Dgt_A3*Dgt_A5*cov35 + 2*Dgt_A3*Dgt_A6*cov36 
				   + 2*Dgt_A5*Dgt_A6*cov56 ;
		gammat_err = sqrt(gammat_err);

	}
	void Calculate_Bp(double m)  //m: keV v_ms:m/s
	{
		double v_ms = v*1000000000.0;
		Bp = (m/Z)*1000.0*v_ms/(Vc_ms*sqrt(Vc_ms*Vc_ms-v_ms*v_ms));
		Bp_true=Bp;         //record Bp calculated using AME
	}
	
	//---------- Nm !=0 ------------
	double Get_T_Nm()
	{
		double Nm = turn_middle;
		return (A1-dA1*0.5+2*A2*Nm+3*A3*Nm*Nm ); //[ps]
	}
	double Get_Terr2_Nm()
	{
		double Nm = turn_middle;
		double DTA1,DTA2,DTA3,DTA6;  // 偏导数 T 对A
		DTA1 = 1;
		DTA2 = 2*Nm;
		DTA3 = 3*Nm*Nm;
		DTA6 = -0.5;
		return (DTA1*DTA1*A1err + DTA2*DTA2*A2err + DTA3*DTA3*A3err + DTA6*DTA6*dA1err
		+ 2*DTA1*DTA2*cov12 + 2*DTA1*DTA3*cov13 + 2*DTA1*DTA6*cov16 + 2*DTA2*DTA3*cov23 + 2*DTA2*DTA6*cov26 + 2*DTA3*DTA6*cov36); //[ps]
	}
	double Get_v_Nm(double L ,double ddT)
	{
		double Nm = turn_middle;
		return L /(dA0 + dA1*Nm + ddT*1000); //[m/ps]
	}
	double Get_verr2_Nm(double L ,double ddT )
	{
		double Nm = turn_middle;
		double Cal_v = L /(dA0 + dA1*Nm + ddT*1000);
		 
		double DvA5,DvA6;  // 偏导数 T 对A
		DvA5 = (-1.0)*Cal_v*Cal_v/L;
		DvA6 = DvA5*Nm;
		return (DvA5*DvA5*dA0err + DvA6*DvA6*dA1err + 2*DvA5*DvA6*cov56);
	}
	double Get_Cov_T_v_Nm(double L ,double ddT)
	{
		double Nm = turn_middle;
		double DTA1,DTA2,DTA3,DTA6;  // 偏导数 T 对A
		DTA1 = 1;
		DTA2 = 2*Nm;
		DTA3 = 3*Nm*Nm;
		DTA6 = -0.5;
		double Cal_v = L /(dA0 + dA1*Nm + ddT*1000);
		double DvA5,DvA6;  // 偏导数 T 对A
		DvA5 = (-1.0)*Cal_v*Cal_v/L;
		DvA6 = DvA5*Nm;

		return (DTA1*DvA5*cov15 + DTA2*DvA5*cov25 + DTA3*DvA5*cov35 + DTA6*DvA5*cov56
			  + DTA1*DvA6*cov16 + DTA2*DvA6*cov26 + DTA3*DvA6*cov36 + DTA6*DvA6*dA1err);

	}

	/*
	void Calculate_gt_v2(double L,double ddT,double m,int i)
	{
		double v_1 = L/( (dA0+dA1)*0.001 + ddT);
		double T_1 = T+(2*A2+3*A3)*0.001 ;  // ns
		double v_ms_1 = v_1*1000000000.0;
		double Bp_1 = (m/Z)*1000.0*v_ms_1/(Vc_ms*sqrt(Vc_ms*Vc_ms-v_ms_1*v_ms_1));
		double C_1 = T_1*v_1;
		//if(i%200==0)cout<<v_1<<" "<<T_1<<" "<<Bp_1<<" "<<C_1<<endl;
		gammat_v2 = sqrt( (Bp_1-Bp)*C_1/(C_1-C)/Bp_1 );
	}
	*/
	void operator=(const ION& ion_ref)
	{
		name = ion_ref.name;
		ion_number = ion_ref.ion_number;
		inject_number = ion_ref.inject_number;
		inject_filename = ion_ref.inject_filename;
		Species = ion_ref.Species;  //from 0 to NSpecies
		A = ion_ref.A;
		Z = ion_ref.Z;
		T = ion_ref.T;    //ns
		T_err = ion_ref.T_err;
		C = ion_ref.C;    //m
		C_err2 = ion_ref.C_err2;    //m
		Bp = ion_ref.Bp;   //Tm
		Bp_err2 = ion_ref.Bp_err2;
		Bp_true = ion_ref.Bp_true;   //Tm
		v = ion_ref.v;    //m/ns
		v_err = ion_ref.v_err;    //m/ns
		v_err2 = ion_ref.v_err2;
		gammat = ion_ref.gammat;
		gammat_err = ion_ref.gammat_err;
		A1 = ion_ref.A1; A2 = ion_ref.A2; A3 = ion_ref.A3; A4 = ion_ref.A4;
		A1err = ion_ref.A1err; A2err = ion_ref.A2err;
		dA0 = ion_ref.dA0; dA0err = ion_ref.dA0err;
		dA1 = ion_ref.dA1; dA1err = ion_ref.dA1err;
		cov12 = ion_ref.cov12; cov15 = ion_ref.cov15; cov16 = ion_ref.cov16;
		cov25 = ion_ref.cov25; cov26 = ion_ref.cov26; cov56 = ion_ref.cov56;
		//======= 2017_58Ni
		turn_middle = ion_ref.turn_middle;       // N 没有平移到中间圈数
		A3err=ion_ref.A3err;				  //引入了 A3 及其误差
		cov13=ion_ref.cov13;cov23=ion_ref.cov23;cov35=ion_ref.cov35;cov36=ion_ref.cov36; //引入了A3 相关协方差
		//======= ISF
		dFr = ion_ref.dFr;//参考核误差影响因子
		dFt = ion_ref.dFt;//目标核误差影响因子
	}
	//Print
	void PrintT(){cout<<"T = "<<T<<endl;}
	void PrintInfo()
	{
	    cout<<"|name= "<<name<<"| "<<fixed<<setprecision(10)<<" |Z = "<<Z<<"| |T = "<<T<<"| |C = "<<C<<"| |v = "<<v<<"| |gammat = "<<gammat<<"| |Bp = "<<Bp<<"|"<<endl;
	}
	void PrintInfo_readin(ofstream& outfile)
	{
	    outfile<<A<<" "<<name<<fixed<<setprecision(10)<<" Z= "<<Z<<" T= "<<T<<" +- "<<T_err
	    <<" A1= "<<A1<<" (var=) "<<A1err<<" A2= "<<A2<<" var= "<<A2err<<" A3= "<<A3<<" A4= "<<A4
	    <<" dA0= "<<dA0<<" var= "<<dA0err<<" dA1= "<<dA1<<" var= "<<dA1err
	    <<" cov12= "<<cov12<<" cov15= "<<cov15<<" cov16= "<<cov16
	    <<" cov25= "<<cov25<<" cov26= "<<cov26<<" cov56= "<<cov56
	    <<" time = "<<time  
	    <<endl;
	}
	void PrintInfo_readin_58Ni(ofstream& outfile)
	{
	    outfile<<A<<" "<<name<<fixed<<setprecision(10)<<" Z= "<<Z<<" T= "<<T<<" +- "<<T_err
	    <<" A1= "<<A1<<" (var=) "<<A1err<<" A2= "<<A2<<" var= "<<A2err
	    <<" A3= "<<A3<<" var= "<<fixed<<setprecision(20)<<A3err
	    <<fixed<<setprecision(10)<<" A4= "<<A4
	    <<" dA0= "<<dA0<<" var= "<<dA0err<<" dA1= "<<dA1<<" var= "<<dA1err
	    <<" cov12= "<<cov12<<" cov15= "<<cov15<<" cov16= "<<cov16
	    <<" cov25= "<<cov25<<" cov26= "<<cov26<<" cov56= "<<cov56 
	    <<" cov13= "<<cov13<<" cov23= "<<cov23<<" cov35= "<<cov35<<" cov36= "<<cov36
	    <<fixed<<setprecision(0)<<" time = "<<time 
	    <<endl;
	}
	//after T ,v already set
	void PrintInfo_complete(ofstream& outfile)
	{
		outfile<<A<<" "<<name<<" "<<Z<<" "<<inject_number<<" "<<inject_filename<<" "
		<<fixed<<setprecision(7)<<" T= "<<T<<" +- "<<T_err<<" C= "<<C<<" +- "<<sqrt(C_err2)
		<<fixed<<setprecision(10)<<" v= "<<v<<" +- "<<v_err<<" Bp= "<<Bp<<" +- "<<sqrt(Bp_err2)<<" "
		<<" cov_C_v= "<<cov_C_v<<" cov_Bp_C= "<<cov_Brou_C<<" gt = "<<gammat<<" +- "<<gammat_err<< endl;
	}
};

class ION_UNKNOWN:public ION
{
public:
	int ID;
	double M_cal;                 //final result
	double M_cal_VE;                 // any version of error weighted MASS_VER>2
	double M_cal_err;             //sigma m_calculated
	double M_cal_err_VE;             // treat each ref-cal as independent 20231023
	
	double Mvq;
	double Mvq_VE;
	double Mvq_err_VE;

	double M_AME;				  // Nuclear Mass AME INFO for comparison 
	double M_AME_err;
	
	double BpSum;				  // sum of Bp from each calculaiton using different ion_ref, to be averaged later
	double BpSum2;				  // sum of Bp*Bp from each calculaiton using different ion_ref, for sigma calculation
	//double Bp_err;				  // statistical result: std deviation sigma of Bp
	double Bp_err2;				  //error square
	double Bp_v2;                 //20230223

	double v_err2;
	double cov_C_v;
	double dFr;//参考核误差影响因子
	double dFt;//目标核误差影响因子

	bool use;					  //whether M_cal of this ions_unknown is effective or not
	int ref_n;                   //how many ref ions used

	void operator=(const ION& ion_ref)   //generate it from known ion
	{
		ID = ion_ref.ID;
		name = ion_ref.name;
		ion_number = ion_ref.ion_number;   //unknown ion does not have serial number in one injection
		inject_number = ion_ref.inject_number;
		inject_filename = ion_ref.inject_filename;//20240624
		Species = ion_ref.Species;  //from 0 to NSpecies
		A = ion_ref.A;
		Z = ion_ref.Z;
		T = ion_ref.T;    // ns
		T_err = ion_ref.T_err;
		C = ion_ref.C;    //m
		C_err2 = ion_ref.C_err2;    //m
		Bp = ion_ref.Bp;   //Tm
		Bp_true = ion_ref.Bp_true;   //Tm
		Bp_err2 = ion_ref.Bp_err2;
		v = ion_ref.v;    //m/ns
		v_err = ion_ref.v_err;    //m/ns
		v_err2 = ion_ref.v_err2;    //(m/ps)^2
		gammat = ion_ref.gammat;
		gammat_err = ion_ref.gammat_err;
		cov_C_v = ion_ref.cov_C_v;
		dFr = ion_ref.dFr;
		dFt = ion_ref.dFt;

		A1 = ion_ref.A1; A2 = ion_ref.A2; A3 = ion_ref.A3; A4 = ion_ref.A4;
		A1err = ion_ref.A1err; A2err = ion_ref.A2err;
		dA0 = ion_ref.dA0; dA0err = ion_ref.dA0err;
		dA1 = ion_ref.dA1; dA1err = ion_ref.dA1err;
		cov12 = ion_ref.cov12; cov15 = ion_ref.cov15; cov16 = ion_ref.cov16;
		cov25 = ion_ref.cov25; cov26 = ion_ref.cov26; cov56 = ion_ref.cov56;
		

		//==2017 58Ni
		turn_middle = ion_ref.turn_middle;
		A3err=ion_ref.A3err;
		cov13=ion_ref.cov13;cov23=ion_ref.cov23;cov35=ion_ref.cov35;cov36=ion_ref.cov36;
	}
	double Calculate_M(double Bp_in)  // ,T*0.001 ->converted to :[ns]
	{
		double v_ns=v;
		//double v_err_ns=v_err/1000000000;
		//cout<<fixed<<setprecision(10)<<"v="<<v<<endl<<"v_ns="<<v_ns<<endl;
		double ss=1.0/v_ns/v_ns - 1/V_c/V_c;
		
		return  Z*u*(Bp_in/10.3642686577806)*sqrt(ss);
		
	}
	void Calculate_M_err()    // mass and mass_err and mvq
	{	
		//!!! 如何从国际单位的数据得到 质量 (实则为能量 单位 keV)
		///// E[keV] = 1000000* Z * Bp[Tm] * Vc[m/ns]*Vc * sqrt( 1/(v[m/ns]*v) - 1/(Vc*Vc) )

		double v_ns=v;
		double v_err_ns=v_err;
		//cout<<fixed<<setprecision(10)<<"v="<<v<<endl<<"v_ns="<<v_ns<<endl;
		double ss=1.0/v_ns/v_ns - 1/V_c/V_c;
		
		//double AA = Z*u/10.3642686577806;
		double AA = Z*V_c*V_c*1000000.0; //20230223 revise
		
		M_cal =  AA*Bp*sqrt(ss);
		Mvq=M_cal/(Z*u);
	//cout<<fixed<<setprecision(10)<<"AA*s = "<<AA*sqrt(ss)<<endl<<" M= "<<M_cal;
		v_err_ns=0; //
	
		//M_cal_err = AA*sqrt( ss*Bp_err2 + pow (Bp*v_err_ns/(v_ns*v_ns*v_ns) , 2) /ss ); 
	///cout<<"$$$$$ ------AA*s = "<<AA*sqrt(ss)<<endl;
	}

	
	void PrintInfo()
	{
	    cout<<" |name= "<<name<<"| "<<fixed<<setprecision(10)<<" |Z = "<<Z<<"| |T = "<<T<<"| |C = "<<C<<"| |v = "<<v<<"| |gammat = "<<gammat<<"| |Bp = "<<Bp<<"|"
	    	<<"| |AME = "<<M_AME<<" +- "<<M_AME_err<<endl;
	}
	void PrintMassCalInfo()
	{
		cout<<fixed<<setprecision(10)<<" |AME = "<<M_AME<<" +- "<<M_AME_err<<" |Mass_calculated = "<<M_cal<<endl
			<<" |###delta mass = "<<M_cal - M_AME<<"|"<<" +- "<<M_cal_err<<endl;
	}

};

class IONSpecies
{

public:
	
	TString name;
	TString name_latex;
	TString Aname;
	TString folder_path;  // ->Print(ionspecies[i].folder_path+"h_C_and_h2_gtC.png");
	int Species;
	int N;
	int N_gammat;
	int N_unknown;
	int A,Z;  double Tz;
	int IAS_ID;   int IAS_n;  int IsomerID;  int Isomer_n; 

	bool HasResult;
	bool MassUnknown;
	bool IsRef;
	
	bool have_Fit_BpC;
	bool have_Fit_mvqC;
	bool have_Fit_mvqC_VE;
	bool gtC_CHOSEN;

	int h_m_Gauss_fit_opt;  // 0:insufficient counts 1: one gauss 2: two gauss 3:....

	double AveT;
	double AveC;
	double AveBp;
	double Avev;
	double Avegammat;
	
	double SumT;
	double SumC;
	double SumBp;
	double Sumv;
	double Sumgammat;
	
	double SigmaT;
	double SigmaC;
	double SigmaBp;
	double Sigmav;
	double Sigmagammat;

	double SkewnessC;

	//20240709 再整理
	double Tmin,Tmax,Cmin,Cmax,vmin,vmax,Bpmin,Bpmax,gtmin,gtmax;
	TH2F* h2_gtC;
	TGraph* gr_gt_C; 
	TH1F* h_C;TH1F* h_T;    TH1F* h_v;TH1F* h_Bp;TH1F* h_gt;

	
	double AME;        //atomic mass excess [keV]
	double AME_err;   //[keV]
	double exc; double exc_err;  // for ISOMER
	double Mass;     // AME total nuclear mass  [keV],  Mass=A*u + AME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
	double Mvq_AME;
	double Mvq_AME_err;
	double HalfLife;  //[ns]
	double HalfLife_err;  //[ns]
	
	//---- v1:average v2:error weighted v3:h1_fit v4: error dependent  --------
	double Mass_cal;    double Mass_cal_err;     //MASS_VER=1,total nuclear mass  [keV]
	double Mass_cal_v2; double Mass_cal_err_v2;  //MASS_VER=2 fit histogram
	double Mass_cal_VE; double Mass_cal_err_VE;  //MASS_VER>=3 err

	double deltaMass;    double deltaMass_err;
	double deltaMass_v2; double deltaMass_err_v2;  
	double deltaMass_VE; double deltaMass_err_VE;   
	
	double MassExcess_cal;
	double MassExcess_cal_v2;
	double MassExcess_cal_VE;//MASS_VER>=4

	double Mvq_cal;   double Mvq_cal_err;
	double Mvq_cal_v2;   double Mvq_cal_err_v2;
	double Mvq_cal_VE;   double Mvq_cal_err_VE;

	double stdDeviation_V1,stdDeviation_VE;
	
	double Mass_cal_err_VE_sca;   // need error bar 对于已有单个离子误差棒的情况下算 scattering error 与average error 进行比较
	double Mass_cal_err_VE_ave;

	//---------------------------------------------------------------------------------------------------
	//0928

	TGraph* gr_BpC ; 
	double k_BpC;
	
	TGraph* gr_lnBpC ;  
	double k_lnBpC;
	TF1* fitfun_pol1_BpC ;
	TF1* fitfun_pol1_lnBpC ;

	//1016
	TGraph* 	  gr_mvqC; 
	TGraphErrors* grerr_mvqC_VE ; 
	TGraph*       gr_dmC ; 
	TGraphErrors* grerr_dmC_VE ; 
	//20240703 

	TGraph* gr_dmdC ;

	TH2F* h2_mvqC;
	
	TF1* f0_mvq_AME;
	
	double k_mvqC;
	TF1* fitfun_pol1_mvqC ;
	double k_mvqC_v2;
	TF1* fitfun_pol1_mvqC_VE ;

	// own gtC curve
	bool Has_own_avegtC;
	//int C_Division_n[100];double avegt   [100];double avegt2  [100];double sigma_gt[100];  // 换成外部函数 BuildAvegtCCurve_each()
	TGraphErrors* grerr_avegtC ;
	int count_avegtC_n;
	//20240702 use own gtC
	TGraph* gr_gtC_own;TGraph* gr_gtC_own_u;TGraph* gr_gtC_own_d;
	TGraph* gr_gtC_shifted_own;TGraph* gr_gtC_shifted_own_u;TGraph* gr_gtC_shifted_own_d;
	double GT_SHIFT_own;
	
	
	//===================20230518 dA0_T
	TGraphErrors* grerr_dA0_T; int grerr_dA0_T_n;
	TGraph* gr_T_C; int gr_T_C_n;
	double k_dA0_T;double k_dA0_T_err;
	double dA0TfitMIN,dA0TfitMAX;
	bool Has_k_dA0_T;
	TF1* ff_dA0_T_pol1;TF1* ff_dA0_T_pol6;
	//________________________________
	//20230526
	TGraphErrors* grerr_T_t; 	int grerr_T_t_n;   // !LOOP
	TGraphErrors* grr_Bp_time; TGraphErrors* grr_Bp_inj;
	TGraphErrors* grr_mvq_time;

	int Cfilter_n;    //!LOOP
	double Cfilter_aveC;double Cfilter_avegt;double Cfilter_aveCstd;double Cfilter_avegtstd;
	double Cfilter_avev;double Cfilter_aveT;double Cfilter_avevstd;double Cfilter_aveTstd;

	//20230811
	TH1F* h_mvq;   
	TH1F* h_dm;  // nuclear mass-AME [keV]
	TH1F* h_mvq_FD; // F-D bin determined after h_dm
	TF1* fitfun_gaus_dm;
	TF1* fitfun_gaus_mvq;
	// fitting function Gauss: 3 parameters: (A, mu, sigma)
	double fitfun_gaus_dm_A;double fitfun_gaus_dm_mu;double fitfun_gaus_dm_sigma;
	double fitfun_gaus_dm_A_err;double fitfun_gaus_dm_mu_err;double fitfun_gaus_dm_sigma_err;

	//---- Tfix
	TH1F* h_Tfix0;   
	TH1F* h_Tfix1;
	TH1F* h_Tfix2;
	TH1F* h_Tfix3; TF1* fitfun_h_Tfix3;   
	std::vector<double> vector_Tfix3_v3;//需要正确的Bp Cfix才能有正确中心值，这里不需要，只算std
	TGraph* gr_TC0; int gr_TC0_n;
	TGraph* gr_TC1; int gr_TC1_n;
	TGraph* gr_TC2; int gr_TC2_n;
	TGraph* gr_TC3; int gr_TC3_n;
	double sigma_Tfix0,sigma_Tfix1,sigma_Tfix2,sigma_Tfix3;
	double sigma_Tfix3_v3;//需要正确的Bp Cfix才能有正确中心值，这里不需要，只算std
	//---- 20250210 m0 from T0
	TH1F* h_m0_v1;
	TH1F* h_m0_v2;
	TH1F* h_m0_v3; std::vector<double> vector_m0_v3;
	TH1F* h_m1;
	double stddev_m0_v1,stddev_m0_v2,stddev_m0_v3;
	int m0_n;


	//20231030 err ana
	TH1F* h_iont_mass_err;
	TH1F* h_each_ref_cal_mass_err;
	TH1F* h_iont_chi_n;

	//20240908 for IQR, F-D bin
	vector<double> vector_dm_v1;
	double F_D_binIwidth ;

	//20240710
	ofstream outfile_mass;

	IONSpecies():N(0),N_gammat(0), N_unknown(0), 
	IAS_ID(0),HasResult(0),MassUnknown(0),have_Fit_BpC(0),have_Fit_mvqC(0),IsRef(0),h_m_Gauss_fit_opt(0),
	Mvq_cal(0),Mvq_cal_err(0),
	AveT(0),AveC(0),Avev(0),AveBp(0),Avegammat(0),
	SumT(0),SumC(0),Sumv(0),SumBp(0),Sumgammat(0),
	SigmaT(0),SigmaC(0),Sigmav(0),SigmaBp(0),Sigmagammat(0),SkewnessC(0),
	AME(0),AME_err(0),HalfLife(0),HalfLife_err(0),
	Mass(0),Mass_cal(0),Mass_cal_err(0),Mass_cal_v2(0),Mass_cal_err_v2(0),stdDeviation_V1(0),stdDeviation_VE(0),
	deltaMass(0),deltaMass_VE(0),deltaMass_err(0),deltaMass_err_v2(0),deltaMass_err_VE(0)
	{
		// 所有 TGraph TH1成员变量的new 操作在主程序 循环之前进行 
		GT_SHIFT_own=0;
		Cfilter_aveC=0;Cfilter_avegt=0;Cfilter_aveCstd=0;Cfilter_avegtstd=0;Cfilter_n=0;
		Cfilter_avev=0;Cfilter_aveT=0;Cfilter_avevstd=0; Cfilter_aveTstd=0;

		Has_k_dA0_T = false;dA0TfitMIN=-999;dA0TfitMAX=999;
		gtC_CHOSEN = false;
		Has_own_avegtC = false;
		m0_n = 0;
		
	};

	void RESET_in_LOOP()
	{
		HasResult = 0;IsRef=0;
		AveT=0;AveC=0;AveBp=0;Avev=0;Avegammat=0;
		SumT=0;SumC=0;SumBp=0;Sumv=0;Sumgammat=0;
		SigmaT=0;SigmaC=0;SigmaBp=0;Sigmav=0;Sigmagammat=0;
		SkewnessC=0;
		//Tmin,Tmax,Cmin,Cmax,vmin,vmax,Bpmin,Bpmax,gtmin,gtmax; 使用前会重置

		Mass_cal=0;    Mass_cal_err=0;     //MASS_VER=1,total nuclear mass  [keV]
		Mass_cal_v2=0; Mass_cal_err_v2=0;  //MASS_VER=2 fit histogram
		Mass_cal_VE=0; Mass_cal_err_VE=0;  //MASS_VER>=3 err
	
		deltaMass=0;    deltaMass_err=0;
		deltaMass_v2=0; deltaMass_err_v2=0;  
		deltaMass_VE=0; deltaMass_err_VE=0;   
		
		MassExcess_cal=0;
		MassExcess_cal_v2=0;
		MassExcess_cal_VE=0;
	
		Mvq_cal=0;   Mvq_cal_err=0;
		Mvq_cal_v2=0;Mvq_cal_err_v2=0;
		Mvq_cal_VE=0;Mvq_cal_err_VE=0;
	
		stdDeviation_V1=0;       stdDeviation_VE=0;
		Mass_cal_err_VE_sca=0;   Mass_cal_err_VE_ave=0;

		count_avegtC_n = 0;
	}
	
	//Set
	void DealWithAname_in(string Aname_in)
	{
		int name_start = 0;
		string str_A;
		for (int i = 0; i < Aname_in.length(); ++i)
		{
			if( 48<=Aname_in[i]&&Aname_in[i]<=57)name_start++;
			else if( (97<=Aname_in[i]&&Aname_in[i]<=122)|| (65<=Aname_in[i]&&Aname_in[i]<=90) )
			{
				str_A = Aname_in.substr(0,name_start);
				A= atoi(str_A.c_str());
				name = Aname_in.substr(name_start);
				Z = convert_name_to_z(name);
				break;
			}
		}
	}
	
	void SetTz()
	{
		Tz= double(A)/2-double(Z);
	}
	
	void SetName_form()
	{
		
		Aname = name.Format("%d",A)+name;	
		name_latex=name.Format("^{%d}",A) + name;
	}
	void SetName_form_exc(TString info_tmp)
	{
		Aname = name.Format("%d",A)+name+"_"+ info_tmp;
		name_latex=name.Format("^{%d}",A) + name;
			
	}
	void SetMass()  // total nuclear mass
	{
		//Mass = A*u+AME; // total Atomic Mass
		double Me=510.9989;
		Mass=A*u + AME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
	}
	double GetMassExcess_cal(double nuclear_mass)  // from total nuclear mass to atomic mass excess
	{
		//Mass = A*u+AME;
		double Me=510.9989;
		return  nuclear_mass - A*u + Z*Me - ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
	}
	void DoSetMass_v2(double fit_mu,double dm_fit,double m_fit_err)
	{
		deltaMass_v2 = dm_fit;
    	Mass_cal_v2  = fit_mu*Z*u;
    	MassExcess_cal_v2 = GetMassExcess_cal(Mass_cal_v2);
    	Mass_cal_err_v2 = m_fit_err;
    	deltaMass_err_v2 = sqrt( pow(Mass_cal_err_v2,2)+pow(AME_err,2)  );
	}
	void SetMvq_AME()
	{
		Mvq_AME = Mass/(u*Z);
	}
	//20230606
	double Calculate_only_return_Bp(double vns, double m_ref)
	{
		return m_ref/Z*0.000001/(V_c*V_c)/sqrt( (1/(vns*vns)) - (1/(V_c*V_c)) );
	}
	double Calculate_Bp_from_L(double L, double Ts, double td)
	{
		double v_ns = L/(Ts+td);  // m/ns
		return Calculate_only_return_Bp(v_ns, Mass);

	}
	//0928
	double Get_k_lnBpC()
	{
		TF1* fitfun_pol1 = new TF1("fitfun_pol1","pol1",4,5);
		gr_lnBpC->Fit(fitfun_pol1,"q");
		k_lnBpC = fitfun_pol1->GetParameter(1);
		delete fitfun_pol1;
		return k_lnBpC;
	}
	double Get_k_BpC()
	{
		TF1*fitfun_pol1 = new TF1("fitfun_pol1","pol1",4,5);
		gr_BpC->Fit(fitfun_pol1,"q");
		k_BpC = fitfun_pol1->GetParameter(1);
		delete fitfun_pol1;
		return k_BpC;
	}
	void Fit_BpC(int opt)
	{
		if(opt==0)
		{
			gr_BpC->Fit(fitfun_pol1_BpC,"q");
			k_BpC = fitfun_pol1_BpC->GetParameter(1);
		}
		else if(opt==1)
		{
			gr_lnBpC->Fit(fitfun_pol1_lnBpC,"q");
			k_lnBpC = fitfun_pol1_lnBpC->GetParameter(1);	
		}
	}
	void Fit_mvqC(int v)
	{
		if(v==1)
		{
		gr_mvqC->Fit(fitfun_pol1_mvqC,"q");			
		k_mvqC = fitfun_pol1_mvqC->GetParameter(1);
		}
		else if(v>1)
		{
		grerr_mvqC_VE->Fit(fitfun_pol1_mvqC_VE,"q");			
		k_mvqC_v2 = fitfun_pol1_mvqC->GetParameter(1);
		}
	}
	void Fit_dA0_T(int poln)
	{
		if(poln==1)
		{
			ff_dA0_T_pol1 = new TF1("ff_dA0_T_pol1","pol1",500,700);
			grerr_dA0_T->Fit(ff_dA0_T_pol1,"q");
			k_dA0_T = ff_dA0_T_pol1->GetParameter(1);
			k_dA0_T_err = ff_dA0_T_pol1->GetParErrors()[1];

		}
		else if(poln==6)
		{
			ff_dA0_T_pol6 = new TF1("ff_dA0_T_pol6","pol6",500,700);
			grerr_dA0_T->Fit(ff_dA0_T_pol6,"q");
		}
		else{cout<<"error! poln of fit dA0 T!!!"<<endl;}
		

	}
	void Show_h_dm_fit_gaus_paras()
	{
		cout<<Aname<<" gaus fit 3 paras :"<<endl;
		cout<<" A = "<<fitfun_gaus_dm_A<<" +- "<<fitfun_gaus_dm_A_err<<" "
			<<" mu = "<<fitfun_gaus_dm_mu<<" +- "<<fitfun_gaus_dm_mu_err<<" "
			<<" sigma = "<<fitfun_gaus_dm_sigma<<" +- "<<fitfun_gaus_dm_sigma_err<<endl;
	}
	void Show_h_dm_fit_gaus_paras(ofstream& outfile)
	{
		outfile<<Aname<<" gaus fit 3 paras :"<<endl;
		outfile<<" A = "<<fitfun_gaus_dm_A<<" +- "<<fitfun_gaus_dm_A_err<<" "
			<<" mu = "<<fitfun_gaus_dm_mu<<" +- "<<fitfun_gaus_dm_mu_err<<" "
			<<" sigma = "<<fitfun_gaus_dm_sigma<<" +- "<<fitfun_gaus_dm_sigma_err<<endl;
	}
	//0930---- double[subregion_n]
	/*
	void ResetGtArrays(int subregion_n)
	{
		if(subregion_n>100){cout<<" error!!! in IONSpecies.ResetGtArrays() subregion_n > 100"<<endl;return;}
		for(int i=0;i<100  ;i++)
		{
			C_Division_n[i]=0;
		    avegt[i]=0;
            avegt2[i]=0;
            sigma_gt[i]=0;
		}
	}
	void Fill_ion_gt(const ION& a ,int subregion_n,double gt_err_upper_limit)
	{
		if(subregion_n>100){cout<<" error!!! in IONSpecies.Fill_ion_gt() subregion_n > 100"<<endl;return;}
		if(a.gammat_err<=gt_err_upper_limit)
		{
			if(a.C_region>=0&&a.C_region<=subregion_n-1)
			{
				C_Division_n[a.C_region]++;
				avegt[a.C_region]     += a.gammat;
				avegt2[a.C_region]    += a.gammat*a.gammat;
				count_avegtC_n ++;
			}
			
		}
	}
	void BuildAvegtCCurve(int subregion_n, const TH1F* h1 )
	{
		for(int i=0;i<subregion_n;i++)
		{
			if(C_Division_n[i]>3)
        	{
        	    avegt[i] /= C_Division_n[i];
        	    avegt2[i] /= C_Division_n[i];
        	    sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;

        	    grerr_avegtC->SetPoint( grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
        	    grerr_avegtC->SetPointError(grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;
        	    
        	    gr_gtC_own->SetPoint( gr_gtC_own->GetN(),h1->GetBinCenter(i+1), avegt[i] );  // 20240702

        	    //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
        	}    
		}
		Has_own_avegtC = true;
	}
	*/
	//20240702
	void Create_gr_gtC_shifted_own(TGraph* gr_gtC_chosen,TGraph* gr_gtC_chosen_u,TGraph* gr_gtC_chosen_d)
	{
		double xtmp,ytmp=0;
		for(int i=0;i<gr_gtC_chosen->GetN() ; i++)
		{
			gr_gtC_chosen->GetPoint(i,xtmp,ytmp);
			gr_gtC_shifted_own->SetPoint(i,xtmp,ytmp+GT_SHIFT_own);
			gr_gtC_chosen_u->GetPoint(i,xtmp,ytmp);
			gr_gtC_shifted_own_u->SetPoint(i,xtmp,ytmp+GT_SHIFT_own);
			gr_gtC_chosen_d->GetPoint(i,xtmp,ytmp);
			gr_gtC_shifted_own_d->SetPoint(i,xtmp,ytmp+GT_SHIFT_own);
		}
	}
	void operator=(const IONSpecies& ion_ref)
	{
		name = ion_ref.name;
		IAS_ID = ion_ref.IAS_ID;
		A = ion_ref.A;
		Z = ion_ref.Z;
	}
	void Record(const ION& a)
	{
		N++;
		SumT+=a.T;SumC+=a.C;Sumv+=a.v;SumBp+=a.Bp;
		if(a.gammat>1.3&&a.gammat<1.4){ Sumgammat+=a.gammat;N_gammat++; }
	}
	void Record_ion_BpC(const ION& a)
	{
		//0928
		gr_BpC->SetPoint(gr_BpC->GetN(),a.C,a.Bp);
		gr_lnBpC->SetPoint(gr_lnBpC->GetN(),log(a.C),log(a.Bp));
	}
	
	//Calculate
	void CalculateAve()
	{
		AveT = SumT / N;AveC =SumC / N;AveBp =SumBp / N;Avev = Sumv / N;
		Avegammat =Sumgammat / N_gammat;
	}
	void SigmaPreAdd(const ION& a)
	{	//
		SigmaT += pow( (a.T-AveT),2 );SigmaC += pow( (a.C-AveC),2 );
		Sigmav += pow( (a.v-Avev),2 );SigmaBp += pow((a.Bp-AveBp),2);
		Sigmagammat += pow( (a.gammat-Avegammat),2 );
	}
	void CalculateSigma()
	{   //use after SigmaPreAdd of all ions is done
		if(N>1)
		{
			SigmaT = sqrt(SigmaT/(N-1)) ;SigmaC =sqrt(SigmaC/(N-1));SigmaBp = sqrt(SigmaBp/(N-1)); 
			Sigmav = sqrt(Sigmav/(N-1));Sigmagammat = sqrt(Sigmagammat/(N_gammat-1));
		}		
	}
	void CalculateSkewness()
	{
		SkewnessC=SkewnessC/(N*SigmaC*SigmaC*SigmaC);
	}
	//record
	
	void ClearRecord()
	{
		N=0;N_gammat=0;
		SumT=0;SumC=0;Sumv=0;SumBp=0;Sumgammat=0;
		AveT = 0;AveC =0;Avev = 0;AveBp =0;Avegammat =0;
		SigmaT =0;SigmaC =0;Sigmav = 0;SigmaBp =0;Sigmagammat = 0;
		
	}
	//Print
	
	void SkewnessPreAdd(const ION& a)
	{
		SkewnessC += pow( (a.C-AveC),3 );
	}
	void PrintInfo()
	{
		cout<<"-------------------------------------------------------------------"<<endl;
	    cout<<"\033[33m|species No. "<<Species<<endl<<"|name= "<<name<<"| "<<"|A= "<<A<<" |Z= "<<Z<<"\033[0m|N = "<<N<<fixed<<setprecision(7)<<endl
	    <<" |AveT = "<<AveT<<" |SigmaT = "<<SigmaT<<"| |Avev = "<<Avev<<"| |Sigmav = "<<Sigmav<<"| |AveC = "<<AveC<<"| |SigmaC = "<<SigmaC<<endl
	    <<"| |AveBp = "<<AveBp<<"| |SigmaBp = "<<SigmaBp<<"| |Avegammat = "<<Avegammat<<"| "<<"| |Sigmagammat = "<<Sigmagammat<<"|"<<endl
	    <<" |AME = "<<AME<<" +- "<<AME_err<<" MassUnknown: "<<MassUnknown<<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
	void PrintInfo(ofstream& outfile)
	{
		
	    outfile<<"|speciesNo. "<<Species<<" |element= "<<name<<" "<<"|A= "<<A<<" |Z= "<<Z<<" |count = "<<N<<fixed<<setprecision(7)
	    <<" |AveT= "<<AveT<<" |SigmaT= "<<SigmaT<<" |Avev= "<<Avev<<" |Sigmav= "<<Sigmav<<" |AveC= "<<AveC<<" |SigmaC= "<<SigmaC
	    <<" |AveBp= "<<AveBp<<" |SigmaBp= "<<SigmaBp<<" |Avegammat= "<<Avegammat<<" |Sigmagammat= "<<Sigmagammat
	    <<fixed<<setprecision(1)<<" |AME_ME= "<<AME<<" |AME_error= "<<AME_err<<" |AME_nuc_Mass= "<<Mass
	    <<fixed<<setprecision(10)<<" |AME_mvq= "<<Mvq_AME    
	    <<" IsMassUnknown: "<<MassUnknown<<" IsRef: "<<IsRef<<" gtC_CHOSEN: "<<gtC_CHOSEN
	    <<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
	void PrintInfo2()
	{
		cout<<"-------------------------------------------------------------------"<<endl;
	    cout<<"|name= "<<name<<"| "<<"|species No. "<<Species<<endl<<"|N = "<<N<<fixed<<setprecision(7)<<" |SigmaT = "<<SigmaT<<"| |Sigmav = "<<Sigmav<<"| |Sigmagammat = "<<Sigmagammat<<"|"<<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
};

class NUBASE_IONSpecies
{

public:
	
	TString name;
	TString name_latex;
	TString Aname;
	TString folder_path;
	TString nuclide_info;
	
	TString Jpi_AME20;
	TString DiscoverYear;
	TString ProductionMethod;   //20230115 
	TString DecayMode_str;   //20230117 

	string str_line;   //保留原始整行信息
	string HalfLife_str;
	string HalfLifeUnit_str;
	string HalfLifeErr_str;
	string br_info_str;
	string AME_str;
	string AME_err_str;

	int Species;

	int N_gammat;
	int N_unknown;
	int A,Z,N;
	int state_i;  //AME20 str[7] 的一位数字 i=0 (gs); i=1,2 (isomers); i=3,4 (levels); i=5 (resonance); i=8,9 (IAS)
                  //i=3,4,5,6 can also indicate isomers (when more than two isomers are presented in a nuclide)
	int DiscoverYear_int;
	int ProductionMethod_id;
	int method_classified_id;	
	bool HasIAS;
	bool HasResult;
	bool MassUnknown;
	bool IsStable;
	bool IsIS;

	
	double Mvq;
	double Mvq_AME;
	double Tz;
	double T_isospin;

	//double Sn,Sp,S2n,S2p,Sn_err,Sp_err,S2n_err,S2p_err;
	//int Sn_info,Sp_info,S2n_info,S2p_info;  // info: 0=ini, 1=exp, 2=# systematic ,3=* no value

	double AME;        //unit kev mass excess
	double AME_err;   //unit kev
	double exc;        //unit kev excited state mass excess
	double exc_err;   //unit kev
	double HalfLife;  //unit ns
	double HalfLife_err;  //unit ns
	double Mass;     // total nuclear mass  [keV]
	double AtomicMass;//total atomic mass [keV]
	double BE;       // nuclear binding energy 
	double AtomicBE; // AME2020 DEFINITION 
	double Mass_cal;
	double MassExcess_cal;
	
	double stdDeviation;

	//vector<int> vector_exc_state;
	int next_nuclide_i; //
	NUBASE_IONSpecies():N(0),N_gammat(0), N_unknown(0), HasIAS(0),HasResult(0),MassUnknown(0),IsStable(0),
	Mvq(0),AME(0),AME_err(0),HalfLife(0),HalfLife_err(0),Mass(0),Mass_cal(0),AtomicMass(0),
	//Sn(0),Sp(0),S2n(0),S2p(0),Sn_err(0),Sp_err(0),S2n_err(0),S2p_err(0),Sn_info(0),Sp_info(0),S2n_info(0),S2p_info(0),
	DiscoverYear_int(0),ProductionMethod_id(0),method_classified_id(0)
	{
		
	};
	
	
	//Set
	void DealWithAname_in(string Aname_in)
	{
		int name_start = 0;
		string str_A;
		for (int i = 0; i < Aname_in.length(); ++i)
		{
			if( 48<=Aname_in[i]&&Aname_in[i]<=57)name_start++;
			else if( (97<=Aname_in[i]&&Aname_in[i]<=122)|| (65<=Aname_in[i]&&Aname_in[i]<=90) )
			{
				str_A = Aname_in.substr(0,name_start);
				A= atoi(str_A.c_str());
				name = Aname_in.substr(name_start);
				Z = convert_name_to_z(name);
				//cout<<"debug in DealWithAname_in "<<A<<" "<<name<<" "<<Z<<endl;
				break;

			}
		}
		SetName_form();

	}
	
	void SetTz()
	{
		Tz= double(A)/2-double(Z);
	}
	void SetIsospinT()
	{
		T_isospin = abs(Tz);
		if(nuclide_info=='i')T_isospin+=1;
		else if(nuclide_info=='j')T_isospin+=2;
	}
	
	void SetName_form()
	{
		name_latex=name.Format("^{%d}",A) + convert_z_to_name(Z);
		Aname = name.Format("%d",A)+name;	
		while(Aname.Index(" ")!=-1){Aname.Remove(Aname.Index(" "),1);}
	}
	void SetMass()  // total nuclear mass
	{
		//Mass = A*u+AME;
		double Me=510.9989;
		Mass=A*u + AME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
		AtomicMass = A*u + AME;

	}
	void GetMassExcess_cal()  // from total nuclear mass to atomic mass excess
	{
		//Mass = A*u+AME;
		double Me=510.9989;
		MassExcess_cal = Mass - A*u + Z*Me - ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;

	}
	void SetMvq_AME()
	{
		Mvq_AME = Mass/(u*Z);
	}
	void SetBindingEnergy() // after setting Nuc Mass
	{
		BE = Z* mass_proton + N* mass_neutron -  Mass;
		AtomicBE = Z*mass_Hydrogen+N*mass_neutron - AtomicMass;
	}
	double GetME_from_BE(double AtomicBE_in)
	{
		return Z*mass_Hydrogen+N*mass_neutron- AtomicBE_in -A*u;
	}
	bool Is_out_of_NucChart(int Z,int N, int* isotope_L, int*isotope_R)
	{
		if(Z<0||N<0){ cout<<"error input Z<0||N<0 in func: If_out_of_NucChart"<<endl; return true;}
		if(N>isotope_R[Z]||N<isotope_L[Z]){return false;}
		else return true;
	}
	//////////////////////////////////////////////////////////////
	void operator=(const NUBASE_IONSpecies& ion_ref)
	{
		name = ion_ref.name;
		HasIAS = ion_ref.HasIAS;
		A = ion_ref.A;
		Z = ion_ref.Z;
	}

	//Print	
	void PrintInfo()
	{
		cout<<"-------------------------------------------------------------------"<<endl;
	    cout<<"\033[33m |species No. "<<Species<<endl<<" |name= "<<name<<"| "<<" |Aname= "<<Aname<<"| "<<" state: "<<nuclide_info
	    <<" |A= "<<A<<" |Z= "<<Z<<"\033[0m |N = "<<N<<fixed<<setprecision(1)<<endl
	    	<<" |AME = "<<AME<<" +- "<<AME_err<<" MassUnknown: "<<MassUnknown<<endl
	    	<<" |nucmass= "<<Mass<<" |AtomicMass= "<<AtomicMass<<" |BE= "<<BE<<endl
	    	<<" |next_nuclide_i = "<<next_nuclide_i<<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
	void PrintInfo(ofstream & outfile)
	{
		outfile<<"-------------------------------------------------------------------"<<endl;
	    outfile<<" |species No. "<<Species<<endl<<" |name= "<<name<<"| "<<" |Aname= "<<Aname<<"| "<<" state: "<<nuclide_info
	    <<" |A= "<<A<<" |Z= "<<Z<<fixed<<setprecision(1)<<" Tz= "<<Tz<<" T_isospin= "<<T_isospin<<endl
	    <<fixed<<setprecision(1)<<" |AME = "<<AME<<" +- "<<AME_err<<" MassUnknown: "<<MassUnknown<<endl
	    <<" |nucmass= "<<Mass<<" |AtomicMass= "<<AtomicMass<<" |BE= "<<BE<<endl
	    <<" |next_nuclide_i = "<<next_nuclide_i<<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
	
};


//////////////////////////////////////////////////////////////////////////////////////////
//----------- 2024 0913 Integrated 整合了NUBASE 2020 文件读入所有 5843 种 核素信息（包括 IAS ij Isomer m,n)

void NUBASE2020_Build(NUBASE_IONSpecies* ISS,TString NUBASE_FILE)
{

bool AMEreadin_check_ON 		=0 ;
bool OUTPUT_check_ON 			=0 ;


gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
TString strtmp;
ifstream infile;
ofstream outfile;

TLatex* lat_n = new TLatex();   
TString lat_text;
lat_n->SetTextColor(kAzure+7);
lat_n->SetTextFont(43);
lat_n->SetTextSize(20);
lat_n->SetTextAlign(11);
for(int i=0;i<ZN_MAX_Z;i++)
	for(int j=0;j<ZN_MAX_N;j++)
		{ZN_AME[i][j]=-1;}
// ============================================== AME DATA FORMAT =========================================
//# Nubase2020 citation: F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, and G. Audi, Chin. Phys. C45, 030001 (2021)
//#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//#                                     Data FORMAT 
//#     column   quantity   format      description
//#       1: 3   AAA           a3       Mass Number (AAA)
//#       5: 8   ZZZi          a4       Atomic Number (ZZZ); i=0 (gs); i=1,2 (isomers); i=3,4 (levels); i=5 (resonance); i=8,9 (IAS)
//#                                     i=3,4,5,6 can also indicate isomers (when more than two isomers are presented in a nuclide)
//#     12: 16   A El          a5       A Element 
//#     17: 17   s             a1       s=m,n (isomers); s=p,q (levels); s=r (reonance); s=i,j (IAS); 
//#                                     s=p,q,r,x can also indicate isomers (when more than two isomers are presented in a nuclide)
//#     19: 31   Mass #     f13.6       Mass Excess in keV (# from systematics)
//#     32: 42   dMass #    f11.6       Mass Excess uncertainty in keV (# from systematics)
//#     43: 54   Exc #      f12.6       Isomer Excitation Energy in keV (# from systematics)
//#     55: 65   dE #       f11.6       Isomer Excitation Energy uncertainty in keV (# from systematics)
//#     66: 67   Orig          a2       Origin of Excitation Energy  
//#     68: 68   Isom.Unc      a1       Isom.Unc = *  (gs and isomer ordering is uncertain) 
//#     69: 69   Isom.Inv      a1       Isom.Inv = &  (the ordering of gs and isomer is reversed compared to ENSDF) 
//#     70: 78   T #         f9.4       Half-life (# from systematics); stbl=stable; p-unst=particle unstable
//#     79: 80   unit T        a2       Half-life unit 
//#     82: 88   dT            a7       Half-life uncertainty 
//#     89:102   Jpi */#/T=    a14      Spin and Parity (* directly measured; # from systematics; T=isospin) 
//#    103:104   Ensdf year    a2       Ensdf update year 
//#    115:118   Discovery     a4       Year of Discovery 
//#    120:209   BR            a90      Decay Modes and their Intensities and Uncertanties in %; IS = Isotopic Abundance in %
//#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//



//================================================ read in AME2020 ========================================================


int str_line_len_min=9999;int str_line_len_max=-1;int str_line_len=0;
string str_Aname_in;
string str_mass_in;
string str_dmass_in;
string str_exc_in;
string str_dE_in;
string str_orig_in,str_Iso_Unc_in,str_Iso_Inv,str_T_in,str_T_Unit_in,str_dT_in,str_Jpi_in;
string str_ensdfyear_in,str_discover_in,str_br_in;
ofstream outfile_str_check;

if(AMEreadin_check_ON)
{
	outfile_str_check.open("OUTPUT//AMEreadin_check.txt");
	cout<<"----NUBASE2020_Build:   ---------- AMEreadin_check_ON ----------------"<<endl;
}
/////////////////////////////////////////
int n=0;
bool str_end_space_remove = 0;

infile.open(NUBASE_FILE);   //NUBASE V4数据
if(!infile){cerr<<"error !! !infile ";exit(1);}
if(OUTPUT_check_ON){ outfile.open("output//AME_data_check.txt");   cout<<"----NUBASE2020_Build:   --------- output_check on ------------"<<endl;}
string str_line;

//nubase_3.mas20.txt , ignore 25 lines at beginning
for(int i=0;i<25;i++){getline(infile,str_line);}
//=========

while(getline(infile ,str_line))
{

	/*
	if(n<100)cout<<str_line.substr(0,8)<<" "<< atoi((str_line.substr(0,3)).c_str())<<" "
		         <<atoi((str_line.substr(4,3)).c_str())<<" "<<atoi((str_line.substr(7,1)).c_str())<<endl;
	*/
	//outfile<<str_line<<endl;
	/// 带上尾部空格 len min=75  max=217

	//去除尾部空格
	str_line.erase(str_line.find_last_not_of(" ") + 1); 
	if(!str_end_space_remove)str_end_space_remove=1;
	
	str_line_len = str_line.length();
	if(str_line_len>str_line_len_max)str_line_len_max = str_line_len;
	if(str_line_len<str_line_len_min)str_line_len_min = str_line_len;
	/// 去除尾部空格 len min=58  max=195

	ISS[n].str_line = str_line;
	ISS[n].Species = n;
    ISS[n].A = atoi((str_line.substr(0,3)).c_str());
    ISS[n].Z = atoi((str_line.substr(4,3)).c_str());
    ISS[n].N = ISS[n].A-ISS[n].Z;
    ISS[n].state_i = atoi((str_line.substr(7,1)).c_str());
    ISS[n].SetTz();

    
    str_Aname_in = str_line.substr(11,5);
    //cout<<"len= "<<str_line.length()<<" Aname_in "<<str_Aname_in<<endl;
    str_mass_in  = str_line.substr(18,13);
    str_dmass_in = str_line.substr(31,11);
    if(str_line_len>54){str_exc_in       = str_line.substr(42,12);
    if(str_line_len>65){str_dE_in        = str_line.substr(54,11);
    if(str_line_len>67){str_orig_in      = str_line.substr(65,2);
    if(str_line_len>68){str_Iso_Unc_in   = str_line.substr(67,1);
    if(str_line_len>69){str_Iso_Inv      = str_line.substr(68,1);
    if(str_line_len>78){str_T_in         = str_line.substr(69,9);  
    if(str_line_len>80){str_T_Unit_in    = str_line.substr(78,2);    
    if(str_line_len>88){str_dT_in        = str_line.substr(81,7);           
    if(str_line_len>102){str_Jpi_in       = str_line.substr(88,14);         
    if(str_line_len>104){str_ensdfyear_in = str_line.substr(102,2);       
    if(str_line_len>118){str_discover_in  = str_line.substr(114,4);         
    if(str_line_len>120){str_br_in        = str_line.substr(119);}  
	}
	//else if(str_line_len>=115){str_discover_in  = str_line.substr(114);}对于后有空格的处理
	else{str_discover_in   = str_line.substr(114);}
	}else{str_ensdfyear_in = str_line.substr(102);}
	}else{str_Jpi_in       = str_line.substr(88);}
	}else{str_dT_in        = str_line.substr(81);}
	}else{str_T_Unit_in    = str_line.substr(78);}
    }else{str_T_in         = str_line.substr(69);}
	}else{str_Iso_Inv      = str_line.substr(68);}
	}else{str_Iso_Unc_in   = str_line.substr(67);}
	}else{str_orig_in      = str_line.substr(65);}
	}else{str_dE_in        = str_line.substr(54);}
	}else{str_exc_in       = str_line.substr(42);}
    
    //if(n>4060)cout<<"now "<<n;
    if(AMEreadin_check_ON)
    {//检查读入的字符串    	
    	outfile_str_check<<str_Aname_in<<" " 
    	<<str_mass_in      <<" "<<" | "
    	<<str_dmass_in     <<" "<<" | "
    	<<str_exc_in       <<" "<<" | "
    	<<str_dE_in        <<" "<<" | "
    	<<str_orig_in      <<" "<<" | "
    	<<str_Iso_Unc_in   <<" "<<" | "
    	<<str_Iso_Inv      <<" "<<" | "
    	<<str_T_in         <<" "<<" | "
    	<<str_T_Unit_in    <<" "<<" | "
    	<<str_dT_in        <<" "<<" | "     
    	<<str_Jpi_in       <<" "<<" | "    
    	<<str_ensdfyear_in <<" "<<" | "  
    	<<str_discover_in  <<" "<<" | "    
    	<<str_br_in        <<" "<<endl;

    }

    ISS[n].name = convert_z_to_name(ISS[n].Z);
    ISS[n].SetName_form();
    while(str_Aname_in.find(" ")!=string::npos)str_Aname_in.erase(str_Aname_in.find(" "),1);
    //ISS[n].DealWithAname_in(str_Aname_in);
    if( !(ISS[n].Aname==str_Aname_in)){cout<<" Aname =|"<<ISS[n].Aname<<"|"<<str_Aname_in<<"|"<<endl;cerr<<"error !! !(ISS[n].Aname==str_Aname_in)";exit(1);}
    
    ISS[n].nuclide_info = str_line.substr(16,1);
    ISS[n].SetIsospinT();//需要用到nuclide_info

    string str_Aname_tmp;
    if(ISS[n].state_i==0)str_Aname_tmp = str_Aname_in + str_line.substr(16,1);
    else          		str_Aname_tmp = str_Aname_in + "_" + str_line.substr(16,1);
    while(str_Aname_tmp.find(" ")!=string::npos)str_Aname_tmp.erase(str_Aname_tmp.find(" "),1);   //去除可能的空格
    ISS[n].Aname = str_Aname_tmp;

    ISS[n].AME_str = str_mass_in;
    ISS[n].AME_err_str = str_dmass_in;
    ISS[n].AME = atof(str_mass_in.c_str());
    ISS[n].AME_err = atof(str_dmass_in.c_str());
    if(str_mass_in.find("#")!=string::npos)
    {
    	//cout<<ISS[n].Aname<<" "<<str_mass_in<<" "<< ISS[n].AME<<endl;
    	ISS[n].MassUnknown=true;
    }
    
    ISS[n].SetMass();// nucmass: Mass, AtomicMass
    ISS[n].SetBindingEnergy();
    ISS[n].exc = atof(str_exc_in.c_str());
    ISS[n].exc_err = atof(str_dE_in.c_str());

    if(str_T_in.find("stbl")!=string::npos ){ISS[n].IsStable = 1;}
    else   									{ISS[n].IsStable = 0;}
    ISS[n].HalfLife_str = str_T_in;
	ISS[n].HalfLifeUnit_str=str_T_Unit_in;
	ISS[n].HalfLifeErr_str = str_dT_in;
    ISS[n].Jpi_AME20 = str_Jpi_in;
    ISS[n].DiscoverYear = str_discover_in;
    ISS[n].DiscoverYear_int = atoi(str_discover_in.c_str());

    ISS[n].br_info_str = str_br_in;
    //if(n%34==0)cout<<ISS[n].Aname<<" "<<str_mass_in<<" "<<atof(str_mass_in.c_str())<<endl;

    
    
    if(ISS[n].state_i==0) //ground state
    {
    	ZN_AME[ISS[n].Z][ISS[n].A-ISS[n].Z] = n;
    }
    
    n++;
}
if(str_end_space_remove)cout<<"----NUBASE2020_Build: str end space removed "<<endl;
else cout<<"----NUBASE2020_Build: str end space not removed "<<endl;

cout<<" len max = "<<str_line_len_max<<" len min = "<<str_line_len_min<<endl;
infile.close();
if(AMEreadin_check_ON)
{
	cout<<" check ZN_AME[][]: ISS[ZN_AME[6][6]].str_line= "<< ISS[ZN_AME[6][6]].str_line<<endl;
}
if(AMEreadin_check_ON)outfile_str_check.close();
//_______________________ infile __________________________________

cout<<"----NUBASE2020_Build completed : n of ISS[i] = "<<n<<endl;

for(int i=0;i<n  ;i++)
{
	if(ISS[i].state_i!=0){ISS[i].next_nuclide_i=-1; continue;}
	int j = i;
	while(ISS[j].A==ISS[i].A &&ISS[j].Z==ISS[i].Z)
	{

		j++;
	}
	ISS[i].next_nuclide_i=j;
}

if(OUTPUT_check_ON)
{
	for(int i=0;i<n  ;i++)
	{
		ISS[i].PrintInfo(outfile);
	}
	
}
/*
//====================== ZN_AME =========================
//方格占据--同位素
ZN_n=0;
bool FIND_L = false;
for(int i=0;i<120;i++)
{
	FIND_L = false;
	for(int j=0;j<200;j++)
	{
		if(ZN_AME[i][j]>-1)
		{
			ZN_n++ ;
			if(!FIND_L){ Isotope_L[i]=j;Isotope_R[i]=j;FIND_L=true; }
			else {Isotope_R[i]=j;}	
		}
	}
	if(OUTPUT_check_ON)outfile<<i<<" "<<convert_z_to_name(i)<<" "<<Isotope_L[i]<<" ~ "<<Isotope_R[i]<<endl;
}
cout<<"ZN_n = "<<ZN_n<<endl;

FIND_L = false;
for(int j=0;j<200;j++)
{
	FIND_L = false;
	for(int i=0;i<120;i++)
	{
		if(ZN_AME[i][j]>-1)
		{
			if(!FIND_L){ Isotone_B[j]=i;Isotone_T[j]=i;FIND_L=true; }
			else {Isotone_T[j]=i;}	
		}
	}
	if(OUTPUT_check_ON)outfile<<"N= "<<j<<" "<<Isotone_B[j]<<" ~ "<<Isotone_T[j]<<endl;
}
//_______________________ZN_________________________________
*/
if(OUTPUT_check_ON)outfile.close();

}//__________________________________________________________


