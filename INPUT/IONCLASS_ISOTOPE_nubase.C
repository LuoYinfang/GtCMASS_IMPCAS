////// IMME 
/////   20230113
////


const Double_t Vc=0.2997924580;      //m/ns
const Double_t Vc_ms=299792458.0;      //m/ns
const Double_t u  = 931494.102417;   //KeV
const Double_t mass_neutron = 939565.420 ;   // u+ 8071.318062(440)
const Double_t mass_proton = 938272.089 ;   // u+ 7288.971064(13) - 511
const Double_t mass_Hydrogen = 938783.0734810; //



#define ZN_MAX_Z 120
#define ZN_MAX_N 200
int ZN_n;
int ZN[ZN_MAX_Z][ZN_MAX_N];  //对应于ISS序号--ground state -- NUBASE没有的是 -1



int Isotope_L[200];
int Isotope_R[200];
int Isotone_B[300];
int Isotone_T[300];

double v_to_gamma(double v)
{
	return pow((1-(v*v/Vc_ms/Vc_ms)),-0.5);
}
int convert_name_to_z(TString name)
{
	int tmp=0;
	while(name.Index(" ")!=-1) { name.Remove(name.Index(" "),1); tmp++;if(tmp>1000){cout<<"while loop!! break"<<endl;break;}}


	if(name.CompareTo("n")==0 ) return 0;
	else if(name.CompareTo("H")==0 ) return 1;
	else if(name.CompareTo("He")==0 ) return 2;
	else if(name.CompareTo("Li")==0 ) return 3;
	else if(name.CompareTo("Be")==0 ) return 4;
	else if(name.CompareTo("B")==0 ) return 5;
	else if(name.CompareTo("C")==0 ) return 6;
	else if(name.CompareTo("N")==0 ) return 7;
	else if(name.CompareTo("O")==0 ) return 8;
	else if(name.CompareTo("F")==0 ) return 9;
    
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
	



	else { cout<<"convert_name_to_z failed!!!  --- name = "<<name<<endl;return -1;}

}
TString convert_z_to_name(int z)
{
	switch(z)
	{
		case 0:return"n" ; break;
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
void Aname_to_pn(string Aname_in, int& p, int& n)
{
	TString name;
	int name_start = 0;
	string str_A;
	for (int i = 0; i < Aname_in.length(); ++i)
	{
		if( 48<=Aname_in[i]&&Aname_in[i]<=57)name_start++;
		else if( (97<=Aname_in[i]&&Aname_in[i]<=122)|| (65<=Aname_in[i]&&Aname_in[i]<=90) )
		{
			str_A = Aname_in.substr(0,name_start);
			
			name = Aname_in.substr(name_start);
			
			p = convert_name_to_z(name);
			n= atoi(str_A.c_str()) - p;
			//cout<<"debug in DealWithAname_in "<<A<<" "<<name<<" "<<Z<<endl;
			break;
		}
	}

}

/////////===================================================================
class ION
{
public:
	TString name;
	int ion_number;					//serial number in this injection
	int inject_number;              //appear in which injection
	int Species;  //from 0 to NSpecies
	int A;
	int Z;
	int C_region;
	double T;    //ps
	double T_err;    //ps
	double TOF;
	double TOF_err;
	double dA1;
	double _2A2_6A3;
	double C;    //m
	double C_err2;
	double Bp_true;   //Tm
	double Bp;   //Tm
	double Bp_err2;   //Tm
	double v;    //m/s
	double v_err;    //m/s
	//double v_ns;    //m/ns
	//double v_err_ns;    //m/ns
	double gammat;
	bool gt_abnormal;

	//Get
	ION():T(0),C(0),Bp(0),v(0),gammat(0),gt_abnormal(0){};

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
	void Calculate_Bp(double m)  //m: keV v:m/s
	{
		Bp = (m/Z)*1000.0*v/(Vc_ms*sqrt(Vc_ms*Vc_ms-v*v));
		Bp_true=Bp;         //record Bp calculated using AME
		
	}
	void operator=(const ION& ion_ref)
	{
		name = ion_ref.name;
		ion_number = ion_ref.ion_number;   
		inject_number = ion_ref.inject_number;
		Species = ion_ref.Species;  //from 0 to NSpecies
		A = ion_ref.A;
		Z = ion_ref.Z;
		T = ion_ref.T;    //ps
		T_err = ion_ref.T_err;
		C = ion_ref.C;    //m
		C_err2 = ion_ref.C_err2;    //m
		Bp= ion_ref.Bp;   //Tm
		Bp_true= ion_ref.Bp_true;   //Tm
		v = ion_ref.v;    //m/s
		v_err = ion_ref.v_err;    //m/s
		//v_ns = ion_ref.v_ns;      // m/ns
		//v_err_ns = ion_ref.v_err_ns;
		gammat = ion_ref.gammat;
	}
	//Print
	void PrintT(){cout<<"T = "<<T<<endl;}
	void PrintInfo()
	{
	    cout<<"|name= "<<name<<"| "<<fixed<<setprecision(10)<<" |Z = "<<Z<<"| |T = "<<T<<"| |C = "<<C<<"| |v = "<<v<<"| |gammat = "<<gammat<<"| |Bp = "<<Bp<<"|"<<endl;
	}
};

class ION_UNKNOWN:public ION
{
public:
	double M_cal;                 //final result
	double M_cal_square;                 //final result
	double M_cal_err;             //sigma m_calculated
	double Mvq;
	//double* M_cals;				  //result from each ref ion
	double M_AME;				  // AME INFO for comparison 
	double M_AME_err;
	double BpSum;				  // sum of Bp from each calculaiton using different ion_ref, to be averaged later
	double BpSum2;				  // sum of Bp*Bp from each calculaiton using different ion_ref, for sigma calculation
	//double Bp_err;				  // statistical result: std deviation sigma of Bp
	double Bp_err2;				  //error square
	bool use;					  //whether M_cal of this ions_unknown is effective or not
	int ref_n;                   //how many ref ions used
	void operator=(const ION& ion_ref)   //generate it from known ion
	{
		name = ion_ref.name;
		ion_number = ion_ref.ion_number;   //unknown ion does not have serial number in one injection
		inject_number = ion_ref.inject_number;
		Species = ion_ref.Species;  //from 0 to NSpecies
		A = ion_ref.A;
		Z = ion_ref.Z;
		T = ion_ref.T;    //ps
		T_err = ion_ref.T_err;
		C = ion_ref.C;    //m
		C_err2 = ion_ref.C_err2;    //m
		Bp= ion_ref.Bp;   //Tm
		Bp_true= ion_ref.Bp_true;   //Tm
		v = ion_ref.v;    //m/s
		v_err = ion_ref.v_err;    //m/s
		//v_ns = ion_ref.v_ns;      // m/ns
		//v_err_ns = ion_ref.v_err_ns;
		gammat = ion_ref.gammat;
	}
	double Calculate_M(double Bp_in)  // ,T*0.001 ->converted to :[ns]
	{
		double v_ns=v/1000000000;
		//double v_err_ns=v_err/1000000000;
		//cout<<fixed<<setprecision(10)<<"v="<<v<<endl<<"v_ns="<<v_ns<<endl;
		double ss=1.0/v_ns/v_ns - 1/Vc/Vc;
		
		return  Z*u*(Bp_in/10.3642686577806)*sqrt(ss);
		
	}
	void Calculate_M_err()    // mass and mass_err and mvq
	{	

		double v_ns=v/1000000000;
		double v_err_ns=v_err/1000000000;
		//cout<<fixed<<setprecision(10)<<"v="<<v<<endl<<"v_ns="<<v_ns<<endl;
		double ss=1.0/v_ns/v_ns - 1/Vc/Vc;
		double AA = Z*u/10.3642686577806;

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
	    <<" |A= "<<A<<" |Z= "<<Z<<"\033[0m |N = "<<N<<fixed<<setprecision(7)<<endl
	    	<<" |AME = "<<AME<<" +- "<<AME_err<<" MassUnknown: "<<MassUnknown<<endl
	    	<<"nucmass= "<<Mass<<"AtomicMass= "<<AtomicMass<<"BE= "<<BE<<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
	void PrintInfo(ofstream & outfile)
	{
		outfile<<"-------------------------------------------------------------------"<<endl;
	    outfile<<" |species No. "<<Species<<endl<<" |name= "<<name<<"| "<<" |Aname= "<<Aname<<"| "<<" state: "<<nuclide_info
	    <<" |A= "<<A<<" |Z= "<<Z<<fixed<<setprecision(1)<<" Tz= "<<Tz<<" T_isospin= "<<T_isospin<<endl
	    <<fixed<<setprecision(7)<<" |AME = "<<AME<<" +- "<<AME_err<<" MassUnknown: "<<MassUnknown<<endl
	    <<"nucmass= "<<Mass<<"AtomicMass= "<<AtomicMass<<"BE= "<<BE<<endl;
		//cout<<"-------------------------------------------------------------------"<<endl;
	}
	
};




void NUBASE2020_Build(NUBASE_IONSpecies* ISS,TString NUBASE_FILE)
{
cout<<" debug "<<endl;

for(int i=0;i<ZN_MAX_Z;i++)
	for(int j=0;j<ZN_MAX_N;j++)
		{ZN[i][j]=-1;}


bool AMEreadin_check_ON 		=1 ;
bool OUTPUT_check_ON 			=1 ;
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
	cout<<"----NUBASE2020_Build:---------- AMEreadin_check_ON ----------------"<<endl;
}
/////////////////////////////////////////
int n=0;
bool str_end_space_remove = 0;

infile.open(NUBASE_FILE);   //NUBASE V4数据
if(!infile){cerr<<"error !! !infile ";exit(1);}
if(OUTPUT_check_ON){ outfile.open("output//AME_data_check.txt");   cout<<"--------- output_check on ------------"<<endl;}
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

    string str_Aname_tmp = str_Aname_in + str_line.substr(16,1);
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

    if(OUTPUT_check_ON)ISS[n].PrintInfo(outfile);
    
    if(ISS[n].state_i==0) //ground state
    {
    	ZN[ISS[n].Z][ISS[n].A-ISS[n].Z] = n;
    }
    
    n++;
}
if(str_end_space_remove)cout<<"----NUBASE2020_Build: str end space removed "<<endl;
else cout<<"----NUBASE2020_Build: str end space not removed "<<endl;

cout<<" len max = "<<str_line_len_max<<" len min = "<<str_line_len_min<<endl;
infile.close();
if(AMEreadin_check_ON)outfile_str_check.close();
//_______________________ infile __________________________________

cout<<"----NUBASE2020_Build:   AME20 lines = n of ISS[i] = "<<n<<endl;


/*
//====================== ZN =========================
//方格占据--同位素
ZN_n=0;
bool FIND_L = false;
for(int i=0;i<120;i++)
{
	FIND_L = false;
	for(int j=0;j<200;j++)
	{
		if(ZN[i][j]>-1)
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
		if(ZN[i][j]>-1)
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


