//--------------------------------IMP CSR IMS DToF GtC ---------------------------------------------------------
//| Institue of Modern Physics, Cooler Storage Ring group, Isochronous Mass Spectrometry, Double Time of Flight|
//|                                            gammat-C method                                                 |
//|                                            Main Function                                                   |
//| author: Luo Yinfang, Xing Yuanming*, Lv Jiahao,                                                            |
//| run by : ROOT 6.26/06,6.24/06  WINDOWS/UBUNBU        https://root.cern                                     |
//|                                                                                          |since 2021 0916| |
//________________________________ version :   v0.97 20250228__________________________________________________//
#include <math.h>
//#include <unistd.h>   // do not include this on WINDOWS
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include"INPUT//GtC_DrawingFormat.h"     //头文件，画图格式 包含各种TGraph TGraphErrors Th1...
#include"INPUT//GtC_IonClass_NUBASE.h"          //头文件，设置了各种类，包括ION离子 IONSpecies 离子种类， NUBASE2020
#include"INPUT//GtC_Funcions.h"
//=====================因为经常要设置，所以尽量放在前面容易看到=========================================
//--------------------- LOOP PARAMETERS ----------------------------  
bool LOOP_ON     =         0   ;   // control all the canvas below  画布开关/循环开关   开启循环了就不画图
//----------------------------------------------------------------------------------------------------
bool L_ddt_correlation_ON =   1;  //最低卡方处 L和ddT近似一个线性关系，打开此开关可以只在此最佳线附近进行扫描，加快扫描速度
//----------------------- L --------------------------------------

//|| 2017_58Ni 18.046 0.1470
//|| 2021_36Ar_SET3 18.057 0.0957
//|| 2021_36Ar_SET2 18.052 0.0760
double L_down  = 18.046; //L循环起始位置 18.046
int L_n     =    1;      //L循环次数 不循环就是1
double dL=0.001; //  步长           
//----------------------- ddT --------------------------------------
double ddT_down= 0.1470; //0.1470 ddT循环起始位置  0.1470
int ddT_n   =    1;      //ddT循环次数 不循环就是1
double dddT=0.0002;//  步长（L_ddt_correlation_ON 开启后无效） //0.0761 //-0.058
bool scan_k_ddtC_ON=0; bool scan_C_inter_ON=0;
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define ions_unknown_n 11   //待计算的粒子数（ions_n读入的全部粒子）
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define MASS_VER       1        //v1 v2: unweighted  v3~ : error weighted 
//Mass Version质量计算方法 ：1.等权;2.等权 质量直方图高斯拟合;3.单次注入等权，结果加权（偏微分a系数底层误差传递;4.粗暴认为单次计算全独立;5.？;6.dF加权算法
#define USE_WEIGHTED   0        //历史遗留 v3
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//_____________________因为经常要设置，所以尽量放在前面容易看到________________________________________

//////////////////////////// for different experiment settings  //////////////////////
//--------- 2017_58Ni,2021_36Ar_SET2,2021_36Ar_SET3, 2021_36Ar_SET1,            ----------------------
    //人为指定本次实验数据。因为在一次程序运行中不会更改所以采用手动注释的方式来更改定义
    //TString THIS_EXP= "2017_58Ni";                              ////用于研发与验证 gtC方法的标准范本数据
    //TString THIS_EXP= "2021_36Ar_SET2";                         ////22Si,23Si
    TString THIS_EXP= "2021_36Ar_SET3";                       ////22Al,26P,27S,31Ar,23Si
    //TString THIS_EXP = "2021_36Ar_SET1";                      ////23Si,22Al,26P,27S,

//!!! NOTE:在程序中搜索关键词进行修改   @EXP@ ---- 实验依赖   ##ARTIFICIAL ----  人为设定  
////////////////////////////  switch and fundamental settings  ///////////////////////////////////
#define ON_UBUNTU  444
#define ON_WINDOWS 555

#define MAX_IONSPECIES 120      // 用于开数组  核素种类数上限
#define MAX_IONS       120000   // 用于开数组  所有离子数目上限 (上限来自：2021 36Ar SET2 有11万多个)
#define MAX_INJECTIONS 35000      // 用于开数组  注入数上限
#define MAX_SCAN_n 5000         // limit of scan loop times          扫描数上限
#define subregion_n 50         // divide C to subregions  γt(C)曲线，将C分成100份
#define dsubregion_n 10         // divide subregion to interval steps   将上面的小区间再分10份
#define ref_n_MAX 30            //  参考核上限
#define IAS_N_MAX 60            //  IAS态上限
#define Isomer_N_MAX 30            //  Isomer上限
#define TzN 6            //指定出现了多少Tz 值 应当替换为 vector All_Tz 
//---------------------- gtC ------------------------------------

#define h1_STEP1_EL_PERCENT 99.5           //轨道C 有效区间，去掉两侧过大过小C散点， 中间留下的有效部分数量占比
#define squeeze_k 5            // 去散点迭代次数 5
#define n_sigma 3               // filtering boundary how many times of stderror Z-score方法去散点： 偏离平均值多少倍标准差的，去掉
#define set_smooth_k 2          // 平滑取左右各 set_smooth_k 个点
#define choose_largeZ_min -114.514           //（已无效）
#define choose_largeZ_max -22.333
#define C_DIVISION_CHOSEN_MIN 10        //构建gtC曲线时，C小区间内gt点的数目最小值，必须超过这个数，这个区间才有平均值
#define Z_division 7
#define use_gtC_type 1      //1: gtC,readin or not || 2:gt0, 3: inj 1-4 4:fit poln 
                            //对应不同的方法得到γt(C)曲线，1.γt(C)曲线是 直接读取的(readin)或者 计算出的(正常运行时使用)!!!!!!  2.γt为固定值，3.γt(C)随注入数变化(现已验证不变)，4.多项式拟合
#define gtC_chosen_gtMIN 1.3   // 要用于 each_gtC 的Fill_ion_gt
#define gtC_chosen_gtMAX 1.4

int CONDITION_gt_Z = 6;    //@EXP 14      //绘制γt  C 二维分布的选项，对 A Z 进行筛选
int CONDITION_gt_A =10;   //注意！后续根据不同实验会修改！！//@EXP 24 18  


double scan_k_ddtC_down=0, scan_k_ddtC_up=0.020, dscan_k_ddtC=0.01; //扫描开关(已无效)
double scan_C_inter_down=128.75, dscan_C_inter=0.005;               //扫描开关(已无效)
double gt0_down=1.384, dgt0=0.001;
int scan_k_ddtC_n    =   1;
int scan_C_inter_n   =   1;
int gt0_n   =   1;

int IsRef_MIN_N =100;          // @EXP 参考核最小数目，计数过小不用这个核做参考核
double IsRef_AME_ERR_MIN= 5;      //@EXP 参考核的AME 误差上限

//------------------------------------ 以下是变量设置 ----------------------------------------------------------------

// choose your running environment which is decisive when making directories
//int THIS_ENVIRONMENT = ON_UBUNTU;            //选择运行环境
int THIS_ENVIRONMENT = ON_WINDOWS;
//------------------------------------- INPUT DATA -------------------------------------
TString INFILENAME_1;                                        //输入文件   已经经过核素识别的数据文件
TString INFILENAME_2; //为本设置准备好的AME20核数据文件
TString INFILENAME_3; //同质异位素相似态（isobaric analog state简称IAS）
TString INFILENAME_4; //同核异能态(Isomer)
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ options @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//==================================  各种 bool 开关  ========================================
//------------  OUTPUT FOLDER  -------------------------
bool MAKE_DIR = 1;                              //创建文件夹
bool Create_ionspecies_folder_ON   = 1;           //创建每个核素的文件夹，分核素储存图片
//------------ INPUT CHECK ----------------------------
bool READ_prepared_EXP_AME         =        1   ; // =1 读取每次实验专门准备的部分AME文件 速度更快   =0 直接从nubase20读取所有的AME 信息
    bool Isomer_IAS_check_ON       =        1   ;//READ_prepared_EXP_AME=0 直接从nubase20读取所有的AME 信息  
bool ReadInIons_check_on           =        0   ; //检查读入： 每一个 ion 的所有信息 输出到文件 主动中断 return
bool SetIonSpecies_check_on        =        0   ; //检查读入： 每种 ionspecies
bool ReadInAME_check_on            =        0   ; //输出读入的AME质量，检查一下

//----------- MASS calculation ----------------------
bool Show_MASS_VER1_ON              =    0;      //在MASS_VER不是1时 仍然显示等权的质量结果 v1
bool Do_mass_err_scatter_ON         =    1;      // scattering error requires MASS_VER>=3, error bar is needed   计算散布误差（如果比加权的平均值误差大，则采用散布误差）
bool Only_draw_mass_VE              =    0  ;    //1: only draw error-weighted mass 


//------ merr 
bool READIN_MERR_ON       =      0;              //读入下面的输出文件
bool Do_iont_error_out_ON =      0;              //输出每个粒子的相关信息以及误差（等权情况下无效）
bool Do_iont_error_in_ON  =      0;              //外部导入算好的误差

//------------- C filtering -------------------------------
bool C_filter_ON    =       1;                    //限制轨道 
double C_filter_min = 128.600, C_filter_max=129.000;//@EXP 在后面还要赋值 原128.58  128.61   129.05

//------------- INJ SKIP ------------------------------
bool INJ_SKIP_ON    = 0;              // use bool INJ_SKIP[] to skip certain injections   跳过某些注入
    TString INJ_SKIP_INFILENAME;
    bool INJ_SKIP_THIS  = 0;              // tmp in loop                                      配合上面使用
    bool INJ_sample_ON  = 0;              // sample certain injection for test                取样（某些注入）测试
bool Do_PD_convergence_test =0;       //20240703 PD_convergence_test || in mass_v2 

//------------- gtC construction --------------------
bool Show_c_squeeze_tmp_ON      =  0   ;
bool Do_gtC_smooth_ON           =  1  ;             //平滑γt(C)曲线
int  smooth_opt                  =  1 ;             //平滑γt(C)曲线的方法： 0 不平滑，1 附近点取加权平均，2 附近点的(1+γt^2*ΔC/C)连乘取开方
bool Show_c_smooth_ON           =  0  ;
double  gtC_ERR_upper_bound  = 0.02     ;   //step1 筛选γt散点， γt误差上限 @EXP 0.02
//+++++++++++++++++++++++++++++++++++++++++
bool Show_C_Division_n_chosen_ON   =        0   ;  // cout ions number in every subregion of C  显示每个C区间的粒子计数    

bool DO_gtC_v2_ON               =0;   //20230420
bool gtC_fit_ON                 =0;
int  gtC_fit_pol                =0;  // poln =0 , poln=1, 2
bool READIN_gtC_chosen_ON       =0;  // in use_gtC_type 1 
// 各个step 的一系列统计直方图 h1 h2 散点图 gr,grr
bool SHOW_gtC_step0             =1;
bool SHOW_gtC_step1             =0;
bool SHOW_gtC_step2             =0;
bool SHOW_gtC_step3             =0;
bool SHOW_gtC_step4             =0;
// 各个step轨道分布
bool SHOW_h1_STEP0              =1;
bool SHOW_h1_STEP1              =1;   // this h1 provide distribution for determining effective C region
bool SHOW_h1_E_step2            =0;
//显示各个阶段处理得到的gtC曲线
bool SHOW_grerr_gtC_chosen_steps_ON =1;   // usually on 应当常驻显示

bool GT_SHIFT_Control_ON        =0;
    double GT_SHIFT             =-0.001;    //0.004
    double GT_SHIFT_Left        = 128.1;         // GT_SHIFT REGION 精确调控适用区间
    double GT_SHIFT_Right       = 129.99;
bool DM_FILTER      = 0;// 
int DM_FILTER_n = 0;

bool Do_gtC_with_time_ON =0;  //按注入 分别处理gtC
int gtC_inj_min_1;int gtC_inj_max_1;
int gtC_inj_min_2;int gtC_inj_max_2;
int gtC_inj_min_3;int gtC_inj_max_3;
int gtC_inj_min_4;int gtC_inj_max_4;

//============== test, trial , verification, other process =============
TH1D* h_test = new TH1D("h_test", "h_test", 900, 0.1, 1.9);       //测试用直方图，各项按需改动 
//----------- time injection analysis -------------------
bool time_divide_filter_ON  =    0;
int  choose_inj_option     =     0;  //injection analysis 4 parts
//------------20230614 Do_dA0_T ---------------------------- 
bool Do_dA0_T_ON             =   0;
    bool Recalculate_after_dA0T_ON = 0;
    bool dA0T_ions_filter_ON       = 0;   //只选取TST窗口内离子
    bool Draw_each_T_C_ON          = 0;      //20230614 各种离子C - T
//------------------------------------
bool Do_gtC_largeZline_ON   =    0;
bool Do_Fit_BpC_ON          =    0;       //Fit BpC for all ionspecies// 如果不做scan C0 会自动开启=true
//----------------------------------
bool TEST_USE_OWN_gtC_ON = 0; //2024 1113 use own gtC-ISF




//==================== SHOW  Draw_ON=============================

//--------------- show gtC
bool c_gtC_chosen_ON =           0;//绘制γt(C)曲线并保存
bool c_h2_gtC_all_ON =           0;
bool c_gtC_Tz_ON            =    0;
bool c_gtC_time_all_ON      =    0;
bool c_avegt_ionspecies_ON  =    0;   // 历史遗留
bool c_avegtC_ionspecies_ON =    0;     // BuildAvegtCCurve_each() 查看每种核各自的gtC
    bool Do_gtC_divide_Z    = 0;   //两条gtC 一个是小于Z_divide 的核 另一个是大于Z_divide 的核
    bool Do_gtC_divide_A    = 1;   //两条gtC 一个是小于A_divide 的核 另一个是大于A_divide 的核
bool c_gtC_v2_ON            =    0;

//---- each -----
bool c_each_ionspecies_ON   =    0;            //绘制每种核素相关信息
    bool show_each_h_C_gtC_ON   =1; bool savefile_each_h_C_ON =1;
    bool show_each_h_T_ON   =0;     bool savefile_each_h_T_ON =1;
    bool show_each_h_mvq_ON =0;     bool outfile_each_h_mvq_ON=1;
//_____ each ______


//------------------ MASS RESULT --------------------------------
bool Show_Mass_Result_ON   =     0;      //显示最终的质量结果
int Show_Mass_Result_opt_n =     3;          //显示设置 v1 v2 VE
int Show_Mass_Result_opt [3];            // control the presentation of each version of mass results  //控制各算法的结果输出
//opt[0]:v1, opt[1] :v2, opt[2]: v3, opt[3]:v3_SCA, opt[4]:v4

//---- mass ver 2 ---------
bool mass_v2_show_each_fit_ON = 0; // show each h_dm fit in MASSVER==2
bool Re_IonIdentify_on = 1;   // EXP 2021_36Ar_SET3                                                     //利用质量结果进行重新识别

//------------ MASS result check ------------------
bool Show_mvqC_each_ON      =    0;   // 显示每一种类离子的 mvq-C v1 散点图 ,各自一个canvas ,存储图片
int  Show_mvqC_each_MIN     =    100;  // 一种核的计数要大于等于这个值，才会有窗口显示，对mvq-C 散点图进行线性拟合
    bool Do_Generate_k_mvqC_ON   =1;  // 对于每种核的mvqC散点图做直线拟合 得到各种核的斜率
bool Show_dmC_each_ON       =    0;
bool DoShow_mvqC_each_h2_ON   =    0;  
bool Show_h_iont_merr_each_ON =  0;
bool Show_h_refcal_merr_each_ON= 0;
bool Show_h_iont_refs_chi_each_ON= 0;
bool Draw_gtC_all_ON       =     0;   // use function Draw_gtC_all
//--------------- with injection time ---------------



//~~~~~~~~~~~~~~~~~~~~~~~~~~ canvases ~~~~~~~~~~~~~~~~~~~~~~~~~~
bool c0_ON =                     0;
bool c_h1_CDistribution_ON =     0;
bool c1_ON =                     0;  // Bp-C
bool c1_Bp_C_inj_ON =            0;  // Bp-C-injection TGraph 2D  Bp-C-注入 三维空间散点图
bool c1_Bp_C_ionspecies_ON =     0;
bool c2_ON =                     1;   // gtC all and chosen 展示初始gtC 散点和最终的gtC 曲线！！
bool c3_h_INJ_ON =               0;  // injection-ions 看每次注入的离子数分布情况
    bool Do_outfile_inj_ions   = 1;  // 将每次注入有多少离子输出到文件
bool c4_Do_analysis_with_time_ON=0;  // vary with time 磁场 周期 等等 随时间 / 注入数 变化 横轴可以用特定格式显示 时间
    bool Generate_each_T_time     =   0;bool Generate_each_mvq_time   =   0;bool Generate_each_Bp_time    =   1;
    bool Generate_each_Bp_inj     =   0;
bool c5_Do_each_T_fix_ON       = 0;  // Isochronous curve
    bool Do_Isochronous_mass_ON =1; // Isochronous curve in mass stddev (keV)
bool c6_CdC_ON                 = 0;
bool c7_ON                     = 0;  // h_T_all
bool c7_b_ON                   = 0;  // h_mvq_all
bool c8_ON                     = 0;  // gt- mass
bool c9_h_Mvq_ON               = 0;  // h_Mvq 历史遗留 在计算过程中填入h_Mvq 
bool c10_ON                    = 0;  // fit BpC
bool c11_ON                    = 0;  // BpC lnBpC
bool c_dmdC_ON                 = 0;  // do gr_dmdc->set point
    bool c_gr_dmdC_ON = 0;bool c_h2_dmdC_ON=0;
bool c_dmISO_ON                 = 0;
bool c_ERRANA_on               = 0;   // 20230419 err ana
bool Do_Draw_2024ERRANA_ON     = 0;

bool c_h_Mvq_2gaus_on          = 0;  //绘制几种核的m/q-C的直方图（测试用）

//________________________________ SHOW  Draw_ON ______________________________
// out of loop L ddT 扫描循环外
bool cc1_ON                 =    0;
bool cc_xn_v1_ON            =    1;
//#######################################################
bool OUTPUT_MASSRESULT            =      1; // always ON  输出最终质量结果文件
//#######################################################

bool outputswitch1                =       false;
bool outputswitch2                =       0   ;
bool OUTPUT_allions               =      false;
bool OUTPUT_singleion             =      0    ;  //outfile_single_ion 
bool OUTPUT_check_abnormal_gammat =      0;      //输出异常γt的粒子进行检查
 
bool OUTFILE_each_ion_complete_ON =      0;
bool OUTFILE_exc_state_ON         =      0;
bool OUTFILE_v2err_gtC_ON         =      0;
bool outfile_each_mass_data_ON    =      0;    // 每种离子的质量信息单独保存txt  LOOPON 时应该关掉！ 以免存储大量txt
bool outfile_ERR_ANA_ON           =      0;    //20231018 errors output file 

//_________________________________  bool 开关 结束 ______________________________________________



//------------------------------------以下是全局变量----------------------------------------------------------------

double L=   0.0,ddT= 0.0 ,scan_k_ddtC=0.0  ,scan_C_inter=0.0 ,  gt0 =0.0 ;
double L0=18.053,ddT0=0.099;

int scan_loop_i=0;                        //这是扫描的第几次循环
int scan_loop_total=0;                    //预计的总循环次数
///////////////////////////////////////////////////////////////==================///////////////////////////////////////////////

//const Double_t V_c=0.2997924580;      //m/ns 定义在 IonClass.h
const Double_t k0  = 10.3642686577806;  //coefficient of MvQ / Bp
//const Double_t u  = 931494.102417;   //keV

TDatime TIMENOW;
///////////////////////=========  fstream  ==========//////////////////////////////// 
ifstream infile;                     //打开文件
ofstream outfile_logs;               //记录文件 日志logs
ofstream outfile;
ofstream outfile1;
ofstream outfileBp;
ofstream outfile_single_ion;
ofstream outfile_Lddt_Xn2;
ofstream outfile_errors;
ofstream outfile_errors_check;       // v2 mass errors check 
ofstream outfile_v2error_gtC;
ofstream outfile_PD;                         //20240703 PD_convergence_test
ofstream outfile_ERR_ANA;
ofstream outfile_C_Division_n;
ofstream outfile_sca_ave;            // 用于输出每种离子 ave/sca error 比较结果
ofstream outfile_test;               //测试用，按需调用
ofstream outfile_VER2_FIT_info;      // 在 mass v2 输出直方图拟合相关信息
ofstream outfile_each_ion_complete; //
ofstream outfile_grerr_txt;         //  保存 TGraphErrors 的信息

///////////////////=========== directory ====================///////////////
TString makedir;                    // 指令字符串，根据时间创建文件夹
TString FILEPATH;                   // 文件路径
TString makedir_m_s; TString FILEPATH_m_s;
TString makedir_c1_saved; TString FILEPATH_c1_saved;
TString makedir_chosen; TString FILEPATH_chosen;
TString makedir_ionspecies;TString FILEPATH_ionspecies;     //每种核素对应的文件夹
TString makedir_k_T; TString FILEPATH_k_T;
TString makedir_DodA0T; TString FILEPATH_DodA0T;
TString filenameadd="_z";

///////////////////=========== filename ====================///////////////
string filename0;// logs 日志文件
string filename1;// output single ion
string filename2;// output species final mass results 最终质量结果

///////////////////=========== root file ====================///////////////
TFile* Tfile_save_hist; // 用于存储直方图到root文件


///////////////////=========== latex strtmp ====================///////////////
TLatex* lat_n = new TLatex();   
TString lat_text;
TString strtmp;                 //read in useless info

TString ThisParaInfo;//   L=... ddT=...
TString gtC_fit_info;          //历史遗留

Int_t NSpecies;                 //species of ions         离子种类数
Int_t NSpecies_IAS;             //species of IAS     IAS态数
Int_t NSpecies_Isomer;             //species of isomers   
Int_t NSpecies_other_exc;             //species of isomers      
int Z_readin[MAX_IONSPECIES];
int A_readin[MAX_IONSPECIES];
int N_readin[MAX_IONSPECIES];
TString name[MAX_IONSPECIES];               //ion names
TString namestr,ResultFileName; //namestr 读入的核素名

///////////////////=========== 质量刻度卡方相关 因为每次循环都要记录 所以是全局的====================///////////////
double L_ddT_Xn2[MAX_SCAN_n]        = { 0 };             //每一次循环都有一个卡方值， [scan_loop_i]
double L_ddT_Xn2_v2[MAX_SCAN_n]     = { 0 };
double L_ddT_Xn2_VE[MAX_SCAN_n]     = { 0 };
double L_ddT_Xn2_sys[MAX_SCAN_n]    = { 0 };             //加上系统误差 add err_sys
double L_ddT_Xn2_v2_sys[MAX_SCAN_n] = { 0 };
double L_ddT_Xn2_VE_sys[MAX_SCAN_n] = { 0 };

TGraph* gr_errsys_chin;             //系统误差对卡方的影响 err_sys ~ chi_n 全局 在Show_Mass_Result 里面使用 global
int gr_errsys_chin_n;

int IonNumber ;                 //IonNumber:occurrence of  this ion in one injection, not continuous!!  
int IonNumber_new;              //set new number to ions in one injection to make sure the number is continuous
int InjectNumber;               //InjectNumber from original data, not continuous!!
int InjectNumber_pre=0,InjectNumber_new=0; //previous input data , new number to make it continuous
int Injection[MAX_INJECTIONS];           //record how many ions in all the injections.                         这次注入有几个核
int Injection_m[MAX_INJECTIONS]={0};     //record how many ions that can have mass results in this injections. 一次注入最后有多少个离子有质量结果
int Injection_subfix[MAX_INJECTIONS];    //injection to subfix , the starting subfix for each injection        每次注入的第一个粒子的序号
int total_injection=0;          //record how many times we had effective injections , starts from 1   总注入数
int ion_subfix=0;               //the starting subfix at certain Injection
int ZN_ID[120][220];            //(p,n)--- ionspecies ID                                              核子的质子数，中子数
int ions_n=0;                   // ions_n  =  -----all ions 总粒子数
int A=0;int Z=0;                // A核子数 Z电荷数

double m_tmp;                   // to save the mass calculated in the loop of varying gammat
double Ave_deltaM_all=0;               //所有核 计算值与AME的差值 的平均值
double Ave_deltaM_all_v_others=0;
double Cmax=-1,Cmin=999,C_ave_all=0.0; //Cmax Cmin 计算得到的C的上下限    C_ave_all C的平均值
double C_Shift_ratio=0;
double gtC_chosen_CMAX=0.0,gtC_chosen_CMIN=99999.0;
double gt_max=-1;                      //γt最大值
double gt_min=999;                     //γt最小值
double T_RANGE_MIN=580;                // T 周期范围 人为选定要展示的范围 
double T_RANGE_MAX=680;
double mvq_RANGE_MIN=1.8;                // mvq 范围 人为选定要展示的范围 
double mvq_RANGE_MAX=2.0;

double lr = 0.01;      //limit_ratio: 求导小区间，正比于参数误差 f'=dy/dx, dx = limit_ratio*para_err

int count_gtC_masserr_condition_1 =0;      //20230707 add gt err->mass v2 err
int count_gtC_masserr_condition_2 =0;
int count_gtC_masserr_condition_1_VE =0;   //20230707 add gt err->mass v2 err
int count_gtC_masserr_condition_2_VE =0;
int count_condition1_onemasserr=0;          int count_condition2_onemasserr=0;
int count_condition1_onemasserr_gtC=0;      int count_condition2_onemasserr_gtC=0;

Double_t T,T_err,v_err,v,C,gammat;// T周期，T_err周期误差，v_err速度误差，v速度，C周长，gammat γt

Double_t gtC[subregion_n];
Double_t gtC2[subregion_n];            //C square divided into sub regions
Double_t Sigma_gtC[subregion_n];
Int_t C_Division_n[subregion_n];       // count for numbers in each subregion
int n_gtC_chosen_C = 0;                //计算用的γt(C)区间个数

int Z_LowEdge=8;

//------------------- 0929 gtC_Tz 后续应当替换为 vector All_Tz
double Tz_all_ions[TzN] = {-0.5,-1.0,-1.5,-2.0,-2.5,-3.0};//此次实验中各种核可能的Tz
double gtC_Tz[TzN][subregion_n];
double gtC2_Tz[TzN][subregion_n];
double sigma_gtC_Tz[TzN][subregion_n];
int C_Division_n_Tz[TzN][subregion_n];
int count_Tz_n[TzN]={0};
//------------------- 0929 gtC_Z
//int A_color[25]={0,0,0,0,0, 0,0,1,2,3, 4,5,6,7,8,  9,kPink+1,kSpring-6,kAzure+7,kOrange-3, kCyan+1,kViolet+6,38,41,46};




int gr_n=0;                      // common count of number of points in a TGraph

bool IsNew=true;                 //flag of new species in finding new species at first      //是否为 新核素
bool ERROR_FLAG[5]={0};

int ii=0 ;  // tmp index
bool INJ_SKIP[MAX_IONS]={0};            //20240703
TString SKIP_FILENAME[5000];           //20240705 
int SKIP_n=0;

TF1 * f_zero;            //水平虚线 y=0
double xtmp=0;double ytmp=0;double yerrtmp=0;

//------------------------- ISF --------------------
long double ISF[subregion_n * dsubregion_n] = { 0 };//Integral Step Factor 步长积分因子（加快运算）
long double ISF_u[subregion_n * dsubregion_n] = { 0 };//γt(C)曲线误差上限
long double ISF_d[subregion_n * dsubregion_n] = { 0 };//下限
/*
long double ISF_1  [subregion_n * dsubregion_n] = { 0 };
long double ISF_u_1[subregion_n * dsubregion_n] = { 0 };
long double ISF_d_1[subregion_n * dsubregion_n] = { 0 };
long double ISF_2  [subregion_n * dsubregion_n] = { 0 };
long double ISF_u_2[subregion_n * dsubregion_n] = { 0 };
long double ISF_d_2[subregion_n * dsubregion_n] = { 0 };
*/
long double ISF_inj1[subregion_n * dsubregion_n] = { 0 }; // 暂时先只有center 不考虑 上下界
long double ISF_inj2[subregion_n * dsubregion_n] = { 0 };
long double ISF_inj3[subregion_n * dsubregion_n] = { 0 };
long double ISF_inj4[subregion_n * dsubregion_n] = { 0 };
//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////         Object      /////////////////////////////////////////////
ION ions[MAX_IONS];  //class ION , all ions in this experiment(result of FastFortran)  定义了一个类class ION（头文件 IONCLASS_TVMASSAr36.h）
double ReadIn_merr[MAX_IONS]={0};
double ReadIn_m_v2[MAX_IONS]={0};

//ION ions_copy[57000];

IONSpecies ionspecies[MAX_IONSPECIES];                //  all ion species in this experiment 本次实验所有的核素种类
IONSpecies IAS_ionspecis[IAS_N_MAX];
IONSpecies Isomer_ionspecis[Isomer_N_MAX];

NUBASE_IONSpecies ISS_AME[5888];   // nubase 2020 : 5843 nuclides
vector<double> All_Tz;
vector<TString>All_Tz_str;
vector<int> All_Tz_COLOR;
int All_Tz_n;



//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////         FUNCTION      /////////////////////////////////////////////

//利用γt(C)曲线，计算未知核的Bρ
double Calculate_Unknown_Bp(ION& ion_ref , ION_UNKNOWN& ion_unknown, TGraph* gr_gammat_C)
{
    
    double nsteps = 100.0;
    double deltaC = (ion_unknown.C - ion_ref.C) / nsteps;
    double Bp=0.0,C=0.0;
    double gt=0.0;
    C = ion_ref.C;
    Bp= ion_ref.Bp;
    for(int i=0;i<nsteps;i++)
    {
        gt=gr_gammat_C->Eval(C);
        Bp = Bp*(1+gt*gt*deltaC/C);
        C = C+deltaC;
    }
    return Bp;
}
////////////////////=============================
////================================================  bp ==============================================///////////////////////

//TGraph 计算目标核的Bρ  
double Calculate_Unknown_Bp_1(double Bpi,double Ci,double Ct, TGraph* gr_gammat_C,double& gammat_ave,int nsteps,int inject_number )
{//MODE 1
    //outfile.open("abnormal_gt.txt",ios::out|ios::app);
    // calculate Bp by changing C step by step, fixed step length,
    //int nsteps = 1000;
    //if(Ci>gtC_chosen_CMAX)Ci = gtC_chosen_CMAX;
    //if(Ci<gtC_chosen_CMIN)Ci = gtC_chosen_CMIN;
    //if(Ci>gtC_chosen_CMAX||Ci<gtC_chosen_CMIN||Ct>gtC_chosen_CMAX||Ct<gtC_chosen_CMIN)outfile<<"Eval("<<Ci<<") = "<<gr_gammat_C->Eval(Ci)<<"  Eval("<<Ct<<") = "<<gr_gammat_C->Eval(Ct)<<endl;

    double deltaC = (Ct - Ci) / nsteps;
    double Bp=0.0,C=0.0;
    double gt=0.0;
    gammat_ave=0.0;
    C = Ci;
    Bp= Bpi;
    
    for(int i=0;i<nsteps;i++)
    {
        
        gt=gr_gammat_C->Eval(C);
        //gt=gr_gammat_C->Eval(C); 
        //if(gt<1.1||gt>1.4) outfile<<"@@@@@@@@@@@   "<<inject_number<<" "<<Ci<<" "<<Ct<<endl;
        Bp = Bp*(1+gt*gt*deltaC/C);
        //Bp = Bp*(1+gr_gammat_C->Eval(C)*deltaC/C);
        C = C+deltaC;
        gammat_ave+=gt; 
    }
    //outfile.close();
    //cout<<"final C check ="<<C- ion_unknown.C<<endl;
    gammat_ave = gammat_ave/nsteps;
    return Bp;

}
//TSpline 计算目标核的Bρ  
double Calculate_Unknown_Bp_1(double Bpi,double Ci,double Ct, TSpline3* sp3_gammat_C,double& gammat_ave,int nsteps,int inject_number )
{// MODE 1
    
    double deltaC = (Ct - Ci) / nsteps;
    double Bp=0.0,C=0.0;
    double gt=0.0;
    gammat_ave=0.0;
    C = Ci;
    Bp= Bpi;
    
    for(int i=0;i<nsteps;i++)
    {
        gt=sp3_gammat_C->Eval(C);
        //gt=gr_gammat_C->Eval(C); 
        //if(gt<1.1||gt>1.4) outfile<<"@@@@@@@@@@@   "<<inject_number<<" "<<Ci<<" "<<Ct<<endl;
        Bp = Bp*(1+gt*gt*deltaC/C);
        //Bp = Bp*(1+gr_gammat_C->Eval(C)*deltaC/C);
        C = C+deltaC;
        gammat_ave+=gt; 
    }
    //outfile.close();
    //cout<<"final C check ="<<C- ion_unknown.C<<endl;
    gammat_ave = gammat_ave/nsteps;
    return Bp;

}
//////////////////

//得到(1 + γt * γt * ΔC/ C)数组，即 ISF 步长积分因子
void Calculate_GtC_ISF(TGraph* gr_gammat_C, long double* cal_ISF)
{    //超出C_MIN C_MAX，直接变为 0
    long double gt = 0.0;
    long double C = 0.0, C_i = 0.0;
    long double C_MIN = gr_gammat_C->GetX()[0], C_MAX = gr_gammat_C->GetX()[gr_gammat_C->GetN() - 1];//已经得到的γt(C)曲线的范围
    long double cal_Cstep = (C_MAX - C_MIN) / (n_gtC_chosen_C - 1.0) / dsubregion_n;//细分后的小区间
    int n = 1;//从1开始，对应 h_locate，超出范围的ISF[]为0
    for (int i = 0; i < n_gtC_chosen_C; i++)
    {
        C_i = gr_gammat_C->GetX()[i];
        for (int j = 0; j < n_gtC_chosen_C * dsubregion_n; j++)
        {
            C = C_i + cal_Cstep * (j - dsubregion_n / 2.0);//很奇怪，如果每次用 上次计算周长加一步 的方法，会有明显的误差(很小)（改成这样后大部分情况下能完全对上，但是有时还是有误差，不可控，不清楚原因）
            gt = gr_gammat_C->Eval(C); //插值来计算gt值
            cal_ISF[n] = 1 + gt * gt * cal_Cstep / C; n++;
            //cal_ISF[n] = C; n++;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  测试
            if ((j + 1) % dsubregion_n == 0)
            {
                if (i < gr_gammat_C->GetN() - 1)
                {
                    if (C + cal_Cstep * (j + 1) > gr_gammat_C->GetX()[i + 1] - dsubregion_n / 2.0 * cal_Cstep - cal_Cstep / 100.0)//采用这种算法，避免γt(C)曲线中间有空点
                    {
                        break;
                    }
                }
                else//最后一个点，gr_gammat_C->GetX()[i + 1]不存在，也算十个
                {
                    cal_ISF[n - 1] = 0;//超范围了，这个点变为0
                    break;
                }
            }
        }
    }
}

//采用等步长的算法，计算目标核Bρ
long double Calculate_Unknown_Bp_dF(double Bpi, double Cr, double Ct, TGraph* gr_gammat_C, long double* ISF, TH1D* h_locate)
{
    //gr_gammat_C γt(C)曲线    long double* cal_ISF 
    //h_locate 寻找每个粒子所处的周长区间用 
    //这个直方图奇怪的范围是为了和ISF区间对齐,使得
    //h_locate->GetBinLowEdge(h_locate->FindBin(C))  //此周长的左边界       =        ISF[h_locate->FindBin(C)]对应的C值

    if (Cr > Ct)//为了正向反向计算得到的值统一
    {
        long double f = Calculate_Unknown_Bp_dF(1.0, Ct, Cr, gr_gammat_C, ISF, h_locate);
        if (f > 0.001)//不为零
        {
            return Bpi / f;
        }
        else
        {
            return 0;
        }
    }
    long double C_MIN = gr_gammat_C->GetX()[0], C_MAX = gr_gammat_C->GetX()[gr_gammat_C->GetN() - 1];
    long double Cr_start = h_locate->GetBinLowEdge(h_locate->FindBin(Cr)), Ct_stop = h_locate->GetBinLowEdge(h_locate->FindBin(Ct));//开始结束所处的 bin 左边界
    long double Bp = 0.0, C = 0.0;//C:从 Cr 到 Ct

    Bp = Bpi;//初始的Bρ（参考核）

    for (int i = h_locate->FindBin(Cr); i < h_locate->FindBin(Ct); i++)
    {
        Bp = Bp * ISF[i];                                                                // 由 γt^2 = (ΔBρ/Bρ) / (ΔC/C) 得
    }
    Bp = Bp / (1 + pow(gr_gammat_C->Eval(Cr), 2) * (Cr - Cr_start) / Cr) * (1 + pow(gr_gammat_C->Eval(Ct), 2) * (Ct - Ct_stop) / Ct);
    //(1 + pow(gr_gammat_C->Eval(Cr), 2) * (Cr - Cr_start) / Cr) 开始与中心值相差的一小段，多走了，去掉
    //(1 + pow(gr_gammat_C->Eval(Ct), 2) * (Ct - Ct_stop ) / Ct) 结束与中心值相差的一小段，少走了，加上
    return Bp;

}

//利用误差因子dF的方法，给出一组参考核的下，目标核的质量结果，和参考核、目标核误差导致的误差 （v6）
void Calculate_Mass_with_err_dF(ION_UNKNOWN& ion_t, ION* ions_ref, int ref_n, TGraph* gr_gammat_C, long double* ISF, TH1D* h_locate, double& mass_result, double& mass_result_err2)
{   // 进来一个目标核 iont 多个参考核的数组 ions_ref, 参考核数量 ref_n, 使用的gr_gammat_C , 记录最终结果的 mass, mass_err
    //M: keV
    double mass_result_each[ref_n_MAX];          //每个核给出的质量结果
    double mass_err2_result_each[ref_n_MAX];      //每个参考核误差导致的误差

    //认为以下误差直接相互独立
    double mass_result_R = 0;                         //参考核数组给出的质量中心值                  (最终结果)

    double mass_err2_R = 0;                           //由于参考核误差导致的误差(方差)
    //double mass_err2_R = 0;                           //由于参考核误差导致的误差(方差)，取以下两者最大的那个        (不对，已无效)
    //double mass_err2_R_av = 0;                        //参考核给出的平均值误差 (前期用作权重因子和)
    //double mass_err2_R_scat = 0;                      //参考核给出的分布误差

    double mass_err2_T = 0;                           //由于目标核误差导致的误差(方差)

    double ref_eff_n = 0;                             //能给出正常质量的参考核数目
    //double test[10];//测试用数组


    for (int i = 0; i < ref_n; i++)
    {
        //  vt = L/(a5t*0.001+ddT);vr = L/(a5r*0.001+ddT);          //两个核的速度  V_c:  m/ns
        mass_result_each[i] = Calculate_Unknown_Bp_dF(ions_ref[i].Bp, ions_ref[i].C, ion_t.C, gr_gammat_C, ISF, h_locate);                               //这一步给出的是Bρ  
        mass_result_each[i] = 1000000.0 * ion_t.Z * V_c * V_c * sqrt((1.0 / (ion_t.v * ion_t.v)) - (1.0 / (V_c * V_c))) * mass_result_each[i];     //每个核给出的中心值
        mass_err2_result_each[i] = ions_ref[i].dFr * mass_result_each[i] * ions_ref[i].dFr * mass_result_each[i];                                                                               //由于参考核误差导致的误差
        if (mass_result_each[i] > 0.1)//超出范围的粒子作为参考核结果为0，只记录正常结果
        {
            mass_err2_R += 1.0 / mass_err2_result_each[i];                    //加权因子
            mass_result_R += mass_result_each[i] / mass_err2_result_each[i];  //加权质量结果
            ref_eff_n++;
        }
    }
    if (ref_eff_n > 0)//有可用的参考核，只有个别参考核超出范围
    {
        mass_result_R = mass_result_R / mass_err2_R;                                   //参考核数组给出的质量中心值                  (最终结果)
        mass_result = mass_result_R;                                                      //传回质量结果
        mass_err2_R = 1 / mass_err2_R;                                              //参考核给出的平均值误差 (方差)
        mass_err2_T = mass_result_R * mass_result_R * ion_t.dFt * ion_t.dFt;              //目标核误差导致的误差(方差)
        mass_result_err2 = mass_err2_R + mass_err2_T;                                     //参考核、目标核误差导致的误差(方差)
        return;
    }
    else
    {
        mass_result = 0;//目标核超出范围，或者全部参考核超出范围，那就没必要继续算了
        mass_result_err2 = 0;
        return;
    }
    // mass = 1000000.0 * ion_t_Z * V_c * V_c * sqrt((1.0 / (vt * vt)) - (1.0 / (V_c * V_c))) * Bpt;
}
//利用误差因子dF的方法，给出一组参考核的下，目标核的质量结果，和参考核、目标核、γt(C)曲线误差导致的误差 （v6）
void Calculate_Mass_with_err_dF_all(ION_UNKNOWN& ion_t, ION* ions_ref, int ref_n, TGraph* gr_gammat_C, TGraph* gr_gammat_C_u, TGraph* gr_gammat_C_d, TH1D* h_locate, double& mass_result, double& mass_result_err)
{
    double mass_result_R_T = 0;      //由一组参考核得到的 目标核 质量结果 
    double mass_err2_R_T = 0;        //参考核、目标核误差导致的误差(方差)

    double mass_err2_gtC = 0;        //由于γt(C)曲线误差导致的误差(方差),认为其与上面独立
    double mass_result_gtC_u = 0;    //取这两个 相对中心值偏大的那个作为γt(C)曲线的误差
    double mass_result_gtC_d = 0;
    double mass_err2_gtC_u = 0;      //这两个变量没什么用，
    double mass_err2_gtC_d = 0;

    Calculate_Mass_with_err_dF(ion_t, ions_ref, ref_n, gr_gammat_C, ISF, h_locate, mass_result_R_T, mass_err2_R_T);
    if (mass_result_R_T > 0.1)       //中心值正常
    {
        Calculate_Mass_with_err_dF(ion_t, ions_ref, ref_n, gr_gammat_C_u, ISF_u, h_locate, mass_result_gtC_u, mass_err2_gtC_u);
        Calculate_Mass_with_err_dF(ion_t, ions_ref, ref_n, gr_gammat_C_d, ISF_d, h_locate, mass_result_gtC_d, mass_err2_gtC_d);
        if ((mass_result_gtC_u - mass_result_R_T) / (mass_result_R_T - mass_result_gtC_d) > 0)          //  变动γt(C)曲线，对质量结果的影响是单向的。
        {
            mass_err2_gtC = 0.5 * abs(mass_result_gtC_u - mass_result_gtC_d);
        }
        else
        {
            mass_err2_gtC = 0.5 * (abs(mass_result_gtC_u - mass_result_R_T) + abs(mass_result_R_T - mass_result_gtC_d));// γt(C)曲线上下限得到的结果都比中心值大（小）
        }
        //cout << fixed << setiosflags(ios::left) << setw(10) << setprecision(9) << mass_err2_gtC << " " << sqrt(mass_err2_R_T) << endl;
        mass_err2_gtC = mass_err2_gtC * mass_err2_gtC;
        mass_result = mass_result_R_T;
        mass_result_err = sqrt(mass_err2_gtC + mass_err2_R_T);
        return;
    }
    else
    {
        mass_result = 0;           //中心值不正常，跳了
        mass_result_err = 0;
        return;
    }
}



////////////////////=============================
// TGraph 计算目标核Bρ（γt 为常数）
double Calculate_Unknown_Bp_gt0(double Bpi,double Ci,double Ct, int nsteps,double gt_constant)
{// MODE 2
    
    //int nsteps = 1000;
    double deltaC = (Ct - Ci) / nsteps;
    
    double Bp=0.0,C=0.0;
    double gt=0.0;
    
    C = Ci;
    Bp= Bpi;
    for(int i=0;i<nsteps;i++)
    {
        //gt=gr_gammat_C->Eval(C)+d_gt;
        gt= gt_constant;
        Bp = Bp*(1+gt*gt*deltaC/C);
        //Bp = Bp*(1+gr_gammat_C->Eval(C)*deltaC/C);
        C = C+deltaC;
       
    }
    //cout<<"final C check ="<<C- ion_unknown.C<<endl;
    
    return Bp;
}

/////////////////////////////////

//在得到γt(C)曲线时，根据轨道不同，计算目标Bρ
double Calculate_Unknown_Bp(double Bpi,double Ci,double Ct,TGraph* gr_gammat_C)
{// MODE 3
    int nsteps = 100;                         //计算次数
    double deltaC = (Ct - Ci) / nsteps;       //每步步长
    double Bp=0.0,C=0.0;                      //目标核
    double gt=0.0;
    C = Ci;                                   //从参考核的 Bρ C 开始
    Bp= Bpi;
    for(int i=0;i<nsteps;i++)
    {
        gt=gr_gammat_C->Eval(C);              //获取该点γt
        Bp = Bp*(1+gt*gt*deltaC/C);           // 由 γt^2 = (dBρ/Bρ)/(dC/C) 得到
        //Bp = Bp*(1+gr_gammat_C->Eval(C)*deltaC/C);
        C = C+deltaC;
    }
    
    return Bp;
}
double Calculate_Unknown_Bp(double Bpi,double Ci,double Ct,TSpline3* gr_gammat_C)
{// MODE 3
    int nsteps = 100;
    double deltaC = (Ct - Ci) / nsteps;
    double Bp=0.0,C=0.0;
    double gt=0.0;
    C = Ci;
    Bp= Bpi;
    for(int i=0;i<nsteps;i++)
    {
        gt=gr_gammat_C->Eval(C);
        Bp = Bp*(1+gt*gt*deltaC/C);
        //Bp = Bp*(1+gr_gammat_C->Eval(C)*deltaC/C);
        C = C+deltaC;
    }
    return Bp;
}
double Calculate_Unknown_Bp_4(double Bpi,double Ci,double Ct, TF1* f_gammat_C,double& gammat_ave,int nsteps,int inject_number )
{//fitfun use gtC type==4

    double deltaC = (Ct - Ci) / nsteps;
    double Bp=0.0,C=0.0;
    double gt=0.0;
    gammat_ave=0.0;
    C = Ci;
    Bp= Bpi;
    
    for(int i=0;i<nsteps;i++)
    {
        gt=f_gammat_C->Eval(C);
        //gt=gr_gammat_C->Eval(C); 
        //if(gt<1.1||gt>1.4) outfile<<"@@@@@@@@@@@   "<<inject_number<<" "<<Ci<<" "<<Ct<<endl;
        Bp = Bp*(1+gt*gt*deltaC/C);
        //Bp = Bp*(1+gr_gammat_C->Eval(C)*deltaC/C);
        C = C+deltaC;
        gammat_ave+=gt; 
    }
    //outfile.close();
    //cout<<"final C check ="<<C- ion_unknown.C<<endl;
    gammat_ave = gammat_ave/nsteps;
    return Bp;

}
//////////___________________________________________ bp __________________________________________________//////

void From_grr_gtC_to_h_locate(TGraphErrors* grr, TH1D* h_l,TString h_name,double s_Cmin,double s_Cmax)
{
    int n_gtC_chosen_C = grr->GetN();
    double gtC_chosen_CMIN = grr->GetX()[0];
    double gtC_chosen_CMAX = grr->GetX()[grr->GetN()-1];
    h_l = new TH1D(h_name, h_name, n_gtC_chosen_C * dsubregion_n, 
                            gtC_chosen_CMIN - 0.5* ((s_Cmax - s_Cmin) / subregion_n), 
                            gtC_chosen_CMAX + 0.5* ((s_Cmax - s_Cmin) / subregion_n));
}

//========================================= gtC curve ===================================================

//利用散点，得到γt(C)曲线
void GetGtC_Curve(TGraphErrors* grerr_in , TGraphErrors* grerr_avegtC, int s_n, TH1F* h1 )
{
    s_n = subregion_n;                   //将 C 分为 subregion_n 个区间
    int* C_Division_n = new int[s_n];    //此周长区间的粒子计数
    double* avegt     = new double[s_n]; //γt平均值
    double* avegt2    = new double[s_n]; //γt平方平均值（用于计算误差）
    double* sigma_gt  = new double[s_n]; //γt标准差
    int C_region=0;                      //此粒子的C所在的区间
    for(int i=0;i<s_n  ;i++) //初始化
    {
        C_Division_n[i]=0;
        avegt[i]=0;
        avegt2[i]=0;
        sigma_gt[i]=0;
    }

    int grerr_in_n =grerr_in->GetN();//总粒子数
    double x,y;
    for(int i=0;i<grerr_in_n;i++)
    { 
        grerr_in->GetPoint(i,x,y);    //得到该粒子的C γt
        if(y>0&&y<50)
        {
            C_region = h1->FindBin(x)-1;
            if(C_region>s_n-1||C_region<0){continue;}
            C_Division_n[C_region]++;         //此区间计数
            avegt[C_region]     += y;
            avegt2[C_region]    += y*y;
        }
    }

    for(int i=0;i<s_n;i++)//设置γt(C)曲线各点
    {
        if(C_Division_n[i]>3)
        {
            avegt[i] /= C_Division_n[i];
            avegt2[i] /= C_Division_n[i];
            sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;
            grerr_avegtC->SetPoint( grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
            grerr_avegtC->SetPointError(grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;
   
            //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
        }    
//cout<<" C_Division_n["<<i<<"] :"<<C_Division_n[i]<<avegt[i]<<" "<<avegt2[i]<<endl;
    }
}
void GetBpC_sca_points( TGraph* gr_BpC,  int inj_min,int inj_max,double gt_min,double gt_max, int Z_min, int Z_max )
{
    for(int i=0;i<ions_n;i++)
    { 
        if(ionspecies[ions[i].Species].MassUnknown){continue;} // 只使用已知核的B\rho
        if(ions[i].inject_number>=inj_min&&ions[i].inject_number<=inj_max ) 
        {
            if(ions[i].gammat>=gt_min&&ions[i].gammat<=gt_max) 
            {
                if(ions[i].Z>=Z_min&&ions[i].Z<=Z_max)
                {
                    gr_BpC->SetPoint(gr_BpC->GetN(), ions[i].C, ions[i].Bp);
                }

            }
        }
        
    }
}
void GetBp_inj_sca_points( TGraph* gr_Bp_inj,  int inj_min,int inj_max,double gt_min,double gt_max, int Z_min, int Z_max )
{
    for(int i=0;i<ions_n;i++)
    { 
        if(ionspecies[ions[i].Species].MassUnknown){continue;} // 只使用已知核的B\rho
        if(ions[i].inject_number>=inj_min&&ions[i].inject_number<=inj_max ) 
        {
            if(ions[i].gammat>=gt_min&&ions[i].gammat<=gt_max) 
            {
                if(ions[i].Z>=Z_min&&ions[i].Z<=Z_max)
                {
                    gr_Bp_inj->SetPoint(gr_Bp_inj->GetN(), ions[i].inject_number, ions[i].Bp);
                }
            }
        }
    }
}
//使用一注入区间的散点，得到γt(C)曲线
void GetGtC_Curve( TGraphErrors* grerr_avegtC,  TH1F* h1,int inj_min,int inj_max,double gt_min,double gt_max, int Z_min, int Z_max,
    int& ions_inj_n,int& ions_gt_n,int& ions_Z_n, int* count, int& count_n, TH1F* h_C )
{
    int s_n = subregion_n;
    int* C_Division_n = new int[s_n];
    double* avegt     = new double[s_n];
    double* avegt2    = new double[s_n];
    double* sigma_gt  = new double[s_n];
    int C_region=0;
    
    ions_inj_n =0;
    ions_gt_n=0;
    ions_Z_n=0;

    for(int i=0;i<s_n  ;i++)
    {
        C_Division_n[i]=0;
        avegt[i]=0;
        avegt2[i]=0;
        sigma_gt[i]=0;
    }
    double x,y;
    for(int i=0;i<ions_n;i++)
    { 
        if(ions[i].inject_number>=inj_min&&ions[i].inject_number<=inj_max ) 
        {
            ions_inj_n++;
            if(ions[i].gammat>=gt_min&&ions[i].gammat<=gt_max) 
            {
                ions_gt_n++;
                if(ions[i].Z>=Z_min&&ions[i].Z<=Z_max)
                {

                ions_Z_n++;
                C_region = h1->FindBin(ions[i].C)-1;
                if(C_region>subregion_n-1||C_region<0){continue;}
                
                C_Division_n[C_region]++;
                avegt[C_region]     += ions[i].gammat;
                avegt2[C_region]    += ions[i].gammat*ions[i].gammat;

                h_C->Fill(ions[i].C);
                }

            }
        }
        
    }

    for(int i=0;i<s_n;i++)
    {
        if(C_Division_n[i]>3)
        {
            avegt[i] /= C_Division_n[i];
            avegt2[i] /= C_Division_n[i];
            sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;
            grerr_avegtC->SetPoint( grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
            grerr_avegtC->SetPointError(grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;
            count[grerr_avegtC->GetN()-1] = C_Division_n[i];
            
            //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
        }    
//cout<<" C_Division_n["<<i<<"] :"<<C_Division_n[i]<<avegt[i]<<" "<<avegt2[i]<<endl;
    }
    count_n = grerr_avegtC->GetN();
    delete[] C_Division_n;
    delete[] avegt       ;
    delete[] avegt2      ;
    delete[] sigma_gt    ;
}
void GetGtC_Curve_Z_divide( TGraphErrors* grerr_avegtC,  TH1F* h1,int inj_min,int inj_max,double gt_min,double gt_max, int Z_min, int Z_max,int& count_n)
{
    //按照 Zmin~Zmax 筛选
    int s_n = subregion_n;
    int* C_Division_n = new int[s_n];
    double* avegt     = new double[s_n];
    double* avegt2    = new double[s_n];
    double* sigma_gt  = new double[s_n];
    int C_region=0;
    for(int i=0;i<s_n  ;i++)
    {
        C_Division_n[i]=0;avegt[i]=0;avegt2[i]=0;sigma_gt[i]=0;
    }
    double x,y;
    for(int i=0;i<ions_n;i++)
    { 
        if(ions[i].inject_number>=inj_min&&ions[i].inject_number<=inj_max ) 
        {  
            if(ions[i].gammat>=gt_min&&ions[i].gammat<=gt_max) 
            {
                if(ions[i].Z>=Z_min&&ions[i].Z<=Z_max)
                {
                count_n++;
                C_region = h1->FindBin(ions[i].C)-1;
                if(C_region>s_n-1||C_region<0){continue;}
                
                C_Division_n[C_region]++;
                avegt[C_region]     += ions[i].gammat;
                avegt2[C_region]    += ions[i].gammat*ions[i].gammat;
                }
            }
        }
    }
    for(int i=0;i<s_n;i++)
    {
        if(C_Division_n[i]>3)
        {
            avegt[i] /= C_Division_n[i];
            avegt2[i] /= C_Division_n[i];
            sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;
            grerr_avegtC->SetPoint( grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
            grerr_avegtC->SetPointError(grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;

            //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
        }    
//cout<<" C_Division_n["<<i<<"] :"<<C_Division_n[i]<<avegt[i]<<" "<<avegt2[i]<<endl;
    }
    
    delete[] C_Division_n;
    delete[] avegt       ;
    delete[] avegt2      ;
    delete[] sigma_gt    ;
}
void GetGtC_Curve_A_divide( TGraphErrors* grerr_avegtC,  TH1F* h1,int inj_min,int inj_max,double gt_min,double gt_max, int A_min, int A_max,int& count_n)
{
    //按照 Amin~Amax 筛选
    int s_n = subregion_n;
    int* C_Division_n = new int[s_n];
    double* avegt     = new double[s_n];
    double* avegt2    = new double[s_n];
    double* sigma_gt  = new double[s_n];
    int C_region=0;
    for(int i=0;i<s_n  ;i++)
    {
        C_Division_n[i]=0;avegt[i]=0;avegt2[i]=0;sigma_gt[i]=0;
    }
    double x,y;
    for(int i=0;i<ions_n;i++)
    { 
        if(ions[i].inject_number>=inj_min&&ions[i].inject_number<=inj_max ) 
        {  
            if(ions[i].gammat>=gt_min&&ions[i].gammat<=gt_max) 
            {
                if(ions[i].A>=A_min&&ions[i].A<=A_max)
                {
                count_n++;
                C_region = h1->FindBin(ions[i].C)-1;
                if(C_region>s_n-1||C_region<0){continue;}
                
                C_Division_n[C_region]++;
                avegt[C_region]     += ions[i].gammat;
                avegt2[C_region]    += ions[i].gammat*ions[i].gammat;
                }
            }
        }
    }
    for(int i=0;i<s_n;i++)
    {
        if(C_Division_n[i]>3)
        {
            avegt[i] /= C_Division_n[i];
            avegt2[i] /= C_Division_n[i];
            sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;
            grerr_avegtC->SetPoint( grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
            grerr_avegtC->SetPointError(grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;

            //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
        }    
//cout<<" C_Division_n["<<i<<"] :"<<C_Division_n[i]<<avegt[i]<<" "<<avegt2[i]<<endl;
    }
    
    delete[] C_Division_n;
    delete[] avegt       ;
    delete[] avegt2      ;
    delete[] sigma_gt    ;
}

void GetGtC_Curve_step12( TGraphErrors* grerr_avegtC,  TH1F* h1,int inj_min,int inj_max,double gt_min,double gt_max, 
    int& ions_inj_n,int& ions_gt_n,int& ions_choose_n, int* count, int& count_n, TH1F* h_C )
{
    int s_n = subregion_n;
    int* C_Division_n = new int[s_n];
    double* avegt     = new double[s_n];
    double* avegt2    = new double[s_n];
    double* sigma_gt  = new double[s_n];
    int C_region=0;
    
    ions_inj_n =0;
    ions_gt_n=0;
    ions_choose_n=0;

    for(int i=0;i<s_n  ;i++)
    {
        C_Division_n[i]=0;
        avegt[i]=0;
        avegt2[i]=0;
        sigma_gt[i]=0;
    }
    double x,y;
    for(int i=0;i<ions_n;i++)
    { 
        if(ions[i].inject_number>=inj_min&&ions[i].inject_number<=inj_max ) 
        {
            ions_inj_n++;
            if(ions[i].gammat>=gt_min&&ions[i].gammat<=gt_max) 
            {
                ions_gt_n++;
                if( ionspecies[ions[i].Species].gtC_CHOSEN == true  && ions[i].gammat_err<=gtC_ERR_upper_bound)
                {   
                    ions_choose_n++;
                    C_region = h1->FindBin(ions[i].C)-1;
                    if(C_region>subregion_n-1||C_region<0){continue;}
                
                    C_Division_n[C_region]++;
                    avegt[C_region]     += ions[i].gammat;
                    avegt2[C_region]    += ions[i].gammat*ions[i].gammat;

                    h_C->Fill(ions[i].C);
                }

            }
        }
        
    }

    for(int i=0;i<s_n;i++)
    {
        if(C_Division_n[i]>3)
        {
            avegt[i] /= C_Division_n[i];
            avegt2[i] /= C_Division_n[i];
            sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;
            grerr_avegtC->SetPoint( grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
            grerr_avegtC->SetPointError(grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;
            count[grerr_avegtC->GetN()-1] = C_Division_n[i];
            //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
        }    
//cout<<" C_Division_n["<<i<<"] :"<<C_Division_n[i]<<avegt[i]<<" "<<avegt2[i]<<endl;
    }
    count_n = grerr_avegtC->GetN();
    delete[] C_Division_n;
    delete[] avegt       ;
    delete[] avegt2      ;
    delete[] sigma_gt    ;
}

//================================= 20240712 steps =====================

//筛选γt误差过大的点，返回剩余的点数，更新后的γt(C)曲线通过avegt，sigma_gt，C_SR_n传递
int BuildAvegtSR_step1(double* avegt,double* sigma_gt,int* C_SR_n,TH1F* h1_E )
{
    int C_bin;      //start from 1 ~ subregion_n
    int C_region=0; //start from 0
    int ions_left=0;
    int count_out_E_step1 =0;
    for(int i=0;i<subregion_n  ;i++){avegt[i]=0;sigma_gt[i]=0;C_SR_n[i]=0;}//初始化 avegt中心值，sigma_gt误差，C_SR_n计数

    for(int i=0;i<ions_n  ;i++)
    {
        if(ions[i].gammat_err<=gtC_ERR_upper_bound)  // condition for step1 筛选γt误差过大的点，保留误差较小的点
        {
            C_bin = h1_E->FindBin(ions[i].C); // normal:1~subregion_n
            if(C_bin<=0||C_bin>=subregion_n+1){count_out_E_step1++;}//超出范围，不要
            else if(C_bin>=1&&C_bin<=subregion_n)//正常，进行计算
            {
                ions_left++; //chose this one
                C_region= C_bin-1;
                C_SR_n[C_region]++;
                avegt[C_region]     += ions[i].gammat;
                sigma_gt[C_region]    += ions[i].gammat*ions[i].gammat;
            }
            else{ cout<<"this ion C= "<<ions[i].C<<" error!! C_bin = "<<C_bin<<endl;}
        }
    }
    cout<<" BuildAvegtSR_step1 count_out_E_step1: "<<count_out_E_step1<<endl;
    // obtain arrays
    for(int i=0;i<subregion_n  ;i++)
    {
        avegt[i]/=C_SR_n[i]; // mean  γt中心值
        sigma_gt[i]/=C_SR_n[i]; // mean of square 
        sigma_gt[i] = sqrt(sigma_gt[i]-avegt[i]*avegt[i] );   // standard_deviation = sqrt( mean_(x^2) - (mean_x)^2 )
        sigma_gt[i]/=sqrt(double(C_SR_n[i]) );            // sigma of mean = std/sqrt(n)   γt中心值的误差
    }
    return ions_left;
}
//未使用
int BuildAvegtSR_step2(double* avegt,double* sigma_gt,int* C_SR_n,TH1F* h1_E )
{
    int C_bin;      //start from 1 ~ subregion_n
    int C_region=0; //start from 0
    int ions_left=0;
    int count_out_E_step1 =0;
    for(int i=0;i<subregion_n  ;i++){avegt[i]=0;sigma_gt[i]=0;C_SR_n[i]=0;}

    for(int i=0;i<ions_n  ;i++)
    {
        if(ions[i].gammat_err<=gtC_ERR_upper_bound && ionspecies[ions[i].Species].gtC_CHOSEN)  // condition for step2
        {
            C_bin = h1_E->FindBin(ions[i].C); // normal:1~subregion_n
            if(C_bin==0||C_bin==subregion_n+1){count_out_E_step1++;}
            else if(C_bin>=1&&C_bin<=subregion_n)
            {
                ions_left++; //chose this one
                C_region= C_bin-1;
                C_SR_n[C_region]++;
                avegt[C_region]     += ions[i].gammat;
                sigma_gt[C_region]    += ions[i].gammat*ions[i].gammat;
            }
            else{ cout<<"this ion C= "<<ions[i].C<<" error!! C_bin = "<<C_bin<<endl;}
        }
    }
    cout<<" BuildAvegtSR_step2 count_out_E_step1:"<<count_out_E_step1<<endl;
    // obtain arrays
    for(int i=0;i<subregion_n  ;i++)
    {
        avegt[i]/=C_SR_n[i]; // mean  
        sigma_gt[i]/=C_SR_n[i]; // mean of square
        sigma_gt[i] = sqrt(sigma_gt[i]-avegt[i]*avegt[i] );   // standard_deviation = sqrt( mean_(x^2) - (mean_x)^2 )
        sigma_gt[i]/=sqrt(double(C_SR_n[i]) );            // sigma of mean = std/sqrt(n)
    }
    return ions_left;
}

//γt(C)随时间的变化 以及对应的 B\rho C 情况
void Do_gtC_with_time(TH1F* h1,TGraphErrors*grerr_avegtC_inj_1,TGraphErrors*grerr_avegtC_inj_2,TGraphErrors*grerr_avegtC_inj_3,TGraphErrors*grerr_avegtC_inj_4 )
{
    cout<<endl<<endl<<" ================ Do_gtC_with_time() on ========================"<<endl;
    bool Z_min_ON=0;
    bool Show_count_ON=0;
    TString condition;
    int ions_inj_n =0;
    int ions_gt_n=0;
    int ions_Z_n=0;
    int count_n=0;
    int count[300]={0};
    double posx,posy;
    TLatex* lat_n = new TLatex();
    TString lat_text;
    lat_n->SetTextSize(0.04);
    lat_n->SetTextFont(42);
    //gStyle->SetPadBottomMargin(0.05);
    //gStyle->SetPadTopMargin(0.01);
    //gStyle->SetPadLeftMargin(0.1);
    double Yaxis_min = 1.335;
    double Yaxis_max = 1.371;
    int inj_color_1 = kRed;
    int inj_color_2 = kAzure;
    int inj_color_3 = kOrange;
    int inj_color_4 = 8;
    //------------------------------1 
    int inj_min_1 = 2000; gtC_inj_min_1 = inj_min_1;
    int inj_max_1 = 4000; gtC_inj_max_1 = inj_max_1;
    double gtC_min_1 = 1.3;
    double gtC_max_1 = 1.4;
    int Z_min_1=10;
    int Z_max_1=30;
    condition = strtmp.Format("inj:%d~%d, gtC:%.2f~%.2f, %d<=Z<=%d",inj_min_1, inj_max_1,gtC_min_1, gtC_max_1, Z_min_1,Z_max_1);
    
    AxisFormat(grerr_avegtC_inj_1,condition ,""," #gamma_{t}",inj_color_1);
    TH1F* h_C_inj_1 = new TH1F("h_C_inj_1","h_C_inj_1",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
    AxisFormat(h_C_inj_1,""," C [m] "," count",1);
    h_C_inj_1->GetXaxis()->SetTitleOffset(-0.5);
    
    //GetGtC_Curve( grerr_avegtC_inj_1,  h1,  inj_min_1, inj_max_1,   gtC_min_1, gtC_max_1 ,Z_min_1,Z_max_1 ,ions_inj_n,ions_gt_n, ions_Z_n,count,count_n,h_C_inj_1);
    GetGtC_Curve_step12( grerr_avegtC_inj_1,h1, inj_min_1,inj_max_1,  gtC_min_1, gtC_max_1 ,ions_inj_n,ions_gt_n, ions_Z_n,count,count_n,h_C_inj_1 );
    cout<<"|| inj1  n ="<<count_n<<endl;
    
    TGraph* gr_BpC_inj_1 = new TGraph();
    AxisFormat(gr_BpC_inj_1,condition ," C(m) "," B#rho (Tm) ",inj_color_1);
    GetBpC_sca_points(gr_BpC_inj_1,inj_min_1, inj_max_1,   gtC_min_1, gtC_max_1 ,Z_min_1,Z_max_1);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_time_1 = new TCanvas("c_gtC_time_1","c_gtC_time_1",1000,1000);
    c_gtC_time_1->Divide(1,2,0.01,0.01);
    c_gtC_time_1->cd(1);
    grerr_avegtC_inj_1->GetXaxis()->SetRangeUser(128.5,129.0);
    grerr_avegtC_inj_1->GetYaxis()->SetRangeUser(Yaxis_min,Yaxis_max);
    grerr_avegtC_inj_1->DrawClone("apl");
    
    int ions_n_1=ions_gt_n;
    lat_n->DrawLatex(128.803,Yaxis_min+0.3*(Yaxis_max-Yaxis_min),strtmp.Format("inj %d~%d: %d",inj_min_1, inj_max_1,ions_inj_n));
    lat_n->DrawLatex(128.803,Yaxis_min+0.2*(Yaxis_max-Yaxis_min),strtmp.Format("then gtC %.2f~%.2f %d",gtC_min_1, gtC_max_1,ions_gt_n ));
    lat_n->DrawLatex(128.803,Yaxis_min+0.1*(Yaxis_max-Yaxis_min),strtmp.Format("then %d<=Z<=%d : %d",Z_min_1,Z_max_1,ions_Z_n));
    ions_n_1=ions_Z_n;
    c_gtC_time_1->cd(2);
    h_C_inj_1->Draw();
    if(Show_count_ON)
    {
        for(int i=0;i<count_n  ;i++)
        {
            grerr_avegtC_inj_1->GetPoint(i,posx,posy);
            lat_text = strtmp.Format("%d",count[i]);
            lat_n->DrawLatex(posx,2,lat_text);
        }
    }
    h_C_inj_1->GetXaxis()->SetRangeUser(128.5,129.0);
    c_gtC_time_1->Print(FILEPATH+"gtC_injection_1.png");
    c_gtC_time_1->Print(FILEPATH+"gtC_injection_1.root");

    //------------------------------2 
    int inj_min_2 = 6000; gtC_inj_min_2 = inj_min_2;
    int inj_max_2 = 8000; gtC_inj_max_2 = inj_max_2;
    double gtC_min_2 = 1.3;
    double gtC_max_2 = 1.4;
    int Z_min_2=10;
    int Z_max_2=30;
    condition = strtmp.Format("inj:%d~%d, gtC:%.2f~%.2f, %d<=Z<=%d",inj_min_2, inj_max_2,gtC_min_2, gtC_max_2, Z_min_2,Z_max_2);
    
    AxisFormat(grerr_avegtC_inj_2,condition ,""," #gamma_{t}",inj_color_2);
    TH1F* h_C_inj_2 = new TH1F("h_C_inj_2","h_C_inj_2",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
    AxisFormat(h_C_inj_2,""," C [m] "," count",1);
    h_C_inj_2->GetXaxis()->SetTitleOffset(-0.5);
    
    //GetGtC_Curve( grerr_avegtC_inj_2,  h1,  inj_min_2, inj_max_2,   gtC_min_2, gtC_max_2 ,Z_min_2,Z_max_2 ,ions_inj_n,ions_gt_n, ions_Z_n,count,count_n,h_C_inj_2);
    GetGtC_Curve_step12( grerr_avegtC_inj_2,h1, inj_min_2,inj_max_2,  gtC_min_2, gtC_max_2 ,ions_inj_n,ions_gt_n, ions_Z_n,count,count_n,h_C_inj_2 );
    
    cout<<"|| inj2  n ="<<count_n<<endl;
    
    TGraph* gr_BpC_inj_2 = new TGraph();
    AxisFormat(gr_BpC_inj_2,condition ," C(m) "," B#rho (Tm) ",inj_color_2);
    GetBpC_sca_points(gr_BpC_inj_2,inj_min_2, inj_max_2,   gtC_min_2, gtC_max_2 ,Z_min_2,Z_max_2);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_time_2 = new TCanvas("c_gtC_time_2","c_gtC_time_2",1000,1000);
    c_gtC_time_2->Divide(1,2,0.01,0.01);
    c_gtC_time_2->cd(1);
    grerr_avegtC_inj_2->GetXaxis()->SetRangeUser(128.5,129.0);
    grerr_avegtC_inj_2->GetYaxis()->SetRangeUser(Yaxis_min,Yaxis_max);
    grerr_avegtC_inj_2->DrawClone("apl");
    
    int ions_n_2=ions_gt_n;
    lat_n->DrawLatex(128.803,Yaxis_min+0.3*(Yaxis_max-Yaxis_min),strtmp.Format("inj %d~%d: %d",inj_min_2, inj_max_2,ions_inj_n));
    lat_n->DrawLatex(128.803,Yaxis_min+0.2*(Yaxis_max-Yaxis_min),strtmp.Format("then gtC %.2f~%.2f %d",gtC_min_2, gtC_max_2,ions_gt_n ));
    lat_n->DrawLatex(128.803,Yaxis_min+0.1*(Yaxis_max-Yaxis_min),strtmp.Format("then %d<=Z<=%d : %d",Z_min_2,Z_max_2,ions_Z_n));
    ions_n_2=ions_Z_n;
    c_gtC_time_2->cd(2);
    h_C_inj_2->Draw();
    if(Show_count_ON)
    {
        for(int i=0;i<count_n  ;i++)
        {
            grerr_avegtC_inj_2->GetPoint(i,posx,posy);
            lat_text = strtmp.Format("%d",count[i]);
            lat_n->DrawLatex(posx,2,lat_text);
        }
    }
    h_C_inj_2->GetXaxis()->SetRangeUser(128.5,129.0);
    c_gtC_time_2->Print(FILEPATH+"gtC_injection_2.png");
    c_gtC_time_2->Print(FILEPATH+"gtC_injection_2.root");
    
    
    //------------------------------3
    int inj_min_3 = 6000; gtC_inj_min_3 = inj_min_3;
    int inj_max_3 = 8000; gtC_inj_max_3 = inj_max_3;
    double gtC_min_3 = 1.3;
    double gtC_max_3 = 1.4;
    int Z_min_3=10;
    int Z_max_3=30;
    condition = strtmp.Format("inj:%d~%d, gtC:%.2f~%.2f, %d<=Z<=%d",inj_min_3, inj_max_3,gtC_min_3, gtC_max_3, Z_min_3,Z_max_3);
   
    AxisFormat(grerr_avegtC_inj_3,condition ,""," #gamma_{t}",inj_color_3);
    TH1F* h_C_inj_3 = new TH1F("h_C_inj_3","h_C_inj_3",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
    AxisFormat(h_C_inj_3,""," C [m] "," count",1);
    h_C_inj_3->GetXaxis()->SetTitleOffset(-0.5);
    GetGtC_Curve( grerr_avegtC_inj_3,  h1,  inj_min_3, inj_max_3,   gtC_min_3, gtC_max_3 ,Z_min_3,Z_max_3 ,ions_inj_n,ions_gt_n, ions_Z_n,count,count_n,h_C_inj_3);
    cout<<"|| inj3  n ="<<count_n<<endl;
    TGraph* gr_BpC_inj_3 = new TGraph();
    AxisFormat(gr_BpC_inj_3,condition ," C(m) "," B#rho (Tm) ",inj_color_3);
    GetBpC_sca_points(gr_BpC_inj_3,inj_min_3, inj_max_3,   gtC_min_3, gtC_max_3 ,Z_min_3,Z_max_3);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_time_3 = new TCanvas("c_gtC_time_3","c_gtC_time_3",1000,1000);
    c_gtC_time_3->Divide(1,2,0.01,0.01);
    c_gtC_time_3->cd(1);
    grerr_avegtC_inj_3->GetXaxis()->SetRangeUser(128.5,129.0);
    grerr_avegtC_inj_3->GetYaxis()->SetRangeUser(Yaxis_min,Yaxis_max);
    grerr_avegtC_inj_3->DrawClone("apl");
    
    int ions_n_3=ions_gt_n;
    lat_n->DrawLatex(128.803,Yaxis_min+0.3*(Yaxis_max-Yaxis_min),strtmp.Format("inj %d~%d: %d",inj_min_3, inj_max_3,ions_inj_n));
    lat_n->DrawLatex(128.803,Yaxis_min+0.2*(Yaxis_max-Yaxis_min),strtmp.Format("then gtC %.2f~%.2f %d",gtC_min_3, gtC_max_3,ions_gt_n ));
    lat_n->DrawLatex(128.803,Yaxis_min+0.1*(Yaxis_max-Yaxis_min),strtmp.Format("then %d<=Z<=%d : %d",Z_min_3,Z_max_3,ions_Z_n));
    ions_n_3=ions_Z_n;
    c_gtC_time_3->cd(2);
    h_C_inj_3->Draw();
    if(Show_count_ON)
    {
        for(int i=0;i<count_n  ;i++)
        {
            grerr_avegtC_inj_3->GetPoint(i,posx,posy);
            lat_text = strtmp.Format("%d",count[i]);
            lat_n->DrawLatex(posx,2,lat_text);
        }
    }
    h_C_inj_3->GetXaxis()->SetRangeUser(128.5,129.0);
    c_gtC_time_3->Print(FILEPATH+"gtC_injection_3.png");
    c_gtC_time_3->Print(FILEPATH+"gtC_injection_3.root");

    //------------------------------4
    int inj_min_4 = 8000;  gtC_inj_min_4 = inj_min_4;
    int inj_max_4 = 10000; gtC_inj_max_4 = inj_max_4;
    double gtC_min_4 = 1.3;
    double gtC_max_4 = 1.4;
    int Z_min_4=10;
    int Z_max_4=30;
    condition = strtmp.Format("inj:%d~%d, gtC:%.2f~%.2f, %d<=Z<=%d",inj_min_4, inj_max_4,gtC_min_4, gtC_max_4, Z_min_4,Z_max_4);
    
    AxisFormat(grerr_avegtC_inj_4,condition ,""," #gamma_{t}",inj_color_4);
    TH1F* h_C_inj_4 = new TH1F("h_C_inj_4","h_C_inj_4",h1->GetNbinsX(),h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
    AxisFormat(h_C_inj_4,""," C [m] "," count",1);
    h_C_inj_4->GetXaxis()->SetTitleOffset(-0.5);
    GetGtC_Curve( grerr_avegtC_inj_4,  h1,  inj_min_4, inj_max_4,   gtC_min_4, gtC_max_4 ,Z_min_4,Z_max_4 ,ions_inj_n,ions_gt_n, ions_Z_n,count,count_n,h_C_inj_4);
    cout<<"|| inj4  n ="<<count_n<<endl;
    TGraph* gr_BpC_inj_4 = new TGraph();
    AxisFormat(gr_BpC_inj_4,condition ," C(m) "," B#rho (Tm) ",inj_color_4);
    GetBpC_sca_points(gr_BpC_inj_4,inj_min_4, inj_max_4,   gtC_min_4, gtC_max_4 ,Z_min_4,Z_max_4);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_time_4 = new TCanvas("c_gtC_time_4","c_gtC_time_4",1000,1000);
    c_gtC_time_4->Divide(1,2,0.01,0.01);
    c_gtC_time_4->cd(1);
    grerr_avegtC_inj_4->GetXaxis()->SetRangeUser(128.5,129.0);
    grerr_avegtC_inj_4->GetYaxis()->SetRangeUser(Yaxis_min,Yaxis_max);
    grerr_avegtC_inj_4->DrawClone("apl");
    
    int ions_n_4=ions_gt_n;
    lat_n->DrawLatex(128.803,Yaxis_min+0.3*(Yaxis_max-Yaxis_min),strtmp.Format("inj %d~%d: %d",inj_min_4, inj_max_4,ions_inj_n));
    lat_n->DrawLatex(128.803,Yaxis_min+0.2*(Yaxis_max-Yaxis_min),strtmp.Format("then gtC %.2f~%.2f %d",gtC_min_4, gtC_max_4,ions_gt_n ));
    lat_n->DrawLatex(128.803,Yaxis_min+0.1*(Yaxis_max-Yaxis_min),strtmp.Format("then %d<=Z<=%d : %d",Z_min_4,Z_max_4,ions_Z_n));
    ions_n_4=ions_Z_n;
    c_gtC_time_4->cd(2);
    h_C_inj_4->Draw();
    if(Show_count_ON)
    {
        for(int i=0;i<count_n  ;i++)
        {
            grerr_avegtC_inj_4->GetPoint(i,posx,posy);
            lat_text = strtmp.Format("%d",count[i]);
            lat_n->DrawLatex(posx,2,lat_text);
        }
    }
    h_C_inj_4->GetXaxis()->SetRangeUser(128.5,129.0);
    c_gtC_time_4->Print(FILEPATH+"gtC_injection_4.png");
    c_gtC_time_4->Print(FILEPATH+"gtC_injection_4.root");

    //------------------------------------- compare
    //gStyle->SetPadBottomMargin(0.1);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_time_all = new TCanvas("c_gtC_time_all","c_gtC_time_all",1800,1000);
    //grerr_avegtC_inj_1->SetTitle(strtmp.Format("inj:%d~%d, gtC:%.2f~%.2f",1, total_injection,gtC_min_1, gtC_max_1) );
    grerr_avegtC_inj_1->GetXaxis()->SetTitle("C [m]");
    grerr_avegtC_inj_1->GetYaxis()->SetRangeUser(Yaxis_min,Yaxis_max);

    grerr_avegtC_inj_1->Draw("apl");
    grerr_avegtC_inj_2->Draw("plsame");
    //grerr_avegtC_inj_3->Draw("plsame");
    //grerr_avegtC_inj_4->Draw("plsame");

    auto legend_gtC_time = new TLegend(0.75,0.10,0.90,0.35); //downright
    // 随注入变化 
    legend_gtC_time->AddEntry(grerr_avegtC_inj_1,strtmp.Format("injection: %d~%d",inj_min_1,inj_max_1),"ple" );
    legend_gtC_time->AddEntry(grerr_avegtC_inj_2,strtmp.Format("injection: %d~%d",inj_min_2,inj_max_2),"ple" );
    //legend_gtC_time->AddEntry(grerr_avegtC_inj_3,strtmp.Format("injection: %d~%d",inj_min_3,inj_max_3),"ple" );
    //legend_gtC_time->AddEntry(grerr_avegtC_inj_4,strtmp.Format("injection: %d~%d",inj_min_4,inj_max_4),"ple" );
    
    /*
    // 筛选条件 逐步展示
    legend_gtC_time->AddEntry(grerr_avegtC_inj_1,strtmp.Format("all ions :%d",ions_n_1),"ple" );
    legend_gtC_time->AddEntry(grerr_avegtC_inj_2,strtmp.Format("Z>=9 :%d",ions_n_2),"ple" );
    legend_gtC_time->AddEntry(grerr_avegtC_inj_3,strtmp.Format("1.3<#gamma_{t}<1.4 :%d",ions_n_3),"ple" );
    legend_gtC_time->AddEntry(grerr_avegtC_inj_4,strtmp.Format("Z>=9 && 1.3<#gamma_{t}<1.4 :%d",ions_n_4),"ple" );
    */

    legend_gtC_time->Draw("same");
    c_gtC_time_all->Print(FILEPATH+"c_gtC_time_all.png");
    c_gtC_time_all->Print(FILEPATH+"c_gtC_time_all.root");
    delete lat_n;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_BpC_inj_all = new TCanvas("c_BpC_inj_all","c_BpC_inj_all",1800,1000);
    gr_BpC_inj_1->SetTitle(" ");
    gr_BpC_inj_1->SetMarkerSize(1);
    gr_BpC_inj_2->SetMarkerSize(1);
    gr_BpC_inj_3->SetMarkerSize(1);
    gr_BpC_inj_4->SetMarkerSize(1);
    gr_BpC_inj_1->Draw("ap");
    gr_BpC_inj_2->Draw("psame");
    //gr_BpC_inj_3->Draw("psame");
    //gr_BpC_inj_4->Draw("psame");
    auto legend_BpC_time = new TLegend(0.75,0.10,0.90,0.35); //downright
    // 随注入变化 
    legend_BpC_time->AddEntry(gr_BpC_inj_1,strtmp.Format("injection: %d~%d",inj_min_1,inj_max_1),"p" );
    legend_BpC_time->AddEntry(gr_BpC_inj_2,strtmp.Format("injection: %d~%d",inj_min_2,inj_max_2),"p" );
    //legend_BpC_time->AddEntry(gr_BpC_inj_3,strtmp.Format("injection: %d~%d",inj_min_3,inj_max_3),"p" );
    //legend_BpC_time->AddEntry(gr_BpC_inj_4,strtmp.Format("injection: %d~%d",inj_min_4,inj_max_4),"p" );
    legend_BpC_time->Draw("same");

    TGraph_to_outfile(FILEPATH+"gr_BpC_inj_1_data.txt",gr_BpC_inj_1,7);
    TGraph_to_outfile(FILEPATH+"gr_BpC_inj_2_data.txt",gr_BpC_inj_2,7);

    //============================ 58Ni 红蓝两段 Bp-inj 散点图 =================================
    TGraph* gr_Bp_inj_1 = new TGraph();
    AxisFormat(gr_Bp_inj_1,condition ," injection "," B#rho (Tm) ",inj_color_1);
    gr_Bp_inj_1->SetMarkerSize(1);
    GetBp_inj_sca_points(gr_Bp_inj_1,inj_min_1, inj_max_1,   gtC_min_1, gtC_max_1 ,Z_min_1,Z_max_1);
    TGraph* gr_Bp_inj_2 = new TGraph();
    AxisFormat(gr_Bp_inj_2,condition ," injection "," B#rho (Tm) ",inj_color_2);
    gr_Bp_inj_2->SetMarkerSize(1);
    GetBp_inj_sca_points(gr_Bp_inj_2,inj_min_2, inj_max_2,   gtC_min_2, gtC_max_2 ,Z_min_2,Z_max_2);
    TGraph* gr_Bp_inj_all = new TGraph();
    AxisFormat(gr_Bp_inj_all,condition ," injection "," B#rho (Tm) ");
    gr_Bp_inj_all->SetMarkerSize(1);
    GetBp_inj_sca_points(gr_Bp_inj_all,0, 9999999,   gtC_min_1, gtC_max_1 ,Z_min_1,Z_max_1);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_Bp_inj_redblue = new TCanvas("c_Bp_inj_redblue","c_Bp_inj_redblue",1800,1000);
    gr_Bp_inj_all->Draw("ap");
    gr_Bp_inj_1->Draw("psame");
    gr_Bp_inj_2->Draw("psame");

    TGraph_to_outfile(FILEPATH+"gr_Bp_inj_all_data.txt",gr_Bp_inj_all);
    TGraph_to_outfile(FILEPATH+"gr_Bp_inj_1_data.txt",gr_Bp_inj_1);
    TGraph_to_outfile(FILEPATH+"gr_Bp_inj_2_data.txt",gr_Bp_inj_2);




    //============================  2 段 比较  用于 36Ar set2=====================================
    /*
    int inj_divide_a1=0;int inj_divide_a2=2261;
    int inj_divide_b1=20000;int inj_divide_b2=22000;
    TGraph* gr_Bp_time_BG = new TGraph();
    AxisFormat(gr_Bp_time_BG,""," time ", " B#rho (Tm)");
    TGraph* gr_Bp_time_a = new TGraph();
    AxisFormat(gr_Bp_time_a,""," time ", " B#rho (Tm)",kRed);
    TGraph* gr_Bp_time_b = new TGraph();
    AxisFormat(gr_Bp_time_b,""," time ", " B#rho (Tm)",kAzure);

    gr_Bp_time_BG->SetMarkerSize(1);
    gr_Bp_time_a->SetMarkerSize(1);
    gr_Bp_time_b->SetMarkerSize(1);
    for(int i=0;i<ions_n  ;i++)
    {
        if(ionspecies[ions[i].Species].Aname=="10C")
        {
            gr_Bp_time_BG->SetPoint(gr_Bp_time_BG->GetN(),ions[i].time, ions[i].Bp);
            if(ions[i].inject_number>=inj_divide_a1&&ions[i].inject_number<=inj_divide_a2)
            {
                gr_Bp_time_a->SetPoint(gr_Bp_time_BG->GetN(),ions[i].time, ions[i].Bp);
            }
            if(ions[i].inject_number>=inj_divide_b1&&ions[i].inject_number<=inj_divide_b2)
            {
                gr_Bp_time_b->SetPoint(gr_Bp_time_BG->GetN(),ions[i].time, ions[i].Bp);
            }
        }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_Bp_time_ab = new TCanvas("c_Bp_time_ab","c_Bp_time_ab",1800,1000);
    gr_Bp_time_BG->Draw("Ap");
    gr_Bp_time_a->Draw("psame");
    gr_Bp_time_b->Draw("psame");
    */

    //============================= graph 做差 ===============================
    /*
    double dd_limit=0.001; //控制点的对齐容忍度
    TGraphErrors* grerr_avegtC_inj_d = new TGraphErrors();
    AxisFormat(grerr_avegtC_inj_d,"" ,"C"," #gamma_{t}");
    
    Get_difference_of_two_TGraphErrors(grerr_avegtC_inj_1,grerr_avegtC_inj_2,grerr_avegtC_inj_d,dd_limit);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_dd = new TCanvas("c_gtC_dd","c_gtC_dd",1800,1000);
    grerr_avegtC_inj_d->Draw("Ap");
    */

    
}











//=============================== gtC -ALL ============================================================
void Draw_gtC_all(int num, int opt) //opt==0 all , opt==1 dA0T chosen
{
    double Cmin=128.70;
    double Cmax=128.85;
    double time_divide = 900000;

TH2F* h2_gammat_C = new TH2F ("h2_gammat_C"+strtmp.Format("_%d",num),"h2_gammat_C"+strtmp.Format("_%d",num), 250,128.5,129.0, 100,1.3,1.4);
AxisFormat(h2_gammat_C,""," C [m] ","#gamma_{t} ");
for(int i=0;i< ions_n ;i++)
{
    if(opt==0)h2_gammat_C->Fill(ions[i].C,ions[i].gammat);  // 不筛选 原始数据所有点
    else if(opt==1)  //Do_dA0_T 用到的
    {
        if(ions[i].Do_dA0_T_flag){h2_gammat_C->Fill(ions[i].C,ions[i].gammat);}
    }
    else if(opt==2)  // cmin~cmax
    {
        if(ions[i].C>Cmin&&ions[i].C<Cmax){h2_gammat_C->Fill(ions[i].C,ions[i].gammat);}
    }
    else if(opt==3)  // 注入时间筛选
    {
        if(ions[i].time<time_divide){h2_gammat_C->Fill(ions[i].C,ions[i].gammat);}
    }
    else if(opt==4)  //种类筛选
    {
        if(ionspecies[ions[i].Species].Aname=="8B"||ionspecies[ions[i].Species].Aname=="13O"){continue;}
        else    {h2_gammat_C->Fill(ions[i].C,ions[i].gammat);}
    }
    else if(opt==5)  // 筛选
    {
        //if(ions[i].Z>=7&&ions[i].time>time_divide)
        if(ions[i].Z<=6&&ions[i].time>time_divide)
        {h2_gammat_C->Fill(ions[i].C,ions[i].gammat);}
    }
    else if(opt==6)  // gtC condition
    {
        if( (!(ions[i].Z<=CONDITION_gt_Z&&ions[i].A<=CONDITION_gt_A) )&&(ionspecies[ions[i].Species].MassUnknown==0) )
        {h2_gammat_C->Fill(ions[i].C,ions[i].gammat);}
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_h2_gammat_C = new TCanvas("c_h2_gammat_C"+strtmp.Format("_%d",num),"c_h2_gammat_C"+strtmp.Format("_%d",num),1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h2_gammat_C->Draw("colz");
}
//_______________________________ gtC -ALL ___________________________________________________________

//=============================== BpC -ALL ============================================================
void Draw_BpC_all(int num, int opt) //opt==0 all , opt==1 dA0T chosen
{
    double Cmin=128.70;
    double Cmax=128.85;
    double time_divide = 900000;

TH2F* h2_lnBp_lnC = new TH2F ("h2_lnBp_lnC"+strtmp.Format("_%d",num),"h2_lnBp_lnC"+strtmp.Format("_%d",num), 250,log(128.5),log(129.0), 100,log(4.7),log(5.0));
AxisFormat(h2_lnBp_lnC,""," C [m] ","#gamma_{t} ");
for(int i=0;i< ions_n ;i++)
{
    if(ionspecies[ions[i].Species].MassUnknown==0)
    {
        h2_lnBp_lnC->Fill(log(ions[i].C),log(ions[i].Bp));
    }
}
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_h2_lnBp_lnC = new TCanvas("c_h2_lnBp_lnC"+strtmp.Format("_%d",num),"c_h2_lnBp_lnC"+strtmp.Format("_%d",num),1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h2_lnBp_lnC->Draw("colz");
}




//======================================= calculate mass ========================================================

double Calculate_Mass(ION_UNKNOWN& ion_t,ION* ions_ref, int ref_n,int j_now, TGraph* gr_gammat_C,int k1,int k2)
{   
    double v_original,T_original,vi_original,Ti_original;          // i :known ref ions
    double v,vi,v_err,vi_err,T,T_err,Ti,Ti_err,mi,mi_err;
    double v_ns,v_err_ns;
    double Ct,Ci,C_err,Ci_err;
    double A,Ave_Bp,s;
    double Bp,Bpi,Bp_i=0.0;            //Bp_i calculated value using each ref
    double mass=0;
    v_original=ion_t.v; 
    T_original=ion_t.T;
    v_err=ion_t.v_err;
    T_err=ion_t.T_err;
    switch(k1)
    {
        case 0:v=v_original      ; T=T_original;       break;
        case 1:v=v_original+v_err; T=T_original;       break;
        case 2:v=v_original-v_err; T=T_original;       break;
        case 3:v=v_original      ; T=T_original+T_err ;break;
        case 4:v=v_original      ; T=T_original-T_err ;break;
        default:cout<<"default!!";break;
    }
    //now v ,T can be used

    v_ns=v/1000000000.0;
    v_err_ns=v_err/1000000000.0;
    Ct=v*T*0.000000000001;   Bp=ion_t.Bp;
    A=ion_t.Z*u/k0;
    s=sqrt(1.0/v_ns/v_ns - 1.0/V_c/V_c);
    
    
    for(int j=0;j<ref_n  ;j++)
    {
        vi_original=ions_ref[j].v;
        Ti_original=ions_ref[j].T;
        vi_err = ions_ref[j].v_err;
        Ti_err = ions_ref[j].T_err;
        if(j==j_now)
        {
            switch(k2)
            {
                case 0:vi=vi_original       ; Ti=Ti_original       ; break;
                case 1:vi=vi_original+vi_err; Ti=Ti_original       ; break;
                case 2:vi=vi_original-vi_err; Ti=Ti_original       ; break;
                case 3:vi=vi_original;        Ti=Ti_original+Ti_err; break;
                case 4:vi=vi_original       ; Ti=Ti_original-Ti_err; break;
                default:cout<<"default!!";break;
            }
        }
        else
        {
            vi=vi_original       ; Ti=Ti_original       ; 
        }
        
        //now vi ,Ti can be used
        Ci=vi*Ti*0.000000000001;
        mi=ionspecies[ions_ref[j].Species].Mass;
        Bpi= (mi/ions_ref[j].Z)*1000.0*vi/(Vc_ms*sqrt(Vc_ms*Vc_ms-vi*vi));
        Bp_i = Calculate_Unknown_Bp(Bpi,Ci,Ct,gr_gammat_C); 
        
        mass+=A * Bp_i * s;
    }
        
    return mass/ref_n;   
    
}

//================================== LEVEL 1 one mass value from 1t+1r =================================
//根据一个参考核计算目标核质量
double Calculate_one_Mass(int ion_t_Z, int ion_r_Z, double a1t,double a2t,double a5t,double a6t,double a1r,double a2r,double a5r,double a6r, 
    double m_ref,TGraph* gr_gammat_C)
{
    // t= target r= reference

    double Tt,vt,Ct,Bpt;
    double Tr,vr,Cr,Bpr;
    double mass=0;
    vt = L/(a5t*0.001+ddT);vr = L/(a5r*0.001+ddT);              //两个核的速度
    Tt = (a1t-a6t*0.5)*0.001;Tr = (a1r-a6r*0.5)*0.001;          //两个核的周期
    Ct = vt*Tt; Cr = vr*Tr;                                     //周长
    Bpr = m_ref*0.000001*vr/(ion_r_Z*V_c*sqrt(V_c*V_c-vr*vr ) );   //参考核的Bρ
    if(use_gtC_type==1)Bpt = Calculate_Unknown_Bp(Bpr,Cr,Ct,gr_gammat_C);          // ！主要流程，根据得到的γt(C)曲线，两核的轨道长度，计算目标核Bρ
    else if(use_gtC_type==2)Bpt = Calculate_Unknown_Bp_gt0(Bpr,Cr, Ct ,100,gt0);
    else if(use_gtC_type==3)Bpt = Calculate_Unknown_Bp(Bpr,Cr,Ct,gr_gammat_C);

//cout<<"ion_t_Z= "<<ion_t_Z<<"ion_r_Z= "<<ion_r_Z
//<<" vt = "<<vt<<" vr= "<<vr<<" Tt= "<<Tt<<" Tr= "<<Tr<<" Ct= "<<Ct<<" Cr ="<<Cr<<endl
//<<" Bpt= "<<Bpt<<" Bpr= "<<Bpr<<endl;     //debug
    mass = 1000000.0*ion_t_Z*V_c*V_c*sqrt( (1.0/(vt*vt)) - (1.0/(V_c*V_c)))*Bpt;
//cout<<" DEBUG : Calculate_one_Mass return: "<<mass<<endl; //debug
    return mass;
}
//利用一个参考核，计算目标核质量和Bρ（Bρ通过 double& this_Bpt 返回）
double Calculate_one_Mass_thisBp(int ion_t_Z, int ion_r_Z, double a1t,double a2t,double a5t,double a6t,double a1r,double a2r,double a5r,double a6r, 
    double m_ref,TGraph* gr_gammat_C,double& this_Bpt)
{
    // t= target r= reference

    double Tt,vt,Ct,Bpt;
    double Tr,vr,Cr,Bpr;
    double mass=0;
    vt = L/(a5t*0.001+ddT);vr = L/(a5r*0.001+ddT);               //两个核的速度
    Tt = (a1t-a6t*0.5)*0.001;Tr = (a1r-a6r*0.5)*0.001;           //两个核的周期
    Ct = vt*Tt; Cr = vr*Tr;                                      //周长
    Bpr = m_ref*0.000001*vr/(ion_r_Z*V_c*sqrt(V_c*V_c-vr*vr ) );    //参考核的Bρ
    
    if(use_gtC_type==1)Bpt = Calculate_Unknown_Bp(Bpr,Cr,Ct,gr_gammat_C);         // ！主要流程，根据得到的γt(C)曲线，两核的轨道长度，计算目标核Bρ
    else if(use_gtC_type==2)Bpt = Calculate_Unknown_Bp_gt0(Bpr,Cr, Ct ,100,gt0);
    else if(use_gtC_type==3)Bpt = Calculate_Unknown_Bp(Bpr,Cr,Ct,gr_gammat_C);

    
    this_Bpt = Bpt;
//cout<<"ion_t_Z= "<<ion_t_Z<<"ion_r_Z= "<<ion_r_Z
//<<" vt = "<<vt<<" vr= "<<vr<<" Tt= "<<Tt<<" Tr= "<<Tr<<" Ct= "<<Ct<<" Cr ="<<Cr<<endl
//<<" Bpt= "<<Bpt<<" Bpr= "<<Bpr<<endl;     //debug
    mass = 1000000.0*ion_t_Z*V_c*V_c*sqrt( (1.0/(vt*vt)) - (1.0/(V_c*V_c)))*Bpt;      //目标核质量
//cout<<" DEBUG : Calculate_one_Mass return: "<<mass<<endl; //debug
    return mass;
}


//==================================LEVEL 2 error of one mass calculation(1t+1r)=================================

//利用一个参考核得到目标核质量，计算此结果误差，采用变动A系数求偏微分的方法
double Calculate_one_Mass_err_usePDCOV(ION_UNKNOWN& ion_t, ION& ion_ref, TGraph* gr_gammat_C, double lr,double mass_center)
{
    // use PD partial derivative and COV covariance
    // one target + one ref --- 9 parameters ---- 21 errors
    double one_mass_err=0;
    int ion_t_Z = ion_t.Z;
    int ion_r_Z = ion_ref.Z;
    double tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref = 0;//tA目标核系数，rA参考核系数。m_ref参考核质量
    double errors[21]={0};   // 一次计算，9个参数，参考与目标各有6项协方差， 9+2*6 =21 项 误差来源 
    double errors_total;
    double tA1err,tA2err,tA5err,tA6err,rA1err,rA2err,rA5err,rA6err,m_ref_err = 0;  // sqrt_variance
    double at_err[4]={0};double ar_err[5]={0};                                     // variance
    double t_cov[6]={0};   double r_cov[6]={0};
    double PD_t[4]={0};  // partial derivative
    double PD_r[5]={0};  // partial derivative of ref ions 最后第五项偏导[4]是来自已知质量
     

    // target 4 paras
    tA1=ion_t.A1;tA2=ion_t.A2;tA5=ion_t.dA0;tA6=ion_t.dA1;                                                     //目标核各项系数
    // target 4 para err-variance
    at_err[0]=ion_t.A1err;at_err[1]=ion_t.A2err;at_err[2]=ion_t.dA0err;at_err[3]=ion_t.dA1err;                 //目标核各项误差的平方
    // target 4 para err 
    tA1err=sqrt(ion_t.A1err);tA2err=sqrt(ion_t.A2err);tA5err=sqrt(ion_t.dA0err);tA6err=sqrt(ion_t.dA1err);     //开方，标准差
    // target 6 para covariance
    t_cov[0] = ion_t.cov12;t_cov[1] = ion_t.cov15;t_cov[2] = ion_t.cov16;                                      //协方差
    t_cov[3] = ion_t.cov25;t_cov[4] = ion_t.cov26;t_cov[5] = ion_t.cov56;
    
    // ref 5 paras
    rA1=ion_ref.A1; rA2=ion_ref.A2; rA5=ion_ref.dA0; rA6=ion_ref.dA1;                                          //参考核各项系数
    m_ref = ionspecies[ion_ref.Species].Mass;                                                                  //AME质量值
    // ref 5 para err-variances
    ar_err[0]=ion_ref.A1err;ar_err[1]=ion_ref.A2err;ar_err[2]=ion_ref.dA0err;ar_err[3]=ion_ref.dA1err;         //参考核各项误差
    ar_err[4]=ionspecies[ion_ref.Species].AME_err*ionspecies[ion_ref.Species].AME_err;                         //质量误差
    
    // ref 5 para err 
    rA1err=sqrt(ion_ref.A1err);rA2err=sqrt(ion_ref.A2err);rA5err=sqrt(ion_ref.dA0err);rA6err=sqrt(ion_ref.dA1err);//开方
    m_ref_err = ionspecies[ion_ref.Species].AME_err;  // sigma not square
    r_cov[0]=ion_ref.cov12;r_cov[1]=ion_ref.cov15;r_cov[2]=ion_ref.cov16;                                      //协方差
    r_cov[3]=ion_ref.cov25;r_cov[4]=ion_ref.cov26;r_cov[5]=ion_ref.cov56;


    //mass_center=Calculate_one_Mass(ion_t_Z,tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    //质量中心值由此函数获取，那么各参数进行一个微小的变化，质量值也会有相应的变化，则 Δm/ΔA 即为各系数的偏导数
    // 
    // calculate partial derivative of 9 parameters                 9个参数的偏导数
    PD_t[0]= (Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1+lr*tA1err, tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*tA1err) ;
    PD_t[1]= (Calculate_one_Mass(ion_t_Z,ion_r_Z,tA1, tA2+lr*tA2err,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*tA2err) ;
    PD_t[2]= (Calculate_one_Mass(ion_t_Z,ion_r_Z,tA1,tA2, tA5+lr*tA5err,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*tA5err) ;
    PD_t[3]= (Calculate_one_Mass(ion_t_Z,ion_r_Z,tA1,tA2,tA5,tA6+lr*tA6err ,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*tA6err) ;

    PD_r[0]= (Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1,tA2,tA5,tA6,  rA1+lr*rA1err,rA2,rA5,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*rA1err) ;
    PD_r[1]= (Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1,tA2,tA5,tA6, rA1, rA2+lr*rA2err,rA5,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*rA2err) ;
    PD_r[2]= (Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1,tA2,tA5,tA6, rA1,rA2, rA5+lr*rA5err,rA6,m_ref,gr_gammat_C) - mass_center)/(lr*rA5err) ;
    PD_r[3]= (Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1,tA2,tA5,tA6, rA1,rA2,rA5, rA6+lr*rA6err,m_ref,gr_gammat_C) - mass_center)/(lr*rA6err) ;
    PD_r[4]= (Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1,tA2,tA5,tA6, rA1,rA2,rA5,rA6,m_ref+lr*m_ref_err,gr_gammat_C) - mass_center)/(lr*m_ref_err) ;


//for(int i=0;i<4;i++)cout<<" PD_t"<<i<< " "<<PD_t[i]<<endl;
//for(int i=0;i<5;i++)cout<<" PD_r "<<i<< " "<<PD_r[i]<<endl;
 //debug

    // errors 21 项 4 : [0~3] :aterr+  6: [4~9] : t_cov + 5 :[10~14]: arerr + 6:[15~20]: r_cov  
    for(int i=0;i<4;i++)errors[i] = PD_t[i]*PD_t[i]*at_err[i];   //4 : [0~3] :aterr  目标核参数误差
    errors[4] = 2*PD_t[0]*PD_t[1]*t_cov[0];                      //目标核参数协方差
    errors[5] = 2*PD_t[0]*PD_t[2]*t_cov[1];
    errors[6] = 2*PD_t[0]*PD_t[3]*t_cov[2];
    errors[7] = 2*PD_t[1]*PD_t[2]*t_cov[3];
    errors[8] = 2*PD_t[1]*PD_t[3]*t_cov[4];
    errors[9] = 2*PD_t[2]*PD_t[3]*t_cov[5];                     //6: [4~9] : t_cov

    for(int i=0;i<5;i++){ errors[10+i] = PD_r[i]*PD_r[i]*ar_err[i] ;}  //5 :[10~14]: arerr  参考核参数误差

    errors[15] = 2*PD_r[0]*PD_r[1]*r_cov[0];                     //参考核参数协方差
    errors[16] = 2*PD_r[0]*PD_r[2]*r_cov[1];
    errors[17] = 2*PD_r[0]*PD_r[3]*r_cov[2];
    errors[18] = 2*PD_r[1]*PD_r[2]*r_cov[3];
    errors[19] = 2*PD_r[1]*PD_r[3]*r_cov[4];
    errors[20] = 2*PD_r[2]*PD_r[3]*r_cov[5];                 //6:[15~20]: r_cov 
    
    errors_total=0;
    for(int i=0;i<21;i++){  errors_total+=errors[i];        }    //总方差
    
//for(int i=0;i<21;i++){ cout<<" errors ["<<i<<"] = "<< errors[i]<<endl;        } //debug
    
    one_mass_err = sqrt(errors_total);                           //开方，标准差
    return one_mass_err;

}

void Calculate_one_Mass_witherr_useDY(ION_UNKNOWN& ion_t, ION& ion_ref, TGraph* gr_gammat_C,TGraph* gr_gammat_C_u,TGraph* gr_gammat_C_d,
    double& onemass, double & onemasserr)
{
    //use DY: error = ErrorFrom_UD_CENTER[ y(x+sigma_x) , y(x-sigma_x) , y(x+0) ]
    //======================================================================
    
    double one_mass_err=0;

    int ion_t_Z = ion_t.Z;
    int ion_r_Z = ion_ref.Z;
    double tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref = 0;
    int error_n = 9;
    double errors[9]={0};   // 一次计算，9个参数，9项 误差来源 
    double error_gtC =0;
    double errors_total;
    double tA1err,tA2err,tA5err,tA6err,rA1err,rA2err,rA5err,rA6err,m_ref_err = 0;  // sqrt_variance
    double at_err[4]={0};double ar_err[5]={0};                                     // variance

    // target 4 paras
    tA1=ion_t.A1;tA2=ion_t.A2;tA5=ion_t.dA0;tA6=ion_t.dA1;
    // target 4 para err-variance
    at_err[0]=ion_t.A1err;at_err[1]=ion_t.A2err;at_err[2]=ion_t.dA0err;at_err[3]=ion_t.dA1err;
    // target 4 para err 
    tA1err=sqrt(ion_t.A1err);tA2err=sqrt(ion_t.A2err);tA5err=sqrt(ion_t.dA0err);tA6err=sqrt(ion_t.dA1err);
    
    // ref 5 paras
    rA1=ion_ref.A1; rA2=ion_ref.A2; rA5=ion_ref.dA0; rA6=ion_ref.dA1;
    m_ref = ionspecies[ion_ref.Species].Mass; 
    // ref 5 para err-variances
    ar_err[0]=ion_ref.A1err;ar_err[1]=ion_ref.A2err;ar_err[2]=ion_ref.dA0err;ar_err[3]=ion_ref.dA1err;
    ar_err[4]=ionspecies[ion_ref.Species].AME_err*ionspecies[ion_ref.Species].AME_err;
    rA1err=sqrt(ion_ref.A1err);rA2err=sqrt(ion_ref.A2err);rA5err=sqrt(ion_ref.dA0err);rA6err=sqrt(ion_ref.dA1err);
    m_ref_err = ionspecies[ion_ref.Species].AME_err;  // sigma not square

    //------------------------------ 上 中 下 三点 定误差
    double mu,md,mc=0;
    mc=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1,tA2,tA5,tA6,  rA1,rA2,rA5,rA6,m_ref,   gr_gammat_C);
    //cout<<" debug!  mc - AME:"<<mc-ion_t.M_AME<<endl;
    
    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1+tA1err, tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1-tA1err, tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    errors[0]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2+tA2err,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2-tA2err,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    errors[1]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5+tA5err,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5-tA5err,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    errors[2]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6+tA6err,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6-tA6err,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
    errors[3]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );
    //---- ref
    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1+rA1err,rA2,rA5,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1-rA1err,rA2,rA5,rA6,m_ref,gr_gammat_C);
    errors[4]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2+rA2err,rA5,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2-rA2err,rA5,rA6,m_ref,gr_gammat_C);
    errors[5]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2,rA5+rA5err,rA6,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2,rA5-rA5err,rA6,m_ref,gr_gammat_C);
    errors[6]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2,rA5,rA6+rA6err,m_ref,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2,rA5,rA6-rA6err,m_ref,gr_gammat_C);
    errors[7]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2,rA5,rA6,m_ref+m_ref_err,gr_gammat_C);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,  rA1,rA2,rA5,rA6,m_ref-m_ref_err,gr_gammat_C);
    errors[8]=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr, count_condition2_onemasserr );

    //----- gtC error
    mu=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,  gr_gammat_C_u);
    md=Calculate_one_Mass(ion_t_Z,ion_r_Z, tA1, tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,  gr_gammat_C_d);
    error_gtC=ErrorFrom_UD_CENTER_count(mu, md, mc, count_condition1_onemasserr_gtC, count_condition2_onemasserr_gtC );
    
    //----- sum 
    errors_total=0;
    for(int i=0;i<error_n;i++){  errors_total+=errors[i]*errors[i];        }
    //gtC error
    errors_total+=error_gtC*error_gtC;
    one_mass_err = sqrt(errors_total);  
    //________________________________________________________
    onemass= mc;
    onemasserr = one_mass_err;

    //cout<<" debug! in Calculate_one_Mass_witherr_useDY "<<mc-ion_t.M_AME<<" +- "<<one_mass_err<<endl;


}

//==================================LEVEL 3 final mass value of one iont (by using all the refs)=================================


double Calculate_iont_Mass_v1(ION_UNKNOWN& ion_t,ION* ions_ref, int ref_n,TGraph* gr_gammat_C,int k, double delta)
{

    int ion_t_Z=ion_t.Z;
    int ion_r_Z=0;
    double tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref = 0;


    double mass_result_each[ref_n_MAX];
    double mass_resulterr_each[ref_n_MAX];
    double mass_result_sum,mass_result_err_sum=0;
    int para_n = 4+5*ref_n;

    
    double* para = new double[para_n]; // max para number = 4+5*ref_n
    for(int i=0;i<para_n;i++)para[i]=0.0;
    // ======= set parameters ============= 
    para[0]=ion_t.A1;para[1]=ion_t.A2;para[2]=ion_t.dA0;para[3]=ion_t.dA1;
    for(int i=0;i<ref_n  ;i++)
    {
        int j=4+i*5;
        para[j+0]=ions_ref[i].A1; para[j+1]=ions_ref[i].A2; para[j+2]=ions_ref[i].dA0; para[j+3]=ions_ref[i].dA1; 
        para[j+4]= ionspecies[ions_ref[i].Species].Mass ;
    }
    
    // k: change one parameter k=1,2...(4+5*ref_n) , k=0:mass center value
    // delta = lr*err
    if(k==0){}
    else if(k>=1&&k<=para_n){ para[k-1]+=delta; }
    else {cout<<"error k out of range para_n!!!!"<<endl;return(-1);}

    //_______________________ change one para _________________________________

    tA1=para[0];tA2=para[1];tA5=para[2];tA6=para[3];
    mass_result_sum=0;mass_result_err_sum=0;
    for(int i=0;i<ref_n  ;i++)
    {
        int j=4+i*5;
        rA1=para[j+0];rA2=para[j+1];rA5=para[j+2];rA6=para[j+3];
        m_ref = para[j+4];
        ion_r_Z = ions_ref[i].Z; 
        mass_result_each[i]=Calculate_one_Mass(ion_t_Z,ion_r_Z,tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C);
        mass_resulterr_each[i] = Calculate_one_Mass_err_usePDCOV(ion_t,ions_ref[i], gr_gammat_C, lr,mass_result_each[i]);
        mass_result_sum += mass_result_each[i]/ (mass_resulterr_each[i]*mass_resulterr_each[i]);
        mass_result_err_sum += 1/ (mass_resulterr_each[i]*mass_resulterr_each[i]);

//cout<<" mass_result_each[i] "<<i<<" "<<mass_result_each[i]<<endl; //debug
//cout<<" mass_resulterr_each[i] " << i << " "<<mass_resulterr_each[i]<<endl; //debug
    }
    // weighted average
//cout<<"return mass_result_sum/mass_result_err_sum = "<<mass_result_sum/mass_result_err_sum<<endl; //debug

    return mass_result_sum/mass_result_err_sum;


}


//在有一组参考核的情况下计算目标核质量
double Calculate_iont_Mass_v2(ION_UNKNOWN& ion_t,ION* ions_ref, double* para,int ref_n,TGraph* gr_gammat_C,int k, double delta,bool use_weighted_ON)
{
    
    int ion_t_Z=ion_t.Z;
    int ion_r_Z=0;
    double tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref = 0;  //tA目标核系数，rA参考核系数。m_ref参考核质量
    double mass_result_each[ref_n_MAX];                //每个参考核给出的结果
    double mass_resulterr_each[ref_n_MAX];             //每个结果的误差
    double mass_result_sum,mass_result_err_sum=0;
    int para_n = 4+5*ref_n;                            //参数个数

    /*
    double* para = new double[para_n]; // max para number = 4+5*ref_n
    for(int i=0;i<para_n;i++)para[i]=0.0;
    // ======= set parameters ============= 
    para[0]=ion_t.A1;para[1]=ion_t.A2;para[2]=ion_t.dA0;para[3]=ion_t.dA1;
    for(int i=0;i<ref_n  ;i++)
    {
        int j=4+i*5;
        para[j+0]=ions_ref[i].A1; para[j+1]=ions_ref[i].A2; para[j+2]=ions_ref[i].dA0; para[j+3]=ions_ref[i].dA1; 
        para[j+4]= ionspecies[ions_ref[i].Species].Mass ;
    }
    */
    // k: change one parameter k=1,2...(4+5*ref_n) , k=0:mass center value !!!
    // delta = lr*err
    if(k==0){}//质量中心值，不做变动
    else if(k>=1&&k<=para_n){ para[k-1]+=delta; }//对此参数进行微小变动，以计算偏导数
    else {cout<<"error k out of range para_n!!!!"<<endl;return(-1);}

    //_______________________ change one para _________________________________

    tA1=para[0];tA2=para[1];tA5=para[2];tA6=para[3];
    mass_result_sum=0;mass_result_err_sum=0;
    for(int i=0;i<ref_n  ;i++)
    {
        int j=4+i*5;
        rA1=para[j+0];rA2=para[j+1];rA5=para[j+2];rA6=para[j+3];
        m_ref = para[j+4];
        ion_r_Z = ions_ref[i].Z; 
        mass_result_each[i]=Calculate_one_Mass_thisBp(ion_t_Z,ion_r_Z,tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref,gr_gammat_C,ion_t.Bp_v2);//利用此核得到的目标核质量，Bρ
//if(k==0)cout<<"ref "<<i<<" Bpt = "<<ion_t.Bp_v2<<" mass of ion_t "<<mass_result_each[i]-ion_t.M_AME <<endl; ///debug
        // weighted
        if(use_weighted_ON)//加权平均？ 认为加权平均是不合理的，各项的权重，也就是误差并不相互独立
        {
            mass_resulterr_each[i] = Calculate_one_Mass_err_usePDCOV(ion_t,ions_ref[i], gr_gammat_C, lr,mass_result_each[i]);
            mass_result_sum += mass_result_each[i]/ (mass_resulterr_each[i]*mass_resulterr_each[i]);                             //加权平均
            mass_result_err_sum += 1/ (mass_resulterr_each[i]*mass_resulterr_each[i]);
            if(MASS_VER>=4&&k==0){ionspecies[ion_t.Species].h_each_ref_cal_mass_err->Fill(mass_resulterr_each[i]); }
        }
        //unweighted
        else{
            mass_result_sum += mass_result_each[i];                                                                              //等权重平均
        }

//cout<<" mass_result_each[i] "<<i<<" "<<mass_result_each[i]<<endl; //debug
//cout<<" mass_resulterr_each[i] " << i << " "<<mass_resulterr_each[i]<<endl; //debug
    }
    
//cout<<"return mass_result_sum/mass_result_err_sum = "<<mass_result_sum/mass_result_err_sum<<endl; //debug
        

    if(use_weighted_ON) // directly set ions_unknown.M_cal_VE and err
    {
        ion_t.M_cal_err_VE = sqrt(1/mass_result_err_sum);  // without gtC ERR
        ion_t.M_cal_VE =   mass_result_sum/mass_result_err_sum;
    }
           
    if(use_weighted_ON)return mass_result_sum/mass_result_err_sum;// weighted average
    else return mass_result_sum/ref_n;     // unweighted average
//_______________________ return_________________________


}
//================== MASS_VER = 5 =================================
void Calculate_iont_Mass_useDY(ION_UNKNOWN& ion_t,ION* ions_ref, int ref_n,TGraph* gr_gammat_C,TGraph* gr_gammat_C_u,TGraph* gr_gammat_C_d,
    double& iont_Mass, double& iont_Mass_err)
{
    double mass_result_each[ref_n_MAX];
    double mass_resulterr_each[ref_n_MAX];
    double mass_result_sum = 0;
    double mass_result_err_sum=0;
    //cout<<" initial "<<mass_result_sum<<endl;
    double chi_n=0;
    for(int i=0;i<ref_n  ;i++)
    {

        Calculate_one_Mass_witherr_useDY(ion_t,ions_ref[i], gr_gammat_C,gr_gammat_C_u,gr_gammat_C_d,  mass_result_each[i],mass_resulterr_each[i]);
        mass_result_sum += mass_result_each[i]/ (mass_resulterr_each[i]*mass_resulterr_each[i]);
        mass_result_err_sum += 1/ (mass_resulterr_each[i]*mass_resulterr_each[i]);
        ionspecies[ion_t.Species].h_each_ref_cal_mass_err->Fill(mass_resulterr_each[i]);
        chi_n+= pow(mass_result_each[i]-ion_t.M_AME , 2) /  ( pow(mass_resulterr_each[i],2) + pow(ion_t.M_AME_err,2) );

    }

    iont_Mass =   mass_result_sum/mass_result_err_sum;
    //cout<<"debug! Calculate_iont_Mass_useDY :"<< iont_Mass-ion_t.M_AME<<endl;
    iont_Mass_err = sqrt(1/mass_result_err_sum);  

    chi_n = sqrt(chi_n/ref_n);
    ionspecies[ion_t.Species].h_iont_chi_n->Fill(chi_n);
}




////==================================LEVEL 4 mass value and error bar of one iont(by using all the refs )=================================


void Calculate_Mass_with_err(ION_UNKNOWN& ion_t,ION* ions_ref, int ref_n,TGraph* gr_gammat_C,TGraph* gr_gammat_C_u,TGraph* gr_gammat_C_d,
    double& mass_result,double& mass_result_err)
{
    // 进来一个目标核 iont 多个参考核的数组 ions_ref, 参考核数量 ref_n, 使用的gtC-curve , 记录最终结果的 mass, mass_err
    int ion_t_Z = ion_t.Z;
    double tA1,tA2,tA5,tA6,rA1,rA2,rA5,rA6,m_ref = 0;
    double mass_result_each[ref_n_MAX];
    int para_n = 4+ref_n*5;   //  the number of all the parameters that transfer errors in mass calculation //误差传递的参数个数
    double *paras = new double[para_n];                      //所有参数，包括目标核的A系数，参考核的A系数，参考核的质量
    for(int i=0;i<para_n;i++){paras[i]=0;}
    
    double *PD = new double[para_n];// partial derivatives    计算结果对各参数的偏微分
    for(int i=0;i<para_n;i++){PD[i]=0;}
    int error_n = 10+11*ref_n;  // error items: ion_t 4+6=10,  5+6=11 per ref ion           4项系数自身误差 + 6项协方差 + 1AME质量误差
    double *errors = new double[error_n];  for(int i=0;i<error_n;i++){errors[i]=0;}

    double errors_ref_ave[11]={0};    // average errors from ref ions
    double error_ref_ave = 0;
    double error_self=0;    

    double errors_total=0;
    double tA1err,tA2err,tA5err,tA6err,rA1err,rA2err,rA5err,rA6err,m_ref_err = 0;  // sqrt_variance  各项标准差 t目标核 r参考核
    double *var_err=new double[error_n];  for(int i=0;i<error_n;i++){var_err[i]=0;}  // variance       方差和协方差                          
    
    double *sigma_err=new double[error_n];for(int i=0;i<error_n;i++){sigma_err[i]=0;} // 各系数的标准差 sigma = sqrt(variance) 有空位！
    

    double ave_ref_dC=0;
    //----------------------------------------------------------------
    bool use_weighted_ON=USE_WEIGHTED;
    //----------------------------------------------------------------

    // target 4 paras
    paras[0]=ion_t.A1; paras[1]=ion_t.A2; paras[2]=ion_t.dA0; paras[3]=ion_t.dA1;
    // target 4 para err-variance
    var_err[0]=ion_t.A1err;var_err[1]=ion_t.A2err;var_err[2]=ion_t.dA0err;var_err[3]=ion_t.dA1err;
    // target 4 para err 
    sigma_err[0]=sqrt(ion_t.A1err);sigma_err[1]=sqrt(ion_t.A2err);sigma_err[2]=sqrt(ion_t.dA0err);sigma_err[3]=sqrt(ion_t.dA1err);
    // target 6 para covariance
    var_err[4] = ion_t.cov12;var_err[5] = ion_t.cov15;var_err[6] = ion_t.cov16;
    var_err[7] = ion_t.cov25;var_err[8] = ion_t.cov26;var_err[9] = ion_t.cov56;
    

    
    for(int i=0;i<ref_n;i++)//存入每个参考核 系数 以及 误差
    {
        // ref 5 paras
        paras[i*5+4]=ions_ref[i].A1; 
        paras[i*5+5]=ions_ref[i].A2; 
        paras[i*5+6]=ions_ref[i].dA0; 
        paras[i*5+7]=ions_ref[i].dA1;
        paras[i*5+8] = ionspecies[ions_ref[i].Species].Mass;
        
        // ref 5 para err-variances
        var_err[i*11+10]=ions_ref[i].A1err; 
        var_err[i*11+11]=ions_ref[i].A2err; 
        var_err[i*11+12]=ions_ref[i].dA0err; 
        var_err[i*11+13]=ions_ref[i].dA1err;
        var_err[i*11+14] = ionspecies[ions_ref[i].Species].AME_err * ionspecies[ions_ref[i].Species].AME_err;

        // ref 5 para err
        for(int j=0;j<5;j++) 
        { sigma_err[i*11+10+j]=sqrt(var_err[i*11+10+j]); }
        
        var_err[i*11+15]=ions_ref[i].cov12;
        var_err[i*11+16]=ions_ref[i].cov15;
        var_err[i*11+17]=ions_ref[i].cov16;
        var_err[i*11+18]=ions_ref[i].cov25; 
        var_err[i*11+19]=ions_ref[i].cov26;
        var_err[i*11+20]=ions_ref[i].cov56;
    }

   
//for(int i=0;i<error_n;i++){cout<<"var_err["<<i<<"]: "<<var_err[i]<<endl;}///debug
    
    //get center value
    double ion_t_mass_center = Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C,0, 0,use_weighted_ON);//  根据这一组参考核，得到了目标核质量

    if(INJ_sample_ON&&Do_PD_convergence_test)outfile_PD<<endl<<" --------- " <<ion_t.A<<" "<<ion_t.Z<<" "<<ion_t.inject_number<<endl;
    double para_delta = 0;
    for(int i=0;i<4;i++)
    {

        para_delta = lr*sigma_err[i];
        PD[i] =  (Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C, i+1, para_delta,use_weighted_ON) - ion_t_mass_center)/para_delta  ;//变动目标核的4个A系数，计算偏导数
        ///20240703 test============================
        if(INJ_sample_ON&&Do_PD_convergence_test)
        {
            outfile_PD<<"para_delta= "<<para_delta<<" PD["<<i<<"]= "<<PD[i]<<endl;
            for(int jj=1;jj<=20  ;jj++)
            {
                double PD_tmp=0; 
                para_delta = (lr+0.05*jj)*sigma_err[i];
                double PD_dy = (Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C, i+1,para_delta,use_weighted_ON) - ion_t_mass_center);//变动参考核的4个A系数和质量，计算偏导数
                PD_tmp =  PD_dy/para_delta  ;
                outfile_PD<<" "<<para_delta<<" "<<PD_dy<<" "<<PD_tmp<<endl;
            }
        }
        ///20240703 test_____________________

    }
    for(int i=0;i<ref_n;i++) 
    {
        for(int j=0;j<5;j++)
        {
            para_delta = lr*sigma_err[i*11+10+j];
            PD[4+i*5+j]=(Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C, i*5+4+j+1,para_delta,use_weighted_ON) - ion_t_mass_center)/para_delta;

            ///20240703 test ====================
            if(INJ_sample_ON&&Do_PD_convergence_test)
            {
                outfile_PD<<"para_delta= "<<para_delta<<" PD["<<4+i*5+j<<"]= "<<PD[4+i*5+j]<<endl;
                for(int jj=1;jj<=20  ;jj++)
                {
                    double PD_tmp=0;  
                    para_delta = (lr+0.05*jj)*sigma_err[i*11+10+j];
                    double PD_dy = (Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C, i*5+4+j+1, para_delta,use_weighted_ON) - ion_t_mass_center);
                    PD_tmp =  PD_dy/para_delta  ;
                    outfile_PD<<" "<<para_delta<<" "<<PD_dy<<" "<<PD_tmp<<endl;
                }
            }
            ///20240703 test_____________________
        }        
    }
//for(int i=0;i<para_n;i++){cout<<"PD["<<i<<"]: "<<PD[i]<<endl;}///debug


    for(int i=0;i<4  ;i++)
    {
        errors[i] = PD[i]*PD[i]*var_err[i];//目标核各项导致的误差
    }
    
    errors[4] = 2*PD[0]*PD[1]*var_err[4];//目标核各项协方差
    errors[5] = 2*PD[0]*PD[2]*var_err[5];
    errors[6] = 2*PD[0]*PD[3]*var_err[6];
    errors[7] = 2*PD[1]*PD[2]*var_err[7];
    errors[8] = 2*PD[1]*PD[3]*var_err[8];
    errors[9] = 2*PD[2]*PD[3]*var_err[9];
    // ref -- 11 errors per ref 
    for(int i=0;i<ref_n  ;i++)
    {
        for(int j=0;j<5  ;j++)  
        {errors[10+11*i+j] = PD[4+i*5+j]*PD[4+i*5+j]*var_err[10+11*i+j];}   // 5 errors from ii per ref 参考核各项导致的误差
    }
    for(int i=0;i<ref_n  ;i++)
    {
        errors[i*11+15+0] = 2*PD[4+i*5+0]*PD[4+i*5+1]*var_err[i*11+15+0];   //参考核各项协方差
        errors[i*11+15+1] = 2*PD[4+i*5+0]*PD[4+i*5+2]*var_err[i*11+15+1];
        errors[i*11+15+2] = 2*PD[4+i*5+0]*PD[4+i*5+3]*var_err[i*11+15+2];
        errors[i*11+15+3] = 2*PD[4+i*5+1]*PD[4+i*5+2]*var_err[i*11+15+3];
        errors[i*11+15+4] = 2*PD[4+i*5+1]*PD[4+i*5+3]*var_err[i*11+15+4];
        errors[i*11+15+5] = 2*PD[4+i*5+2]*PD[4+i*5+3]*var_err[i*11+15+5];   // 6 errors from ij per ref
    }
    //___________ errors completed _______________________

    // erros total excluding err_gtC
    errors_total=0;
    for(int i=0;i<error_n;i++){errors_total+=errors[i];}                    //各参数  导致的总误差
//for(int i=0;i<error_n;i++){cout<<"errors["<<i<<"]: "<<errors[i]<<endl;} ///debug

    error_self=0; 
    for(int i=0;i<10;i++){error_self+=errors[i];}// the first 10 items are the errors from target ion itself
    error_self=sqrt(error_self);
    for(int i=0;i<ref_n  ;i++)
    {
        for(int j=0;j<11  ;j++)
        {
            errors_ref_ave[j] +=errors[10+11*i+j];
        }
    }
    for(int j=0;j<11  ;j++){errors_ref_ave[j]/=ref_n;}
    error_ref_ave=0;    
    for(int j=0;j<11  ;j++){error_ref_ave+=errors_ref_ave[j];}
    error_ref_ave/=11;
    error_ref_ave=sqrt(error_ref_ave);


    //add gtC error  20230706
    double error_gtC = 0; //error_gtC 由于γt(C)曲线的误差所导致的 质量结果误差
    double m_iont_gtC_u,m_iont_gtC_d=0;
    m_iont_gtC_u = Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C_u,0, 0,use_weighted_ON);//  根据这一组参考核，和γt(C)曲线上限，得到目标核质量
    m_iont_gtC_d = Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C_d,0, 0,use_weighted_ON);//  根据这一组参考核，和γt(C)曲线下限，得到目标核质量
    //error_gtC= 0.5*abs(m_iont_gtC_u - m_iont_gtC_d);
    if( (m_iont_gtC_u- ion_t_mass_center)/( ion_t_mass_center-m_iont_gtC_d) >0 )
    {   
        count_gtC_masserr_condition_1++;
        error_gtC= 0.5*abs(m_iont_gtC_u - m_iont_gtC_d);
    }
    else
    {   
        count_gtC_masserr_condition_2++;
        error_gtC= 0.5* ( abs(m_iont_gtC_u- ion_t_mass_center)+abs( ion_t_mass_center-m_iont_gtC_d) );
    }

    ///######################### get final one mass and err #######################################
    mass_result = ion_t_mass_center;
    mass_result_err = sqrt(errors_total);    //由于参数误差导致的质量误差

    if(OUTFILE_v2err_gtC_ON)outfile_v2error_gtC<<error_gtC<<" / "<<mass_result_err<<" = "<<error_gtC/mass_result_err
    <<" | "<<(mass_result_err+error_gtC);//debug 20230706
    
    //------ add err_gtC to final error
    mass_result_err =  sqrt( mass_result_err*mass_result_err+error_gtC*error_gtC);                    //参数误差+γt(C)曲线误差
    ///############################################################################################
    if(OUTFILE_v2err_gtC_ON)outfile_v2error_gtC<<" | "<<mass_result_err<<endl;

    if(outfile_ERR_ANA_ON)
    {
    //----------------- OUTPUT FILE ERR_ANA
    outfile_ERR_ANA<<ion_t.A<<ion_t.name<<" "<<ion_t.inject_number<<" "<<ion_t.ion_number<<" "<<ref_n;
    //for(int i=0;i<error_n;i++){outfile_ERR_ANA<<errors[i]<<" ";}
    //--fix ouput items 20231020
    outfile_ERR_ANA<<" | "<<mass_result_err<<" | "<<error_self<<" "<<error_ref_ave<<" "<<error_gtC<<" | ";
    
    // output 10 self+ 11(ave ref) note that covariance can be negative
    for(int j=0;j<10  ;j++){outfile_ERR_ANA<<errors[j]<<" ";}
     outfile_ERR_ANA<<" || ";
    for(int j=0;j<11  ;j++){outfile_ERR_ANA<<errors_ref_ave[j]<<" ";}
    outfile_ERR_ANA<<endl;

    }
//cout<<" Calculate_Mass_with_err: mass = "<<mass_result<<" +- "<<mass_result_err<<" "<<ion_t.A<<ion_t.name<<endl; ///debug

}

////================================== 4_v4 =================================
void Calculate_Mass_with_err_v4(ION_UNKNOWN& ion_t,ION* ions_ref, int ref_n,TGraph* gr_gammat_C,TGraph* gr_gammat_C_u,TGraph* gr_gammat_C_d,
    double& mass_result,double& mass_result_err)
{
    int para_n = 4+ref_n*5;
    double *paras = new double[para_n];   
    for(int i=0;i<para_n;i++){paras[i]=0;}
    // target 4 paras
    paras[0]=ion_t.A1; paras[1]=ion_t.A2; paras[2]=ion_t.dA0; paras[3]=ion_t.dA1;
    for(int i=0;i<ref_n;i++)
    {
        // ref 5 paras
        paras[i*5+4]=ions_ref[i].A1; 
        paras[i*5+5]=ions_ref[i].A2; 
        paras[i*5+6]=ions_ref[i].dA0; 
        paras[i*5+7]=ions_ref[i].dA1;
        paras[i*5+8] = ionspecies[ions_ref[i].Species].Mass;
    }
    double ion_t_mass_center = Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C,0, 0, true); // use weighted to set ion_t.M_cal_err_VE
    //ion_t.M_cal_err_VE already set in Calculate_iont_Mass_v2
    double err_without_gtC = ion_t.M_cal_err_VE;

    // add err_gtC
    double error_gtC = 0;
    double m_iont_gtC_u,m_iont_gtC_d=0;
    //!!!!! this will change ion_t.M_cal_err_VE and M_cal_VE!!!
    m_iont_gtC_u = Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C_u,0, 0,true);
    m_iont_gtC_d = Calculate_iont_Mass_v2(ion_t,ions_ref, paras,ref_n,gr_gammat_C_d,0, 0,true);
    //error_gtC= 0.5*abs(m_iont_gtC_u - m_iont_gtC_d);
    if( (m_iont_gtC_u- ion_t_mass_center)/( ion_t_mass_center-m_iont_gtC_d) >0 )
    {   
        count_gtC_masserr_condition_1_VE++;
        error_gtC= 0.5*abs(m_iont_gtC_u - m_iont_gtC_d);
    }
    else
    {   
        count_gtC_masserr_condition_2_VE++;
        error_gtC= 0.5* ( abs(m_iont_gtC_u- ion_t_mass_center)+abs( ion_t_mass_center-m_iont_gtC_d) );
    }
    ion_t.M_cal_err_VE = sqrt(pow(err_without_gtC,2)+pow(error_gtC,2) );
    ion_t.M_cal_VE = ion_t_mass_center;
    mass_result = ion_t_mass_center;
    mass_result_err = ion_t.M_cal_err_VE;
    //============== outfile
    if(outfile_ERR_ANA_ON)
    {
        outfile_ERR_ANA<<ion_t.A<<ion_t.name<<" "<<ion_t.inject_number<<" "<<ion_t.ion_number<<" "<<ref_n
        <<" | "<<mass_result_err<<" | "<<err_without_gtC<<" "<<error_gtC<<" | "<<endl;   
    }

    if(MASS_VER>=4)ionspecies[ion_t.Species].h_each_ref_cal_mass_err->Fill(mass_result_err);
}



//____________________________________  Calculate Mass  _____________________________________

//////////////========================= smooth algorithm =========================

//取 ±set_smooth_k 范围内的点，平滑γt(C)曲线
void Smooth_DW(TGraph* gr_in, TGraph* gr_out, int kn, int opt)
{//opt=0,不平滑；1 附近点取加权平均；2 附近点的(1+γt^2*ΔC/C)连乘取开方

    //TString test_name = "smooth2.txt";
    //ofstream outfile_smooth_test;
    //outfile_smooth_test.open(test_name);

    int n = gr_in->GetN();
    double B[1000] = { 0 };
    double x[1000] = { 0 };
    double y[1000] = { 0 };
    double dx;
    int j_min, j_max;
    int k;

    if (n > 1000) { cout << "!!! n >1000 !!!" << endl; return; }
    //Smooth_distance_weighted  --  2*kn+1 points weighted average
    double tmpsum = 0;
    cout << "smooth_opt = " << opt << "  " << endl;
    for (int i = 0; i < n; i++)
    {
        gr_in->GetPoint(i, x[i], y[i]);//先存入数组中
        //cout<<i<<" "<<x[i]<<" "<<y[i]<<endl;
    }
    if (opt == 0)
    {
        k = 0;
    }
    else if (opt == 1 || opt == 2)
    {
        k = kn;
    }
    else
    {
        cout << "!!! the smooth_opt is wrong. !!!" << endl;
        return;
    }

    for (int i = 0; i < n; i++)
    {
        if (opt == 2)
        {
            B[i] = 1;//连乘，因此初始化为1
        }
        else
        {
            B[i] = 0;
        }
        tmpsum = 0;
        j_min = i - k;
        j_max = i + k;
        if (j_min < 0)j_min = 0;
        if (j_max > n - 1)j_max = n - 1;
        if (opt == 0 || opt == 1)
        {
            for (int j = j_min; j <= j_max; j++)
            {
                B[i] += pow(0.5, abs(x[i] - x[j])) * y[j];//0.5^(dx)*y
                tmpsum += pow(0.5, abs(x[i] - x[j]));//0.5^(dx) 各项权重
            }
            B[i] /= tmpsum;
        }
        else if (opt == 2)
        {
            if (i == 0)//避免数组 x[i-1]越界
            {
                dx = x[i + 1] - x[i];
            }
            else
            {
                dx = x[i] - x[i - 1];
            }
            for (int j = j_min; j <= j_max; j++)
            {
                B[i] = B[i] * pow((1 + pow(y[j], 2) * dx / x[j]), pow(0.5, abs(i - j)));//      (1+γt^2*ΔC/C)^(0.5^(dx))
                tmpsum += pow(0.5, abs(i - j));                                                      //      0.5^(dx) 各项权重
                //cout << "i=" << i << "  j=" << j << "  dx=" << dx << "  B[i]=" << B[i] << "  tmpsum=" << tmpsum << endl;
            }
            B[i] = pow((pow(B[i], 1 / tmpsum) - 1) * x[i] / dx, 0.5);//回到γt
            //cout << "i=" << i <<"  B[i]=" << B[i] << "  tmpsum=" << tmpsum << endl;
        }
        //cout<<B[i]<<endl;
        gr_out->SetPoint(i, x[i], B[i]);
        //cout << i << "  " << x[i] << "  " << B[i] << endl;
    }
    //outfile_smooth_test.close();
    return;
}
//平滑误差
void Smooth_DW_err(TGraphErrors* grerr_in, TGraphErrors* grerr_out, int kn, int opt)
{
    //TString test_name = "smooth_err1.txt";
    //ofstream outfile_smooth_test;
    //outfile_smooth_test.open(test_name);


    int n = grerr_in->GetN();
    double B[1000] = { 0 };
    double Berr[1000] = { 0 };
    double x[1000] = { 0 };
    double y[1000] = { 0 };
    double yerr[1000] = { 0 };
    double dx;
    int j_min, j_max;
    int k;
    if (n > 1000) { cout << "!!! n >1000 !!!" << endl; return; }
    if (opt == 0)
    {
        k = 0;
    }
    else if (opt == 1 || opt == 2)
    {
        k = kn;
    }
    else
    {
        cout << "!!! the smooth_opt is wrong. !!!" << endl;
        return;
    }
    //Smooth_distance_weighted  --  2*kn+1 points weighted average
    double tmpsum = 0;
    for (int i = 0; i < n; i++)
    {
        grerr_in->GetPoint(i, x[i], y[i]);
        yerr[i] = grerr_in->GetErrorY(i);
        //cout<<i<<" "<<x[i]<<" "<<y[i]<<endl;
    }
    for (int i = 0; i < n; i++)
    {
        if (opt == 2)
        {
            B[i] = 1;//连乘，因此初始化为1
        }
        else
        {
            B[i] = 0;
        }
        Berr[i] = 0;
        tmpsum = 0;
        j_min = i - k;
        j_max = i + k;
        if (j_min < 0)j_min = 0;
        if (j_max > n - 1)j_max = n - 1;
        if (opt == 0 || opt == 1)
        {
            for (int j = j_min; j <= j_max; j++)
            {
                B[i] += pow(0.5, abs(x[i] - x[j])) * y[j];
                Berr[i] += pow(0.5, 2 * abs(x[i] - x[j])) * pow(yerr[j], 2);
                tmpsum += pow(0.5, abs(x[i] - x[j]));
            }
            B[i] /= tmpsum;
            Berr[i] = sqrt(Berr[i]) / tmpsum;
        }
        else if (opt == 2)
        {
            if (i == 0)//避免数组 x[i-1]越界
            {
                dx = x[i + 1] - x[i];
            }
            else
            {
                dx = x[i] - x[i - 1];
            }
            for (int j = j_min; j <= j_max; j++)
            {
                B[i] = B[i] * pow((1 + pow(y[j], 2) * dx / x[j]), pow(0.5, abs(i - j)));//      (1+γt^2*ΔC/C)^(0.5^(dx))
                Berr[i] += pow(pow(0.5, abs(i - j)) * y[j] / x[j] / (1 + pow(y[j], 2) * dx / x[j]), 2) * pow(yerr[j], 2);//由误差传递公式得
                tmpsum += pow(0.5, abs(i - j));                                         //      0.5^(dx) 各项权重
            }
            B[i] = pow((pow(B[i], 1 / tmpsum) - 1) * x[i] / dx, 0.5);                         //回到γt
            Berr[i] = Berr[i] / pow(tmpsum * y[i] / x[i] / (1 + pow(y[i], 2) * dx / x[i]), 2);//此处为方差
            Berr[i] = sqrt(Berr[i]);                                                          //开方，得到标准差
        }
        //cout<<B[i]<<endl;
        grerr_out->SetPoint(i, x[i], B[i]);
        grerr_out->SetPointError(i, 0, Berr[i]);
        //outfile_smooth_test << i << "  " << x[i] << "  " << B[i] << "  " << Berr[i] << endl;
    }
    return;
}

//void Smooth_DW(TGraph*gr_in,TGraph*gr_out,int kn)
//{
//    int n = gr_in->GetN();
//    double B[1000]={0};
//    double x[1000]={0};
//    double y[1000]={0};
//    int j_min,j_max;
//    if(n>1000) {cout<<"!!! n >1000 !!!"<<endl;return;}
//    //Smooth_distance_weighted  --  2*kn+1 points weighted average
//    double tmpsum=0;
//    for(int i=0;i<n  ;i++)
//    {
//        gr_in->GetPoint(i,x[i],y[i]);
//        //cout<<i<<" "<<x[i]<<" "<<y[i]<<endl;
//    }
//    for(int i=0;i<n  ;i++)
//    {
//        B[i]=0;
//        tmpsum=0;
//        j_min = i-kn;
//        j_max = i+kn;
//        if(j_min<0)j_min=0;
//        if(j_max>n-1)j_max=n-1;
//        for(int j = j_min;j<=j_max;j++)
//        {
//            
//            B[i]+=pow(0.5,abs(x[i]-x[j]))*y[j];
//            tmpsum+=pow(0.5,abs(x[i]-x[j]));
//        }
//        B[i] /= tmpsum;
//        //cout<<B[i]<<endl;
//        gr_out->SetPoint(i,x[i],B[i]);
//    }
//    return;
//}
//void Smooth_DW_err(TGraphErrors*grerr_in,TGraphErrors*grerr_out,int kn)
//{
//    int n = grerr_in->GetN();
//    double B[1000]={0};
//    double Berr[1000]={0};
//    double x[1000]={0};
//    double y[1000]={0};
//    double yerr[1000]={0};
//    int j_min,j_max;
//    if(n>1000) {cout<<"!!! n >1000 !!!"<<endl;return;}
//    //Smooth_distance_weighted  --  2*kn+1 points weighted average
//    double tmpsum=0;
//    for(int i=0;i<n  ;i++)
//    {
//        grerr_in->GetPoint(i,x[i],y[i]);
//        yerr[i] =  grerr_in->GetErrorY(i);
//
//        //cout<<i<<" "<<x[i]<<" "<<y[i]<<endl;
//    }
//    for(int i=0;i<n  ;i++)
//    {
//        B[i]=0;
//        Berr[i]=0;
//        tmpsum=0;
//        j_min = i-kn;
//        j_max = i+kn;
//        if(j_min<0)j_min=0;
//        if(j_max>n-1)j_max=n-1;
//        for(int j = j_min;j<=j_max;j++)
//        {
//            B[i]+=pow(0.5,abs(x[i]-x[j]))*y[j];
//            Berr[i]+= pow(0.5,2*abs(x[i]-x[j]))*pow(yerr[j],2);
//            tmpsum+=pow(0.5,abs(x[i]-x[j]));
//        }
//        B[i] /= tmpsum;
//        Berr[i] = sqrt(Berr[i]) / tmpsum;
//        //cout<<B[i]<<endl;
//        grerr_out->SetPoint(i,x[i],B[i]);
//        grerr_out->SetPointError(i,0,Berr[i]);
//    }
//    return;
//}

//================ 20230616

// 筛掉问题太大，超过n_sigma的异常点，重新计算 gtC 数组，
void squeeze_gtC(double* gtC,double* gtC_sigma, int* C_n, int opt_n,TH1F* h1, TGraph* gr)
{
    // tmp 临时数组 先进行初始化    opt_n 选取几倍σ误差范围内
    double gtC_i [1000];
    double gtC_sigma_i [1000];
    int  C_n_i [1000];
    if(subregion_n>1000){cout<<" subregion_n>1000 in squeeze_gtC !! array overflow!! return "; return;}    
    for(int i=0;i<subregion_n;i++){gtC_i[i]=0;gtC_sigma_i[i]=0;C_n_i[i]=0;}
    for(int i=0;i<subregion_n;i++){gtC_sigma[i]*=sqrt(C_n[i]-1);}
    // std dev
    int C_region_i=0;//bin索引，对应数组
    int gr_n=0;
    //根据传进来的gtC 筛选范围内的散点
    for(int i=0;i<ions_n;i++)
    {
        C_region_i = h1->FindBin(ions[i].C)-1; //0:超出下界 1~subregion_n ,
        if( (C_region_i<0) || C_region_i >(subregion_n-1) ){ continue;}


        //if( (abs(ions[i].gammat-gtC[C_region_i]) < opt_n*gtC_sigma[C_region_i]) &&ions[i].Z>=choose_largeZ_min &&ions[i].Z<=choose_largeZ_max)
        //if((abs(ions[i].gammat-gtC[C_region_i]) < opt_n*gtC_sigma[C_region_i]) &&(!(ions[i].Z<=CONDITION_gt_Z&&ions[i].A<=CONDITION_gt_A) )&&(ionspecies[ions[i].Species].MassUnknown==0) )
        if(( abs(ions[i].gammat-gtC[C_region_i]) < opt_n*gtC_sigma[C_region_i]) &&ionspecies[ions[i].Species].gtC_CHOSEN )
        {//筛选规则：首先，满足之前选取核的要求（重荷且质量已知）；其次，此核的γt在 opt_n倍σ误差范围内（偏差太大的就扔掉了）
            gtC_i[C_region_i]+=ions[i].gammat;
            gtC_sigma_i[C_region_i]+=ions[i].gammat*ions[i].gammat;
            C_n_i[C_region_i]++;
            gr->SetPoint(gr_n++,ions[i].C,ions[i].gammat);
            //gr zone 记录筛选后的散点
        }
    }
    //用新的范围内散点生成gtC 数组
    for(int i=0;i<subregion_n  ;i++)
    {
        //if(C_n_i[i]>C_DIVISION_CHOSEN_MIN)
        {
            gtC_i[i]/=C_n_i[i];
            gtC_sigma_i[i]/=C_n_i[i];
            gtC_sigma_i[i]=sqrt(gtC_sigma_i[i]-gtC_i[i]*gtC_i[i]);
            gtC_sigma_i[i] /= sqrt(C_n_i[i] -1);
            //刷新 gtC 数组 refresh 筛掉异常点，重新计算的结果
            gtC[i]=gtC_i[i];
            gtC_sigma[i]=gtC_sigma_i[i];
            C_n[i] = C_n_i[i];
        }   
    }   
    
}
//___________________ 20230616
//20221030 ===============integration of gtC largeZline ================
void Do_gtC_largeZline(ION* ions, int ions_n, TH1F* h1)
{
    Double_t gtCZ[Z_division][subregion_n];
    Double_t gtC2Z[Z_division][subregion_n];  //C square divided into sub regions
    Double_t Sigma_gtCZ[Z_division][subregion_n];
    Int_t C_Division_nZ[Z_division][subregion_n];       // count for numbers in each subregion
    int Z_large_equal[Z_division]={3,4,5,6,8,9};   // >= //  Z_division =6
    int Selected_zn[Z_division];
    for(int i=0;i<Z_division  ;i++){Selected_zn[i]=0;}
    for(int j=0;j<Z_division;j++)
    {   for(int i=0;i<subregion_n  ;i++){    gtCZ[j][i]=0.0;gtC2Z[j][i]=0.0;Sigma_gtCZ[j][i]=0.0;C_Division_nZ[j][i]=0;}  }
    

    TGraph* gr_gtC_largeZline[6][3];
    TGraph* gr_gtC_largeZ[6];
    TGraphErrors* grerr_avegtC_largeZ[6];
    for(int i=0;i<Z_division;i++)
    {
        gr_gtC_largeZ[i]=new TGraph();
        grerr_avegtC_largeZ[i]=new TGraphErrors();
        for(int j=0;j<3;j++)
            {gr_gtC_largeZline[i][j]=new TGraph();}
    } 


    for(int i=0;i<ions_n  ;i++)
    {
        for(int j=0;j<Z_division  ;j++)
        {
            if(ions[i].Z>=Z_large_equal[j] && (ions[i].gammat >1.3&&ions[i].gammat <1.4))
            {
                if( ions[i].C_region<0|| ions[i].C_region>subregion_n-1){continue;}
                gr_gtC_largeZ[j]->SetPoint(Selected_zn[j]++,ions[i].C,ions[i].gammat);
                gtCZ[j][ions[i].C_region]+=ions[i].gammat;
                gtC2Z[j][ ions[i].C_region ] += ions[i].gammat*ions[i].gammat;        
                C_Division_nZ[j][ ions[i].C_region]++;
            }
        }
    }

    //================   gr_gtC_largeZline  =====================
    int gr_n=0;
    for(int i=0;i<Z_division  ;i++)
    {
        gr_n=0;
        for(int j=0;j<subregion_n  ;j++)
        {
            if(C_Division_nZ[i][j]>C_DIVISION_CHOSEN_MIN)
            {
            gtCZ[i][j] /= C_Division_nZ[i][j];
            gtC2Z[i][j] /= C_Division_nZ[i][j];
            Sigma_gtCZ[i][j] = sqrt (gtC2Z[i][j]-gtCZ[i][j]*gtCZ[i][j]) ;
            //cout<<i<<"  "<<j<<"  "<<gtCZ[i][j]<<endl;
    
            gr_gtC_largeZline[i][0]->SetPoint(gr_n,h1->GetBinCenter(j+1), gtCZ[i][j] );
            gr_gtC_largeZline[i][1]->SetPoint(gr_n,h1->GetBinCenter(j+1), gtCZ[i][j]+n_sigma* Sigma_gtCZ[i][j] );
            gr_gtC_largeZline[i][2]->SetPoint(gr_n,h1->GetBinCenter(j+1), gtCZ[i][j]-n_sigma* Sigma_gtCZ[i][j] );
            grerr_avegtC_largeZ[i]->SetPoint(gr_n,h1->GetBinCenter(j+1),gtCZ[i][j]);
            grerr_avegtC_largeZ[i]->SetPointError(gr_n,0,Sigma_gtCZ[i][j]/sqrt(C_Division_nZ[i][j]) );
            gr_n++;
            }
                    
        }
        AxisFormat(gr_gtC_largeZline[i][0],strtmp.Format("Z>=%d",Z_large_equal[i]),"","",2+i);gr_gtC_largeZline[i][0]->SetMarkerSize(1);
        AxisFormat(gr_gtC_largeZline[i][1],strtmp.Format("Z>=%d +2sigma",Z_large_equal[i]),"","",2+i);gr_gtC_largeZline[i][1]->SetMarkerSize(1);
        AxisFormat(gr_gtC_largeZline[i][2],strtmp.Format("Z>=%d-2sigma",Z_large_equal[i]),"","",2+i);gr_gtC_largeZline[i][2]->SetMarkerSize(1);
        AxisFormat(grerr_avegtC_largeZ[i],strtmp.Format("Z>=%d",Z_large_equal[i]),"","",2+i);
        grerr_avegtC_largeZ[i]->SetLineWidth(3);
    } 
    
    bool c_gtC_largeZ_ON =1;
    if(c_gtC_largeZ_ON)
    {//================================c_gtC_largeZ_ON===============================================
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c_gtC_largeZ = new TCanvas("c_gtC_largeZ","c_gtC_largeZ  all ",1000,500);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    c_gtC_largeZ->cd(1);
    auto mg_gtC_largeZ=new TMultiGraph();
    for(int i=0;i<6  ;i++)
    {
        AxisFormat(gr_gtC_largeZ[i],strtmp.Format("#gamma_{t} ions selected Z>=%d",Z_large_equal[i]),"C [m]","#gamma_{t}",1);
        gr_gtC_largeZ[i]->SetMarkerSize(1);
    }
    gr_gtC_largeZ[0]->SetMarkerColor(kRed);
    gr_gtC_largeZ[1]->SetMarkerColor(kAzure);
    gr_gtC_largeZ[2]->SetMarkerColor(kOrange);
    gr_gtC_largeZ[3]->SetMarkerColor(kCyan);
    gr_gtC_largeZ[4]->SetMarkerColor(kViolet);
    gr_gtC_largeZ[0]->SetLineColor(kRed);
    gr_gtC_largeZ[1]->SetLineColor(kAzure);
    gr_gtC_largeZ[2]->SetLineColor(kOrange);
    gr_gtC_largeZ[3]->SetLineColor(kCyan);
    gr_gtC_largeZ[4]->SetLineColor(kViolet);
    gr_gtC_largeZ[5]->SetMarkerColor(1);
    for(int i=0;i<6  ;i++)
    {
        mg_gtC_largeZ->Add(gr_gtC_largeZ[i]);
    }
    AxisFormat(mg_gtC_largeZ,"#gamma_{t} for large Z","C [m]","Ave #gamma_{t}");
    
    
    mg_gtC_largeZ->DrawClone("ap");
    c_gtC_largeZ->BuildLegend();
    lat_n->SetTextFont(43);
    lat_n->SetTextSize(30);
    for(int i=0;i<6  ;i++){ lat_n->DrawLatex(128.9,1.35-i*0.005,strtmp.Format("Z>=%d,ions selected: %d",Z_large_equal[i],Selected_zn[i]));}
    
    lat_n->SetTextSize(40);
    }//_____________________________________c_gtC_largeZ_ON________________________________________________
}
//_____________________________________________________       



//2022================  Generate_BpC_dispersion ================================
void Generate_BpC_dispersion(TString c_name, IONSpecies* ionspecies,int NSpecies,int opt,bool c_BpC_dispersion_ON,
    bool c_BpC_dd_ON,double& C_intersection )
{
    cout<<"--------------Generate_BpC_dispersion ON---------------------------"<<endl;
    int have_fit_n =0;
    double tmp[100]={0};
    double C_start =128 , C_end = 129.5;
    double lnC_start = log(C_start);    double lnC_end = log(C_end);
    int steps = 1000;
    double C_one_step = (C_end -C_start)/steps; double lnC_one_step = (lnC_end -lnC_start)/steps;
    
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ionspecies[i].N>10){have_fit_n++;}
    }
    cout<<endl<<" if(ionspecies[i].N>10) have fit n ="<<have_fit_n<<endl;

    bool chosen[100]={0};
    int chosen_n=0;
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ionspecies[i].Aname=="10C"||ionspecies[i].Aname=="8B"||ionspecies[i].Aname=="13O"||ionspecies[i].Aname=="11C")
            {chosen[i]=1;}
        chosen_n+=chosen[i];
    }
    cout<<" chosen_n =  "<<chosen_n<<endl;
    

    TGraph* gr_BpC_dispersion = new TGraph();   
    TGraph* gr_lnBpC_dispersion = new TGraph();
    AxisFormat(gr_BpC_dispersion,"gr_BpC_dispersion","C","variance of Bpc");
    AxisFormat(gr_lnBpC_dispersion,"gr_lnBpC_dispersion","lnC","variance of lnBpc");
    gr_BpC_dispersion->SetMarkerSize(1.5);
    gr_lnBpC_dispersion->SetMarkerSize(1.5);
    
    TGraph* gr_BpC_dd[51];
    //int gr_BpC_dd_n[100]
    for(int i=0;i<chosen_n  ;i++)
    {   
        gr_BpC_dd[i] = new TGraph();      
        AxisFormat(gr_BpC_dd[i],"","C","#Delta Bp",1+i%9);
        gr_BpC_dd[i]->SetMarkerSize(1.5);
        gr_BpC_dd[i]->SetMarkerStyle(20+i/10);
    }
    
    //TGraph* gr_BpC_ddd[50];
    //TGraph* gr_BpC_dd = new TGraph();
    //AxisFormat(gr_BpC_dd,"","C","#Delta Bp");
      //  gr_BpC_dd->SetMarkerSize(1.5);
    //int gr_BpC_dispersion_n=0;  
    for(int j=0;j<steps; j++ ) //128--129
    {
        int tmp_n=0;
        for(int i=0;i<NSpecies  ;i++)
        {
            if(chosen[i])
            {
                gr_BpC_dd[tmp_n]->SetTitle(ionspecies[i].Aname);
                if(opt==0){ tmp[tmp_n]=ionspecies[i].fitfun_pol1_BpC->Eval(C_start+j*C_one_step);}
                else if(opt==1){tmp[tmp_n]=ionspecies[i].fitfun_pol1_lnBpC->Eval(lnC_start+j*lnC_one_step);}
                tmp_n++;
            }
        }
        if(tmp_n!=chosen_n)cout<<"error! tmp_n!=chosen_n"<<endl;
        //if(j==584){for(int i=0;i<have_fit_n  ;i++){cout<<fixed<<setprecision(9)<<tmp[i]<<endl;}}

        for(int i=0;i<chosen_n  ;i++)
        {
            gr_BpC_dd[i]->SetPoint(j,C_start+j*C_one_step,tmp[i]-tmp[0]);
        }
        //gr_BpC_dd->SetPoint(j,C_start+j*C_one_step,tmp[have_fit_n/2]-tmp[0]);
        

        if(opt==0) {gr_BpC_dispersion->SetPoint(j , C_start+j*C_one_step, variance(tmp,chosen_n));}
        else if(opt==1) {gr_lnBpC_dispersion->SetPoint(j , lnC_start+j*lnC_one_step, variance(tmp,chosen_n));}

    }
    
    if(c_BpC_dispersion_ON)
    {//===========================================================================
    TCanvas *c_BpC_dispersion = new TCanvas(c_name,"BpC_dispersion",1000,500);    
    if(opt==0) {gr_BpC_dispersion->Draw("Ap");}
    else if(opt==1) {gr_lnBpC_dispersion->Draw("Ap");}
    }//____________________________________________________________________________

    gr_BpC_dispersion->Sort(&TGraph::CompareY);
    double x_out,y_out;
    gr_BpC_dispersion->GetPoint(0, x_out, y_out);
    cout<<"------ lowest point : ("<<x_out<<" "<<fixed<<setprecision(15)<<y_out<<setprecision(3)<<" )"<<endl;
    C_intersection = x_out;

    if(c_BpC_dd_ON)
    {//===========================================================================
    TCanvas *c_BpC_dd = new TCanvas(c_name+"_2","#Delta Bp_C",1000,500);
    
    gr_BpC_dd[0]->Draw("apl");
    for(int i=0;i<chosen_n;i++)
    {
        gr_BpC_dd[i]->Draw("samep");
    }

    auto legend_c_BpC_dd = new TLegend(0.80,0.10,0.95,0.35); //downright
        //legend_c_BpC_dd->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"C");
    for(int i=0;i<chosen_n  ;i++)
    {
        legend_c_BpC_dd->AddEntry(gr_BpC_dd[i],gr_BpC_dd[i]->GetTitle(),"p" );
    }
    legend_c_BpC_dd->Draw("same");
    }//____________________________________________________________________________
    
}



//======================== 20230531 mvq_C each species
void Show_mvqC_each(IONSpecies* ISS)
{

    TGraphErrors* grerr_kT = new TGraphErrors();
    AxisFormat(grerr_kT,ThisParaInfo,"T[ns]","k_mvqC");

    if(MASS_VER==1)
    {
        TCanvas* c_mvqC_each[MAX_IONSPECIES];
        for(int i=0;i<NSpecies;i++)
        {
            if(ISS[i].gr_mvqC->GetN()<=Show_mvqC_each_MIN){continue;}
            c_mvqC_each[i] =  new TCanvas("c_mvqC_"+ISS[i].Aname ,"c_mvqC_"+ISS[i].Aname ,1600,800);
            ISS[i].gr_mvqC->Draw("ap");  
            ISS[i].f0_mvq_AME->Draw("same");
            c_mvqC_each[i]->Print(FILEPATH+ ISS[i].Aname+"_mvqC.png");
            //cout<<"ISS[i].gr_mvqC_n= "<<ISS[i].gr_mvqC->GetN()<<endl;
    
            if(MASS_VER>2)delete c_mvqC_each[i] ;
        }
        if(Do_Generate_k_mvqC_ON)
        {
            cout<<"--------------Generate  k_mvqC ON in Show_mvqC_each for v1---------------------------"<<endl;
            grerr_kT->SetTitle(ThisParaInfo+"_v1");
            for(int i=0;i<NSpecies  ;i++)
            {
                if(ISS[i].gr_mvqC->GetN()>Show_mvqC_each_MIN)
                {
                    
                    ISS[i].Fit_mvqC(1);  // get ISS[i].k_mvqC
                    ISS[i].have_Fit_mvqC = true;
                    grerr_kT->SetPoint(grerr_kT->GetN(), ISS[i].AveT, ISS[i].fitfun_pol1_mvqC->GetParameter(1) );
                    grerr_kT->SetPointError(grerr_kT->GetN()-1, 0, ISS[i].fitfun_pol1_mvqC->GetParErrors()[1]);
                    
                }
            }
        }

    }
    else if(MASS_VER>=3)
    {
        TCanvas* c_mvqC_each_VE[MAX_IONSPECIES];
        for(int i=0;i<NSpecies;i++)
        {
            if(ISS[i].grerr_mvqC_VE->GetN()<=Show_mvqC_each_MIN){continue;}
            c_mvqC_each_VE[i] =  new TCanvas("c_mvqC_VE_"+ISS[i].Aname ,"c_mvqC_VE_"+ISS[i].Aname ,1600,800);
            ISS[i].grerr_mvqC_VE->Draw("ap");  
            ISS[i].f0_mvq_AME->Draw("same");
            c_mvqC_each_VE[i]->Print(FILEPATH+ ISS[i].Aname+"_mvqC_VE.png");
            //cout<<"ISS[i].gr_mvqC_n= "<<ISS[i].gr_mvqC_n<<endl;
            
            //选择只留下哪些核的窗口
            //if(ISS[i].Aname=="17Ne"||ISS[i].Aname=="13N"||ISS[i].Aname=="11C"||ISS[i].Aname=="25Si"){}
            //else    { delete c_mvqC_each[i] ;}
            
        }
        if(Do_Generate_k_mvqC_ON)
        {
            cout<<"--------------Generate  k_mvqC ON in Show_mvqC_each for vE---------------------------"<<endl;
            grerr_kT->SetTitle(ThisParaInfo+"_VE");
            for(int i=0;i<NSpecies  ;i++)
            {
                if(ISS[i].grerr_mvqC_VE->GetN()>Show_mvqC_each_MIN)
                {
                    ISS[i].Fit_mvqC(2);
                    ISS[i].have_Fit_mvqC_VE = true;
                    grerr_kT->SetPoint(grerr_kT->GetN(), ISS[i].AveT, ISS[i].fitfun_pol1_mvqC_VE->GetParameter(1) );
                    grerr_kT->SetPointError(grerr_kT->GetN()-1, 0, ISS[i].fitfun_pol1_mvqC_VE->GetParErrors()[1]);
                }
            }  
        }
    }
    else{}

    
    if(Do_Generate_k_mvqC_ON)
    {
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TCanvas *c_kT = new TCanvas("c_kT","c_kT",1000,500);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        grerr_kT->SetMarkerSize(2);
        grerr_kT->GetYaxis()->SetRangeUser(-0.001,0.001);
        grerr_kT->DrawClone("ap");
        //grerr_kT->Print();
       
        lat_n->SetTextAngle(90);
        lat_n->SetTextFont(43);
        lat_n->SetTextColor(kAzure+7);
        lat_n->SetTextSize(25);
        for(int i=0;i<NSpecies  ;i++)
        {
            if(MASS_VER==1)
            if(ISS[i].have_Fit_mvqC)
            {
                lat_text = ISS[i].Aname+ strtmp.Format(" : %d ",ISS[i].N);
                lat_n->DrawLatex(ISS[i].AveT,ISS[i].fitfun_pol1_mvqC->GetParameter(1),lat_text);
            }

            if(MASS_VER>=3)
            if(ISS[i].have_Fit_mvqC_VE)
            {
                lat_text = ISS[i].Aname+ strtmp.Format(" : %d ",ISS[i].N);
                lat_n->DrawLatex(ISS[i].AveT,ISS[i].fitfun_pol1_mvqC_VE->GetParameter(1),lat_text);
            }
        }
        lat_n->SetTextAngle(0);

        c_kT->Print(FILEPATH_k_T  +"kT_"+strtmp.Format("%d.png",scan_loop_i) );
        if(LOOP_ON)delete c_kT;
    }
    
}
void Show_dmC_each(IONSpecies* ISS)
{
    TCanvas* c_dmC_each[MAX_IONSPECIES];
    for(int i=0;i<NSpecies;i++)
    {
        c_dmC_each[i] =  new TCanvas("c_dmC_"+ISS[i].Aname ,"c_dmC_"+ISS[i].Aname ,1600,800);
        ISS[i].gr_dmC->DrawClone("ap");  
        f_zero->Draw("same");
        c_dmC_each[i]->Print(FILEPATH+ ISS[i].Aname+"_dmC.png");
        //cout<<"ISS[i] gr_dmC n= "<<ISS[i].gr_dmC->GetN()<<endl;

        if(MASS_VER>2)delete c_dmC_each[i] ;
    }
    if(MASS_VER>2)
    {
        TCanvas* c_dmC_each_VE[MAX_IONSPECIES];
        for(int i=0;i<NSpecies;i++)
        {
            c_dmC_each_VE[i] =  new TCanvas("c_dmC_VE_"+ISS[i].Aname ,"c_dmC_VE_"+ISS[i].Aname ,1600,800);
            ISS[i].grerr_dmC_VE->DrawClone("ap");  
            f_zero->Draw("same");
            c_dmC_each_VE[i]->Print(FILEPATH+ ISS[i].Aname+"_dmC_VE.png");
            //cout<<"ISS[i] grerr_dmC_VE_n= "<<ISS[i].grerr_dmC_VE->GetN()<<endl;
            if(ISS[i].Aname=="17Ne"||ISS[i].Aname=="13N"||ISS[i].Aname=="11C"||ISS[i].Aname=="25Si"){}
            else    { delete c_dmC_each_VE[i] ;}
        }
    }
}

void Show_mvqC_each_h2(IONSpecies* ISS)
{ 
    
    TCanvas* c_mvqC_each_h2[MAX_IONSPECIES];
    
    for(int i=0;i<NSpecies;i++)
    {
        c_mvqC_each_h2[i] =  new TCanvas("c_mvqC_h2_"+ISS[i].Aname ,"c_mvqC_h2_"+ISS[i].Aname ,1600,800);
        ISS[i].h2_mvqC->DrawClone("colz");  
        c_mvqC_each_h2[i]->Print(FILEPATH+ ISS[i].Aname+"_h2mvqC.png");
        //delete ISS[i].h2_mvqC;

    }
    
    /*    //预先把所有种类平均mvq 输出 ， 用于在前端定义h2
    ofstream outfile_mvq;
    outfile_mvq.open("36Ar_22Si_each_avemvq.txt");
    for(int i=0;i<NSpecies;i++)
    {
        outfile_mvq<<ISS[i].Mvq<<endl;
    }
    outfile_mvq.close();
    */
}

//--------------------------------------------
void Show_h_iont_merr_each(IONSpecies* ISS)
{
    TCanvas* c_iontmerr_each[MAX_IONSPECIES];
    for(int i=0;i<NSpecies;i++)
    {
        c_iontmerr_each[i] =  new TCanvas("c_h_iont_merr_"+ISS[i].Aname ,"c_h_iont_merr_"+ISS[i].Aname ,1600,800);
        ISS[i].h_iont_mass_err->Draw();  
        
        c_iontmerr_each[i]->Print(FILEPATH+ ISS[i].Aname+"_h_iont_mass_err.png");
        if(ISS[i].Aname=="17Ne"||ISS[i].Aname=="13N"||ISS[i].Aname=="11C"||ISS[i].Aname=="25Si"){}
        else    { delete c_iontmerr_each[i] ;}
    }
    
}
void Show_h_refcal_merr_each(IONSpecies* ISS)
{
    TCanvas* c_refcalmerr_each[MAX_IONSPECIES];
    for(int i=0;i<NSpecies;i++)
    {
        c_refcalmerr_each[i] =  new TCanvas("c_h_refcal_merr_"+ISS[i].Aname ,"c_h_refcal_merr_"+ISS[i].Aname ,1600,800);
        ISS[i].h_each_ref_cal_mass_err->Draw();  
        
        c_refcalmerr_each[i]->Print(FILEPATH+ ISS[i].Aname+"_h_refcal_merr.png");
        if(ISS[i].Aname=="17Ne"||ISS[i].Aname=="13N"||ISS[i].Aname=="11C"||ISS[i].Aname=="25Si"){}
        else    { delete c_refcalmerr_each[i] ;}
    }
}
void Show_h_iont_refs_chi_each(IONSpecies* ISS)
{
    TCanvas* c_h_iont_chi_each[MAX_IONSPECIES];
    for(int i=0;i<NSpecies;i++)
    {
        c_h_iont_chi_each[i] =  new TCanvas("c_h_iont_chi_each_"+ISS[i].Aname ,"c_h_iont_chi_each_"+ISS[i].Aname ,1600,800);
        ISS[i].h_iont_chi_n->Draw();  
        
        c_h_iont_chi_each[i]->Print(FILEPATH+ ISS[i].Aname+"_h_iont_chi_each.png");
        if(ISS[i].Aname=="17Ne"||ISS[i].Aname=="13N"||ISS[i].Aname=="11C"||ISS[i].Aname=="25Si"){}
        else    { delete c_h_iont_chi_each[i] ;}
    }
}
//_________________________   each species



// ======================== 20230518 所有种类 dA0-T 散点
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FIT_TS(double C, double L, double gt, double k1, double k2, double T)
{
    double tmp=0; //       gamma^2 / (gamma^2-gt^2)
    double v = C/T;  // m/ns
    double gm = 1/sqrt(1-v*v/V_c/V_c);
    tmp = gm*gm/(gm*gm - gt*gt);
    return L/C*(1-tmp) + k1*tmp - k2*v*tmp;
}
//所有种类 dA0-T 散点
void Do_dA0_T(IONSpecies* ISS, double& k_L,double& k_C,double& kgt,double& k1,double& k2,double& ktd,double& CHI, double Ci,double dCi, int kk)
{
    bool Draw_k_dA0_T_ON        = 0;
    bool Draw_each_grerr_dA0T   = 0;
    bool Draw_Ts_T_ON           = 0;
    bool Draw_res_ON            = 0;
    bool Draw_Bp_dispersion_ON  = 0;
    bool Draw_chi_n_ON          = 0;
    bool Draw_k_dTsT_fit_ON     = 0;
    bool Show_each_grerr_dA0T_n = 0;//看一下选取到的每种离子有多少个
    int dA0_T_FIT_POINTS_MIN=30;

    double k_Li=0; 
    double dk_L = 0.001;

    double k_LC; // L/C
    // choose points to be fitted in the relation dTs/dT~T
    int xN = 0;
    double xT[20]={0};
    double ydTsT[20]={0};
    double ydTsT_fit[20]={0};
    double res_ydTsT_fit[20]={0}; //  ((fit-y)
    double ydTsTerr[20]={0};

    double time_divide = 900000;
    double tmpx,tmpy;
    TString lat_text_name[30];
    double Ave_Ts[30],Ave_T[30];  // 非连续序号
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].grerr_dA0_T_n = 0;
        ISS[i].grerr_dA0_T  = new TGraphErrors();
        AxisFormat(ISS[i].grerr_dA0_T,ISS[i].Aname+"  dA0_T"," T[ns] ", "dA0[ns]");
        ISS[i].grerr_dA0_T->SetMarkerSize(1.0);
        ISS[i].grerr_dA0_T->SetLineWidth(1.0);
    }
    
    // predefine fitting ionspecies and its range
    /*
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Aname =="10C") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= 88.24;ISS[i].dA0TfitMAX= 88.38;  continue;  }
        if(ISS[i].Aname =="11C") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 93.22;  continue;  }
        if(ISS[i].Aname =="14O") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 90.74;  continue;  }
        if(ISS[i].Aname =="12N") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 89.75;  continue;  }
        if(ISS[i].Aname =="9C")  {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  } 
        if(ISS[i].Aname =="17Ne"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 89.3;   continue;  }
        if(ISS[i].Aname =="13N") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 93.9;   continue;  }
        if(ISS[i].Aname =="15O") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 94.4;   continue;  }
        if(ISS[i].Aname =="18Ne"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 92.2;   continue;  }
        if(ISS[i].Aname =="21Mg"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
        if(ISS[i].Aname =="20Mg"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 88.35;  continue;  }
        if(ISS[i].Aname =="23Si"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= 87.55;ISS[i].dA0TfitMAX= 87.67; continue;  }
        if(ISS[i].Aname =="7Be") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
        if(ISS[i].Aname =="22Al"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
        if(ISS[i].Aname =="17F") {ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
        if(ISS[i].Aname =="24Si"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
        if(ISS[i].Aname =="20Na"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
        if(ISS[i].Aname =="22Si"){ISS[i].Has_k_dA0_T=1;ISS[i].dA0TfitMIN= -999;ISS[i].dA0TfitMAX= 999;    continue;  }
    }
    */
    // select point to set TGraph
    for(int j=0;j<ions_n  ;j++){ions[j].Do_dA0_T_flag=0;}
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].time>time_divide ){continue;} //divide by time-- T is changed in the experiment
        //dA0 -- ps
        if( ions[j].C>Ci-dCi && ions[j].C<Ci+dCi)
            {
            ions[j].Do_dA0_T_flag=1;
            ISS[ions[j].Species].grerr_dA0_T->SetPoint(ISS[ions[j].Species].grerr_dA0_T_n , ions[j].T, ions[j].dA0*0.001);
            ISS[ions[j].Species].grerr_dA0_T->SetPointError(ISS[ions[j].Species].grerr_dA0_T_n , 0, sqrt(ions[j].dA0err)*0.001);
            ISS[ions[j].Species].grerr_dA0_T_n++;
            }
    }
    // choose available ionspecies 
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].grerr_dA0_T_n>dA0_T_FIT_POINTS_MIN)
        {
            ISS[i].Has_k_dA0_T=1;
            if(Show_each_grerr_dA0T_n)cout<<endl<<ISS[i].Aname<<" "<<ISS[i].AveT<<" grerr_dA0_T_n = "<<ISS[i].grerr_dA0_T_n<<endl;
            //!!!!!! 7Be abnormal
            if(ISS[i].Aname=="7Be"){ISS[i].Has_k_dA0_T=0;}
        }
    }
    //return;/////////####debug
    // do fit
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Has_k_dA0_T)
        {
            ISS[i].Fit_dA0_T(1);
            
            //cout<<Info_fitfun_pol1(ISS[i].ff_dA0_T_pol1)<<endl;
        }
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas*c_dA0_T_ISS[30];
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Has_k_dA0_T==0)continue;

        if(Draw_each_grerr_dA0T)
        {
        c_dA0_T_ISS[i] = new TCanvas("c_dA0_T_"+ISS[i].Aname,"c_dA0_T_"+ISS[i].Aname  ,1000,500);
        ISS[i].grerr_dA0_T->DrawClone("ap");

        }

        //c_dA0_T_ISS[i]->Print(ISS[i].folder_path+"dA0_T.png");
        //c_dA0_T_ISS[i]->Print(FILEPATH+ ISS[i].Aname+"_dA0_T.png");
    }

    //k_dA0_T
    
    double tmpxT=0;
    TGraphErrors* grerr_k_dA0_T = new TGraphErrors();
    int grerr_k_dA0_T_n=0;
    AxisFormat(grerr_k_dA0_T,"","AveT [ns]","d(dA0)/dT");
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Has_k_dA0_T)
        {
            tmpxT = ISS[i].grerr_dA0_T->GetMean(1);
            grerr_k_dA0_T->SetPoint(grerr_k_dA0_T_n,tmpxT,ISS[i].k_dA0_T);
            //ISS[i].grerr_dA0_T->GetMean(1);

            grerr_k_dA0_T->SetPointError(grerr_k_dA0_T_n,0,ISS[i].k_dA0_T_err);
            lat_text_name[grerr_k_dA0_T_n] = ISS[i].name_latex;
            grerr_k_dA0_T_n++;
            // choose points to be fitted
            if(tmpxT>632.0)
            {
                xT[xN] = tmpxT;
                ydTsT[xN] = ISS[i].k_dA0_T;
                ydTsTerr[xN] = ISS[i].k_dA0_T_err;
                xN++;
            }
        }
    }        

    TGraphErrors* grerr_k_dA0_T_data = new TGraphErrors(xN,xT,ydTsT,0,ydTsTerr);
    //grerr_k_dA0_T_data->Print();
    AxisFormat(grerr_k_dA0_T_data,"","AveT [ns]","d(dA0)/dT");

    if(Draw_k_dA0_T_ON)
    {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_k_dA0_T = new TCanvas("c_k_dA0_T","c_k_dA0_T"  ,1200,600);
    grerr_k_dA0_T->Draw("Ap");
    lat_n->SetTextColor(kAzure);
    lat_n->SetTextFont(43);
    lat_n->SetTextSize(20);
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Has_k_dA0_T)
        {
            lat_text = ISS[i].name_latex;
            lat_n->DrawLatex(ISS[i].grerr_dA0_T->GetMean(1),ISS[i].k_dA0_T,lat_text);
        }   
    }
    }//Draw_k_dA0_T_ON

    // line Ts-T
    TF1* ff_TST = new TF1("ff_TST","pol1",500,700);
    TGraphErrors* grerr_Ts_T = new TGraphErrors();
    int grerr_Ts_T_n=0;
    AxisFormat(grerr_Ts_T,"","AveT [ns]","Ts [ns]");
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Has_k_dA0_T)
        {
            grerr_Ts_T->SetPoint(grerr_Ts_T_n,ISS[i].grerr_dA0_T->GetMean(1),ISS[i].grerr_dA0_T->GetMean(2));
            Ave_Ts[i] = ISS[i].grerr_dA0_T->GetMean(2);
            Ave_T[i] = ISS[i].grerr_dA0_T->GetMean(1);
            grerr_Ts_T->SetPointError(grerr_Ts_T_n,
                ISS[i].grerr_dA0_T->GetRMS(1)/sqrt(ISS[i].grerr_dA0_T->GetN()),
                ISS[i].grerr_dA0_T->GetRMS(2)/sqrt(ISS[i].grerr_dA0_T->GetN()) );
            grerr_Ts_T_n++;
        }   
    }

    grerr_Ts_T->Fit(ff_TST,"q");
    cout<<Info_fitfun_pol1(ff_TST)<<endl;
    //########## get td ####################
    ktd = -1.0*ff_TST->GetParameter(0);
    //########## get L/C ####################
    k_LC = ff_TST->GetParameter(1);

    if(Draw_Ts_T_ON)
    {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_Ts_T = new TCanvas("c_Ts_T","c_Ts_T"  ,1200,600);
    grerr_Ts_T->Draw("ap");
    
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ISS[i].Has_k_dA0_T)
        {
            lat_text = ISS[i].name_latex;
            lat_n->DrawLatex(ISS[i].grerr_dA0_T->GetMean(1),ISS[i].grerr_dA0_T->GetMean(2),lat_text);
        }   
    }
    }//Draw_Ts_T_ON

    TGraphErrors* grerr_res2 = new TGraphErrors();
    int grerr_res2_n=0;
    AxisFormat(grerr_res2,"","AveT [ns]","residue #Delta^{2}");
    grerr_res2->SetMarkerColor(2);
    grerr_res2->SetLineColor(2);
    Get_residue(grerr_Ts_T,ff_TST,grerr_res2);
    /*
    TGraphErrors* grerr_res = new TGraphErrors();
    int grerr_res_n=0;
    AxisFormat(grerr_res,"","AveT [ns]","residue #Delta^{2}");
    for(int i=0;i<NSpecies  ;i++)
    {if(ISS[i].Has_k_dA0_T)
        {grerr_res->SetPoint(grerr_res_n,ISS[i].grerr_dA0_T->GetMean(1),   pow( ISS[i].grerr_dA0_T->GetMean(2)-ff_TST->Eval(ISS[i].grerr_dA0_T->GetMean(1)) ,2 )  );
        grerr_res_n++;
        }   }
    */
    if(Draw_res_ON)
    {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_res = new TCanvas("c_res","c_res"  ,1200,600);
    //grerr_res->Draw("ap");
    grerr_res2->Draw("ap");
    for(int i=0;i<grerr_res2->GetN()  ;i++)
    {
        grerr_res2->GetPoint(i,tmpx,tmpy);
        lat_n->DrawLatex(tmpx,tmpy,lat_text_name[i]);       
    }
    }//Draw_res_ON


    // scan L to make all species have the same Bp
    int k_Li_n=0;
    //k_Li = 18.033;
    TGraph* gr_Bp = new TGraph();
    AxisFormat(gr_Bp,"gr_Bp for each ionspecies"," L ","Bp StdDev");

    for(k_Li_n = 0;k_Li_n<100;k_Li_n++)
    {
        k_Li = 18.0+k_Li_n*dk_L;
        //construct Bp array
        double Bp_each[29]={0};
        int Bp_each_n=0;
        for(int i=0;i<NSpecies  ;i++)
        {
            if(ISS[i].Has_k_dA0_T) //choose ionspecies to calculate Bp dispersion
            {
                //cout<<ISS[i].Aname<<" "<<ISS[i].Calculate_Bp_from_L(k_Li,Ave_Ts[i],ktd)<<endl;
                Bp_each[Bp_each_n] = ISS[i].Calculate_Bp_from_L(k_Li,Ave_Ts[i],ktd);
                Bp_each_n++;
            }
        }
        //cout<<k_Li_n<<" "<<k_Li<<" "<<sqrt(variance(Bp_each,Bp_each_n))<<endl;
        gr_Bp->SetPoint(k_Li_n,k_Li, sqrt(variance(Bp_each,Bp_each_n)) );
    }
    if(Draw_Bp_dispersion_ON)
    {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_Bp_dispersion = new TCanvas("c_Bp_dispersion","c_Bp_dispersion"  ,1200,600);
    gr_Bp->Draw("ap");
    }//Draw_Bp_dispersion_ON

    gr_Bp->Sort(&TGraph::CompareY);
    double Bp_variance_min;
    //############### get L #####################
    gr_Bp->GetPoint(0, k_L, Bp_variance_min);
    cout<<"------Bp variance  lowest point : ("<<fixed<<setprecision(4)<<k_L<<" "<<Bp_variance_min<<" )"<<endl;
    //############### get C #####################
    k_C = k_L / k_LC;
    cout<<" L/C ="<<k_LC<<" L= "<<k_L<<" C= "<<k_C<<endl;

    // now we have : k_L,k_LC,k_C,

    //scan 3 para
    TGraph* gr_chi_n_count = new TGraph();
    AxisFormat(gr_chi_n_count,""," count "," #chi_{n}");
    gr_chi_n_count->SetMarkerSize(1);
    TGraph* gr_chi_n_record_abs = new TGraph();
    int count=0;
    int kgti,k1i,k2i=0;
    double kgt_scan,k1_scan,k2_scan;
    double kgt_down= 1.350  ; double dkgt = 0.001;int kgt_n = 50;
    double k1_down = -0.1  ; double dk1 = 0.001;  int k1_n = 200; 
    double k2_down = -0.2  ;  double dk2 = 0.002;  int k2_n = 200; 
    double chi_n =0; // sum( (fit-y)^2/yerr^2 ) /n
    if(kk==199)
    {
    cout<<"------Begin 3 k para scan:"<<endl
    <<" para1: gammat ----   ||"<<kgt_down<<" ~ "<<kgt_down+dkgt*kgt_n<<" || "<<kgt_n<<" x "<<dkgt<<endl
    <<" para2: k1 dL/dC ---- ||"<<k1_down<<" ~ "<<k1_down+dk1*k1_n<<" || "<<k1_n<<" x "<<dk1<<endl
    <<" para3: k2 dtd/dC ----|| "<<k2_down<<" ~ "<<k2_down+dk2*k2_n<<" || "<<k2_n<<" x "<<dk2<<endl
    <<"____________"<<endl;
    }

    for(kgti=0;kgti<kgt_n  ;kgti++)
    {
        kgt_scan = kgt_down+kgti*dkgt;
        for(k1i=0;k1i<k1_n;k1i++)
        {
            k1_scan = k1_down + k1i*dk1;
            for(k2i=0;k2i<k2_n;k2i++)
            {
                k2_scan = k2_down + k2i*dk2;
                
                chi_n=0;
                for(int j = 0;j<xN; j++)
                {
                    ydTsT_fit[j] = FIT_TS(k_C, k_L, kgt_scan, k1_scan,k2_scan,xT[j]); 
                    //(double C, double L, double gt, double k1, double k2, double T)
                    chi_n+= pow( ((ydTsT_fit[j]- ydTsT[j])/ydTsTerr[j] ),2);
                }
                chi_n/=xN;
                chi_n = sqrt(chi_n);
                gr_chi_n_count->SetPoint(count,count,chi_n);
                gr_chi_n_record_abs->SetPoint(count,count,abs(chi_n-1));
                count++;
                //if(count%200000==0)cout<<" k scan now dealing with:  "<<count<<endl;
        
            }
        }
    }
    if(Draw_chi_n_ON)
    {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_chi_n_count = new TCanvas("c_chi_n_count","c_chi_n_count"  ,1200,600);
    gr_chi_n_count->Draw("Ap");
    }//Draw_chi_n_ON
    
    //gr_chi_n_count->Sort(&TGraph::CompareY);
    gr_chi_n_record_abs->Sort(&TGraph::CompareY);
    double count_best,chi_n_best;
    //######################## get chi n best
    gr_chi_n_record_abs->GetPoint(0,count_best,chi_n_best);
    chi_n_best += 1;
    cout<<" best at count : "<<count_best<<" chi_n_best = "<<chi_n_best<<endl;
    CHI = chi_n_best;

    TGraphErrors* grerr_k_dA0_T_fit = new TGraphErrors(xN,xT,0,ydTsT,ydTsTerr);
    AxisFormat(grerr_k_dA0_T_fit,"","AveT [ns]","fit: d(dA0)/dT");
    grerr_k_dA0_T_fit->SetMarkerColor(kRed);

    count=0;
    for(kgti=0;kgti<kgt_n  ;kgti++)
    {
        kgt_scan = kgt_down+kgti*dkgt;
        for(k1i=0;k1i<k1_n;k1i++)
        {
            k1_scan = k1_down + k1i*dk1;
            for(k2i=0;k2i<k2_n;k2i++)
            {
                k2_scan = k2_down + k2i*dk2;

                if(count == count_best)
                {
                    cout<<" ( gt = "<<kgt_scan<<" , k1 = "<<k1_scan<<" , k2 = "<<k2_scan<<" ) "<<endl;
                    //#########################
                    kgt=kgt_scan;
                    k1=k1_scan;
                    k2=k2_scan;
                    for(int j = 0;j<xN; j++)
                    {
                        ydTsT_fit[j] = FIT_TS(k_C, k_L, kgt_scan, k1_scan,k2_scan,xT[j]);
                        grerr_k_dA0_T_fit->SetPoint(j,xT[j],ydTsT_fit[j]);

                    }
                }
                count++;
            }
        }
    }
    
    if(Draw_k_dTsT_fit_ON)
    {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_k_dTsT_T = new TCanvas("c_k_dTsT_T","c_k_dTsT_T"  ,1200,600);
    grerr_k_dA0_T_data->Draw("ap");
    grerr_k_dA0_T_fit->Draw("psame");
    c_k_dTsT_T->Print(FILEPATH_DodA0T+ strtmp.Format("k_dA0_T_fit_%d.png",kk) );
    if(kk>=1)delete c_k_dTsT_T;
    }//Draw_k_dTsT_fit_ON
    

    cout<<" C ="<<k_C<<", L = "<<k_L<<" ddT= "<<ktd<<", gt= "<<kgt<<", dLdC= "<<k1<<", ddTdC= "<<k2<<endl;

    for(int i=0;i<NSpecies  ;i++)
    {
        delete ISS[i].grerr_dA0_T;
    }

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//====================== 20230616 recalculate v C gt Bp
void Recalculate_after_dA0T()
{
    TGraph* gr_CC_from_dA0T_readin = new TGraph();
    AxisFormat(gr_CC_from_dA0T_readin,"","C[m]","C NEW[m]");
    TGraph_from_infile("INPUT//dA0T//gr_CC_from_dA0T_save.txt",gr_CC_from_dA0T_readin);
    TGraph* gr_LC_from_dA0T_readin = new TGraph();
    AxisFormat(gr_LC_from_dA0T_readin,"","C[m]","L[m]");
    TGraph_from_infile("INPUT//dA0T//gr_LC_from_dA0T_save.txt",gr_LC_from_dA0T_readin);
    //TCanvas* c_check_LC_readin = new TCanvas("c_check_LC_readin","c_check_LC_readin",1000,500);
    //gr_LC_from_dA0T_readin->Draw("ap");
    
    TGraph* gr_ddtC_from_dA0T_readin = new TGraph();
    AxisFormat(gr_ddtC_from_dA0T_readin,"","C [m]","#Deltat_{d} [ns]");
    TGraph_from_infile("INPUT//dA0T//gr_ddtC_from_dA0T_save.txt",gr_ddtC_from_dA0T_readin);
    //TCanvas* c_check_ddtC_readin = new TCanvas("c_check_ddtC_readin","c_check_ddtC_readin",1000,500);
    //gr_ddtC_from_dA0T_readin->Draw("ap");

    /*
    TGraph* gr_LC_2= new TGraph();
    AxisFormat(gr_LC_2,"","C[m]","L[m]");
    TGraph* gr_ddtC_2 = new TGraph();
    AxisFormat(gr_ddtC_2,"","C [m]","#Deltat_{d} [ns]");
    TGraph* gr_LC_3= new TGraph();
    AxisFormat(gr_LC_3,"","C[m]","L[m]");
    TGraph* gr_ddtC_3 = new TGraph();
    AxisFormat(gr_ddtC_3,"","C [m]","#Deltat_{d} [ns]");
    Smooth_DW(gr_LC_from_dA0T_readin,gr_LC_2,5);
    Smooth_DW(gr_ddtC_from_dA0T_readin,gr_ddtC_2,5);
    Smooth_DW(gr_LC_2,gr_LC_3,5);
    Smooth_DW(gr_ddtC_2,gr_ddtC_3,5);
    Draw_one_TGraph(gr_LC_2,"c_gr_LC_2");
    Draw_one_TGraph(gr_ddtC_2,"c_gr_ddtC_2");
    Draw_one_TGraph(gr_LC_3,"c_gr_LC_3");
    Draw_one_TGraph(gr_ddtC_3,"c_gr_ddtC_3");
    */
    TGraph* gr_LC_i[20];
    TGraph* gr_ddtC_i[20];
    
    for(int i=0;i<20  ;i++)
    {
        gr_LC_i[i] = new TGraph();
        gr_ddtC_i[i] = new TGraph();
    }
    int smooth_k= 3;
    Smooth_DW(gr_LC_from_dA0T_readin,gr_LC_i[0],5, smooth_opt);
    Smooth_DW(gr_ddtC_from_dA0T_readin,gr_ddtC_i[0],5, smooth_opt);
    for(int i=0;i<smooth_k-1  ;i++)
    {
        Smooth_DW(gr_LC_i[i],gr_LC_i[i+1],5, smooth_opt);
        Smooth_DW(gr_ddtC_i[i],gr_ddtC_i[i+1],5, smooth_opt);
    }
    AxisFormat(gr_LC_i[smooth_k-1],"","C[m]","L[m]");
    AxisFormat(gr_ddtC_i[smooth_k-1],"","C[m]","#Deltat_{d} [ns]");

    Draw_one_TGraph(gr_LC_i[smooth_k-1],"c_gr_LC_smooth");
    Draw_one_TGraph(gr_ddtC_i[smooth_k-1],"c_gr_ddtC_smooth");
    

    TGraph* gr_new_L_C = new TGraph();
    AxisFormat(gr_new_L_C,"","C [m]","L [m]");

    TGraph* gr_new_ddt_C = new TGraph();
    AxisFormat(gr_new_ddt_C,"","C [m]","ddt [ns]");
    TGraph* gr_new_C_C = new TGraph();
    AxisFormat(gr_new_C_C,"","C [m]","new C [m]");
    TH1F* h_v_change_ratio = new TH1F("h_v_change_ratio","h_v_change_ratio",1000,-0.1,0.1);
    C_ave_all=0;
    double C_change=0;
    double ddT_new=0;
    double L_new=0;
    double Cmin=999;Cmax=-1;
    for(int i=0;i<ions_n;i++)
    {
        //C_change =gr_CC_from_dA0T_readin->Eval(ions[i].C);
        C_change = ions[i].C; 
        
        //ddT_new = gr_ddtC_from_dA0T_readin->Eval(C_change);
        //L_new = gr_LC_from_dA0T_readin->Eval(C_change); 
        ddT_new = gr_ddtC_i[smooth_k-1]->Eval(C_change);
        L_new = gr_LC_i[smooth_k-1]->Eval(C_change);    
        v = L_new/( (ions[i].dA0/1000.0) + ddT_new);  //[m/ns]
        h_v_change_ratio->Fill(  (v - ions[i].v) /ions[i].v); 
    
        v_err = v*(v/L)/1000.0*sqrt(ions[i].dA0err);  // v_err = v*v/L * sigma_dA0
        ions[i].v=v;  // ions.v [m/ns]
        ions[i].v_err = v_err;
        C = v*ions[i].T;  //[m]     
        gr_new_C_C->SetPoint(i,ions[i].C,C);
        ions[i].C=C;
        ions[i].C_err2 = sigmaU_XY2(ions[i].C,ions[i].v,ions[i].v_err,ions[i].T,ions[i].T_err);
        C_ave_all+=ions[i].C;
        if(ions[i].C>Cmax)Cmax=ions[i].C;
        if(ions[i].C<Cmin)Cmin=ions[i].C;
    
        gr_new_L_C->SetPoint(i,C,L_new);
        gr_new_ddt_C->SetPoint(i,C,ddT_new);

        //###################################################################################################
        //gammat= sqrt( (1/(1 - pow(v/V_c,2))) / (1 - ions[i].A2*2/ions[i].T*((ions[i].dA0/1000.0)+ddT)/ions[i].dA1 ) );
        ions[i].Calculate_gt(ddT_new);
        
        //###################################################################################################
        //ions[i].gammat=gammat;
        gammat = ions[i].gammat;

        //############ get Bp of all ions
        ions[i].Calculate_Bp(ionspecies[ions[i].Species].Mass);   //AME mass  of this species
    
    }
    //TCanvas* c_h_v_change_ratio = new TCanvas("c_h_v_change_ratio","c_h_v_change_ratio",1000,500);
    //AxisFormat(h_v_change_ratio,"", "v_change_raio","count");
    //h_v_change_ratio->Draw();
    

    //if(Draw_each_T_C_ON)Draw_each_T_C(ionspecies, 1,1);
    Draw_gtC_all(2,0);
    Draw_one_TGraph(gr_new_L_C ,"c_new_L_C");
    Draw_one_TGraph(gr_new_ddt_C ,"c_new_ddt_C");
    Draw_one_TGraph(gr_new_C_C ,"c_new_C_C");
}

// __________________________ 20230518 所有种类 dA0-T 散点

//=============== 20230614 C-T 散点 =========================
void Draw_each_T_C(IONSpecies*ISS, double T_move,double sT_move)
{
    double time_divide = 900000;
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].gr_T_C_n = 0;
        ISS[i].gr_T_C  = new TGraph();
        AxisFormat(ISS[i].gr_T_C,ISS[i].Aname+"  T_C "," T[ns] ", "C[m]");
        ISS[i].gr_T_C->SetMarkerSize(1.0);
        ISS[i].gr_T_C->SetLineWidth(1.0);
    }

    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].time>time_divide ){continue;} //divide by time-- T is changed in the experiment
        //dA0 -- ps
        
        //if( (ions[j].dA0*0.001)>ISS[ions[j].Species].dA0TfitMIN &&( (ions[j].dA0*0.001)<ISS[ions[j].Species].dA0TfitMAX) )
        //    if(    (ions[j].T>ISS[ions[j].Species].AveT*T_move-ISS[ions[j].Species].SigmaT*sT_move)
        //        && (ions[j].T<ISS[ions[j].Species].AveT*T_move+ISS[ions[j].Species].SigmaT*sT_move)        )
            {
                ISS[ions[j].Species].gr_T_C->SetPoint(ISS[ions[j].Species].gr_T_C_n , ions[j].T, ions[j].C);
                ISS[ions[j].Species].gr_T_C_n++;
            }
    }
    TCanvas*c_T_C_ISS[30];
    for(int i=0;i<NSpecies  ;i++)
    {
        //if(ISS[i].Has_k_dA0_T==0)continue;

        c_T_C_ISS[i] = new TCanvas("c_T_C"+ISS[i].Aname,"c_T_C"+ISS[i].Aname  ,1000,500);
        ISS[i].gr_T_C->Draw("ap");

        //c_T_C_ISS[i]->Print(ISS[i].folder_path+"C_T.png");
        //c_T_C_ISS[i]->Print(FILEPATH+ ISS[i].Aname+"C_T.png");
    }
}
void Draw_each_T_C(IONSpecies*ISS)
{
    double time_divide = 900000;
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].gr_T_C_n = 0;
        ISS[i].gr_T_C  = new TGraph();
        AxisFormat(ISS[i].gr_T_C,ISS[i].Aname+"  T_C "," T[ns] ", "C[m]");
        ISS[i].gr_T_C->SetMarkerSize(1.0);
        ISS[i].gr_T_C->SetLineWidth(1.0);
    }

    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].time>time_divide ){continue;} //divide by time-- T is changed in the experiment
        //dA0 -- ps
        //if(ions[j].Do_dA0_T_flag==1)
        {
            ISS[ions[j].Species].gr_T_C->SetPoint(ISS[ions[j].Species].gr_T_C_n , ions[j].T, ions[j].C);
            //ISS[ions[j].Species].gr_T_C->SetPointError(ISS[ions[j].Species].gr_T_C_n , 0, sqrt(ions[j].dA0err)*0.001);
            ISS[ions[j].Species].gr_T_C_n++;
        }
    }
    TCanvas*c_T_C_ISS[30];
    for(int i=0;i<NSpecies  ;i++)
    {
        c_T_C_ISS[i] = new TCanvas("c_T_C"+ISS[i].Aname,"c_T_C"+ISS[i].Aname  ,1000,500);
        ISS[i].gr_T_C->DrawClone("ap");
        delete ISS[i].gr_T_C;
        //c_T_C_ISS[i]->Print(ISS[i].folder_path+"C_T.png");
        //c_T_C_ISS[i]->Print(FILEPATH+ ISS[i].Aname+"C_T.png");
    }
}
void Show_gr2d_LddtC(string filename)
{

    ifstream infile;
    infile.open(filename);
    if(!infile){cerr<<"infile failed"<<endl;exit(1);}
    TGraph2D* gr2d_L_ddt_C = new TGraph2D();
    TGraph* gr_C_L = new TGraph();
    TGraph* gr_C_ddt = new TGraph();
    AxisFormat(gr2d_L_ddt_C,"","L [m]","ddT [ns]","C [m]");

    AxisFormat(gr_C_L,"","L [m]","C [m]");
    gr_C_L->SetMarkerSize(2);
    AxisFormat(gr_C_ddt,"","ddT [ns]","C [m]");
    gr_C_ddt->SetMarkerSize(2);
    int n=0;
    double L,ddt,C;
    while(infile>>L>>ddt>>C)
    {
         gr2d_L_ddt_C->SetPoint(n,L,ddt,C);
         gr_C_L->SetPoint(n,L,C);
         gr_C_ddt->SetPoint(n,ddt,C);
         n++;
    }
    infile.close();
    TCanvas* c_Show_gr2d_L_ddt_C = new TCanvas("c_Show_gr2d_L_ddt_C","c_Show_gr2d_L_ddt_C",1200,600); 
    gr2d_L_ddt_C->DrawClone("pcol");
    
    TCanvas* c_Show_gr_C_L = new TCanvas("c_Show_gr_C_L","c_Show_gr_C_L",1200,600);
    gr_C_L->Draw("Ap"); 
    TCanvas* c_Show_C_ddt = new TCanvas("c_Show_C_ddt","c_Show_C_ddt",1200,600); 
    gr_C_ddt->Draw("Ap");
    
}

//=========== 20230908 Tfix =============
double Get_Tfix1(ION& ions_i,double C0,double gt0)
{
    double Ti = ions_i.T;
    double Ci = ions_i.C ;
    double gamma = ions_i.Calculate_only_return_gamma();
    return Ti+(1- pow(gt0,2)/pow(gamma,2)) * (C0-Ci)/Ci*Ti;
}
double Get_Tfix2(ION& ions_i,double C0,TGraph* gr_gtC)
{
    double dC=0.001;
    double Ti = ions_i.T;
    double Ci = ions_i.C ;
    double gamma = ions_i.Calculate_only_return_gamma();
    double C,T,gtC=0;
    int n_step = abs(C0-Ci)/dC+1;   
    dC = (C0-Ci)/abs(C0-Ci)*dC;//可正可负
    C = Ci;
    T = Ti;
    //cout<<"debug! n_step= "<<n_step<<endl;
    for(int i=0;i<n_step  ;i++)  
    {
        gtC = gr_gtC->Eval(C+dC/2);
        T=T+(1- pow(gtC,2)/pow(gamma,2)) * dC/C*T;
        C+=dC;
    }
    return T;
}
double Get_Tfix3(ION& ions_i,double T_ave_in)
{
    double Ti = ions_i.T;
    double mi = ions_i.mvq_v1*ions_i.Z*u;  // mass result v1
    //double T_ave = ionspecies[ions_i.Species].AveT ;
    double T_ave = T_ave_in ;
    double m_ave = ionspecies[ions_i.Species].Mass_cal ;

    double gamma = ions_i.Calculate_only_return_gamma();
    return T_ave+(1/ pow (gamma,2) ) * (mi-m_ave)/m_ave*T_ave;
}
//======== 20250210 
double Get_m_from_T0_v1(ION& ions_i,double T_ave_in)
{
    //从原始周期按照 deltaT/T = (1/gamma^2) * delta m/m 得到 m0 分散很大
    double Ti = ions_i.T;
    
    //double T_ave = ionspecies[ions_i.Species].AveT ;
    double T_ave = T_ave_in ;
    double m_ave = ionspecies[ions_i.Species].Mass_cal ;

    double gamma = ions_i.Calculate_only_return_gamma();
    return m_ave*(1+ gamma*gamma*(Ti-T_ave)/T_ave);
}
double Get_m_from_T0_v2(ION& ions_i, double Bp_fix)
{
    //!!! 如何从国际单位的数据得到 质量 (实则为能量 单位 keV)
    ///// E[keV] = 1000000* Z * Bp[Tm] * Vc[m/ns]*Vc * sqrt( 1/(v[m/ns]*v) - 1/(Vc*Vc) )

    double v_ns=ions_i.v;
    //cout<<fixed<<setprecision(10)<<"v="<<v<<endl<<"v_ns="<<v_ns<<endl;
    double ss=1.0/v_ns/v_ns - 1/V_c/V_c;
    
    double AA = ions_i.Z*V_c*V_c*1000000.0; 
    
    return  AA*Bp_fix*sqrt(ss);
    
}
double Get_m_from_T0_v3(ION& ions_i,double Bp_fix,double C_fix)
{
    //!!! 如何从国际单位的数据得到 质量 (实则为能量 单位 keV)
    ///// E[keV] = 1000000* Z * Bp[Tm] * Vc[m/ns]*Vc * sqrt( 1/(v[m/ns]*v) - 1/(Vc*Vc) )

    double v_ns=C_fix/ions_i.T;
    //cout<<fixed<<setprecision(10)<<"v="<<v<<endl<<"v_ns="<<v_ns<<endl;
    double ss=1.0/v_ns/v_ns - 1/V_c/V_c;
    
    double AA = ions_i.Z*V_c*V_c*1000000.0; 
    
    return  AA*Bp_fix*sqrt(ss);
    
}
double Get_Tfix_from_m_v3(ION& ions_i,double Bp_fix,double C_fix)
{
    //!!! 如何从国际单位的数据得到 质量 (实则为能量 单位 keV)
    ///// E[keV] = 1000000* Z * Bp[Tm] * Vc[m/ns]*Vc * sqrt( 1/(v[m/ns]*v) - 1/(Vc*Vc) )

    double NucMass_tmp = ions_i.mvq_v1*ions_i.Z*u;
    double v_ns_tmp=0.3;
    double AA = ions_i.Z*V_c*V_c*1000000.0; 
    double ss = pow(NucMass_tmp/AA/Bp_fix,2);
    v_ns_tmp = sqrt( 1 / (ss +1/V_c/V_c) );
    return  C_fix/v_ns_tmp;
    
}

void Build_gr_TC0(IONSpecies* ISS)
{
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].gr_TC0_n = 0;
        ISS[i].gr_TC0  = new TGraph();
        AxisFormat(ISS[i].gr_TC0,ISS[i].Aname+"_TC0 ","C[m]", "T[ns]");
        ISS[i].gr_TC0->SetMarkerSize(1.0);
        ISS[i].gr_TC0->SetLineWidth(1.0);
    }
    for(int j=0;j<ions_n  ;j++)
    {
        //dA0 -- ps
        if(ions[j].C>C_filter_min&&ions[j].C<C_filter_max)
        //if( (ions[j].inject_number>0&&ions[j].inject_number<2000)||(ions[j].inject_number>23175&&ions[j].inject_number<25175)||SET2_PART2 )     
        //if( (ions[j].time>time_divide_1&&ions[j].time<time_divide_2)||SET2_PART2 )     
        //20230912  按照 SET2 all data, 前后各取2000 次注入               
        {
            ISS[ions[j].Species].gr_TC0->SetPoint(ISS[ions[j].Species].gr_TC0_n ,  ions[j].C,ions[j].T);
            ISS[ions[j].Species].gr_TC0_n++;
            ISS[ions[j].Species].h_Tfix0->Fill(ions[j].T);
        }
    } 
}
void Build_gr_TC1(IONSpecies* ISS)
{
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].gr_TC1_n = 0;
        ISS[i].gr_TC1  = new TGraph();
        AxisFormat(ISS[i].gr_TC1,ISS[i].Aname+"_TC1 ","C[m]", "T[ns]",kRed);
        ISS[i].gr_TC1->SetMarkerSize(1.0);
        ISS[i].gr_TC1->SetLineWidth(1.0);
    }
    for(int j=0;j<ions_n  ;j++)
    {                      
        if(ions[j].C>C_filter_min&&ions[j].C<C_filter_max)
        {
        //if((ions[j].inject_number>0&&ions[j].inject_number<2000)||(ions[j].inject_number>23175&&ions[j].inject_number<25175)||SET2_PART2 )
        //if( (ions[j].time>time_divide_1&&ions[j].time<time_divide_2)||SET2_PART2 )     
        {
            ISS[ions[j].Species].gr_TC1->SetPoint(ISS[ions[j].Species].gr_TC1_n ,  ions[j].C,Get_Tfix1(ions[j],128.808,1.375) );
            ISS[ions[j].Species].gr_TC1_n++;
            ISS[ions[j].Species].h_Tfix1->Fill(Get_Tfix1(ions[j],128.808,1.375));
        }
        }
    } 
}
void Build_gr_TC2(IONSpecies* ISS,TGraph* gr_gtC_chosen)
{
    double Tfix2=0;
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].gr_TC2_n = 0;
        ISS[i].gr_TC2  = new TGraph();
        AxisFormat(ISS[i].gr_TC2,ISS[i].Aname+"_TC2 ","C[m]", "T[ns]",kRed);
        ISS[i].gr_TC2->SetMarkerSize(1.0);
        ISS[i].gr_TC2->SetLineWidth(1.0);
    }
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].C>C_filter_min&&ions[j].C<C_filter_max)
        {
        //if((ions[j].inject_number>0&&ions[j].inject_number<2000)||(ions[j].inject_number>23175&&ions[j].inject_number<25175)||SET2_PART2 )
        //if( (ions[j].time>time_divide_1&&ions[j].time<time_divide_2)||SET2_PART2 )     
        {
            Tfix2 = Get_Tfix2(ions[j],128.808,gr_gtC_chosen);
            ISS[ions[j].Species].gr_TC2->SetPoint(ISS[ions[j].Species].gr_TC2_n ,  ions[j].C, Tfix2 );
            ISS[ions[j].Species].gr_TC2_n++;
            ISS[ions[j].Species].h_Tfix2->Fill(Tfix2);
        }
        }
    } 
}
void Build_gr_TC3(IONSpecies* ISS)
{   //依靠质量结果进行Tfix3修正
    double Tfix3=0;
    double Tfix3_v3=0;
    double T_ave=0;
    for(int i=0;i<NSpecies  ;i++)
    {
        ISS[i].gr_TC3_n = 0;
        ISS[i].gr_TC3  = new TGraph();
        AxisFormat(ISS[i].gr_TC3,ISS[i].Aname+"_TC3 ","C[m]", "T[ns]",kRed);
        ISS[i].gr_TC3->SetMarkerSize(1.0);
        ISS[i].gr_TC3->SetLineWidth(1.0);
    }
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].C>C_filter_min&&ions[j].C<C_filter_max)
        {
        //if((ions[j].inject_number>0&&ions[j].inject_number<2000)||(ions[j].inject_number>23175&&ions[j].inject_number<25175)||SET2_PART2 )
        //if( (ions[j].time>time_divide_1&&ions[j].time<time_divide_2)||SET2_PART2 )     
        {
            //不是所有的离子都有质量结果
            T_ave = ISS[ions[j].Species].h_Tfix0->GetMean(); // 需要先 Build_gr_TC0
            if(ions[j].mvq_v1 >0)
            {
                Tfix3 = Get_Tfix3(ions[j],T_ave);
                Tfix3_v3 = Get_Tfix_from_m_v3(ions[j],4.82,128.808); //随便一个Bpfix  Cfix 反正只看标准差不关心中心值
                ISS[ions[j].Species].gr_TC3->SetPoint(ISS[ions[j].Species].gr_TC3_n ,  ions[j].C, Tfix3 );
                ISS[ions[j].Species].gr_TC3_n++;
                ISS[ions[j].Species].h_Tfix3->Fill(Tfix3);

                ISS[ions[j].Species].vector_Tfix3_v3.push_back(Tfix3_v3);
            }
        }
        }
    } 

}
//==== 20250210
void Build_h_m0(IONSpecies* ISS)
{   //依靠原始T 来得到一个 初始 m0  
    double m0_v1=0;   //利用 deltaT/T = (1/gamma^2) * delta m/m
    double m0_v2=0; // 固定 Bp v2不靠谱,不使用
    double m0_v3=0;//固定 Bp 和C  --- v3和v1 是等价的
    double T_ave=0;

    double Bp_fix =0; double C_fix =0;
    if(THIS_EXP=="2017_58Ni"){Bp_fix=5.4758;C_fix=128.86;}
    else {Bp_fix=4.82;C_fix=128.808;}
    //else{cout<<"error!! in Build_h_m0 NO THIS_EXP"<<endl<<endl<<endl; return;}
    
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].C>C_filter_min&&ions[j].C<C_filter_max)
        //不是所有的离子都有质量结果
        if(ions[j].mvq_v1>0)  // 必须有最终质量，这样对比的是这些有最终质量的核，前后 m0 和 m1 的分散
        {
            T_ave = ISS[ions[j].Species].h_Tfix0->GetMean();// 需要先 Build_gr_TC0
            m0_v1 = Get_m_from_T0_v1(ions[j],T_ave);    
            ISS[ions[j].Species].h_m0_v1->Fill(m0_v1);
            //m0_v2 = Get_m_from_T0_v2(ions[j],Bp_fix);    
            //ISS[ions[j].Species].h_m0_v2->Fill(m0_v2);

            m0_v3 = Get_m_from_T0_v3(ions[j],Bp_fix,C_fix);    
            ISS[ions[j].Species].h_m0_v3->Fill(m0_v3);  //不一定落在直方图范围
            ISS[ions[j].Species].vector_m0_v3.push_back(m0_v3);

            ISS[ions[j].Species].h_m1->Fill(ions[j].mvq_v1*ions[j].Z*u);
            

            ISS[ions[j].Species].m0_n++;
        }
    } 

}


//此次注入是在实验开始(2021.10.24)后多久
double filename_to_time_second(string filename)
{
    int start_date = 0;
    int start_hour = 0;
    int pos = filename.find("_");
    string str = filename;// 将注入文件的名字拆成 年 月 日 时 分 秒
    string year_str = str.substr(pos-8,4);
    string date_str = str.substr(pos-4,4);
    string hour_str = str.substr(pos+1,2);
    string minute_str = str.substr(pos+3,2);
    string second_str = str.substr(pos+5,2);
    int year = atoi( (str.substr(pos-8,4)).c_str() );
    int date = atoi( date_str.c_str() );
    int hour = atoi( hour_str.c_str() );
    int minute = atoi( minute_str.c_str() );
    int second = atoi( second_str.c_str() );
    //cout<<str.substr(pos-8,4)<<endl<<str.substr(pos-4,4)<<endl<<hour_str<<endl<<minute_str<<endl<<second_str<<endl;
    //cout<<"int :"<<endl<<year<<endl<<date<<endl;

    //2021 36Ar 以 1024 0点为起点  1101-> 1032 针对这次实验的日期！！ 
    if (THIS_EXP == "2021_36Ar_SET1")//SET1 PART1 从 1019 23点左右开始
    {
        if (date > 1031)date = date - 69;  //@EXP@
        start_date = 1019;
        start_hour = 23;
    }
    if(THIS_EXP=="2021_36Ar_SET2")
    {
        if(date>1031)date = date-69;  //@EXP@ 1024
        start_date = 1024;
        start_hour = 0;
    }
    if(THIS_EXP=="2021_36Ar_SET3")  // 1106 08 ~
    {
        
        start_date = 1106;
        start_hour = 8;
    }
    //2017 58 Ni 这里以 1227 20点为起点  1227_2000~~~1231_0859 
    if(THIS_EXP=="2017_58Ni")
    {
        start_date = 1227;
        start_hour = 20;
    }
    
    double time  = (date-start_date)*24*3600+( (hour-start_hour)*3600.0+minute*60+second);
    return time;

}
TString which_slash(int n=1)
{
    TString result ="";

    if(THIS_ENVIRONMENT==ON_WINDOWS){for(int i=0;i<n;i++)result+= "\\";}
    else if(THIS_ENVIRONMENT==ON_UBUNTU){for(int i=0;i<n;i++)result+= "/";}

    return result;
}
//给每种识别的核素都创建一个对应的文件夹
void CreateFolder(IONSpecies* ionspecies, int N,TString makedir_ionspecies)
{
    TString makedir_each;
    for(int i=0;i<N;i++)
    {
        if(THIS_ENVIRONMENT==ON_WINDOWS)makedir_each = makedir_ionspecies+"\\\\"+ionspecies[i].Aname;
        else if(THIS_ENVIRONMENT==ON_UBUNTU)makedir_each = makedir_ionspecies+"//"+ionspecies[i].Aname;

        system(makedir_each);
        ionspecies[i].folder_path=FILEPATH_ionspecies+ionspecies[i].Aname+which_slash(2);
    }
}


    

//======================= show mass ===================================

//计算卡方 注：这里返回的是 Xn 而不是 Xn^2
double ChiSquare_sqrt(TGraphErrors* grerr )
{
    int n = grerr->GetN();   
    if(n<1){return 0;}
    double tmpx,tmpy;
    double y,yerr;
    double sum=0;
    for(int i=0;i<n;i++)
    {
       grerr->GetPoint(i,tmpx,y);
       yerr=grerr->GetErrorY(i);
       sum+=pow(y/yerr,2);
       //cout<<"##### "<<y<<" "<<yerr<<" "<<sum<<endl;
    }
    //return sqrt(sum/ (n-1) );// 0917 更改为n-1
    return sqrt(sum/ (n) );// 
}
int GetTzColor(double Tz)
{
    if     (Tz==-3.0)  {return kRed        ;}
    else if(Tz==-2.5)  {return kViolet     ;}
    else if(Tz==-2.0)  {return kOrange-3    ;}
    else if(Tz==-1.5)  {return kTeal+4     ;}
    else if(Tz==-1.0)  {return kAzure+7  ;}
    else if(Tz==-0.5)  {return kRed-2      ;}
    else {cout<<"\033[7m\033[31m in GetTzColor()  \033[0m "<<" Tz ="<<Tz <<"\033[7m\033[31m out of range !! \033[0m "<<endl;return 1;}
}
int GetTzColor_v2(double Tz)
{
    if     (Tz==-3.0)  {return kRed        ;}
    else if(Tz==-2.5)  {return kViolet     ;}
    else if(Tz==-2.0)  {return kOrange-3    ;}
    else if(Tz==-1.5)  {return kAzure     ;}
    else if(Tz==-1.0)  {return kRed-2  ;}
    else if(Tz==-0.5)  {return kBlack      ;}
    else {cout<<"\033[7m\033[31m in GetTzColor()  \033[0m "<<" Tz ="<<Tz <<"\033[7m\033[31m out of range !! \033[0m "<<endl;return 1;}
}
//只画某个Tz的核的质量结果
void Build_grerr_mass_Tz(IONSpecies* ISS,TGraphErrors* grerr,int opt,int plot_min,double Tz0,bool PLOTMIN_control_OFF,bool known_unknown_control_off,bool known_unknown_switch)
{
    int grerr_n=0;
    for(int i=0;i<NSpecies;i++)
    {
    //未知和已知的阴影都画出
    //ISS[i].PrintInfo(i);
        // 筛选条件生成 grerr AME shadow ，
        if(!ISS[i].HasResult)continue;

        // grerr1 所有核 限制N
        if( 
        ( ((ISS[i].IsRef||!known_unknown_switch)&&(!ISS[i].IsRef||known_unknown_switch) )||known_unknown_control_off) 
        && ((ISS[i].N_unknown>plot_min)||PLOTMIN_control_OFF) 
        ) 
        {
            if(ISS[i].Tz==Tz0)
            {
                if(opt==0)
                {
                    grerr->SetPoint(grerr_n, ISS[i].AveT ,ISS[i].deltaMass );
                    grerr->SetPointError(grerr_n,0,ISS[i].Mass_cal_err);
                    grerr_n++;
                }
                else if(opt==1)
                {
                    grerr->SetPoint(grerr_n, ISS[i].AveT ,ISS[i].deltaMass_v2 );
                    grerr->SetPointError(grerr_n,0,ISS[i].Mass_cal_err_v2);
                    grerr_n++;
                }
                
                else if(opt==2)
                {
                    grerr->SetPoint(grerr_n, ISS[i].AveT ,ISS[i].deltaMass_VE );
                    grerr->SetPointError(grerr_n,0,ISS[i].Mass_cal_err_VE);
                    grerr_n++;
                }
            }
        }
    }
    grerr->SetMarkerColor(GetTzColor(Tz0));
}

//========================== ISS own gtC ===========================================
//ionspecies[ts].grerr_avegtC ionspecies[ts].gr_gtC_own
void BuildAvegtCCurve_each(const TH1F* h1, double gt_err_upper_limit )
{
    int C_Division_n[200];
    double avegt   [200];
    double avegt2  [200];
    double sigma_gt[200];
    for(int i=0;i<subregion_n;i++){C_Division_n[i]=0;avegt[i]=0;avegt2[i]=0;sigma_gt[i]=0;}

    if(subregion_n>200){cout<<" error!!! in  BuildAvegtCCurve_each() subregion_n > 200"<<endl;return;}
    
    for(int ts=0;ts<NSpecies  ;ts++) // 每一个种类 逐个构建 own gtC
    {
        if(ionspecies[ts].N<50){continue;}  // 这种离子的计数太少不能生成自己的gtC curve
        //===RESET===
        for(int i=0;i<subregion_n;i++){C_Division_n[i]=0;avegt[i]=0;avegt2[i]=0;sigma_gt[i]=0;}
        ionspecies[ts].count_avegtC_n = 0;
        //----  fill ions.gt of this species
        for(int j=0;j<ions_n  ;j++)
        {
            if(ions[j].Species==ionspecies[ts].Species)
            {
                if(ions[j].gammat_err<=gt_err_upper_limit)
                {
                    if(ions[j].C_region>=0&&ions[j].C_region<=subregion_n-1)
                    {
                        C_Division_n[ions[j].C_region]++;
                        avegt[ions[j].C_region]     += ions[j].gammat;
                        avegt2[ions[j].C_region]    += ions[j].gammat*ions[j].gammat;
                        ionspecies[ts].count_avegtC_n ++;
                    }
                }
            }
        }
        // build gt arrays, own gtC curve
        for(int i=0;i<subregion_n;i++)
        {
            if(C_Division_n[i]>3)
            {
                avegt[i] /= C_Division_n[i];
                avegt2[i] /= C_Division_n[i];
                sigma_gt[i] = sqrt (avegt2[i]-avegt[i]*avegt[i]) ;
                ionspecies[ts].grerr_avegtC->SetPoint( ionspecies[ts].grerr_avegtC->GetN(),h1->GetBinCenter(i+1), avegt[i] );
                ionspecies[ts].grerr_avegtC->SetPointError(ionspecies[ts].grerr_avegtC->GetN()-1, 0, sigma_gt[i]/sqrt( C_Division_n[i]-1) ) ;
                
                ionspecies[ts].gr_gtC_own->SetPoint( ionspecies[ts].gr_gtC_own->GetN(),h1->GetBinCenter(i+1), avegt[i] );  // 20240702
                //cout<<avegt[i]<<" "<<avegt2[i] <<" "<<sigma_gt[i]<<endl;
            }    
        }
        Grerr_sigma_to_gr(ionspecies[ts].grerr_avegtC,ionspecies[ts].gr_gtC_own_u,1);
        Grerr_sigma_to_gr(ionspecies[ts].grerr_avegtC,ionspecies[ts].gr_gtC_own_d,2);
        ionspecies[ts].Has_own_avegtC = true;
        
    }
    
    
}

//============================================== Show_Mass_Result ==================================

//显示质量结果
void Show_Mass_Result(IONSpecies* ISS,int NSpecies,int ID,TString thispara,bool DRAW_ON,bool to_save,double Save_YAxis_LOW,double Save_YAxis_UP, bool err_sys_ON, double err_sys,int* opt,int optn,
    bool PLOTMIN_control_OFF,bool known_unknown_control_off,bool known_unknown_switch)
{
    //version: v1 v2 VE
    //这个函数首先制作不同质量版本的grerr，计算卡方，可输出文件 。 同时，在DRAW_ON 的情况下进行结果画图，可保存

    TLatex* lat_n_v1 = new TLatex();    TString lat_text_v1;   int COLOR_v1 = kBlack;
    TLatex* lat_n_v2 = new TLatex(); TString lat_text_v2;int COLOR_v2 = kAzure-7;
    TLatex* lat_n_VE = new TLatex(); TString lat_text_VE;int COLOR_VE = kAzure;
    TString str_TITLE="";
    
    int  SHADOW_VERSION      =    2  ; //AME 阴影样式 =1： 连接误差棒端点  =2： 矩形盒子
    bool COUT_c1_grerr1      =    0  ; // 输出每个点的值
    bool ShowCount           =    1  ; //在每个点的旁边显示计数
    bool c1_ShowNucleiInfo   =    1  ; //显示信息
    int c1_point_color_opt   =    2  ; //20231211 单独建立TGraph  =1点按照Tz分颜色, =2, RTO: R= reference/ T= Target/ O=others
        bool Only_one_Tz_ON  =    0  ;  //只显示某一种Tz 系列
        double Only_this_Tz    =  -1 ;  //显示哪个Tz 系列
    bool c1_Show_detail_info =    1  ; //显示详细信息
    if(LOOP_ON)c1_Show_detail_info =    false  ;

    if(known_unknown_control_off)c1_point_color_opt   =    2  ; //RTO 只适用于 所有核都画出来
    else c1_point_color_opt   =    1  ;  // 如果只画 REF / NONEREF 则根据Tz color
    
    
    int plot_min =5;                       //用于 shadow addAME   最小计数限制
    //bool PLOTMIN_control_OFF        = 0;  // =1 OFF 不对最小计数筛选      =0 ON 筛选最小计数
    //bool known_unknown_control_off  = 0;  // =1 control off 输出所有      =0 已知或未知筛选
    //bool known_unknown_switch       = 1;  // =1 只输出已知 不输出未知  ， =0  只输出未知 不输出已知
                                        // 在known_unknown_control_off=0 的前提下， 
    // known_unknown_switch =1-- MassUnknow 输出已知 不输出未知   known_unknown_switch =0--- 输出未知 不输出已知
    //根据写法 ( ((!ionspecies[i].MassUnknown||!known_unknown_switch)&&(ionspecies[i].MassUnknown||known_unknown_switch) )||known_unknown_control_off) 
                                    //f = ( (!a||!b)&&(a||b) )||c,      c=1,f=1 //// b=0,f=a// b=1, f=!a

    TString XTitle= "Revolution T [ns]";                               //坐标轴标题
    TString YTitle= "m_{EXP} - m_{AME} [keV]";
    //-------------------- the basic graph at bottom ---------------------
    int grerr_AME_shadow_n = 0;
    TGraphErrors *grerr_AME_shadow = new TGraphErrors();               //阴影为其AME质量误差区间
    AxisFormat(grerr_AME_shadow,"AME2020 errors",XTitle,YTitle,1);
    grerr_AME_shadow->SetFillColor(kGray);
    grerr_AME_shadow->SetFillStyle(3001);                              //阴影填充风格，在此修改
    //------------------result--------------------------
    // ----------mass v1 
    TGraphErrors *grerr1=new TGraphErrors(); 
    int grerr1_n=0;
    AxisFormat(grerr1,"",XTitle,YTitle,COLOR_v1);
    grerr1->SetLineWidth(2);grerr1->SetMarkerSize(2.5);    
    // ----------mass v2  h1
    TGraphErrors *grerr2 =new TGraphErrors(); 
    int grerr2_n=0;
    AxisFormat(grerr2,"",XTitle,YTitle,COLOR_v2);
    grerr2->SetLineWidth(2);grerr2->SetMarkerSize(2.5);
    //-----------mass VE
    TGraphErrors *grerrE = new TGraphErrors(); 
    int grerrE_n =0;
    AxisFormat(grerrE, "", XTitle, YTitle, COLOR_VE);
    grerrE->SetLineWidth(2);grerrE->SetMarkerSize(2.5);

    //--------------------add AME err for chi calculation not for drawing------------------------
    TGraphErrors *grerr1_add_AMEERR=new TGraphErrors();    
    int grerr1_add_AMEERR_n=0;
    TGraphErrors *grerr2_add_AMEERR=new TGraphErrors();
    int grerr2_add_AMEERR_n=0;
    TGraphErrors *grerrE_add_AMEERR = new TGraphErrors();
    int grerrE_add_AMEERR_n=0;
    
    //---------------------results add errsys-----------------------
    TGraphErrors *grerr1_sys=new TGraphErrors(); 
    int grerr1_sys_n=0;
    AxisFormat(grerr1_sys,"",XTitle,YTitle,COLOR_v1);
    grerr1_sys->SetLineWidth(2);
    grerr1_sys->SetMarkerSize(2.5);
    TGraphErrors *grerr1_add_AMEERR_sys=new TGraphErrors();
    int grerr1_add_AMEERR_sys_n=0;

    TGraphErrors *grerr2_sys=new TGraphErrors(); 
    int grerr2_sys_n=0;
    AxisFormat(grerr2_sys,"",XTitle,YTitle,COLOR_v2);
    grerr2_sys->SetLineWidth(2);
    grerr2_sys->SetMarkerSize(2.5);
    TGraphErrors *grerr2_add_AMEERR_sys=new TGraphErrors();
    int grerr2_add_AMEERR_sys_n=0;

    TGraphErrors *grerrE_sys=new TGraphErrors(); 
    int grerrE_sys_n=0;
    AxisFormat(grerrE_sys,"",XTitle,YTitle,COLOR_VE);
    grerrE_sys->SetLineWidth(2);
    grerrE_sys->SetMarkerSize(2.5);
    TGraphErrors *grerrE_add_AMEERR_sys=new TGraphErrors();
    int grerrE_add_AMEERR_sys_n=0;



    //--------------------- grerr Tz -----------------------
    TGraphErrors **grerr1_Tz=new TGraphErrors*[All_Tz_n];
    TGraphErrors **grerr2_Tz=new TGraphErrors*[All_Tz_n];
    TGraphErrors **grerrE_Tz=new TGraphErrors*[All_Tz_n];
    
    for(int j=0;j<All_Tz_n;j++){ grerr1_Tz[j]=new TGraphErrors();AxisFormat(grerr1_Tz[j],"",XTitle,YTitle,GetTzColor(All_Tz[j]));grerr1_Tz[j]->SetLineWidth(2);grerr1_Tz[j]->SetMarkerSize(2.5);}
    for(int j=0;j<All_Tz_n;j++){ grerr2_Tz[j]=new TGraphErrors();AxisFormat(grerr2_Tz[j],"",XTitle,YTitle,GetTzColor(All_Tz[j]));grerr2_Tz[j]->SetLineWidth(2);grerr2_Tz[j]->SetMarkerSize(2.5);}
    for(int j=0;j<All_Tz_n;j++){ grerrE_Tz[j]=new TGraphErrors();AxisFormat(grerrE_Tz[j],"",XTitle,YTitle,GetTzColor(All_Tz[j]));grerrE_Tz[j]->SetLineWidth(2);grerrE_Tz[j]->SetMarkerSize(2.5);}
    
    
    //-----------------------------------------------------------------------------------------------------------
    for(int i=0;i<NSpecies;i++) 
    {
    //未知和已知的阴影都画出
    //ISS[i].PrintInfo(i);
        // 筛选条件生成 grerr AME shadow ，
        if(!ISS[i].HasResult)continue;//此核没有结果，跳过
        if(                                                      //正常情况
        ( ((ISS[i].IsRef||!known_unknown_switch)&&(!ISS[i].IsRef||known_unknown_switch) )||known_unknown_control_off) //第一个 (a+!b)(!a+b) = a!a + ab + !b!a + !bb,  a!a 就是事件a成立的同时，事件a不成立，当然是不可能的
        && ((ISS[i].N_unknown>plot_min)||PLOTMIN_control_OFF)                                                                 //所以a!a=0，!bb=0，原式 = ab + !a!b，含义为ab均成立或ab均不成立
        ) 
        { //条件解释 1.(ISS[i].N_unknown>plot_min)||PLOTMIN_control_OFF  首先粒子数足够，或者粒子数不做限制，满足则下一项
          //         2.known_unknown_control_off                 是否打开了输出所有核开关，没有打开则进行下一项
          //         3.((ISS[i].IsRef||!known_unknown_switch)&&(!ISS[i].IsRef||known_unknown_switch))  该核为参考核，输出已知核，或 该核不为参考核，输出未知核
          //    ISS[i].PrintInfo(i);     绘制此核的    AME质量误差
            grerr_AME_shadow->SetPoint(grerr_AME_shadow_n,ISS[i].AveT,0);
            grerr_AME_shadow->SetPointError(grerr_AME_shadow_n,0.5,ISS[i].AME_err);
            grerr_AME_shadow_n++;
        }

        // grerr1 所有核 限制N
        if( 
        ( ((ISS[i].IsRef||!known_unknown_switch)&&(!ISS[i].IsRef||known_unknown_switch) )||known_unknown_control_off) 
        && ((ISS[i].N_unknown>plot_min)||PLOTMIN_control_OFF) 
        ) 
        {
            grerr1->SetPoint(grerr1_n, ISS[i].AveT ,ISS[i].deltaMass );  //绘制此核的  质量测量值与AME的差值
            grerr1->SetPointError(grerr1_n,0,ISS[i].Mass_cal_err);
            grerr1_n++;
            grerr1_sys->SetPoint(grerr1_sys_n, ISS[i].AveT ,ISS[i].deltaMass );
            grerr1_sys->SetPointError(grerr1_sys_n,0,sqrt(pow(ISS[i].deltaMass_err,2)+pow(1*err_sys,2)) );//err_sys系统误差，一般为0
            grerr1_sys_n++;
            //------- SET Tz grerr
            for(int j=0;j<All_Tz_n;j++)
                {Build_grerr_mass_Tz(ISS,grerr1_Tz[j],0,plot_min,All_Tz[j],PLOTMIN_control_OFF,known_unknown_control_off,known_unknown_switch);}
            
            //v2
            if(opt[1]) //fit 
            {
            grerr2->SetPoint(grerr2_n, ISS[i].AveT ,ISS[i].deltaMass_v2 );
            grerr2->SetPointError(grerr2_n,0,ISS[i].Mass_cal_err_v2);
            grerr2_n++;
            grerr2_sys->SetPoint(grerr2_n, ISS[i].AveT ,ISS[i].deltaMass_v2 );
            grerr2_sys->SetPointError(grerr2_n, 0,sqrt(pow(ISS[i].deltaMass_err_v2,2)+pow(ISS[i].Z*err_sys,2)) );
            grerr2_sys_n++;
            //------- SET Tz grerr
            for(int j=0;j<All_Tz_n;j++)
                {Build_grerr_mass_Tz(ISS,grerr2_Tz[j],1,plot_min,All_Tz[j],PLOTMIN_control_OFF,known_unknown_control_off,known_unknown_switch);}

            }
            //ve
            if(opt[2])
            {
            grerrE->SetPoint(grerrE_n, ISS[i].AveT ,ISS[i].deltaMass_VE );
            grerrE->SetPointError(grerrE_n,0,ISS[i].Mass_cal_err_VE);
            grerrE_n++;
            grerrE_sys->SetPoint(grerrE_n, ISS[i].AveT ,ISS[i].deltaMass_VE );
            grerrE_sys->SetPointError(grerrE_n, 0,sqrt(pow(ISS[i].deltaMass_err_VE,2)+pow(ISS[i].Z*err_sys,2)) );
            grerrE_sys_n++;
            for(int j=0;j<All_Tz_n;j++)
                {Build_grerr_mass_Tz(ISS,grerrE_Tz[j],2,plot_min,All_Tz[j],PLOTMIN_control_OFF,known_unknown_control_off,known_unknown_switch);}
            }
            
            
        }
        //卡方计算一定使用已知核
        if( ISS[i].IsRef && ((ISS[i].N_unknown>plot_min)||PLOTMIN_control_OFF) ) 
        {
            // v1 always on
            {
                grerr1_add_AMEERR->SetPoint(grerr1_add_AMEERR_n, ISS[i].AveT ,ISS[i].deltaMass ); 
                grerr1_add_AMEERR->SetPointError(grerr1_add_AMEERR_n,0,    ISS[i].deltaMass_err   );
                grerr1_add_AMEERR_n++;
                grerr1_add_AMEERR_sys->SetPoint(grerr1_add_AMEERR_sys_n, ISS[i].AveT ,ISS[i].deltaMass);
                grerr1_add_AMEERR_sys->SetPointError(grerr1_add_AMEERR_sys_n,0,sqrt ( pow (ISS[i].deltaMass_err,2)  + pow (1*err_sys,2)) );//误差项  (σ_cal^2 + σ_AME^2)^0.5
                grerr1_add_AMEERR_sys_n++;
            }
            // v2 fit
            if(opt[1])
            {
                grerr2_add_AMEERR->SetPoint(grerr2_add_AMEERR_n, ISS[i].AveT ,ISS[i].deltaMass_v2 ); 
                grerr2_add_AMEERR->SetPointError(grerr2_add_AMEERR_n,0,   ISS[i].deltaMass_err_v2 );
                grerr2_add_AMEERR_n++;
                grerr2_add_AMEERR_sys->SetPoint(grerr2_add_AMEERR_sys_n, ISS[i].AveT ,ISS[i].deltaMass_v2 ); 
                grerr2_add_AMEERR_sys->SetPointError(grerr2_add_AMEERR_sys_n,0,sqrt ( pow (ISS[i].deltaMass_err_v2,2) + pow (ISS[i].Z*err_sys,2)) );
                grerr2_add_AMEERR_sys_n++;
            }
            // VE
            if(opt[2])
            {
                grerrE_add_AMEERR->SetPoint(grerrE_add_AMEERR_n, ISS[i].AveT ,ISS[i].deltaMass_VE ); 
                grerrE_add_AMEERR->SetPointError(grerrE_add_AMEERR_n,0, ISS[i].deltaMass_err_VE   );
                grerrE_add_AMEERR_n++;
                grerrE_add_AMEERR_sys->SetPoint(grerrE_add_AMEERR_sys_n, ISS[i].AveT ,ISS[i].deltaMass_VE ); 
                grerrE_add_AMEERR_sys->SetPointError(grerrE_add_AMEERR_sys_n,0,sqrt ( pow (ISS[i].deltaMass_err_VE,2) + + pow (ISS[i].Z*err_sys,2)) );
                grerrE_add_AMEERR_sys_n++;
            }
        }
    }
    

    outfile_logs<<" ============ Show_Mass_Result() No. "<<ID<<" ========================"<<endl;
    outfile_logs<<"PLOTMIN_control_OFF = "<<PLOTMIN_control_OFF<<" plot_min = "<<plot_min
    <<" known_unknown_control_off: "<<known_unknown_control_off
    <<" known_unknown_switch:"<< known_unknown_switch<<"  // =1 only show known , =0  only show unknown "<<endl;

    outfile_logs<<"grerr1_n = "<<grerr1_n<<" grerr1_add_AMEERR_n= "<<grerr1_add_AMEERR_n<<endl;
    if(opt[1])outfile_logs<<"grerr2_n = "<<grerr2_n<<" grerr2_add_AMEERR_n= "<<grerr2_add_AMEERR_n<<endl;
    if(opt[2])outfile_logs<<"grerrE_n = "<<grerrE_n<<" grerrE_add_AMEERR_n= "<<grerrE_add_AMEERR_n<<endl;
    


    grerr_AME_shadow->Sort(&TGraph::CompareX);//T order 按照x值把点排序
    grerr1->Sort(&TGraph::CompareX);//T order
    grerr1_add_AMEERR->Sort(&TGraph::CompareX);//T order
    grerr1_add_AMEERR_sys->Sort(&TGraph::CompareX);//T order
    //v2
    if(opt[1])
    {   
    grerr2->Sort(&TGraph::CompareX);//T order
    grerr2_add_AMEERR->Sort(&TGraph::CompareX);//T order
    grerr2_add_AMEERR_sys->Sort(&TGraph::CompareX);//T order
    }
    //VE
    if(opt[2])
    {   
    grerrE->Sort(&TGraph::CompareX);//T order
    grerrE_add_AMEERR->Sort(&TGraph::CompareX);//T order
    grerrE_add_AMEERR_sys->Sort(&TGraph::CompareX);//T order
    }

    //Get Chisquare here
    if(!err_sys_ON)
    {
        //不加入系统误差
        L_ddT_Xn2[scan_loop_i] = ChiSquare_sqrt(grerr1_add_AMEERR);
        L_ddT_Xn2_v2[scan_loop_i] = ChiSquare_sqrt(grerr2_add_AMEERR);
        L_ddT_Xn2_VE[scan_loop_i] = ChiSquare_sqrt(grerrE_add_AMEERR);
    }
    else 
    {
        //加入系统误差
        L_ddT_Xn2_sys[scan_loop_i] = ChiSquare_sqrt(grerr1_add_AMEERR_sys);
        L_ddT_Xn2_v2_sys[scan_loop_i] = ChiSquare_sqrt(grerr2_add_AMEERR_sys);
        L_ddT_Xn2_VE_sys[scan_loop_i] = ChiSquare_sqrt(grerrE_add_AMEERR_sys);

        //构建 gr_errsys_chin : 卡方随着加入的系统误差 err_sys 的变化
        if(Only_draw_mass_VE)
        {   
            gr_errsys_chin->SetPoint(gr_errsys_chin_n++, err_sys,ChiSquare_sqrt(grerrE_add_AMEERR_sys) );
        }
        else 
        {   
            if(MASS_VER==1)gr_errsys_chin->SetPoint(gr_errsys_chin_n++, err_sys,ChiSquare_sqrt(grerr1_add_AMEERR_sys) );
            if(MASS_VER==2)gr_errsys_chin->SetPoint(gr_errsys_chin_n++, err_sys,ChiSquare_sqrt(grerr2_add_AMEERR_sys) );
        } 
    }
    
    //Show Chi_n in the title of grerr_AME_shadow---only when showing ref-ions since the chi is calculated only for ref ions
    // ============ title ========================
    str_TITLE = thispara+strtmp.Format(" #chi_{n} = ");
    if(known_unknown_control_off==0&&known_unknown_switch==1)
    {
        if(opt[0]){ if(err_sys_ON) str_TITLE+= strtmp.Format(" %f ",ChiSquare_sqrt(grerr1_add_AMEERR_sys)) ;
                    else           str_TITLE+= strtmp.Format(" %f ",ChiSquare_sqrt(grerr1_add_AMEERR));     }
        if(opt[1]){ if(err_sys_ON) str_TITLE+= strtmp.Format("| %f ",ChiSquare_sqrt(grerr2_add_AMEERR_sys)) ;
                    else           str_TITLE+= strtmp.Format("| %f ",ChiSquare_sqrt(grerr2_add_AMEERR));     }
        if(opt[2]){ if(err_sys_ON) str_TITLE+= strtmp.Format("| %f ",ChiSquare_sqrt(grerrE_add_AMEERR));
                    else           str_TITLE+= strtmp.Format("| %f ",ChiSquare_sqrt(grerrE_add_AMEERR_sys));     }
    }
    grerr_AME_shadow->SetTitle(str_TITLE);

    //================= outfile_Lddt_Xn2 卡方记录文件输出 ========================
    if(use_gtC_type==2)  // gt0
    {
        outfile_Lddt_Xn2<<L<<" "<<ddT<<" "<<gt0<<" "
                         <<L_ddT_Xn2[scan_loop_i]<<" "<<L_ddT_Xn2_v2[scan_loop_i]<<" "<<L_ddT_Xn2_VE[scan_loop_i]<<endl;
    }
    else if(use_gtC_type==1||use_gtC_type==3)
    {
        if(gtC_fit_ON)
        {
            outfile_Lddt_Xn2<<L<<" "<<ddT<<" "<<gtC_fit_info<<" "
            <<L_ddT_Xn2[scan_loop_i]<<" "<<L_ddT_Xn2_v2[scan_loop_i]<<" "<<L_ddT_Xn2_VE[scan_loop_i]<<endl;
        }
        else
        {   
            outfile_Lddt_Xn2<<L<<" "<<ddT<<" "
            <<L_ddT_Xn2[scan_loop_i]<<" "<<L_ddT_Xn2_v2[scan_loop_i]<<" "<<L_ddT_Xn2_VE[scan_loop_i]<<endl;
        }
    }
    else{ cout<<" wrong  use_gtC_type !! found in Show_Mass_Result()"<<endl;}


    //新增 根据 参考核 目标核 其它核 3分类的颜色显示 由于ROOT 一个TGraph 里面所有点只能是一种颜色，只能把不同颜色的点分别加到不同TGraph
    //-------------------- grerr 3 parts: r=reference t=targets(mass unknown) o=others(mass known but not ref)
    int RTO[MAX_IONSPECIES]={0}; // r = 1, t=2, o=3
    for(int i=0;i<NSpecies;i++)
    {
        if(ISS[i].MassUnknown==1)
        {    RTO[i]=2;}  // target
        else
        {
            if(ISS[i].IsRef){RTO[i]=1;} // ref
            else {RTO[i]=3;}    // others
        }
    }
    
    int COLOR_R = 1;int COLOR_T = kRed;int COLOR_O = kBlue;  // color of r=reference t=targets(mass unknown) o=others(mass known but not ref)
    
    TGraphErrors *grerr1_R=new TGraphErrors();
    TGraphErrors *grerr1_T=new TGraphErrors();
    TGraphErrors *grerr1_O=new TGraphErrors(); 
    AxisFormat(grerr1_R,"",XTitle,YTitle,COLOR_R);
    AxisFormat(grerr1_T,"",XTitle,YTitle,COLOR_T);
    AxisFormat(grerr1_O,"",XTitle,YTitle,COLOR_O);
    grerr1_R->SetLineWidth(2);grerr1_R->SetMarkerSize(2.5);grerr1_R->SetMarkerStyle(20);    
    grerr1_T->SetLineWidth(2);grerr1_T->SetMarkerSize(2.5);grerr1_T->SetMarkerStyle(22);
    grerr1_O->SetLineWidth(2);grerr1_O->SetMarkerSize(2.5);grerr1_O->SetMarkerStyle(21);
    // ----------mass v2  h1
    TGraphErrors *grerr2_R=new TGraphErrors();
    TGraphErrors *grerr2_T=new TGraphErrors();
    TGraphErrors *grerr2_O=new TGraphErrors(); 
    AxisFormat(grerr2_R,"",XTitle,YTitle,COLOR_R);
    AxisFormat(grerr2_T,"",XTitle,YTitle,COLOR_T);
    AxisFormat(grerr2_O,"",XTitle,YTitle,COLOR_O);
    grerr2_R->SetLineWidth(2);grerr2_R->SetMarkerSize(2.5);grerr2_R->SetMarkerStyle(20);    
    grerr2_T->SetLineWidth(2);grerr2_T->SetMarkerSize(2.5);grerr2_T->SetMarkerStyle(22);
    grerr2_O->SetLineWidth(2);grerr2_O->SetMarkerSize(2.5);grerr2_O->SetMarkerStyle(21);
    //-----------mass VE
    TGraphErrors *grerrE_R=new TGraphErrors();
    TGraphErrors *grerrE_T=new TGraphErrors();
    TGraphErrors *grerrE_O=new TGraphErrors(); 
    AxisFormat(grerrE_R,"",XTitle,YTitle,COLOR_R);
    AxisFormat(grerrE_T,"",XTitle,YTitle,COLOR_T);
    AxisFormat(grerrE_O,"",XTitle,YTitle,COLOR_O);
    grerrE_R->SetLineWidth(2);grerrE_R->SetMarkerSize(2.5);grerrE_R->SetMarkerStyle(20);    
    grerrE_T->SetLineWidth(2);grerrE_T->SetMarkerSize(2.5);grerrE_T->SetMarkerStyle(22);
    grerrE_O->SetLineWidth(2);grerrE_O->SetMarkerSize(2.5);grerrE_O->SetMarkerStyle(21);

    for(int i=0;i<NSpecies;i++) 
    {
        if(!ISS[i].HasResult)continue;//此核没有结果，跳过
        if(RTO[i]==1)  // R
        {
            grerr1_R->SetPoint     (grerr1_R->GetN(), ISS[i].AveT ,ISS[i].deltaMass );  
            grerr1_R->SetPointError(grerr1_R->GetN()-1,0,ISS[i].Mass_cal_err);
            grerr2_R->SetPoint     (grerr2_R->GetN(), ISS[i].AveT ,ISS[i].deltaMass_v2 );
            grerr2_R->SetPointError(grerr2_R->GetN()-1,0,ISS[i].Mass_cal_err_v2);
            grerrE_R->SetPoint     (grerrE_R->GetN(), ISS[i].AveT ,ISS[i].deltaMass_VE );
            grerrE_R->SetPointError(grerrE_R->GetN()-1,0,ISS[i].Mass_cal_err_VE);
        }
        else if(RTO[i]==2)  //T
        {
            grerr1_T->SetPoint     (grerr1_T->GetN(), ISS[i].AveT ,ISS[i].deltaMass );  
            grerr1_T->SetPointError(grerr1_T->GetN()-1,0,ISS[i].Mass_cal_err);
            grerr2_T->SetPoint     (grerr2_T->GetN(), ISS[i].AveT ,ISS[i].deltaMass_v2 );
            grerr2_T->SetPointError(grerr2_T->GetN()-1,0,ISS[i].Mass_cal_err_v2);
            grerrE_T->SetPoint     (grerrE_T->GetN(), ISS[i].AveT ,ISS[i].deltaMass_VE );
            grerrE_T->SetPointError(grerrE_T->GetN()-1,0,ISS[i].Mass_cal_err_VE);
        }
        else if(RTO[i]==3) //O
        {
            grerr1_O->SetPoint     (grerr1_O->GetN(), ISS[i].AveT ,ISS[i].deltaMass );  
            grerr1_O->SetPointError(grerr1_O->GetN()-1,0,ISS[i].Mass_cal_err);
            grerr2_O->SetPoint     (grerr2_O->GetN(), ISS[i].AveT ,ISS[i].deltaMass_v2 );
            grerr2_O->SetPointError(grerr2_O->GetN()-1,0,ISS[i].Mass_cal_err_v2);
            grerrE_O->SetPoint     (grerrE_O->GetN(), ISS[i].AveT ,ISS[i].deltaMass_VE );
            grerrE_O->SetPointError(grerrE_O->GetN()-1,0,ISS[i].Mass_cal_err_VE);
        }
    }

//-----------------------------------------------------------------------------------------------------------------------
if(DRAW_ON)
{


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auto c1_AMEcompare_T_shadow = new TCanvas(strtmp.Format("c1_AMEcompare_T_shadow_%d",ID),strtmp.Format("c1_AMEcompare_T_shadow_%d",ID),2000,1000);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//================================ Draw TGraphErrors ==================================
c1_AMEcompare_T_shadow->SetGrid(1);
//auto mg1 = new TMultiGraph();


//###################### AME shadow #############################
// drawing shadow needs points to be set in order
if(SHADOW_VERSION==1)grerr_AME_shadow->Draw("a3");   //"3"   A filled area is drawn through the end points of the vertical error bars.
//grerr_AME_shadow->Draw("a2");   //"2"   Error rectangles are drawn.
else if(SHADOW_VERSION==2)grerr_AME_shadow->Draw("a5");    //like option "2". In addition the contour line around the boxes is drawn.
else grerr_AME_shadow->Draw("a3");


// ------- point color option 1=Tz_color 
if(c1_point_color_opt==1)
{
    for(int j=0;j<All_Tz_n;j++)
    {
        // show only one Tz series
        if(All_Tz[j]==Only_this_Tz||!Only_one_Tz_ON){}
        else {continue;}

        if(opt[0]){grerr1_Tz[j]->Draw("samep"); }
        if(opt[1]){grerr2_Tz[j]->Draw("samep"); }
        if(opt[2]){grerrE_Tz[j]->Draw("samep"); }
    }
}
//------- point color option 2=ref/target/others
else if(c1_point_color_opt==2)
{
    if(opt[0]){grerr1_R->Draw("samep");grerr1_T->Draw("samep");grerr1_O->Draw("samep"); }
    if(opt[1]){grerr2_R->Draw("samep");grerr2_T->Draw("samep");grerr2_O->Draw("samep"); }
    if(opt[2]){grerrE_R->Draw("samep");grerrE_T->Draw("samep");grerrE_O->Draw("samep"); }
}
else
{
    if(opt[0])
    {   
        if(err_sys_ON) grerr1_sys->Draw("samep");
        else           grerr1->Draw("samep");
    }
    if(opt[1])
    {   
        if(err_sys_ON) grerr2_sys->Draw("samep");
        else           grerr2->Draw("samep");
    }
    if(opt[2])
    {   
        if(err_sys_ON) grerrE_sys->Draw("samep");
        else           grerrE->Draw("samep");
    }
}

grerr_AME_shadow->GetYaxis()->SetRangeUser(Save_YAxis_LOW,Save_YAxis_UP);

//________________________________ Draw TGraphErrors ______________________________

///----------------- info ------------------------

TLegend* legend1= new TLegend(0.15,0.75,0.3,0.9);
//legend1->SetHeader(" results compare to AME ","C"); // option "C" allows to center the header
    
legend1->AddEntry(grerr_AME_shadow,"AME2020 errors shadow","F");
if(err_sys_ON)legend1->AddEntry(grerr1_add_AMEERR,strtmp.Format("#sigma_{sys} = %.2f [keV]",err_sys),"P");




if(c1_point_color_opt==1)
{
    for(int j=0;j<All_Tz_n;j++)
    {legend1->AddEntry(grerr1_Tz[j],strtmp.Format("Tz = %.1f",All_Tz[j]),"PE");}
}
if(c1_point_color_opt==2)
{
    legend1->AddEntry(grerr1_R," reference ","PE");
    legend1->AddEntry(grerr1_T," target ","PE");
    legend1->AddEntry(grerr1_O," others ","PE");
}
else
{
    if(opt[0])legend1->AddEntry(grerr1,"mass result v1","PE"); //标签
    if(opt[1])legend1->AddEntry(grerr2,"mass v2 fit histograms result ","PE");
    if(opt[2])legend1->AddEntry(grerrE,"VE error weighted results ","PE");
}
legend1->Draw();

if(COUT_c1_grerr1)
{   
    cout<<"================ c1 grerr1 ============"<<endl;
    for(int i=0;i<NSpecies;i++){cout<<ISS[i].A<<" "<< ISS[i].name <<" "<<ISS[i].N_unknown <<" "<<ISS[i].deltaMass <<" "<<ISS[i].Mass_cal_err<<endl;}
}
//================= c1 show nuclei info ======================
lat_n_v1->SetTextColor(COLOR_v1);lat_n_v2->SetTextColor(COLOR_v2);lat_n_VE->SetTextColor(COLOR_VE);
lat_n_v1->SetTextFont(43);lat_n_v2->SetTextFont(43);lat_n_VE->SetTextFont(43);
lat_n_v1->SetTextSize(20);lat_n_v2->SetTextSize(20);lat_n_VE->SetTextSize(20);
lat_n_v1->SetTextAlign(12);lat_n_v2->SetTextAlign(12);  lat_n_VE->SetTextAlign(12);
double x_pos_shift = 0.001*ionspecies[0].AveT;  // for T~600 ns , 0.001*T~0.6ns
double y_pos_ratio = 5;
if(c1_ShowNucleiInfo)
{
    
    for(int i=0;i<NSpecies;i++)
    {    
        if( 
            ( ((ISS[i].IsRef||!known_unknown_switch)&&(!ISS[i].IsRef||known_unknown_switch) )||known_unknown_control_off) 
            && ((ISS[i].N_unknown>plot_min)||PLOTMIN_control_OFF) 
          ) 
          
        {
            // Tz info filter
            if(ISS[i].Tz==Only_this_Tz||!Only_one_Tz_ON){}
            else {continue;}
            
            if(ShowCount)
            {
                lat_text    = ISS[i].name_latex+ lat_text.Format(": %d",ISS[i].N_unknown) ;
                lat_text_v2 = ISS[i].name_latex+ lat_text.Format(": %d",ISS[i].N_unknown) ; 
                lat_text_VE = ISS[i].name_latex+ lat_text.Format(": %d",ISS[i].N_unknown) ; 
            }
            else 
            {
                lat_text    = ISS[i].name_latex ;       
                lat_text_v2 = ISS[i].name_latex ;         
                lat_text_VE = ISS[i].name_latex ;         
            }
    
            if(c1_Show_detail_info)
            {
                if(err_sys_ON)
                {
                    if(opt[0])lat_text    ="#splitline{"+lat_text+"}"   +"{"+strtmp.Format("%.1f(%.1f)",ISS[i].deltaMass,sqrt(pow(ISS[i].Mass_cal_err,2)+pow(1*err_sys,2)) ) + "}";
                    if(opt[1])lat_text_v2 ="#splitline{"+lat_text_v2+"}"+"{"+strtmp.Format("%.1f(%.1f)",ISS[i].deltaMass_v2,ISS[i].Mass_cal_err_v2+ISS[i].Z*err_sys) + "}";
                    if(opt[2])lat_text_VE ="#splitline{"+lat_text_VE+"}"+"{"+strtmp.Format("%.1f(%.1f)",ISS[i].deltaMass_VE,ISS[i].Mass_cal_err_VE+ISS[i].Z*err_sys) + "}";
                }
                else
                {
                    if(opt[0])lat_text    ="#splitline{"+lat_text+"}"   +"{"+strtmp.Format("%.1f(%.1f)",ISS[i].deltaMass,ISS[i].Mass_cal_err) + "}";
                    if(opt[1])lat_text_v2 ="#splitline{"+lat_text_v2+"}"+"{"+strtmp.Format("%.1f(%.1f)",ISS[i].deltaMass_v2,ISS[i].Mass_cal_err_v2) + "}";
                    if(opt[2])lat_text_VE ="#splitline{"+lat_text_VE+"}"+"{"+strtmp.Format("%.1f(%.1f)",ISS[i].deltaMass_VE,ISS[i].Mass_cal_err_VE) + "}";
                }
            }
            // text color
            if(c1_point_color_opt==1)
            {   
                lat_n_v1->SetTextColor(GetTzColor(ISS[i].Tz));lat_n_v2->SetTextColor(GetTzColor(ISS[i].Tz));lat_n_VE->SetTextColor(GetTzColor(ISS[i].Tz));
            }
            else if(c1_point_color_opt==2)
            {
                if(RTO[i]==1){lat_n_v1->SetTextColor(COLOR_R);lat_n_v2->SetTextColor(COLOR_R);lat_n_VE->SetTextColor(COLOR_R);}
                if(RTO[i]==2){lat_n_v1->SetTextColor(COLOR_T);lat_n_v2->SetTextColor(COLOR_T);lat_n_VE->SetTextColor(COLOR_T);}
                if(RTO[i]==3){lat_n_v1->SetTextColor(COLOR_O);lat_n_v2->SetTextColor(COLOR_O);lat_n_VE->SetTextColor(COLOR_O);}
            }
            else 
            {   
                lat_n_v1->SetTextColor(COLOR_v1);lat_n_v2->SetTextColor(COLOR_v2);lat_n_VE->SetTextColor(COLOR_VE);
            }
            
            //-------------------- lat_n -> Draw
            if(c1_Show_detail_info)
            {
                double xx,yy=0;  
                if(THIS_EXP=="2017_58Ni"&&ID==0&&!LOOP_ON)
                {
                    xx=ISS[i].AveT;yy=ISS[i].deltaMass+abs(ISS[i].deltaMass)/(ISS[i].deltaMass)*y_pos_ratio;
                    if(ISS[i].Aname=="25Si"){xx=ISS[i].AveT-x_pos_shift;yy=-15;}
                    if(ISS[i].Aname=="18Ne"){xx=ISS[i].AveT;yy=-15;}
                    if(ISS[i].Aname=="29S"){xx=ISS[i].AveT+x_pos_shift;yy=-15;}
                    if(ISS[i].Aname=="20Na"){xx=ISS[i].AveT;yy=10;}
                    if(ISS[i].Aname=="22Mg"){xx=ISS[i].AveT;yy=-10;}
                    if(ISS[i].Aname=="11C"){xx=ISS[i].AveT+2;yy=15;}
                    if(ISS[i].Aname=="26Si"){xx=ISS[i].AveT;yy=-10;}
                    if(ISS[i].Aname=="13N"){xx=ISS[i].AveT;yy=8;}
                    if(ISS[i].Aname=="28P"){xx=ISS[i].AveT;yy=10;}
                    if(ISS[i].Aname=="30S"){xx=ISS[i].AveT;yy=15;}
                    if(ISS[i].Aname=="15O"){xx=ISS[i].AveT;yy=-5;}
                    if(ISS[i].Aname=="32Cl"){xx=ISS[i].AveT;yy=10;}
                    if(ISS[i].Aname=="34Ar"){xx=ISS[i].AveT;yy=-5;}
                    if(ISS[i].Aname=="17F"){xx=ISS[i].AveT;yy=5;}
                    if(ISS[i].Aname=="53Ni"){xx=ISS[i].AveT+x_pos_shift;yy=20;}
                    if(ISS[i].Aname=="36K"){xx=ISS[i].AveT;yy=10;}
                    if(ISS[i].Aname=="38Ca"){xx=ISS[i].AveT;yy=25;}
                    if(ISS[i].Aname=="19Ne"){xx=ISS[i].AveT;yy=-5;}
                    if(ISS[i].Aname=="40Sc"){xx=ISS[i].AveT;yy=20;}
                    if(ISS[i].Aname=="42Ti"){xx=ISS[i].AveT;yy=-15;}
                    if(ISS[i].Aname=="21Na"){xx=ISS[i].AveT;yy=-10;}
                    if(ISS[i].Aname=="46Cr"){xx=ISS[i].AveT;yy=-25;}
                    if(ISS[i].Aname=="23Mg"){xx=ISS[i].AveT;yy=-20;}
                    if(ISS[i].Aname=="48Mn"){xx=ISS[i].AveT;yy=20;}
                    if(ISS[i].Aname=="50Fe"){xx=ISS[i].AveT;yy=-7;}
                    if(ISS[i].Aname=="25Al"){xx=ISS[i].AveT;yy=-11;}
                    if(ISS[i].Aname=="27Si"){xx=ISS[i].AveT;yy=15;}
                    if(ISS[i].Aname=="56Cu"){xx=ISS[i].AveT;yy=13;}
                    if(ISS[i].Aname=="29P"){xx=ISS[i].AveT;yy=27;}
                    if(ISS[i].Aname=="31S"){xx=ISS[i].AveT;yy=23;}
                    if(ISS[i].Aname=="33Cl"){xx=ISS[i].AveT;yy=20;}
                    if(ISS[i].Aname=="35Ar"){xx=ISS[i].AveT;yy=10;}
                    if(ISS[i].Aname=="37K"){xx=ISS[i].AveT;yy=-10;}
                    if(ISS[i].Aname=="39Ca"){xx=ISS[i].AveT;yy=15;}
                    if(ISS[i].Aname=="41Sc"){xx=ISS[i].AveT;yy=20;}
                    if(ISS[i].Aname=="47Cr"){xx=ISS[i].AveT+2;yy=10;}
                    if(ISS[i].Aname=="49Mn"){xx=ISS[i].AveT+2;yy=5;}
                    if(ISS[i].Aname=="51Fe"){xx=ISS[i].AveT+2;yy=-10;}
                    if(ISS[i].Aname=="55Ni"){xx=ISS[i].AveT+2;yy=-15;}

                    if(opt[0])lat_n_v1   ->DrawLatex(xx, yy,     lat_text);
                    if(opt[1])lat_n_v2->DrawLatex(xx, yy,     lat_text_v2);
                    if(opt[2])lat_n_VE->DrawLatex(xx, yy,     lat_text_VE);
              
                }
                else
                {
                    if(opt[0])lat_n_v1   ->DrawLatex(ISS[i].AveT+x_pos_shift, ISS[i].deltaMass+abs(ISS[i].deltaMass)/(ISS[i].deltaMass)*y_pos_ratio,     lat_text);
                    if(opt[1])lat_n_v2->DrawLatex(ISS[i].AveT+x_pos_shift, ISS[i].deltaMass_v2 +abs(ISS[i].deltaMass_v2)/(ISS[i].deltaMass_v2)*y_pos_ratio,     lat_text_v2);
                    if(opt[2])lat_n_VE->DrawLatex(ISS[i].AveT+x_pos_shift, ISS[i].deltaMass_VE +abs(ISS[i].deltaMass_VE)/(ISS[i].deltaMass_VE)*y_pos_ratio,     lat_text_VE);
                }
            }
            else 
            {
                if(opt[0])lat_n_v1   ->DrawLatex(ISS[i].AveT, ISS[i].deltaMass+abs(ISS[i].deltaMass)/(ISS[i].deltaMass)*y_pos_ratio,     lat_text   );
                if(opt[1])lat_n_v2->DrawLatex(ISS[i].AveT, ISS[i].deltaMass_v2+abs(ISS[i].deltaMass_v2)/(ISS[i].deltaMass_v2)*y_pos_ratio,  lat_text_v2);
                if(opt[2])lat_n_VE->DrawLatex(ISS[i].AveT, ISS[i].deltaMass_VE+abs(ISS[i].deltaMass_VE)/(ISS[i].deltaMass_VE)*y_pos_ratio,  lat_text_VE);
                
            }
        }
    }// for each species
    
    
}//________________ c1 show nuclei info ____________________

if(to_save) 
{
    //c1_AMEcompare_T_shadow->Print(FILEPATH_c1_saved +"mass_result_"+strtmp.Format("%d%d%d_%d_%d.png",PLOTMIN_control_OFF,known_unknown_control_off,known_unknown_switch,scan_loop_i,err_sys_ON));
    c1_AMEcompare_T_shadow->Print(FILEPATH_c1_saved +"mass_result_"+strtmp.Format("loop_%d-ID_%d.png",scan_loop_i,ID));
    //c1_AMEcompare_T_shadow->Print(FILEPATH_c1_savedw +"mass_result_"+strtmp.Format("%d.root",scan_loop_i));
}
if(LOOP_ON)delete c1_AMEcompare_T_shadow;

if(LOOP_ON)
{
    delete grerr_AME_shadow;
    delete grerr1;
    delete grerr2 ;
    delete grerrE;
    delete grerr1_add_AMEERR;
    delete grerr2_add_AMEERR;
    delete grerrE_add_AMEERR;
    delete grerr1_sys;
    delete grerr1_add_AMEERR_sys;
    delete grerr2_sys;
    delete grerr2_add_AMEERR_sys;
    delete grerrE_sys;
    delete grerrE_add_AMEERR_sys;
    for(int j=0;j<All_Tz_n;j++){delete grerr1_Tz[j];}
    for(int j=0;j<All_Tz_n;j++){delete grerr2_Tz[j];}
    for(int j=0;j<All_Tz_n;j++){delete grerrE_Tz[j];} 
}


}// if Draw ON


}//_____________________________________ show mass __________________________________________________


//跳过某些注入
void Do_Set_INJ_SKIP()
{
    cout<<" Do_Set_INJ_SKIP ON--------------------------- "<<endl;
    ifstream infile;
    //INJ_SKIP_INFILENAME="INPUT//INJ_SKIP.txt";
    //INJ_SKIP_INFILENAME="INPUT//INJ_SKIP_17Ne_20240704.txt";
    //INJ_SKIP_INFILENAME="INPUT//INJ_SKIP_17Ne_13N_11C.txt";
    INJ_SKIP_INFILENAME="INPUT//INJ_SKIP_20Mg.txt";
    cout<<" -----INJ_SKIP_INFILENAME= "<<INJ_SKIP_INFILENAME<<endl;
    infile.open(INJ_SKIP_INFILENAME);  // large_dm 表格的格式
    int inj_in;
    TString skip_filename_in;
    SKIP_n=0;
    while(infile>>strtmp>>strtmp>>strtmp>>inj_in>>skip_filename_in)
    {
        INJ_SKIP[inj_in] = true;
        //cout<<" debug inj_in = "<<inj_in<<endl;
        SKIP_FILENAME[SKIP_n] = skip_filename_in;
        //cout<<" debug skip_filename_in = "<<skip_filename_in<<endl;
        SKIP_n++;
    }        
    infile.close();
    if(SKIP_n>=MAX_INJECTIONS){cout<<" !!!!! "<<endl<<endl<<endl<<" SKIP_n>=MAX_INJECTIONS "<<endl;}
    cout<<" SKIP_n = "<<SKIP_n<<endl;
            
}

///////////////________________________________FUNCTION_______________________________////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TH1F* h2 = new TH1F("h2","gammat_C",10000,1.0,2.0); // gammat distribution
TH1F* h_Mvq = new TH1F("h_Mvq","h_Mvq",50000,1.5,1.97);
TH1F* h_Mvq_AME = new TH1F("h_Mvq_AME","h_Mvq_AME",50000,1.5,1.97);
TH1F* h_Mvq_14O21Mg_ReIdentify = new TH1F("h_Mvq_14O_21Mg", "h_Mvq_14O_21Mg", 1000, 1.75, 1.751);//名字按需改动
double Mvq_test_min = 2, Mvq_test_max = 1;



//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~




//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=






//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
    //# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    //                          ########################################
//                                        #                #
//                                        #                #
//                                        #                #      
//                                        #                #                                  
//                                       #                 #
//                                       #                 #
//                            ###############################################            
//                                      #                  #   
//                                     #                   #   
//                                    #                    #   
//                                    #                    #    
//                                   #                     #    
//                                   #                     #    
//                                  #                      #      
//                                        
//
//                                    #                    #
//                                   #                   #
//                                  #                  #          #  
//                                 #                 #             #                             
//                           ##################    ##################                                                                           
//                               #        #                          #         
//                               #      #                                     
//                                #     #         ###################                           
//                                 #   #          #                 #            
//                                  # #            #               #            
//                                   ##            #               #           
//                                   # #           #               #            
//                                  #   #           ################                            
//                                                                           
//                                                                           
//# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void GtCMASS()
{ 
// START 
gStyle->SetOptStat("nemr");
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
gStyle->SetEndErrorSize(6);
gStyle->SetLabelFont(132,"xyz");   // 全局 xyz三个轴上标签字体 默认 62 =   10*labelfontnumber + precision  
gStyle->SetTitleFont(132,"xyz");

lat_n->SetTextColor(1);
lat_n->SetTextFont(43);
lat_n->SetTextSize(20);

int DATE=TIMENOW.GetDate();
int TIME=TIMENOW.GetTime();

if(THIS_ENVIRONMENT==ON_UBUNTU) //在UBUNTU上运行，生成所需指令字符串 
{
    makedir=makedir.Format("mkdir -p OUTPUT//%d_%ddir",DATE,TIME); //后面那些是子文件夹的生成指令
    FILEPATH=FILEPATH.Format("OUTPUT//%d_%ddir//",DATE,TIME);  
    makedir_m_s =  makedir+"//m_s"; FILEPATH_m_s = FILEPATH+"m_s//";
    makedir_ionspecies =makedir+"//IONSPECIES"; FILEPATH_ionspecies = FILEPATH+"IONSPECIES//";
    makedir_c1_saved = makedir+"//c1_saved"; FILEPATH_c1_saved = FILEPATH+"c1_saved//";
    makedir_DodA0T = makedir+"//DodA0T"; FILEPATH_DodA0T = FILEPATH+"DodA0T//";
    makedir_k_T = makedir+"//k_T";FILEPATH_k_T = FILEPATH+"k_T//";
    makedir_chosen = makedir+"//chosen"; FILEPATH_chosen = FILEPATH+"chosen//";
    
}
else if(THIS_ENVIRONMENT==ON_WINDOWS)//在WINDOWS上运行
{        
    makedir=makedir.Format("mkdir  OUTPUT\\\\%d_%ddir",DATE,TIME); 
    FILEPATH=FILEPATH.Format("OUTPUT\\\\%d_%ddir\\\\",DATE,TIME); 
    makedir_m_s =  makedir+"\\\\m_s"; FILEPATH_m_s = FILEPATH+"m_s\\\\";
    makedir_ionspecies =makedir+"\\\\IONSPECIES"; FILEPATH_ionspecies = FILEPATH+"IONSPECIES\\\\";
    makedir_c1_saved = makedir+"\\\\c1_saved"; FILEPATH_c1_saved = FILEPATH+"c1_saved\\\\";
    makedir_DodA0T = makedir+"\\\\DodA0T"; FILEPATH_DodA0T = FILEPATH+"DodA0T\\\\";
    makedir_k_T = makedir+"\\\\k_T"; FILEPATH_k_T = FILEPATH+"k_T\\\\";
    makedir_chosen = makedir+"\\\\chosen"; FILEPATH_chosen = FILEPATH+"chosen\\\\";
}

if(MAKE_DIR)//是否创建文件夹
{//依次生成
    system(makedir);system(makedir_m_s);system(makedir_ionspecies);system(makedir_c1_saved);system(makedir_chosen);
    system(makedir_k_T);system(makedir_DodA0T);
}
if(MAKE_DIR)cout<<"\033[32m||makedir = "<<"TRUE"<<" ||\033[0m"<<endl;//   \033 转义序列开始  [32m 接下来的字为绿色 \033[0m 转义序列结束
else        cout<<"\033[32m||makedir = "<<"FALSE"<<" ||\033[0m"<<endl;

filename0=FILEPATH+"logs.txt";//日志文件
outfile_logs.open(filename0);

//----- 20250302 TFile
//Tfile_save_hist = new TFile(FILEPATH+"save_hist.root", "RECREATE");

//outfile1.open("LargeDeviation.txt");

NSpecies=0;      //start from 0
NSpecies_IAS=0;  ////number of IAS 
NSpecies_Isomer=0;  ////number of Isomer
NSpecies_other_exc =0; 

int count_use_massknown_but_notREF =0;//20240626
int Selected_n=0;
int gtC_n=0;
int calculated_unknown_ions_n =0;
int calculated_unknown_ions_n_v_others =0;
TGraph* gr_gtC_selected=new TGraph();

f_zero = new TF1("f_zero","pol1",-10000,10000);
f_zero->SetParameters(0,0);
f_zero->SetLineWidth(3);
f_zero->SetLineStyle(9);

// 20230419 error ana 
TGraphErrors* grerr_ERRYX = new TGraphErrors();
grerr_ERRYX->SetMarkerSize(1);
int grerr_ERRYX_n = 0;
TH2F* h2_ERRYX = new TH2F("","",500,0,0.1,100,-0.0002,0.0002);

//initialization
for(int i=0;i<120  ;i++)
    for(int j=0;j<200;j++)
    {  ZN_ID[i][j]=-1;  }
for(int i=0;i<MAX_IONSPECIES;i++)
{
    Z_readin[i]=0;
    A_readin[i]=0;
    N_readin[i]=0;
}

double A1_in,A1err_in,A2_in,A2err_in,A3_in,A4_in,dA0_in,dA0err_in,dA1_in,dA1err_in;//t = A0 + A1*N + A2*N^2 + A3*N^3
double cov12_in,cov15_in,cov16_in,cov25_in,cov26_in,cov56_in;//各项之间的协方差
int ion_readin,inj_readin,ion_number_new,inj_new,inj_pre;//ion_readin此次注入的粒子数 inj_readin读取的是第几次注入 ion_number_new这次注入的第几个核 inj_new这一次是第几次注入 inj_pre记录上一次读取的是第几次注入
string filename_in;

//-------------------------------------- 2017 58iNi -----------------------------------
//----58Ni .log input format 20231107
double TOF,TOF_err,Delta_t_ToF,DeltaTTmp,gammat_in;
double CovMatrix_in[36];
//double time_start_last_up,time_start_last_down;
double A0up_in,A0errup_in,A1up_in,A1errup_in,A2up_in,A2errup_in,A3up_in,A3errup_in;
double A0down_in,A0errdown_in,A1down_in,A1errdown_in,A2down_in,A2errdown_in,A3down_in,A3errdown_in;
double X2up_in,X2down_in,AmpAve_up,AmpAve_down;
int TotalFitPoint_up_in,TotalFitPoint_down_in,TurnFirstUP_in,TurnFirstDOWN_in,TurnLastUP_in,TurnLastDOWN_in;
double N_m; // turn_middle

//----58Ni covmatrix[100] input format 20211203 
double N_middle_in;
double V0u0u,V0u1u,V0u0d,V0u1d,V0u2,V0u3;
double V1u1u,V1u0d,V1u1d,V1u2,V1u3;
double V0d0d,V0d1d,V0d2,V0d3;
double V1d1d,V1d2,V1d3;
double V22,V23;
double V33;
//-------------------------------------- THIS_EXP -----------------------------------
//.....ooo00000oooo..........ooo00000oooo..........ooo00000oooo..........ooo00000oooo..........ooo00000oooo.....
if(THIS_EXP=="2017_58Ni")
{
    //INFILENAME_1 = "INPUT//2017_58Ni//58Ni_gtCalculator_outfile_20231121.txt";      // data input 
    //INFILENAME_1 = "INPUT//2017_58Ni//58Ni_gtCalculator_outfile_20240917.txt";      // 20240917 增加 V23 V33 小数点后到16位
    INFILENAME_1 = "INPUT//2017_58Ni//58Ni_gtCalculator_outfile_20240917-del4936.txt";      // 20240917 增加 V23 V33 小数点后到16位
    
    INFILENAME_2 = "INPUT//2017_58Ni//58Ni_IONSPECIES_AME.txt";                //为本设置准备好的AME20核数据文件
    INFILENAME_4 = "INPUT//2017_58Ni//58Ni_IONSPECIES_ISOMER.txt";         //Isomer 
}
else if (THIS_EXP == "2021_36Ar_SET2")
{
    INFILENAME_1 = "INPUT//2021_36Ar_SET2//IdentifyResult_20230224_RUN2-b_part2.txt";
    //INFILENAME_1 = "INPUT//2021_36Ar_SET2//IdentifyResult_20230224_RUN2-b.txt";
    //INFILENAME_1 = "INPUT//2021_36Ar_SET2//IdentifyResult_20230224_RUN2-b_part1.txt";
    //INFILENAME_1 = "INPUT//2021_36Ar_SET2//IdentifyResult_20230922_RUN3-B-FORM_RUN2_PART2.txt"; //20230922 res_tof

    INFILENAME_2 = "INPUT//2021_36Ar_SET2//36Ar_IONSPECIES_AME.txt";//为本设置准备好的AME20核数据文件
    INFILENAME_3 = "INPUT//2021_36Ar_SET2//36Ar_IONSPECIES_IAS.txt";//同质异位素相似态（isobaric analog state简称IAS）
}
else if(THIS_EXP=="2021_36Ar_SET3")
{
    //INFILENAME_1 = "INPUT//2021_36Ar_SET3//IdentifyResult_20230824_174150-del-2inj.txt";  // identical to 20240626-XYM-Analysis_28S_SET3.txt
    INFILENAME_1 = "INPUT//2021_36Ar_SET3//IdentifyResult_20240815_T_region_620_670-del5.txt";  // 

    //INFILENAME_1 = "INPUT//2021_36Ar_SET3//IdentifyResult_20240815_T_region_620_670-s.txt";  // 
    
    //INFILENAME_1 = "INPUT//2021_36Ar_SET3//IdentifyResult_20240815_T_region_620_670_p1.txt";  // 
    //INFILENAME_1 = "INPUT//2021_36Ar_SET3//IdentifyResult_20240815_T_region_620_670_p2.txt";  // 
    //INFILENAME_1 = "INPUT//2021_36Ar_SET3//IdentifyResult_20240815_T_region_620_670_p3.txt";  // 
    
    //INFILENAME_2 = "INPUT//2021_36Ar_SET3//36Ar_SET3_IONSPECIES_AME.txt";                //为本设置准备好的AME20核数据文件
    INFILENAME_2 = "INPUT//2021_36Ar_SET3//36Ar_SET3_IONSPECIES_AME_620.txt";                //全周期
    
    INFILENAME_3 = "INPUT//2021_36Ar_SET3//36Ar_SET3_IONSPECIES_IAS_ISOMER.txt";         //同质异位素相似态（isobaric analog state简称IAS）
}
else if (THIS_EXP == "2021_36Ar_SET1")
{
    //INFILENAME_1 = "INPUT//2021_36Ar_SET1//IdentifyResult_20250119_part1.txt";  //已经分好段了
    INFILENAME_1 = "INPUT//2021_36Ar_SET1//IdentifyResult_20250119_part2.txt";

    INFILENAME_2 = "INPUT//2021_36Ar_SET1//36Ar_23Si_IONSPECIES_AME.txt";                //为本设置准备好的AME20核数据文件
    INFILENAME_3 = "INPUT//2021_36Ar_SET1//36Ar_23Si_IONSPECIES_IAS_ISOMER.txt";         //同质异位素相似态（isobaric analog state简称IAS）
}

else
{
    cout<<"\033[7m\033[31m  FATAL ERROR ! ! ! ! THIS_EXP is wrong !!!  \033[0m "<<endl;
    cout<<" failed to load input files . return."<<endl;
    return;
}

//_____________________________________________________________________________________________

cout<<  "========================     THIS EXP IS :  \033[7m\033[33m"<<THIS_EXP<<" \033[0m       =============================="<<endl;
outfile_logs<<"======================== THIS EXP IS :  "<<THIS_EXP<<"  =============================="<<endl;
//|| 2017_58Ni 18.046 0.1470
//|| 2021_36Ar_SET3 18.057 0.0957
//|| 2021_36Ar_SET2 18.052 0.0760
if(THIS_EXP=="2017_58Ni")
{    
    L_down = 18.046;
    ddT_down = 0.1470;

    T_RANGE_MIN = 610;   T_RANGE_MAX = 645;
    mvq_RANGE_MIN = 1.77; mvq_RANGE_MAX = 1.97;
    CONDITION_gt_Z= 6;
    CONDITION_gt_A= 18;
    gtC_ERR_upper_bound = 0.02;
    
}

if(THIS_EXP=="2021_36Ar_SET3")
{
    L_down = 18.057;
    ddT_down = 0.0957;


    T_RANGE_MIN = 630;   T_RANGE_MAX = 670;
    mvq_RANGE_MIN = 1.61; mvq_RANGE_MAX = 1.87;
    CONDITION_gt_Z= 6;
    CONDITION_gt_A= 14;
    C_filter_ON = 1;
    C_filter_min =128.63;
    C_filter_max =129.00;
        gtC_ERR_upper_bound = 0.02;

}
if(THIS_EXP=="2021_36Ar_SET2")
{
    L_down = 18.052;
    ddT_down = 0.0760;
    T_RANGE_MIN = 590;   T_RANGE_MAX = 680;
    mvq_RANGE_MIN = 1.50; mvq_RANGE_MAX = 2.20;
    CONDITION_gt_Z= 8;
    CONDITION_gt_A= 15;
    C_filter_ON = 1;
    C_filter_min =128.60;
    C_filter_max =129.05;
        gtC_ERR_upper_bound = 0.02;

}
if(THIS_EXP=="2021_36Ar_SET1")
{
    T_RANGE_MIN = 630;   T_RANGE_MAX = 670;
    mvq_RANGE_MIN = 1.61; mvq_RANGE_MAX = 1.87;
    CONDITION_gt_Z= 6;
    CONDITION_gt_A= 14;
    //C_filter_ON = 1;
    //C_filter_min =128.63;
    //C_filter_max =129.00;
}
//===========================  READ IN ALL IONS =============================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inj_pre = -1;
ion_number_new = 0; inj_new = 0;
int ion_number_max=0;//单次注入的最多粒子数
int input_lines = 0;

ofstream outfile_ion_readin_check;
if(ReadInIons_check_on)outfile_ion_readin_check.open("OUTPUT//ions_readin_check.txt");
// get T for  all ions
NSpecies = 0;
int i_in=-1;//读入的序列数（从0开始，用于数组）

infile.open(INFILENAME_1);
if(!infile){cerr<<"infile failed at all ions readin";exit(1);}
outfile_logs<<"read in all ions from file: "<<INFILENAME_1<<endl;
//----------------- 读入数据 2021_36Ar_SET3 ------------------
//----------------- 读入数据 2021_36Ar_SET2 ------------------
if(THIS_EXP == "2021_36Ar_SET3" || THIS_EXP == "2021_36Ar_SET2")
{
while(infile>>namestr>>A>>Z>>T>>T_err>>ion_readin>>inj_readin
            >>A1_in>>A1err_in>>A2_in>>A2err_in>>A3_in>>A4_in
            >>dA0_in>>dA0err_in>>dA1_in>>dA1err_in
            >>cov12_in>>cov15_in>>cov16_in>>cov25_in>>cov26_in>>cov56_in
            >>filename_in)//正常读入
{    //使用ION类，依次保存每个粒子信息    
    //revolution time T_in(ns)  dA0 input(ns) then -- [ps]   
    // all err unit is [ps*ps]   2023 Feb.
    i_in++;
    ions_n++;
    ions[i_in].ID=i_in;
    ions[i_in].Z=Z;
    ions[i_in].A=A;
    ions[i_in].name=namestr;
    ions[i_in].T=T;           // [ns]
    ions[i_in].T_err = T_err;      // [ns]
    
    ions[i_in].A1=A1_in;ions[i_in].A1err=A1err_in;ions[i_in].A2=A2_in;ions[i_in].A2err=A2err_in;//各项误差
    ions[i_in].A3=A3_in;ions[i_in].A4=A4_in;
    ions[i_in].dA0=dA0_in;ions[i_in].dA0err=dA0err_in;ions[i_in].dA1=dA1_in;ions[i_in].dA1err=dA1err_in;
    ions[i_in].dA0=dA0_in*1000;    // [ps]
    ions[i_in].cov12=cov12_in;ions[i_in].cov15=cov15_in;ions[i_in].cov16=cov16_in;
    ions[i_in].cov25=cov25_in;ions[i_in].cov26=cov26_in;ions[i_in].cov56=cov56_in;
    ions[i_in].time = filename_to_time_second(filename_in);//此次注入是在实验开始(2021.10.24)后多久
    ions[i_in].inject_filename = filename_in;  

    //不要求 ion_readin 和 inj_readin 是连续的， 重新设定离子数和注入数
    if(inj_readin==inj_pre)//与之前的粒子是同一次注入
    {
        ions[i_in].inject_number=inj_new;
        ion_number_new++;
        ions[i_in].ion_number=ion_number_new;//是这次注入的第几个核
    }
    else//新注入
    {
        inj_new++;
        ions[i_in].inject_number=inj_new;
        ion_number_new=1;//是这次注入的第1个核
        ions[i_in].ion_number=ion_number_new;
    }
    inj_pre = inj_readin;
    //now sort the ions read in 
    IsNew=true;  //initialization   //这是一个新核素吗？
    for(int j=0;j<NSpecies ;j++)
    {
        if(Z_readin[j]==Z&&A_readin[j]==A )//不是
        {//already has this name,not a new kind
            ions[i_in].Species=j;          //给这个粒子对应核素的ID
            N_readin[j]++;                 //统计该核素数量
            IsNew = false;
            break;            
        }
    }
    //if(NSpecies==0){NSpecies++;name[0]=namestr;}  //for the first ion, it didnt enter the loop above
    if(IsNew)                               //是新的
    {
        name[NSpecies]=namestr; //记录此核素名字 //record newone in No.NSpecies before NSpecies++
        Z_readin[NSpecies]=Z;    
        A_readin[NSpecies]=A;  
        N_readin[NSpecies]=1;                    
        ions[i_in].Species=NSpecies;      
        NSpecies++;             //总种类数+1
    }

    //injection info record
    Injection[inj_new] = ion_number_new;  //这次注入的第几个核 //subfix of Injection starts from 1. Finally Injection array will save the max IonNumber in one injection
    if(ion_number_new>ion_number_max)ion_number_max=ion_number_new;
    //if(i_in>MAX_IONS)break;
} //___________while infile__________

infile.close();//全部读取完毕
}

//----------------- 读入数据 2021_36Ar_SET1 ------------------
if (THIS_EXP == "2021_36Ar_SET1")
{
    while (infile >> namestr >> A >> Z >> T >> T_err >> dA0_in >> dA0err_in >> ion_readin >> inj_readin >> A2_in >> dA1_in >> filename_in)//正常读入
    {    //使用ION类，依次保存每个粒子信息    
        //revolution time T_in(ns)  dA0 input(ns) then -- [ps]   
        // all err unit is [ps*ps]   2023 Feb.
        i_in++;
        ions_n++;
        ions[i_in].ID = i_in;
        ions[i_in].Z = Z;
        ions[i_in].A = A;
        ions[i_in].name = namestr;
        ions[i_in].T = T;           // [ns]
        ions[i_in].T_err = T_err;      // [ns]

        ions[i_in].dA0 = dA0_in * 1000;//[ps]
        ions[i_in].dA0err = dA0err_in * 1000;//[ns]
        ions[i_in].A2 = A2_in;
        ions[i_in].dA1 = dA1_in;
        ions[i_in].A1 = T * 1000 + dA1_in * 0.5;//[ps]



        ions[i_in].A1err = 0; ions[i_in].A2err = 0;//各项误差
        ions[i_in].A3 = 0; ions[i_in].A4 = 0;
        ions[i_in].dA1err = 0;
        ions[i_in].cov12 = 0; ions[i_in].cov15 = 0; ions[i_in].cov16 = 0;
        ions[i_in].cov25 = 0; ions[i_in].cov26 = 0; ions[i_in].cov56 = 0;
        ions[i_in].time = filename_to_time_second(filename_in);//此次注入是在实验开始(2021.10.24)后多久
        ions[i_in].inject_filename = filename_in;

        //不要求 ion_readin 和 inj_readin 是连续的， 重新设定离子数和注入数
        if (inj_readin == inj_pre)//与之前的粒子是同一次注入
        {
            ions[i_in].inject_number = inj_new;
            ion_number_new++;
            ions[i_in].ion_number = ion_number_new;//是这次注入的第几个核
        }
        else//新注入
        {
            inj_new++;
            ions[i_in].inject_number = inj_new;
            ion_number_new = 1;//是这次注入的第1个核
            ions[i_in].ion_number = ion_number_new;
        }
        inj_pre = inj_readin;
        //now sort the ions read in 
        IsNew = true;  //initialization   //这是一个新核素吗？
        for (int j = 0; j < NSpecies; j++)
        {
            if (Z_readin[j] == Z && A_readin[j] == A)//不是
            {//already has this name,not a new kind
                ions[i_in].Species = j;          //给这个粒子对应核素的ID
                N_readin[j]++;                 //统计该核素数量
                IsNew = false;
                break;
            }
        }
        //if(NSpecies==0){NSpecies++;name[0]=namestr;}  //for the first ion, it didnt enter the loop above
        if (IsNew)                               //是新的
        {
            name[NSpecies] = namestr; //记录此核素名字 //record newone in No.NSpecies before NSpecies++
            Z_readin[NSpecies] = Z;
            A_readin[NSpecies] = A;
            N_readin[NSpecies] = 1;
            ions[i_in].Species = NSpecies;
            NSpecies++;             //总种类数+1
        }

        //injection info record
        Injection[inj_new] = ion_number_new;  //这次注入的第几个核 //subfix of Injection starts from 1. Finally Injection array will save the max IonNumber in one injection
        if (ion_number_new > ion_number_max)ion_number_max = ion_number_new;
        //if(i_in>MAX_IONS)break;
    } //___________while infile__________

    infile.close();//全部读取完毕
}


//----------------- 读入数据 2017_58Ni ------------------
if(THIS_EXP == "2017_58Ni")
{
    while(infile>>Z>>A
            >>T>>T_err>>TOF>>TOF_err>>Delta_t_ToF>>DeltaTTmp>>gammat_in
            >>A0up_in>>A0errup_in>>A1up_in>>A1errup_in>>A2up_in>>A2errup_in>>A3up_in>>A3errup_in
            >>A0down_in>>A0errdown_in>>A1down_in>>A1errdown_in
            >>ion_readin>>inj_readin
            //>>X2up_in>>X2down_in
            >>N_middle_in
            //>>TurnFirstUP_in>>TurnFirstDOWN_in>>TurnLastUP_in>>TurnLastDOWN_in
            >>filename_in  
            // input CovMatrix[] terms v20211203
            >>V0u0u>>V0u1u>>V0u0d>>V0u1d>>V0u2>>V0u3
            >>V1u1u>>V1u0d>>V1u1d>>V1u2>>V1u3
            >>V0d0d>>V0d1d>>V0d2>>V0d3
            >>V1d1d>>V1d2>>V1d3
            >>V22>>V23
            >>V33
                  )    
{    
    //for(int i=0;i<36  ;i++){ infile>>CovMatrix_in[i];  }   // input CovMatrix.root [36]  
    //gtCalculator RevT [ps] TOF [ps]
    //revolution time T_in(ns)  dA0 input(ns) then -- [ps]   
    // all err unit is [ps*ps] A1err A2err --variance  2023 Feb.
    i_in++;
    ions_n++;   
    ions[i_in].ID=i_in;
    ions[i_in].Z=Z;
    ions[i_in].A=A;
    namestr = convert_z_to_name(Z);
    ions[i_in].name=namestr;
    
    ions[i_in].turn_middle = N_middle_in;
    N_m = ions[i_in].turn_middle;
    ions[i_in].T=( (A1up_in+A1down_in)*0.5+2*A2up_in*N_m+3*A3up_in*N_m*N_m )*0.001 ;           // [ns] 
    ions[i_in].T_err = sqrt(0.25*pow(A1errup_in,2)+0.25*pow(A1errdown_in,2) -0.5*V1u1d)*0.001;      // [ns]
//if(i_in%1000==0)cout<<" RevT_in = "<<0.001*T<<" +- "<<0.001*T_err <<" T_cal= "<<ions[i_in].T<<" +- "<<ions[i_in].T_err<<" diff= "<<0.001*T-ions[i_in].T<<endl;///debug
    
    ions[i_in].A1=A1down_in;      ions[i_in].A1err=A1errdown_in*A1errdown_in ;  //[ps]  [ps*ps]
    ions[i_in].A2=A2up_in;        ions[i_in].A2err=A2errup_in*A2errup_in;       //[ps]  [ps*ps]
    ions[i_in].A3=A3up_in;        ions[i_in].A3err=A3errup_in*A3errup_in;    //V33=0 can not be used 

    ions[i_in].dA0=A0down_in - A0up_in ;
    ions[i_in].dA1=A1down_in - A1up_in ;

    ions[i_in].dA0err= pow( A0errup_in ,2)+pow( A0errdown_in ,2)-2*V0u0d ;  // [ps] [ps*ps]
    ions[i_in].dA1err= pow( A1errup_in ,2)+pow( A1errdown_in ,2)-2*V1u1d ;  // [ps] [ps*ps]
    ions[i_in].cov12=V1d2;
    ions[i_in].cov15=V0d1d - V0u1d ;
    ions[i_in].cov16=V1d1d - V1u1d;
    ions[i_in].cov25=V0d2 - V0u2 ;
    ions[i_in].cov26=V1d2 - V1u2 ;
    ions[i_in].cov56=V0u1u - V0u1d - V1u0d + V0d1d ;
    ions[i_in].cov13=V1d3;  
    ions[i_in].cov23=V23;
    ions[i_in].cov35=V0d3- V0u3;
    ions[i_in].cov36=V1d3 - V1u3;
    ions[i_in].inject_filename = filename_in; 
    ions[i_in].time = filename_to_time_second(filename_in);
    //___________________________ 20231120 format __________________________

    if(inj_readin==inj_pre)
    {
        ions[i_in].inject_number=inj_new;
        ion_number_new++;
        ions[i_in].ion_number=ion_number_new;
    }
    else
    {
        inj_new++;
        ions[i_in].inject_number=inj_new;
        ion_number_new=1;
        ions[i_in].ion_number=ion_number_new;
    }
    inj_pre = inj_readin;
    //now sort the ions read in 
    IsNew=true;  //initialization 
    for(int j=0;j<NSpecies ;j++)
    {
        if(Z_readin[j]==Z&&A_readin[j]==A )
        {//already has this name,not a new kind
            ions[i_in].Species=j;
            N_readin[j]++;
            IsNew = false;
            break;            
        }
    }
    //if(NSpecies==0){NSpecies++;name[0]=namestr;}  //for the first ion, it didnt enter the loop above
    if(IsNew)
    {
        if((Z<0||Z>120)||(A<0||A>300)){cout<<endl<<endl<<" !!!!!!fatal error returned! illeagal ! Z="<<Z<<" A= "<<A<<endl<<endl;return;}
        name[NSpecies]=namestr;  //record newone in No.NSpecies before NSpecies++
        Z_readin[NSpecies]=Z;    
        A_readin[NSpecies]=A;  
        N_readin[NSpecies]=1;                    
        ions[i_in].Species=NSpecies;      
        NSpecies++;
    }

    //injection info record
    Injection[inj_new] = ion_number_new;  //subfix of Injection starts from 1. Finally Injection array will save the max IonNumber in one injection
    if(ion_number_new>ion_number_max)ion_number_max=ion_number_new;
    //if(i_in>120000)break;
} //___________while infile__________

infile.close();
}

if(ReadInIons_check_on)//检查
{
    for(int i=0;i<ions_n  ;i++)
    {
        ions[i].PrintInfo_readin_58Ni(outfile_ion_readin_check);//输出读取到的信息
    }
    outfile_ion_readin_check.close();

    cout<<"------------ ReadInIons_check_on :"<<endl;
    for(int i=0;i<NSpecies  ;i++)
    {
        cout<<i<<" Z_readin["<<i<<"]: "<<Z_readin[i]<<" A_readin["<<i<<"]: "<<A_readin[i]<<" N_readin["<<i<<"]: "<<N_readin[i]<<endl;
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
total_injection =  ions[ions_n-1].inject_number;           //the last value is the max InjectNumber. starts from 1 !!!

cout<<" \033[33m------------------------------------  \033[0m "<<endl;
cout<<"| input lines = "<<i_in+1<<endl;
cout<<"| NSpecies = "<<NSpecies <<endl;  //  starts from 1
cout<<"| ions_n = "<<ions_n <<endl;      //  starts from 1
cout<<"| total_injection = "<<total_injection <<endl;  // starts from 1
cout<<"| max ion number in one injection = "<<ion_number_max <<endl;
cout<<" \033[33m------------------------------------  \033[0m "<<endl;
//========================================================================================================
if(ReadInIons_check_on){cout<<"ReadInIons_check_on ! return after input !"<<endl; return;}  //######debug

// ====================================  Set ionspecies   =====================================================
// 
for(int i=0;i<NSpecies  ;i++)//设置各种核素的信息
{
    ionspecies[i].name=name[i];
    ionspecies[i].Z=Z_readin[i];
    ionspecies[i].A=A_readin[i];
    ionspecies[i].SetTz();
    ionspecies[i].SetName_form();
    ionspecies[i].Species=i;//设置该核素的ID
    ionspecies[i].AME = 114514;//先进行初始化
    ZN_ID[ionspecies[i].Z][ionspecies[i].A-ionspecies[i].Z] = i;//在ZN网格中设置ID
    
    if(SetIonSpecies_check_on)//输出检查一下
    {

    if(i==0)cout<<"============= set ionspecies check ==================="<<endl;
    cout<<i<<"   "<<ionspecies[i].name<<" "<<ionspecies[i].name_latex<<" "<<ionspecies[i].A<<" "<<ionspecies[i].Z<<" "<<ionspecies[i].Tz 
    <<" ZN_ID ("<<Z_readin[i]<<","<<A_readin[i]-Z_readin[i]<<") = "<<ZN_ID[Z_readin[i]][A_readin[i]-Z_readin[i]]<<endl;
    if(i==NSpecies-1)cout<<"_____________ set ionspecies check ___________________"<<endl;

    }

}
//20221019=========== ionspecies mkdir
if(Create_ionspecies_folder_ON)CreateFolder(ionspecies,NSpecies, makedir_ionspecies);

//20231211=========== ionspecies Tz

for(int i=0;i<NSpecies  ;i++)
{
    All_Tz.push_back(ionspecies[i].Tz);//每种核素的Tz
}
sort(All_Tz.begin(), All_Tz.end());    // 排序
All_Tz.erase(unique(All_Tz.begin(), All_Tz.end()) , All_Tz.end() ); //去掉重复值
for(const auto& i:All_Tz){cout<<" Tz: "<<i<<" ";}
All_Tz_n = All_Tz.size();
cout<<endl<<" --- All_Tz_n = "<<All_Tz_n<<endl;

for(const auto& i:All_Tz)
{
    if(int(i)!=i)
    {
        All_Tz_str.push_back(strtmp.Format("%d/2", int(2*i) ) );
    }
    else
    {
        All_Tz_str.push_back(strtmp.Format("%d", int(i) ) );
    }
}
for(const auto& i:All_Tz_str){cout<<" Tz_str: "<<i<<" ";}
cout<<endl;
for(const auto& i:All_Tz){All_Tz_COLOR.push_back(GetTzColor(i));} // initialization of Tz color

// ____________________________________ Set ionspecies   _________________________________________________________

//--------------- @EXP@ 2021 36Ar SET3 手动添加 21Mg
if(THIS_EXP=="2021_36Ar_SET3"&&Re_IonIdentify_on)
{
    bool SET3_21Mg_SET = 0;
    for (int i = 0; i < NSpecies; i++){if(ionspecies[i].Aname=="21Mg"){SET3_21Mg_SET=1; cout<<"   21Mg already read in  "<<endl; break;} }
    if(!SET3_21Mg_SET)
    {
        ionspecies[NSpecies].name="Mg";
        ionspecies[NSpecies].Z=12;
        ionspecies[NSpecies].A=21;
        ionspecies[NSpecies].SetTz();
        ionspecies[NSpecies].SetName_form();
        ionspecies[NSpecies].Species=NSpecies;//设置该核素的ID
        ZN_ID[12][9] = NSpecies;//在ZN网格中设置ID
        
        ionspecies[NSpecies].AME=10903.9;                    //质量中心值 ##ARTIFICIAL
        ionspecies[NSpecies].AME_err=0.8;            //质量误差
        ionspecies[NSpecies].MassUnknown =0;  //是否已知
        ionspecies[NSpecies].SetMass();         //use AME mass excess: IONSpecies.Mass(nuclear) = A*u + AME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
        ionspecies[NSpecies].SetMvq_AME();      //Mvq_AME = Mass/(u*Z);  nuclear mvq
        
    
        cout<<"----- build ionspecies 21Mg by hand :"<<endl;
        cout<<NSpecies<<"   "<<ionspecies[NSpecies].Aname<<" "<<ionspecies[NSpecies].name_latex<<" "<<ionspecies[NSpecies].A<<" "<<ionspecies[NSpecies].Z<<" "<<ionspecies[NSpecies].Tz 
        <<" ZN_ID ("<<12<<","<<9<<") = "<<ZN_ID[12][9]<<endl;
        cout<<ionspecies[NSpecies].AME<<" +- "<<ionspecies[NSpecies].AME_err<<" "<<ionspecies[NSpecies].Mass<<endl;
    
        NSpecies++;
        cout<<"after 21Mg built: NSpecies = "<<NSpecies<<endl;
    }
}//if(THIS_EXP=="2021_36Ar_SET3"&&Re_IonIdentify_on)

// ================================  READ IN AME data  for ref ions =====================================================

TString name_append;
double AME_tmp,AME_err_tmp;//AME_tmp：质量过剩，AME_err_tmp：误差
double exc_tmp,exc_err_tmp;
int mass_unknown_in = 0;//质量是否未知 1未知 0已知
int AME_all_n=0;//读取AME的总核素数
int AME_appeared_n=0;//这次试验中，找到了这些核素中的几种核素

if(READ_prepared_EXP_AME) //读取每次实验专门准备的AME 信息， 只包括这次实验出现的核
{
    //读取参考核文件
    
    i_in=0;
    cout<<" \033[33m------------------------------------  \033[0m "<<endl;
    cout<<" read in AME data from: "<<INFILENAME_2<< endl;
    
    infile.open(INFILENAME_2);
    while(infile>>A>>Z>>strtmp>>AME_tmp>>AME_err_tmp>>mass_unknown_in)//读入数据
    {
        AME_all_n++;
        if(AME_all_n>100)break;//太多了
    
        int j = ZN_ID[Z][A-Z];//之前统计了初始数据的核素种类，现在看看和AME文件能不能对上
        if(j<0)
        {
            cout<<A<<strtmp<<" "<<Z<<" not appeared in the data "<<endl;//这个种类的核在本次实验中没发现
            continue;
        }  // info available but does not appear in present IonSpecies
        AME_appeared_n++;
        ionspecies[j].AME=AME_tmp;                    //质量中心值
        ionspecies[j].AME_err=AME_err_tmp;            //质量误差
        ionspecies[j].MassUnknown = mass_unknown_in;  //是否已知
        ionspecies[j].SetMass();         //use AME mass excess: IONSpecies.Mass(nuclear) = A*u + AME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
        ionspecies[j].SetMvq_AME();      //Mvq_AME = Mass/(u*Z);  nuclear mvq
        
        if(ReadInAME_check_on)//输出检查一下
        {
        cout<<AME_appeared_n<<" ZN_ID["<<Z<<"]["<<A-Z<<"] = "<< j<<" "<<ionspecies[j].name<<" "<<ionspecies[j].AME<<" +- "<<ionspecies[j].AME_err<<" "<<ionspecies[j].Mass<<endl;
        }          
        //20230720 
        
        //if(A==14&&Z==8){ionspecies[j].MassUnknown=1;cout<<endl<<" Adjust 14O  to be MassUnknown "<<endl;}
        //if(A==21&&Z==12){ionspecies[j].MassUnknown=1;cout<<endl<<" Adjust 21Mg to be MassUnknown"<<endl;}
    }
    
    infile.close();
    int NSpecies_test = 0;
    for (int i = 0; i < NSpecies; i++)//现在看一下，AME文件是否缺少某些核素
    {
        if (ionspecies[i].AME ==114514)//没有成功存进来
        {
            cout << "\"" << INFILENAME_2 << "\"" << " is missing " << ionspecies[i].A << ionspecies[i].name << " information." << endl;//缺少目标核素AME质量，无法计算
            NSpecies_test++;
        }
    }
    if (NSpecies_test != 0) { cout<<" fatal error! NSpecies_test>=1 !"<<endl;return; }
    
    cout<<endl<<" read in AME data lines : "<<AME_all_n <<" appeared ionspecies : "<<AME_appeared_n<<endl;
    // _______________________________  READ IN AME data  for ref ions ________________________________________________________

    
    
    // ==================READ IN AME mass  for  IAS  ================================================
    //-------------------------------------------------------------------------------------
    if(THIS_EXP=="2021_36Ar_SET3" || "2021_36Ar_SET2")
    {
        
        NSpecies_IAS = 0;
        //Initialization IAS=-1
        for(int i=0;i<NSpecies  ;i++){ionspecies[i].IAS_n=0;}
        infile.open(INFILENAME_3);                        //读入IAS态               
        i_in=0;
        while(infile>>A>>Z>>namestr>>name_append>>AME_tmp>>AME_err_tmp)
        {
            IAS_ionspecis[NSpecies_IAS].A = A;
            IAS_ionspecis[NSpecies_IAS].Z = Z;
            IAS_ionspecis[NSpecies_IAS].name = namestr;
            IAS_ionspecis[NSpecies_IAS].name += ( "_"+name_append);
            IAS_ionspecis[NSpecies_IAS].AME = AME_tmp;
            IAS_ionspecis[NSpecies_IAS].AME_err = AME_err_tmp;
        
            IAS_ionspecis[NSpecies_IAS].SetMass();    // total mass keV
            IAS_ionspecis[NSpecies_IAS].SetMvq_AME();
            IAS_ionspecis[NSpecies_IAS].SetName_form();
            if(name_append=="m"||name_append=="n")
            {
                for(int i=0;i< NSpecies ;i++)
                {
                    if(IAS_ionspecis[NSpecies_IAS].A ==ionspecies[i].A &&IAS_ionspecis[NSpecies_IAS].Z ==ionspecies[i].Z)
                    {
                        ionspecies[i].Isomer_n = 1;
                    }
                }
            }
           
            NSpecies_IAS++;
            if(NSpecies_IAS>=IAS_N_MAX){cout<<" \033[7m\033[31m  FATAL ERROR\033[0m"<<"NSpecies_IAS > "<<IAS_N_MAX<<endl<<endl;return; }
        }
        infile.close();

        cout<<" \033[33m------------------------------------  \033[0m "<<endl;
        cout<<" NSpecies_IAS = "<<NSpecies_IAS<<endl;
        cout<<" \033[33m------------------------------------  \033[0m "<<endl;
       
    }
    //-------------------------------------------------------------------------------------

    // ==================READ IN AME mass  for  IAS  ================================================
    //-------------------------------------------------------------------------------------
    if (THIS_EXP == "2021_36Ar_SET1")
    {

        NSpecies_IAS = 0;
        //Initialization IAS=-1
        for (int i = 0; i < NSpecies; i++) { ionspecies[i].IAS_n = 0; }

        infile.open(INFILENAME_3);                        //读入IAS态               
        i_in = 0;
        while (infile >> A >> Z >> namestr >> name_append >> AME_tmp >> AME_err_tmp)
        {
            IAS_ionspecis[NSpecies_IAS].A = A;
            IAS_ionspecis[NSpecies_IAS].Z = Z;
            IAS_ionspecis[NSpecies_IAS].name = namestr;
            IAS_ionspecis[NSpecies_IAS].name += ("_" + name_append);
            IAS_ionspecis[NSpecies_IAS].AME = AME_tmp;
            IAS_ionspecis[NSpecies_IAS].AME_err = AME_err_tmp;

            IAS_ionspecis[NSpecies_IAS].SetMass();    // total mass keV
            IAS_ionspecis[NSpecies_IAS].SetMvq_AME();
            IAS_ionspecis[NSpecies_IAS].SetName_form();
            if (name_append == "m" || name_append == "n")
            {
                for (int i = 0; i < NSpecies; i++)
                {
                    if (IAS_ionspecis[NSpecies_IAS].A == ionspecies[i].A && IAS_ionspecis[NSpecies_IAS].Z == ionspecies[i].Z)
                    {
                        ionspecies[i].Isomer_n = 1;
                    }
                }
            }

            NSpecies_IAS++;
            if (NSpecies_IAS >= IAS_N_MAX) { cout << " \033[7m\033[31m  FATAL ERROR\033[0m" << "NSpecies_IAS > " << IAS_N_MAX << endl << endl; return; }
        }
        infile.close();

        cout << " \033[33m------------------------------------  \033[0m " << endl;
        cout << " NSpecies_IAS = " << NSpecies_IAS << endl;
        cout << " \033[33m------------------------------------  \033[0m " << endl;

    }

    // ==================READ IN AME mass  for  Isomer  ================================================
    if(THIS_EXP=="2017_58Ni")
    {
       
        NSpecies_Isomer = 0;
        for(int i=0;i<NSpecies  ;i++){ionspecies[i].Isomer_n=0;}
        infile.open(INFILENAME_4);                        //读入IAS态               
        i_in=0;
        while(infile>>A>>Z>>namestr>>name_append>>AME_tmp>>AME_err_tmp>>strtmp>>exc_tmp>>exc_err_tmp)
        {
            Isomer_ionspecis[NSpecies_Isomer].A = A;
            Isomer_ionspecis[NSpecies_Isomer].Z = Z;
            Isomer_ionspecis[NSpecies_Isomer].name = namestr;
            Isomer_ionspecis[NSpecies_Isomer].name += ( "_"+name_append);
            Isomer_ionspecis[NSpecies_Isomer].AME = AME_tmp;
            Isomer_ionspecis[NSpecies_Isomer].AME_err = AME_err_tmp;
            Isomer_ionspecis[NSpecies_Isomer].exc = exc_tmp;
            Isomer_ionspecis[NSpecies_Isomer].exc_err = exc_err_tmp;
        
            Isomer_ionspecis[NSpecies_Isomer].SetMass();    // total mass keV
            Isomer_ionspecis[NSpecies_Isomer].SetMvq_AME();
            Isomer_ionspecis[NSpecies_Isomer].SetName_form();
            
            for(int i=0;i< NSpecies ;i++)
            {
                if(Isomer_ionspecis[NSpecies_Isomer].A ==ionspecies[i].A &&Isomer_ionspecis[NSpecies_Isomer].Z ==ionspecies[i].Z)
                {
                    ionspecies[i].Isomer_n ++;
                }
            }
            NSpecies_Isomer++;
            if(NSpecies_Isomer>=Isomer_N_MAX){cout<<" \033[7m\033[31m  FATAL ERROR\033[0m"<<"NSpecies_Isomer > "<<Isomer_N_MAX<<endl<<endl;return; }

        }
        infile.close();
        
        cout<<" \033[33m------------------------------------  \033[0m "<<endl;
        cout<<" NSpecies_Isomer = "<<NSpecies_Isomer<<endl;
        cout<<" \033[33m------------------------------------  \033[0m "<<endl;
        if(NSpecies_Isomer>Isomer_N_MAX)cout<<"!!!!!!!!!!!!!!!!!"<<endl<<"NSpecies_Isomer > "<<Isomer_N_MAX<<endl<<endl<<endl;
        
        for(int i=0;i<NSpecies_Isomer  ;i++)
        {
            cout<<Isomer_ionspecis[i].name<<" "<<Isomer_ionspecis[i].A<<" "<<Isomer_ionspecis[i].Z<<fixed<<setprecision(1)
            <<" isomerAME mass excess = "<<Isomer_ionspecis[i].AME<<" +- "<<Isomer_ionspecis[i].AME_err<<" Mass= "<<Isomer_ionspecis[i].Mass
            <<" exc = "<<Isomer_ionspecis[i].exc<<" +- "<<Isomer_ionspecis[i].exc_err<<endl;
        }
        
    }
}

// 读取 nubase 2020
else 
{
    cout<<"\033[33m------ Begin  NUBASE2020_Build -------------------------  \033[0m "<<endl;
    NUBASE2020_Build(ISS_AME,"INPUT//nubase_4.mas20.txt");
    cout<<"\033[33m------   NUBASE2020_Build  completed -------------------------  \033[0m "<<endl;

    ofstream outfile_exc_state;
    if(OUTFILE_exc_state_ON)outfile_exc_state.open("OUTPUT//exc_state_INFO.txt");

    NSpecies_Isomer = 0;
    NSpecies_IAS =0 ;
    for(int i=0;i<NSpecies  ;i++)
        {ionspecies[i].Isomer_n=0; ionspecies[i].IAS_n=0;}
    for(int i=0;i<NSpecies  ;i++)
    {
        //======== 设置 AME mass ======
        int AME_i=ZN_AME[ionspecies[i].Z][ionspecies[i].A-ionspecies[i].Z];
        if(AME_i<0){cout<<" \033[7m\033[31m  FATAL ERROR\033[0m AME_i<0 "<<endl; return;}
        ionspecies[i].AME=ISS_AME[AME_i].AME;                    //质量中心值
        ionspecies[i].AME_err=ISS_AME[AME_i].AME_err;            //质量误差
        ionspecies[i].MassUnknown = ISS_AME[AME_i].MassUnknown;  //是否已知
        ionspecies[i].SetMass();         //use AME mass excess: IONSpecies.Mass(nuclear) = A*u + AME - Z*Me + ( 14.4381*pow(Z,2.39) + 1.55468*0.000001*pow(Z,5.35) )/1000.0;
        ionspecies[i].SetMvq_AME();      //Mvq_AME = Mass/(u*Z);  nuclear mvq
        
        //======== 设置 AME Isomer IAS======    
        for(int j=AME_i+1;j<ISS_AME[AME_i].next_nuclide_i  ;j++)  //遍历ISS_AME 中从基态之后开始的，后面紧跟着的几行激发态核素
        {
            if(OUTFILE_exc_state_ON)outfile_exc_state<<ISS_AME[j].str_line<<endl;
            TString info_tmp = ISS_AME[j].nuclide_info;
            if(info_tmp =="i"||info_tmp =="j")
            {
                IAS_ionspecis[NSpecies_IAS].A = ISS_AME[j].A;
                IAS_ionspecis[NSpecies_IAS].Z = ISS_AME[j].Z;
                IAS_ionspecis[NSpecies_IAS].name    = ISS_AME[j].name;
                IAS_ionspecis[NSpecies_IAS].AME     = ISS_AME[j].AME;
                IAS_ionspecis[NSpecies_IAS].AME_err = ISS_AME[j].AME_err;
            
                IAS_ionspecis[NSpecies_IAS].SetMass();    // total mass keV
                IAS_ionspecis[NSpecies_IAS].SetMvq_AME();
                IAS_ionspecis[NSpecies_IAS].SetName_form_exc(info_tmp);
               
                ionspecies[i].IAS_n++;
                if(Isomer_IAS_check_ON)cout<<" find IAS: "<<IAS_ionspecis[NSpecies_IAS].Aname<<" "<<IAS_ionspecis[NSpecies_IAS].AME<<" +- "<<IAS_ionspecis[NSpecies_IAS].AME_err<<endl;
                NSpecies_IAS++;
                if(NSpecies_IAS>=IAS_N_MAX){cout<<" \033[7m\033[31m  FATAL ERROR\033[0m"<<"NSpecies_IAS > "<<IAS_N_MAX<<endl<<endl;return; }
            }
            else if(info_tmp =="m"||info_tmp =="n")
            {
                Isomer_ionspecis[NSpecies_Isomer].A = ISS_AME[j].A;
                Isomer_ionspecis[NSpecies_Isomer].Z = ISS_AME[j].Z;
                Isomer_ionspecis[NSpecies_Isomer].name    = ISS_AME[j].name;
                Isomer_ionspecis[NSpecies_Isomer].AME     = ISS_AME[j].AME;
                Isomer_ionspecis[NSpecies_Isomer].AME_err = ISS_AME[j].AME_err;
                Isomer_ionspecis[NSpecies_Isomer].exc     = ISS_AME[j].exc;
                Isomer_ionspecis[NSpecies_Isomer].exc_err = ISS_AME[j].exc_err;
            
                Isomer_ionspecis[NSpecies_Isomer].SetMass();    // total mass keV
                Isomer_ionspecis[NSpecies_Isomer].SetMvq_AME();
                Isomer_ionspecis[NSpecies_Isomer].SetName_form_exc(info_tmp);
                
                ionspecies[i].Isomer_n ++;
                if(Isomer_IAS_check_ON)cout<<" find Isomer: "<<Isomer_ionspecis[NSpecies_Isomer].Aname<<" "<<Isomer_ionspecis[NSpecies_Isomer].AME<<" +- "<<Isomer_ionspecis[NSpecies_Isomer].AME_err<<endl;
                NSpecies_Isomer++;
                if(NSpecies_Isomer>=Isomer_N_MAX){cout<<" \033[7m\033[31m  FATAL ERROR\033[0m"<<"NSpecies_Isomer > "<<Isomer_N_MAX<<endl<<endl;return; }
            }
            else 
            {
                NSpecies_other_exc++;
                if(Isomer_IAS_check_ON)cout<<" \033[7m\033[31m  Be Careful ! ! ! ! nuclide_info is   \033[0m "<<info_tmp<<" at: "<<endl<<ISS_AME[j].str_line<<endl;
            }
        }
    }
    cout<<" find : NSpecies_Isomer = "<<NSpecies_Isomer<<endl
    <<" find : NSpecies_IAS = "<<NSpecies_IAS<<endl
    <<" find : NSpecies_other_exc = "<<NSpecies_other_exc<<endl;


    if(OUTFILE_exc_state_ON)
    {
        cout<<"---------- outfile_exc_state saved -------------"<<endl;
        outfile_exc_state.close();
    }
    
}// else读取 nubase 2020

//============  20230720 IsRef 设置参考核 ====================
if(THIS_EXP=="2017_58Ni"){IsRef_AME_ERR_MIN=5.0;}
if(THIS_EXP=="2021_36Ar_SET3"){IsRef_AME_ERR_MIN=50;}
if(THIS_EXP=="2021_36Ar_SET2"){IsRef_AME_ERR_MIN= 999 ;}
for(int i=0;i<NSpecies  ;i++)
{
    //通用规则
    //AME质量误差较小，质量已知，没有m、n态，（还有粒子数足够多，在后面）
    if(ionspecies[i].AME_err<IsRef_AME_ERR_MIN && !ionspecies[i].MassUnknown&& ionspecies[i].Isomer_n==0)
    {
        ionspecies[i].IsRef=1;
        //cout<<"debug!! IsRef:(not final choice) "<<ionspecies[i].Aname<<endl;
    }
    else {ionspecies[i].IsRef=0;}

    //@EXP@ 每次实验的单独设置
    if (THIS_EXP == "2021_36Ar_SET1")
    {
        //暂未确认
        if (ionspecies[i].Aname == "20Mg") { ionspecies[i].IsRef = 0; }//##ARTIFICIAL
        if (ionspecies[i].Aname == "21Mg") { ionspecies[i].IsRef = 0; }//##ARTIFICIAL
        if (ionspecies[i].Aname == "14O") { ionspecies[i].IsRef = 0; }//##ARTIFICIAL
        if (ionspecies[i].Aname == "22Mg") { ionspecies[i].IsRef = 0; }//##ARTIFICIAL

        //if (ionspecies[i].Aname == "11C") { ionspecies[i].IsRef = 1; }
    }
    if(THIS_EXP=="2021_36Ar_SET3")
    {
        if(ionspecies[i].Aname=="14O") {ionspecies[i].IsRef=0;}//##ARTIFICIAL
        if(ionspecies[i].Aname=="11C") {ionspecies[i].IsRef=1;}
        if(ionspecies[i].Aname=="10C") {ionspecies[i].IsRef=1;}
        if(ionspecies[i].Aname=="20Mg") {ionspecies[i].IsRef=1;}
        if(ionspecies[i].Aname=="13O") {ionspecies[i].IsRef=1;}
    }
    if(THIS_EXP=="2017_58Ni")
    {
        if(ionspecies[i].Aname=="14O") {ionspecies[i].IsRef=0;}//##ARTIFICIAL
        if(ionspecies[i].Aname=="49Fe") {ionspecies[i].IsRef=0;}//##ARTIFICIAL
    }
    
}

//=========================       output ionspecies info      ================================================
outfile.open("OUTPUT//OUTPUT_IonSpeciesInfo.txt");//检验读入，输出核素的相关信息
for(int i=0;i<NSpecies ;i++)
{
    outfile<<setw(5)<<ionspecies[i].Aname<<" A= "<<ionspecies[i].A<<" Z= "<<ionspecies[i].Z<<" Tz= "<<ionspecies[i].Tz
    <<fixed<<setprecision(2)
    <<" AMENucMass= "<<setw(10)<<ionspecies[i].Mass<<" AME ME= "<<ionspecies[i].AME<<" +- "<<ionspecies[i].AME_err
    <<" IAS_n: "<<ionspecies[i].IAS_n<<" Isomer_n: "<<ionspecies[i].Isomer_n<<endl;// 
}
outfile<<" ----------------- IAS ------------------"<<endl;
for(int i=0;i<NSpecies_IAS ;i++)
{
    outfile<<setw(5)<<IAS_ionspecis[i].Aname<<" A= "<<IAS_ionspecis[i].A<<" Z= "<<IAS_ionspecis[i].Z<<" Tz= "<<IAS_ionspecis[i].Tz
    <<fixed<<setprecision(2)
    <<" AMENucMass= "<<setw(10)<<IAS_ionspecis[i].Mass<<" AME ME= "<<IAS_ionspecis[i].AME<<" +- "<<IAS_ionspecis[i].AME_err
    <<fixed<<setprecision(8)<<" Mvq_AME= "<<IAS_ionspecis[i].Mvq_AME
    <<endl;
}

outfile<<" ----------------- Isomer ------------------"<<endl;
for(int i=0;i<NSpecies_Isomer ;i++)
{
    outfile<<setw(5)<<Isomer_ionspecis[i].Aname<<" A= "<<Isomer_ionspecis[i].A<<" Z= "<<Isomer_ionspecis[i].Z<<" Tz= "<<Isomer_ionspecis[i].Tz
    <<fixed<<setprecision(2)
    <<" AMENucMass= "<<setw(10)<<Isomer_ionspecis[i].Mass<<" AME ME= "<<Isomer_ionspecis[i].AME<<" +- "<<Isomer_ionspecis[i].AME_err
    <<fixed<<setprecision(8)<<" Mvq_AME= "<<Isomer_ionspecis[i].Mvq_AME
    <<fixed<<setprecision(2)<<" Exc_Energy= "<<Isomer_ionspecis[i].exc<<" +- "<<Isomer_ionspecis[i].exc_err
    <<endl;
}

outfile.close();
//==========================================================================================================================================

//cout<<" debug_____________________ check read in _____________"<<endl<<"  return after output ionspecies info"<<endl; return;


//=========================       read in ion merr      ================================================

ifstream infile_merr;
if(READIN_MERR_ON)//读入每个粒子的相关信息以及误差（事先输出好的）
{
infile_merr.open("INPUT//iont_errors_1_0629.txt");
int ions_unknown_ID_in=0;
while(infile_merr>>ions_unknown_ID_in>>strtmp>>strtmp>>strtmp>>strtmp>>ReadIn_m_v2[ions_unknown_ID_in]>>ReadIn_merr[ions_unknown_ID_in])
{}
infile_merr.close();

}
//_________________________       read in ion merr      _____________________________//首先弄好所需的TGraph




TGraphErrors* grerr0= new TGraphErrors();      //SigmaT vs AveT for all ionspecies
int grerr0_n=0;
TGraph* gr1= new TGraph();      //Bp C

//1021-dmdC

TGraph* gr_dmdC = new TGraph();  // each calculation of mass from one ref
TH2F* h2_dmdC = new TH2F("h2_dmdC","h2_dmdC",200,-0.5,0.5, 200,-1000,1000);

int gr_dmdCv2_n=0;
TGraph* gr_dmdCv2 = new TGraph();
int grerr_dmC_VE_all_n=0;
TGraphErrors* grerr_dmC_VE_all = new TGraphErrors();

TGraph* gr_dmISO = new TGraph();
TGraph* gr_TerrISO = new TGraph();
TGraph* gr_verrISO = new TGraph();
TH2F* h2_dmvqISO = new TH2F("h2_dmvqISO","h2_dmvqISO",400,-0.1,0.2,  200,-0.03*0.001,0.03*0.001);

TGraph* gr_gtC_tmp= new TGraph(); //γt C 的散点分布
int gr_gtC_tmp_n = 0;
TF1* fitfun_gr_gtC_all = new TF1("fitfun_gr_gtC_all","pol12",127,129);

TGraph* gr_dm_dA0err_all = new TGraph();
int gr_dm_dA0err_all_n=0;
AxisFormat(gr_dm_dA0err_all,"","","dm[keV]");gr_dm_dA0err_all->SetMarkerSize(1.0);
///////////////////////////////////////////////////////////////////////

int ions_this_injection;            //number of ions in this injection 本次注入的粒子数
int ref_n=0;                        //how many ref ions will be used 用到的参考核数
int j0 =0;                          //the starting subfix of ref ions 这一次注入的第一个粒子的序号
int j_n=0;                          //the subfix of ref ions that will be used in the new array 这一次注入中，各个参考核的序号
double Bp_err2_1=0.0;               //sigma_Bpi square calculated using one ref, produced by AME_err
double gammat_ave,Bp_ave_ref;       //gammat_ave 使用此参考核计算时，计算过程的平均γt  Bp_ave_ref 本次注入其他参考核的平均Bρ
double Ct,Ci,Ct_err,Ci_err,Bpi,vi,vi_err;                       //Ct :unknown target ion 未知核 Ci,Bpi: ref ion 参考核
double ref_m,ref_m_err;             //mass of ref ions 参考核质量，误差

int isomer_this_injection=0;
int isomer_appeared=1;                      //isomer species appeared
////////////////////////////////////////
ION_UNKNOWN ions_unknown;                //ion that is now being calculated  当前正在计算质量的离子
ION_UNKNOWN ions_t;                     
ION ions_r[ref_n_MAX];                          //array of reference ions for calculating the mass of one target ion   对应于一个目标核的参考离子数组  
//////////////////////////////////////

int grerr_n=0;

//------------------------------------------------------

int k_tmp=0;
double gt_inject_tmp=0;
TF1* f_gt_ave = new TF1 ("f_gt_ave","[0]",0,total_injection);


//======================== declare gtC steps ==================================
//20240711
int gtC_STEP_NOW = 0;
int gtC_step0_n  =0;   //original all input gt                                                                                    原始γt数据
int gtC_step1_n  =0;   // step1: gt_error filtering: filter out those ions whose gt_err are too large                             筛去γt误差过大的
int gtC_step2_n  =0;   // step2: species selection : choose certain species(heavy ions)                                           种类挑选（重核）
int gtC_step3_n  =0;   // step3: n_sigma squeeze   : filter out ions that is n_sigma away from central value of this subregion    筛去γt中心值偏离γt(C)曲线过远的
int gtC_step4_n  =0;   // step4: smooth            : weighted average according to the distance from one point                    平滑

double gtC_chosen[subregion_n];              //各周长区间的γt平均值（等权）
double gtC_chosen_v2[subregion_n];           //各周长区间的γt平均值（加权平均）
double gtC_chosen_v2_err[subregion_n];       //各周长区间的γt平均值的标准差（加权平均）
double gtC2_chosen[ subregion_n] ;           //各周长区间的γt平方平均值
double gtC_chosen_err[ subregion_n] ;      //各周长区间的γt平均值的误差（等权）=标准差/(sqrt(N-1))
int    C_Division_n_chosen[subregion_n];     //此周长区间的参考核数
//20240712 SR=subregion
double avegtSR_step1[subregion_n];  
double avegtSR_sigma_step1[subregion_n];
int C_n_SR_step1[subregion_n];  


//--------------------- gtC step 画图展示参数 ---------------------------------------------------
double gtC_step_window_Cmin = 128.5;double gtC_step_window_Cmax = 129.0; double gtC_step_bin_width_C = 0.002;
double gtC_step_window_gtmin = 1.2;double gtC_step_window_gtmax = 1.4;  double gtC_step_bin_width_gt = 0.001;
if(THIS_EXP=="2017_58Ni")
{
    gtC_step_window_Cmin = 128.7;
    gtC_step_window_Cmax = 129.0;
    gtC_step_bin_width_C = 0.002;
    gtC_step_window_gtmin = 1.3;
    gtC_step_window_gtmax = 1.4; 
    gtC_step_bin_width_gt = 0.001;
}
else if(THIS_EXP=="2021_36Ar_SET2")
{
    gtC_step_window_Cmin = 128.5;
    gtC_step_window_Cmax = 129.0;
    gtC_step_bin_width_C = 0.002;
    gtC_step_window_gtmin = 1.3;
    gtC_step_window_gtmax = 1.4; 
    gtC_step_bin_width_gt = 0.001;
}
else if(THIS_EXP=="2021_36Ar_SET3")
{
    gtC_step_window_Cmin = 128.5;
    gtC_step_window_Cmax = 129.0;
    gtC_step_bin_width_C = 0.002;
    gtC_step_window_gtmin = 1.3;
    gtC_step_window_gtmax = 1.4; 
    gtC_step_bin_width_gt = 0.001;
}
else
{
    gtC_step_window_Cmin = 128.5;
    gtC_step_window_Cmax = 129.0;
    gtC_step_bin_width_C = 0.002;
    gtC_step_window_gtmin = 1.2;
    gtC_step_window_gtmax = 1.4; 
    gtC_step_bin_width_gt = 0.001;
}
//step0
int gr_gt_gterr_step0_n=0;
TGraph* gr_gt_gterr_step0;
int gr_gtC_step0_n=0;
TGraph* gr_gtC_step0;   
TH1F* h1_gt_all_step0 ;
TH1F* h1_gterr_all_step0 ;
TH2F* h2_gtC_all_step0 ;
TH2F* h2_gt_gterr_step0 ;

//step1
int gr_gt_gterr_step1_n=1;
TGraph* gr_gt_gterr_step1;
int gr_gtC_step1_n=1;
TGraph* gr_gtC_step1;   
TH1F* h1_gt_all_step1 ;
TH1F* h1_gterr_all_step1 ;
TH2F* h2_gtC_all_step1 ;
TH2F* h2_gt_gterr_step1 ;

//step2
int gr_gt_gterr_step2_n=1;
TGraph* gr_gt_gterr_step2;
int gr_gtC_step2_n=1;
TGraph* gr_gtC_step2;   
TH1F* h1_gt_all_step2 ;
TH1F* h1_gterr_all_step2 ;
TH2F* h2_gtC_all_step2 ;
TH2F* h2_gt_gterr_step2 ;

//________________________ declare gtC step ______________________


//1017 L-ddT-intersectionC
TGraph2D*gr2d_LddtC = new TGraph2D();
// 
TGraph2D*gr2d_Lddt_Xn_v1 = new TGraph2D();

double C_intersection = 0;
// ==== save L-ddt-Chi2

if(!LOOP_ON) {L_n=1;ddT_n=1;scan_k_ddtC_n=1;scan_C_inter_n=1;gt0_n=1;}
scan_loop_total =  L_n*ddT_n*scan_k_ddtC_n*scan_C_inter_n*gt0_n;
if(scan_loop_total>MAX_SCAN_n){cout<<"\033[7m\033[31m  WARNING! scan_loop_total>MAX_SCAN_n !!!  \033[0m "<<endl;}
/*L_ddT_Xn2 = new double[scan_loop_total]; for(int i=0;i<scan_loop_total;i++){L_ddT_Xn2[i]=0;}
L_ddT_Xn2_v2 = new double[scan_loop_total];for(int i=0;i<scan_loop_total;i++){L_ddT_Xn2_v2[i]=0;}
L_ddT_Xn2_sys = new double[scan_loop_total]; for(int i=0;i<scan_loop_total;i++){L_ddT_Xn2_sys[i]=0;}
L_ddT_Xn2_v2_sys = new double[scan_loop_total];for(int i=0;i<scan_loop_total;i++){L_ddT_Xn2_v2_sys[i]=0;}
*/
outfile_Lddt_Xn2.open(FILEPATH+"L_ddT_Xn2.txt");//记录 L ddt 的扫描结果，取Xn最接近1的

//记录各种参数
outfile_logs<<"L_n = "<<L_n<<" ddT_n = "<<ddT_n<<endl<<"scan_k_ddtC_"<<endl
<<" L= "<<L_down<<" "<<dL<<endl
<<"ddT= "<<ddT_down<<" "<<dddT<<endl<<"scan_k_ddtC = "<<scan_k_ddtC<<" "<<scan_k_ddtC_up<<" "<<dscan_k_ddtC<<endl;
outfile_logs<<" use_gtC_type = "<<use_gtC_type<<endl;

ofstream outfile_large_dm;
outfile_large_dm.open(FILEPATH+"large_mass_deviation.txt");
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//================================== ionspecies[i] new =======================================
/////   ionspecies 里面的 TGraph TH1F 在每次循环都需要刷新 但是如果在循环内每次都new 的话， 会占用大量内存
//// ---- 为了节省内存开销，处理方法是： 在循环外 new， 然后每次循环 重置。 一次new，重复利用。
/// reset 重置方法：  h1->Reset() ,  TGraph 使用 RemovePoint() 从后往前删除
for(int i=0;i<NSpecies  ;i++)
{
    // BpC & lnBpC
    ionspecies[i].gr_BpC = new TGraph();    
    AxisFormat(ionspecies[i].gr_BpC, ionspecies[i].Aname+"gr_BpC","C [m]","Bp [Tm]",1+i%9);
    ionspecies[i].gr_BpC->SetMarkerStyle(20+int(i/9));
    ionspecies[i].gr_BpC->SetMarkerColor(1+i%9);
    ionspecies[i].gr_BpC->SetMarkerSize(1);

    ionspecies[i].gr_lnBpC = new TGraph();
    AxisFormat(ionspecies[i].gr_lnBpC,ionspecies[i].Aname+"gr_lnBpC","lnC [m]"," lnBp [Tm]",1+i%9);
    ionspecies[i].gr_lnBpC->SetMarkerStyle(20+int(i/9));
    ionspecies[i].gr_lnBpC->SetMarkerColor(1+i%9);
    ionspecies[i].gr_lnBpC->SetMarkerSize(1);

    ionspecies[i].fitfun_pol1_BpC= new TF1("fitfun_pol1_BpC","pol1",127,129);
    ionspecies[i].fitfun_pol1_lnBpC= new TF1("fitfun_pol1_lnBpC","pol1",4,5);
    ionspecies[i].fitfun_pol1_BpC->SetLineColor(1+i%9);
    ionspecies[i].fitfun_pol1_lnBpC->SetLineColor(1+i%9);
    ionspecies[i].fitfun_pol1_BpC->GetYaxis()->SetTitle("Bp [Tm]");

    //mvq_C
    ionspecies[i].gr_mvqC = new TGraph();    
    AxisFormat(ionspecies[i].gr_mvqC, ionspecies[i].Aname+" gr_mvqC","C [m]","m/q",1);
    ionspecies[i].gr_mvqC->SetMarkerStyle(20+int(i/9));
    ionspecies[i].gr_mvqC->SetMarkerColor(1+i%4);
    ionspecies[i].gr_mvqC->SetMarkerSize(1);
    //20230814 dm_C
    ionspecies[i].gr_dmC = new TGraph();    
    AxisFormat(ionspecies[i].gr_dmC, ionspecies[i].Aname+" gr_dmC","C [m]","#Deltam=M_{exp}-M_{AME}",1);
    ionspecies[i].gr_dmC->SetMarkerStyle(20+int(i/9));
    ionspecies[i].gr_dmC->SetMarkerColor(1+i%4);
    ionspecies[i].gr_dmC->SetMarkerSize(1);
    
    
    // new TH1F h_dm after F-D bin width is determined 
    //ionspecies[i].h_dm = new TH1F(ionspecies[i].Aname+" h_dm",ionspecies[i].Aname+" h_dm ", 100,-1000,1000);   //bin width =20 keV
    //AxisFormat(ionspecies[i].h_dm,ionspecies[i].Aname+" h_dm ","#Deltam=M_{exp}-M_{AME}", "Counts");
    

    ionspecies[i].f0_mvq_AME = new TF1("f0_mvq_AME","pol1",0,130);
    ionspecies[i].f0_mvq_AME->SetParameters(ionspecies[i].Mvq_AME,0);
    ionspecies[i].f0_mvq_AME->SetLineWidth(3);
    ionspecies[i].f0_mvq_AME->SetLineStyle(9);

    ionspecies[i].grerr_mvqC_VE = new TGraphErrors();    
    AxisFormat(ionspecies[i].grerr_mvqC_VE, ionspecies[i].Aname+" grerr_mvqC_VE","C [m]","m/q",1);
    ionspecies[i].grerr_mvqC_VE->SetMarkerStyle(24+int(i/9)%4 );
    ionspecies[i].grerr_mvqC_VE->SetMarkerColor(1+i%4);
    ionspecies[i].grerr_mvqC_VE->SetMarkerSize(1);
    ionspecies[i].grerr_mvqC_VE->SetLineWidth(1);

    ionspecies[i].grerr_dmC_VE = new TGraphErrors();    
    AxisFormat(ionspecies[i].grerr_dmC_VE, ionspecies[i].Aname+" grerr_dmC_VE","C [m]","#Deltam=M_{exp}-M_{AME}",1);
    ionspecies[i].grerr_dmC_VE->SetMarkerStyle(24+int(i/9)%4 );
    ionspecies[i].grerr_dmC_VE->SetMarkerColor(1+i%4);
    ionspecies[i].grerr_dmC_VE->SetMarkerSize(1);
    ionspecies[i].grerr_dmC_VE->SetLineWidth(1);
    

    ionspecies[i].fitfun_pol1_mvqC= new TF1("fitfun_pol1_mvqC","pol1",127,129);
    ionspecies[i].fitfun_pol1_mvqC->SetLineColor(1+i%9);
    ionspecies[i].fitfun_pol1_mvqC->SetLineWidth(6);
    ionspecies[i].fitfun_pol1_mvqC_VE= new TF1("fitfun_pol1_mvqC_VE","pol1",127,129);
    ionspecies[i].fitfun_pol1_mvqC_VE->SetLineColor(1+i%9);
    ionspecies[i].fitfun_pol1_mvqC_VE->SetLineWidth(6);


    ionspecies[i].grerr_avegtC= new TGraphErrors();
    ionspecies[i].gr_gtC_own = new TGraph();ionspecies[i].gr_gtC_own_u = new TGraph();ionspecies[i].gr_gtC_own_d = new TGraph(); 
    ionspecies[i].gr_gtC_shifted_own = new TGraph();ionspecies[i].gr_gtC_shifted_own_u = new TGraph();ionspecies[i].gr_gtC_shifted_own_d = new TGraph(); 
    //--0929 gtC for each ionspecies
    AxisFormat(ionspecies[i].grerr_avegtC,ionspecies[i].Aname,"C [m]"," #gamma_{t}(ave) ",my_root_color_10[i%10] );
    AxisFormat(ionspecies[i].gr_gtC_own,ionspecies[i].Aname,"C [m]"," #gamma_{t}(ave) ",my_root_color_10[i%10] );
    AxisFormat(ionspecies[i].gr_gtC_shifted_own,ionspecies[i].Aname,"C [m]"," #gamma_{t}(ave) ",my_root_color_10[i%10] ); //20240702
    ionspecies[i].grerr_avegtC->SetMarkerStyle(20+int(i/8));
    ionspecies[i].grerr_avegtC->SetMarkerSize(2);
    ionspecies[i].gr_gtC_shifted_own->SetMarkerStyle(20+int(i/8));
    ionspecies[i].gr_gtC_own->SetMarkerStyle(20+int(i/8));
    ionspecies[i].gr_gtC_own_u->SetMarkerStyle(20+int(i/8));
    ionspecies[i].gr_gtC_own_d->SetMarkerStyle(20+int(i/8));
    ionspecies[i].gr_gtC_shifted_own->SetMarkerSize(2);
    ionspecies[i].gr_gtC_own_u->SetMarkerSize(2);
    ionspecies[i].gr_gtC_own_d->SetMarkerSize(2);
    //ionspecies[i].ResetGtArrays(subregion_n);


    //---- gr_gt_C
    ionspecies[i].gr_gt_C  = new TGraph();
    AxisFormat(ionspecies[i].gr_gt_C,ionspecies[i].Aname+"  gr_gt_C"," C[m] ", "#gamma_{t}");
    ionspecies[i].gr_gt_C->SetMarkerSize(1.0);
    ionspecies[i].gr_gt_C->SetMarkerStyle(24);

    //h2_gtC   C bin width
    ionspecies[i].h2_gtC  = new TH2F(ionspecies[i].Aname+"  h2_gtC",ionspecies[i].Aname+"  h2_gtC",250,128.5,129.0, 200,1.28,1.48);
    AxisFormat(ionspecies[i].h2_gtC,ionspecies[i].Aname+"  h2_gtC"," C[m] ", "#gamma_{t}");
    //h_C
    ionspecies[i].h_C = new TH1F(ionspecies[i].Aname+" h_C",ionspecies[i].Aname+" h_C", 50,128.5,129.0); // 250 bin
    AxisFormat(ionspecies[i].h_C,ionspecies[i].Aname+" h_C","C [m]", "Counts");

    //h_T   T[ns] 的分bin： Tmin-0.1~ Tmax+0.1  bin width固定为 20ps
    ionspecies[i].h_T = new TH1F(ionspecies[i].Aname+" h_T",ionspecies[i].Aname+" h_T", 
        int((ionspecies[i].Tmax-ionspecies[i].Tmin+0.2)/0.02),ionspecies[i].Tmin-0.1,ionspecies[i].Tmax+0.1);
    AxisFormat(ionspecies[i].h_T,ionspecies[i].Aname+" h_T","T [ns]", "Counts");

    

    //20240902 histogram mvq_C dm_C
    // bin width: assume A~10 ~ 10^7 keV, one bin for 5 keV, then bin width = 5*mvq(~1.5)*10^-7 ~ 5*10^-7
    ionspecies[i].h_mvq = new TH1F(ionspecies[i].Aname+" h_mvq",ionspecies[i].Aname+" h_mvq ", 400,ionspecies[i].Mvq_AME-0.0001,ionspecies[i].Mvq_AME+0.0001);
    AxisFormat(ionspecies[i].h_mvq,ionspecies[i].Aname+" h_mvq ","#frac{m}{q} [u/e]", "Counts");


    //20231030 ERR_ANA          
    ionspecies[i].h_each_ref_cal_mass_err = new TH1F( "h_each_ref_cal_Merr_"+ionspecies[i].Aname,"h_each_ref_cal_Merr_"+ionspecies[i].Aname,500,0,1000 ); 
    AxisFormat(ionspecies[i].h_each_ref_cal_mass_err,"h_each_ref_cal_Merr_"+ionspecies[i].Aname," Mass Error [keV]","Counts");
    ionspecies[i].h_iont_mass_err = new TH1F( "h_iont_Merr_"+ionspecies[i].Aname,"h_iont_Merr_"+ionspecies[i].Aname,1000,0,2000 );
    AxisFormat(ionspecies[i].h_iont_mass_err,"h_iont_mass_err_"+ionspecies[i].Aname," Mass Error [keV]","Counts");
    ionspecies[i].h_iont_chi_n = new TH1F( "h_iont_chi_n_"+ionspecies[i].Aname,"h_iont_chi_n_"+ionspecies[i].Aname,200,0,10 );
    AxisFormat(ionspecies[i].h_iont_chi_n,"h_iont_chi_n_"+ionspecies[i].Aname," #chi_{N} ","Counts");
}


///////////////////////////////////////////////////////////// L , ddT,.. scan /////////////////////////////////////////////////////////////////////////////////

//=================================================================================================================================================


if(LOOP_ON)cout<<"\033[32m||LOOP ON "<<" ||\033[0m"<<endl;
else        cout<<"\033[32m||LOOP OFF "<<" ||\033[0m"<<endl;

if(!scan_k_ddtC_ON)scan_k_ddtC_n=1;   //不扫描则循环数就是1

if(use_gtC_type!=2)gt0_n=1;           //不扫描则循环数就是1

//===================================================== for LOOP =======================================================================================
for(int L_i=0;L_i<L_n;L_i++)                                                             //变动 L
    for(int ddT_i=0;ddT_i<ddT_n;ddT_i++)                                                 //变动 ddt
        for(int scan_k_ddtC_i=0;scan_k_ddtC_i<scan_k_ddtC_n;scan_k_ddtC_i++)
            for(int scan_C_inter_i = 0; scan_C_inter_i<scan_C_inter_n;scan_C_inter_i++)
                for(int gt0_i=0;gt0_i<gt0_n;gt0_i++)
{//====================L-ddT-scan for loops begin    


    L=L_down+L_i*dL;

    if(L_ddt_correlation_ON&&LOOP_ON)//最低卡方处 L和ddT近似一个线性关系，打开此开关可以  只在此线附近进行扫描，加快扫描速度
    {
        if(THIS_EXP=="2021_36Ar_SET2")
        {
            //ddT=5.0*(L-18.038)-0.001 -0.001+ddT_i*dddT;  // 22Si part1
            //ddT=5.0*(L-18.038)+0.007 -0.001+ddT_i*dddT;  //22Si part2
            //ddT=5.0*(L-18.034)-0.0146 -0.0020+ddT_i*dddT;//20230224 RUN2- <scan_alpha>
            //ddT=5.0*(L-18.040)+0.016 -0.0010+ddT_i*dddT;//20230224 RUN2- <scan_beta>
    
            //ddT=5.0*(L-18.052)+0.075 -0.0010+ddT_i*dddT;//20230224 RUN2- <scan_beta>
            ddT=5.0*(L-18.052)+0.072-0.0010+ddT_i*dddT;//20230224 RUN2- <scan_beta>
        }
        if(THIS_EXP=="2017_58Ni")
        {
            //ddT=5.0*(L-18.044)+0.1375+( ddT_i - ddT_n/2.0 )*dddT;
            ddT=5.0*(L-18.046)+0.1470+( ddT_i - ddT_n/2.0 )*dddT;
        }
        if(THIS_EXP=="2021_36Ar_SET3")
        {
            ddT=5.0*(L-18.057)+0.0963+( ddT_i - ddT_n/2.0 )*dddT;
        }
    }
    else { ddT=ddT_down+ddT_i*dddT;}

    if(scan_k_ddtC_ON)scan_k_ddtC= scan_k_ddtC_down+scan_k_ddtC_i*dscan_k_ddtC;

    scan_C_inter = scan_C_inter_down+scan_C_inter_i*dscan_C_inter;

    if(use_gtC_type==2)gt0 = gt0_down+gt0_i*dgt0;


scan_loop_i++;
if(!LOOP_ON&&scan_loop_i>1){cout<<endl<<"DRAW ON  while  scan_loop_i>1 !!"<<endl;return;}//循环时不能画图！退出

ThisParaInfo=ThisParaInfo.Format("L=%.4f (m), #Deltat_{d}=%.4f (ns)",L,ddT);//   L=... ddT=...


cout<<endl<<endl;
cout<<" \033[33m===============================================================================  \033[0m "<<endl
    <<"       ||| loop at \033[41m"<<scan_loop_i<< " / ( "<<scan_loop_total<<" )\033[0m"                      //循环到多少了
    <<fixed<<setprecision(4)
    <<" \033[32mL  = "<<L<<"\033[0m  "<<int(double(L_i+1)/double(L_n)*100)<<" % "                             //L到百分之多少了
    <<" \033[32mddT= "<<ddT<<"\033[0m  "<<int(double(ddT_i+1)/double(ddT_n)*100)<<"% ||| "<<endl;             //ddT到百分之多少了

    if(scan_k_ddtC_ON)cout<<" scan_k_ddtC = "<<scan_k_ddtC<<"  "<<int(double(scan_k_ddtC_i+1)/double(scan_k_ddtC_n)*100)<<" % "<<endl;//同上，输出这两个扫描到多少了
    if(scan_C_inter_ON)cout<<" scan_C_inter = "<<scan_C_inter<<"  "<<int(double(scan_C_inter_i+1)/double(scan_C_inter_n)*100)<<" % "<<endl;
    if(use_gtC_type==2)cout<<" gt_constant = "<<gt0<<"  "<<int(double(gt0_i+1)/double(gt0_n)*100)<<" % "<<endl;

    cout<<"ions_n = "<<ions_n<<endl<<"total_injection "<<total_injection<<endl
    <<" \033[33m===============================================================================  \033[0m "<<endl;

outfile_logs<<endl
    <<"loop at "<<scan_loop_i<< " / ( "<<L_n*ddT_n*scan_k_ddtC_n<<" )"                                          //日志文件里也存一份
    <<fixed<<setprecision(4)
    <<" L= "<<L<<"  "<<int(double(L_i+1)/double(L_n)*100)<<" % "
    <<" ddT= "<<ddT<<"  "<<int(double(ddT_i+1)/double(ddT_n)*100)<<"% "<<endl;

    if(scan_k_ddtC_ON)outfile_logs<<" scan_k_ddtC = "<<scan_k_ddtC<<"  "<<int(double(scan_k_ddtC_i+1)/double(scan_k_ddtC_n)*100)<<" % "<<endl;
    if(scan_C_inter_ON)outfile_logs<<" scan_C_inter = "<<scan_C_inter<<"  "<<int(double(scan_C_inter_i+1)/double(scan_C_inter_n)*100)<<" % "<<endl;
    if(use_gtC_type==2)outfile_logs<<" gt_constant = "<<gt0<<"  "<<int(double(gt0_i+1)/double(gt0_n)*100)<<" % "<<endl;
    
    outfile_logs<<"ions_n = "<<ions_n<<endl<<"total_injection "<<total_injection<<endl
    <<" ==========   "<<endl;

// =====================================  initialization  =====================================================

count_use_massknown_but_notREF =0;
count_gtC_masserr_condition_1 =0;
count_gtC_masserr_condition_2 =0;
for(int i=0;i<ions_n  ;i++)
{
    ions[i].m_VE=0;
    ions[i].m_VE_err=0;
    ions[i].mvq_v1=0;
}

//============================ ionspecis[i].new =================================
for(int i=0;i<NSpecies  ;i++)//每次循环，重置每个离子种类，相关的内容 
{

    ionspecies[i].Mass_cal=0;            ionspecies[i].Mass_cal_err=0;     //MASS_VER=1,total nuclear mass  [keV]
    ionspecies[i].Mass_cal_v2=0;         ionspecies[i].Mass_cal_err_v2=0;  //MASS_VER=2 fit histogram
    ionspecies[i].Mass_cal_VE=0;         ionspecies[i].Mass_cal_err_VE=0;  //MASS_VER>=3 err
    ionspecies[i].Mass_cal_err_VE_ave=0; ionspecies[i].Mass_cal_err_VE_sca=0;
    ionspecies[i].deltaMass=0;           ionspecies[i].deltaMass_err=0;
    ionspecies[i].deltaMass_v2=0;        ionspecies[i].deltaMass_err_v2=0;  
    ionspecies[i].deltaMass_VE=0;        ionspecies[i].deltaMass_err_VE=0;   
    ionspecies[i].MassExcess_cal=0;      ionspecies[i].MassExcess_cal_v2=0;ionspecies[i].MassExcess_cal_VE=0;
    ionspecies[i].stdDeviation_V1=0;     ionspecies[i].stdDeviation_VE=0;

    ionspecies[i].Mvq_cal=0;      ionspecies[i].Mvq_cal_err=0;
    ionspecies[i].Mvq_cal_v2=0;   ionspecies[i].Mvq_cal_err_v2=0;
    ionspecies[i].Mvq_cal_VE=0;   ionspecies[i].Mvq_cal_err_VE=0; 

    ionspecies[i].h_m_Gauss_fit_opt=0;

    ionspecies[i].vector_dm_v1.clear();

    //20240710
    if(outfile_each_mass_data_ON&&!LOOP_ON)
    {
        ionspecies[i].outfile_mass.open(FILEPATH+ionspecies[i].Aname+"_mass_data.txt");
    }    

}
Selected_n=0;


for(int i=0;i<subregion_n  ;i++)//将各种用到的数组初始化
{    
    gtC_chosen[i]=0.0;gtC2_chosen[i]=0.0;gtC_chosen_err[i]=0.0;C_Division_n_chosen[i]=0;
    gtC_chosen_v2[i]=0.0;gtC_chosen_v2_err[i]=0;
    avegtSR_step1[i]=0;avegtSR_sigma_step1[i]=0;C_n_SR_step1[i]=0;
}

for(int i=0;i<subregion_n ;i++){gtC[i]=0.0;gtC2[i]=0.0;Sigma_gtC[i]=0.0;C_Division_n[i]=0;}
//----0929
for(int j=0;j<TzN;j++)
    for(int i=0;i<subregion_n  ;i++)
        {    gtC_Tz[j][i]=0.0;gtC2_Tz[j][i]=0.0;sigma_gtC_Tz[j][i]=0.0;C_Division_n_Tz[j][i]=0;}
  

for (int i = 0; i < subregion_n * dsubregion_n; i++)
{
    ISF[i] = 0;    ISF_u[i] = 0;    ISF_d[i] = 0;
}

//                             这几项每次循环时需要进行初始化
Cmin = 999; Cmax = -1;
C_ave_all = 0.0;
gtC_chosen_CMAX = 0.0, gtC_chosen_CMIN = 99999.0;


h_Mvq->Reset();            //Reset清除直方图内容 在循环外NEW  在！LOOP_ON 里面才会使用，其实不用Reset
h_Mvq_AME->Reset();        // 在循环外NEW  在！LOOP_ON 里面才会使用，其实不用Reset
h_Mvq_14O21Mg_ReIdentify->Reset();

TGraph* gr_gtC_chosen=new TGraph();                    //γt(C)曲线，计算用
TGraphErrors* grerr_gtC_chosen_v2=new TGraphErrors();  // 20230707 gt err added 加权得出的γt(C)？
TGraphErrors* grerr_gtC_with_err=new TGraphErrors();   // 20230707 gt err added
int grerr_gtC_with_err_n=0;

//-----------gtC step0
if(!LOOP_ON&&SHOW_gtC_step0)//输出γt(C)曲线每步的变化
{
    gr_gt_gterr_step0_n=0;
    gr_gt_gterr_step0=new TGraph();  
    AxisFormat(gr_gt_gterr_step0,"gr_gt_gterr_step0 "+ThisParaInfo,"#gamma_{t}","#gamma_{t} error",kAzure+0);
    gr_gt_gterr_step0->SetMarkerSize(1);
    gr_gt_gterr_step0->SetMarkerStyle(24);
    
    gr_gtC_step0_n = 0;
    gr_gtC_step0 = new TGraph();
    AxisFormat(gr_gtC_step0,"gr_gtC_step0 "+ThisParaInfo," C [m] ","#gamma_{t}",kAzure+0);
    gr_gtC_step0->SetMarkerSize(1);
    gr_gtC_step0->SetMarkerStyle(24);
    
    h1_gt_all_step0 = new TH1F("h1_gt_all_step0","h1_gt_all_step0",3000,0,3);
    AxisFormat(h1_gt_all_step0,"h1_gt_all_step0","#gamma_{t}"," count ", kAzure+0);
    h1_gterr_all_step0 = new TH1F("h1_gterr_all_step0","h1_gterr_all_step0",3000,0,3);
    AxisFormat(h1_gterr_all_step0,"h1_gterr_all_step0","#gamma_{t} error"," count ",kAzure+0);
    h2_gtC_all_step0 = new TH2F ("h2_gtC_all_step0","h2_gtC_all_step0", 
        (gtC_step_window_Cmax-gtC_step_window_Cmin)/gtC_step_bin_width_C ,gtC_step_window_Cmin,gtC_step_window_Cmax, 
        (gtC_step_window_gtmax-gtC_step_window_gtmin)/gtC_step_bin_width_gt ,gtC_step_window_gtmin,gtC_step_window_gtmax);
    AxisFormat(h2_gtC_all_step0,"h2_gtC_all_step0"," C [m] ","#gamma_{t}(C)",kAzure+0);
    //AxisFormat(h2_gtC_all_step0,"h2_gtC_all_step0"," ","",kAzure+0);  //故意不要标题 节省空间
    //h2_gtC_all_step0->GetXaxis()->SetLabelSize(0.04);
    //h2_gtC_all_step0->GetYaxis()->SetLabelSize(0.04);

    h2_gt_gterr_step0 = new TH2F ("h2_gt_gterr_step0","h2_gt_gterr_step0", 
        (1.4-1.2)/gtC_step_bin_width_gt ,1.2,1.4, 
        200,0,0.1);
    AxisFormat(h2_gt_gterr_step0,"h2_gt_gterr_step0"," #gamma_{t} ","#gamma_{t} error",kAzure+0);

}

//-----------gtC step1
if(!LOOP_ON&&SHOW_gtC_step1)
{
    gr_gt_gterr_step1_n=0;
    gr_gt_gterr_step1=new TGraph();  
    AxisFormat(gr_gt_gterr_step1,"gr_gt_gterr_step1 "+ThisParaInfo,"#gamma_{t}","#gamma_{t} error",kAzure+1);
    gr_gt_gterr_step1->SetMarkerSize(1);
    gr_gt_gterr_step1->SetMarkerStyle(24);
    
    gr_gtC_step1_n = 0;
    gr_gtC_step1 = new TGraph();
    AxisFormat(gr_gtC_step1,"gr_gtC_step1 "+ThisParaInfo," C [m] ","#gamma_{t}",kAzure+1);
    gr_gtC_step1->SetMarkerSize(1);
    gr_gtC_step1->SetMarkerStyle(24);
    
    h1_gt_all_step1 = new TH1F("h1_gt_all_step1","h1_gt_all_step1",3000,0,3);
    AxisFormat(h1_gt_all_step1,"h1_gt_all_step1","#gamma_{t}"," count ", kAzure+1);
    h1_gterr_all_step1 = new TH1F("h1_gterr_all_step1","h1_gterr_all_step1",3000,0,3);
    AxisFormat(h1_gterr_all_step1,"h1_gterr_all_step1","#gamma_{t} error"," count ",kAzure+1);
    
    h2_gtC_all_step1 = new TH2F ("h2_gtC_all_step1","h2_gtC_all_step1", 
        (gtC_step_window_Cmax-gtC_step_window_Cmin)/gtC_step_bin_width_C ,gtC_step_window_Cmin,gtC_step_window_Cmax, 
        (gtC_step_window_gtmax-gtC_step_window_gtmin)/gtC_step_bin_width_gt ,gtC_step_window_gtmin,gtC_step_window_gtmax);
    //AxisFormat(h2_gtC_all_step1,"h2_gtC_all_step1"," C [m] ","#gamma_{t}(C)",kAzure+1);
    AxisFormat(h2_gtC_all_step1,"h2_gtC_all_step1"," ","",kAzure+1); //故意不要标题 节省空间
    h2_gtC_all_step1->GetXaxis()->SetLabelSize(0.08);
    h2_gtC_all_step1->GetYaxis()->SetLabelSize(0.08);
    
    h2_gt_gterr_step1 = new TH2F ("h2_gt_gterr_step1","h2_gt_gterr_step1", 
        (gtC_step_window_gtmax-gtC_step_window_gtmin)/gtC_step_bin_width_gt ,gtC_step_window_gtmin,gtC_step_window_gtmax, 
        100,0,gtC_ERR_upper_bound);
    AxisFormat(h2_gt_gterr_step1,"h2_gt_gterr_step1"," #gamma_{t} ","#gamma_{t} error",kAzure+1);
}

if(!LOOP_ON&&SHOW_gtC_step2)
{
    //-----------gtC step2
    gr_gt_gterr_step2_n=0;
    gr_gt_gterr_step2=new TGraph();  
    AxisFormat(gr_gt_gterr_step2,"gr_gt_gterr_step2 "+ThisParaInfo,"#gamma_{t}","#gamma_{t} error",kAzure+2);
    gr_gt_gterr_step2->SetMarkerSize(1);
    gr_gt_gterr_step2->SetMarkerStyle(24);
    
    gr_gtC_step2_n = 0;
    gr_gtC_step2 = new TGraph();
    AxisFormat(gr_gtC_step2,"gr_gtC_step2 "+ThisParaInfo," C [m] ","#gamma_{t}",kAzure+2);
    gr_gtC_step2->SetMarkerSize(1);
    gr_gtC_step2->SetMarkerStyle(24);
    
    h1_gt_all_step2 = new TH1F("h1_gt_all_step2","h1_gt_all_step2",3000,0,3);
    AxisFormat(h1_gt_all_step2,"h1_gt_all_step2","#gamma_{t}"," count ", kAzure+2);
    h1_gterr_all_step2 = new TH1F("h1_gterr_all_step2","h1_gterr_all_step2",100,0,gtC_ERR_upper_bound);
    AxisFormat(h1_gterr_all_step2,"h1_gterr_all_step2","#gamma_{t} error"," count ",kAzure+2);
    h2_gtC_all_step2 = new TH2F ("h2_gtC_all_step2","h2_gtC_all_step2", 
        (gtC_step_window_Cmax-gtC_step_window_Cmin)/gtC_step_bin_width_C ,gtC_step_window_Cmin,gtC_step_window_Cmax, 
        (gtC_step_window_gtmax-gtC_step_window_gtmin)/gtC_step_bin_width_gt ,gtC_step_window_gtmin,gtC_step_window_gtmax);
    AxisFormat(h2_gtC_all_step2,"h2_gtC_all_step2"," C [m] ","#gamma_{t}(C)",kAzure+1);
    h2_gt_gterr_step2 = new TH2F ("h2_gt_gterr_step2","h2_gt_gterr_step2", 
        (gtC_step_window_gtmax-gtC_step_window_gtmin)/gtC_step_bin_width_gt ,gtC_step_window_gtmin,gtC_step_window_gtmax, 
        100,0,gtC_ERR_upper_bound);
    AxisFormat(h2_gt_gterr_step2,"h2_gt_gterr_step2"," #gamma_{t} ","#gamma_{t} error",kAzure+2);
}


TGraphErrors* grerr_gtC_chosen=new TGraphErrors();
TGraphErrors* grerr_gtC_all=new TGraphErrors();int grerr_gtC_all_n=0;   //γt(C)散点
TGraphErrors* grerr_gtC_s=new TGraphErrors();int grerr_gtC_s_n=0;       //select
TGraphErrors* grerr_gtC_ave=new TGraphErrors();int grerr_gtC_ave_n=0;   //得到的γt平均值
TGraphErrors* grerr_gtC_v2=new TGraphErrors(); int grerr_gtC_v2_n = 0;

TH2F* h2_gtC_v2 = new TH2F ("h2_gtC_v2","h2_gtC_v2", 250,128.5,129.0, 100,1.3,1.4);
TF1* fitfun_gtC = new TF1("fitfun_gtC","pol1",0,140);
TF1* fitfun_gtC_s = new TF1("fitfun_gtC_s","pol1",0,140);
TF1* fitfun_gtC0_s = new TF1("fitfun_gtC0_s","pol0",0,140);

AxisFormat(grerr_gtC_chosen,"grerr_chosen"+ThisParaInfo,"C [m]","#gamma_{t}",4);
AxisFormat(gr_gtC_chosen,"gr_chosen"+ThisParaInfo,"C [m]","#gamma_{t}",4);
AxisFormat(grerr_gtC_chosen_v2,"grerr_chosen_v2"+ThisParaInfo,"C [m]","#gamma_{t}",kAzure);
AxisFormat(grerr_gtC_with_err,"grerr_gtC_with_err"+ThisParaInfo,"C [m]","#gamma_{t}",kAzure);
grerr_gtC_with_err->SetMarkerSize(1);

AxisFormat(grerr_gtC_all,"all "+ThisParaInfo,"C [m]","#gamma_{t}",3);
AxisFormat(grerr_gtC_s,"all_Selected "+ThisParaInfo,"C [m]","#gamma_{t}",3);
grerr_gtC_all->SetMarkerSize(0.4);
grerr_gtC_s->SetMarkerSize(1);
AxisFormat(grerr_gtC_ave,"grerr_gtC_ave","C [m]","#gamma_{t}",6);


//_____0928
if(DoShow_mvqC_each_h2_ON&&!LOOP_ON)//对每种核 分别设置 h2_mvq_C_的区间（未开启）
{
    ifstream infile_mvq;
    infile_mvq.open("INPUT//36Ar_22Si_each_avemvq_part2.txt");
    if(!infile_mvq)cout<<"infile_mvq failed to open  INPUT//36Ar_22Si_each_avemvq_part2.txt"<<endl<<endl;
    double mvq0;
    int mvqin_i=0;
    while(infile_mvq>>mvq0)
    {
        ionspecies[mvqin_i].h2_mvqC = new TH2F("h2_mvqC_"+ionspecies[mvqin_i].Aname,"h2_mvqC_"+ionspecies[mvqin_i].Aname,200,128.5,129.0, 500, mvq0-0.0003,mvq0+0.0003);
        AxisFormat(ionspecies[mvqin_i].h2_mvqC,"h2_mvqC_"+ionspecies[mvqin_i].Aname,"C [m]", "#frac{m}{q} [u/e]");

        mvqin_i++;
    } 
}
if(Do_iont_error_out_ON&&MASS_VER>2)//非等权算法，可以给出单粒子的误差，因此可以输出单粒子的误差文件
{ outfile_errors.open(FILEPATH+strtmp.Format("iont_errors_%d.txt",scan_loop_i ) );
outfile_errors_check.open(FILEPATH+strtmp.Format("iont_errors_check_%d.txt",scan_loop_i ) );
}
//______________________________________ initialization ________________________________________________  初始化结束



//================================================== obtain v c  gammat 
outfile_C_Division_n.open(FILEPATH+strtmp.Format("C_Division_n_%d.txt",scan_loop_i) );

if(OUTPUT_check_abnormal_gammat)outfile.open("OUTPUT//check_gammat_abnormal.txt");//可以检查异常的γt
if(OUTFILE_each_ion_complete_ON&&!LOOP_ON)
{
    cout<<"  ---------------------  OUTFILE_each_ion_complete_ON --------------------"<<endl; 
    outfile_each_ion_complete.open("OUTPUT//each_ion_complete.txt");
}

//新增部分，计算 T v Bρ C误差平方以及协方差，最终计算每个核的影响因子    cal计算值  Sigma2方差 Cov协方差

    //思路：从 Cr实际值 -> Cr中心值 ----------------------> Ct中心值 -> Ct实际值
//误差影响：            dFr           γt(C)曲线的误差               dFt
//这里默认 Cr < Ct，反过来也差不多

double cal_T = 0, sigma2_T = 0, cal_v = 0, sigma2_v = 0, cov_T_v = 0, cal_gamma = 0;//T周期 v速度 cal计算出的中心值 sigma2方差 cov协方差
double cal_Brou = 0, sigma2_Brou = 0, cal_C = 0, sigma2_C = 0, cov_Brou_C = 0;//Brou Bρ   C周长
double cov_C_v = 0;
//2021 36Ar 的数据 对应于Nm =0 
for(int i=0;i<ions_n;i++)//计算每个粒子
{
    //cal_T     = ions[i].A1 - 0.5*ions[i].dA1;                        //ps Nm =0  
    //if(i%2000==0)cout<<cal_T<<" test-- ";
    cal_T     = ions[i].Get_T_Nm();   //Nm !=0 
    //if(i%2000==0)cout<<cal_T<<" -- | ";
    //sigma2_T  = ions[i].A1err + 0.25*ions[i].dA1err - ions[i].cov16; // T方差 ps^2
    //if(i%2000==0)cout<<sigma2_T<<" test-- ";
    sigma2_T  = ions[i].Get_Terr2_Nm(); //Nm !=0 
    //if(i%2000==0)cout<<sigma2_T<<" -- | ";

    //cal_v     = L / (ions[i].dA0 + ddT*1000);                        // m/ps Nm =0 
    //if(i%2000==0)cout<<cal_v<<" test-- ";
    cal_v     = ions[i].Get_v_Nm(L,ddT); //Nm !=0 
    //if(i%2000==0)cout<<cal_v<<" -- | ";

    //sigma2_v  = L * L / ( pow((ions[i].dA0 + ddT*1000),4) ) * ions[i].dA0err;// v 方差(m/ps)^2 Nm =0 
    //if(i%2000==0)cout<<sigma2_v<<" test-- ";
    sigma2_v   = ions[i].Get_verr2_Nm(L,ddT); //Nm !=0 
    //if(i%2000==0)cout<<sigma2_v<<" -- | ";


    cal_gamma = 1 / sqrt( 1 - pow((cal_v*1000),2)/pow(V_c,2) );
    
    //cov_T_v   = cal_v*cal_v/L * (-ions[i].cov15 + 0.5*ions[i].cov56);//ps * m/ps //Nm =0
    //if(i%2000==0)cout<<cov_T_v<<" test-- "; 
    cov_T_v   = ions[i].Get_Cov_T_v_Nm(L,ddT); //Nm !=0 
    //if(i%2000==0)cout<<cov_T_v<<" test-- "<<endl; 

    cal_Brou = (ionspecies[ions[i].Species].Mass / ions[i].Z) * cal_v * cal_gamma / V_c / V_c / 1000;//Tm   //M: keV  cal_v:m/ps  V_c:m/ns    
    sigma2_Brou = pow((cal_Brou / ionspecies[ions[i].Species].Mass), 2) * pow(ionspecies[ions[i].Species].AME_err, 2) + cal_Brou * cal_Brou * pow(cal_gamma, 4) / cal_v / cal_v * sigma2_v;//Tm^2 这里考虑了质量的误差
    cal_C = cal_v * cal_T;//m
    sigma2_C = cal_v * cal_v * sigma2_T + cal_T * cal_T * sigma2_v + 2 * cal_C * cov_T_v;//m^2
    
    cov_Brou_C = cal_Brou * cal_gamma*cal_gamma*( cal_T/cal_v*sigma2_v + cov_T_v );//Tm * m
    
    cov_C_v = cal_T * sigma2_v + cal_v * cov_T_v;//m*m/ps

    // set ions :
    ions[i].T=cal_T/1000.0;  // ions.v [ns]
    ions[i].T_err = sqrt(sigma2_T)/1000;  //[ns]
    ions[i].v=cal_v*1000.0;  // ions.v [m/ns]
    ions[i].v_err = sqrt(sigma2_v)*1000;  //[m/ns]
    ions[i].v_err2 = sigma2_v;  //[ m/ps*m/ps]
    ions[i].C= cal_C;
    ions[i].C_err2 = sigma2_C;
    ions[i].Bp = cal_Brou;
    ions[i].Bp_err2 = sigma2_Brou;
    ions[i].C_err2 = sigma2_C;
    ions[i].cov_Brou_C = cov_Brou_C;
    ions[i].cov_C_v = cov_C_v;

    C_ave_all+=ions[i].C;
    if(ions[i].C>Cmax)Cmax=ions[i].C;
    if(ions[i].C<Cmin)Cmin=ions[i].C;

    //############################ 单个离子 gt #######################################################################

    //ions[i].Calculate_gt_from_a_with_err(L,ddT);
    ions[i].Calculate_gt_from_a_with_err_58Ni_Nm(L,ddT,i);
    //if(i%500==0){cout<<"debug: "<<ions[i].gammat_err<<endl;}
    
    //############################ 单个离子 gt #######################################################################
    
    gammat = ions[i].gammat;
    if(gammat>gt_max) gt_max=gammat;
    if(gammat<gt_min) gt_min=gammat;
    ions[i].gt_abnormal=false;
    //------------- abnormal gt ???历史遗留 2024待检查--------------------------------
    if( (1 - ions[i].A2*2/ions[i].T*((ions[i].dA0/1000.0)+ddT)/ions[i].dA1 ) <0 )
    { 
        ions[i].gt_abnormal=true;
        if(OUTPUT_check_abnormal_gammat)
        {
            outfile<<ions[i].gammat<<"  "<<ions[i].A2*2/ions[i].T*(ions[i].dA0+ddT)/ions[i].dA1
            <<" v = "<<ions[i].v<<" +- "<<ions[i].v_err<<" C = "<<ions[i].C
            <<"  2A2 = "<<ions[i].A2*2<<" dA0= "<<ions[i].dA0<<" dA0err^2 "<<ions[i].dA0err<<" dA1= "<<ions[i].dA1<<endl;
        }
    }
    
    //########################## get Bp of all ions ###############################
    //ions[i].Calculate_Bp(ionspecies[ions[i].Species].Mass);   //AME mass  of this species
    //gt v2
    //ions[i].Calculate_gt_v2(L,ddT,ionspecies[ions[i].Species].Mass,i);
    //if(i%200==0)cout<<i<<" "<<ions[i].A<<ions[i].name<<" "<<ions[i].gammat<<" "<<ions[i].gammat_v2<<" "<<ions[i].gammat_v2-ions[i].gammat<<endl;///debug
    
    if(OUTFILE_each_ion_complete_ON&&!LOOP_ON){ ions[i].PrintInfo_complete(outfile_each_ion_complete);}

}
if(OUTPUT_check_abnormal_gammat)outfile.close();
if(OUTFILE_each_ion_complete_ON&&!LOOP_ON){ outfile_each_ion_complete.close();}

C_ave_all /= ions_n;
cout<<fixed<<setprecision(6)<<"\033[35m=== for all ions:Cmin="<<Cmin<<"======== Cmax= "<<Cmax<<"====C_ave_all : "<<C_ave_all<<"\033[0m"<<endl;
cout<<"\033[34m gt_min = "<< gt_min<<" , gt_max= "<< gt_max<<"\033[0m"<<endl;





//==========================  complete ionspecies   ========================================
for(int i=0;i<ions_n  ;i++)
{
    ionspecies[ions[i].Species].Record(ions[i]);   //对 T C Bρ γt 总粒子数进行统计
    //0928 add gr_Bp_C gr_lnBpC
    //if(Do_Fit_BpC_ON)ionspecies[ions[i].Species].Record_ion_BpC(ions[i]);   

}
for(int i=0;i<NSpecies  ;i++)
{
    ionspecies[i].CalculateAve(); //determine average T C v Bp gammat 计算平均值
}

//with ave known, calculate sigma T C v Bp gammat
for(int i=0;i<ions_n ;i++)
{
    ionspecies[ions[i].Species].SigmaPreAdd(ions[i]);          //计算T v γ γt标准差                            这一步统计，下一步计算
    ionspecies[ions[i].Species].SkewnessPreAdd(ions[i]);       //推测 可能是计算 C 的 偏度系数（为什么只算C--- 曾经尝试过 历史遗留 保留这种偏度算法功能）
}
for(int i=0;i<NSpecies  ;i++)
{
    ionspecies[i].CalculateSigma();
    ionspecies[i].CalculateSkewness();
}
//___________________ complete ionspecies ____________________________

//-------------------------- make Injection_subfix[]
// work in the order of ions_n, use Injection to determine the subfix that will be used within one injection
Injection_subfix[0]=0;
for(int i=0;i<ions_n  ;i++)
{
    if(ions[i].ion_number==1){Injection_subfix[ions[i].inject_number] = i;}//这一次注入的第一个粒子，记录下来这个粒子的 序号   Injection_subfix[i]:每次注入的第一个离子的序号 
}

//.....ooo00000oooo..........ooo00000oooo..........ooo00000oooo..........ooo00000oooo..........ooo00000oooo.....
if(OUTPUT_allions)outfile.open("OUTPUT//IonInfoOutput.txt");//输出所有粒子
for(int i=0;i<ions_n  ;i++)
{
    if(OUTPUT_allions)
    {
        outfile<<ions[i].name<<" "<<ions[i].Z<<" T= "<<fixed<<setprecision(10)<<ions[i].T<<" +- "<<ions[i].T_err
    <<" v= "<<fixed<<setprecision(6)<<ions[i].v<<" +- "<<ions[i].v_err
    <<" C= "<<fixed<<setprecision(10)<<ions[i].C<<" +- "<<sqrt(ions[i].C_err2)<<" gt= "<<ions[i].gammat
    <<" Bp= "<<ions[i].Bp<<" "<<ions[i].ion_number<<" "<<ions[i].inject_number<<endl;
    }
}    
if(OUTPUT_allions)outfile.close();


//=========================       output ionspecies info in loop     ================================================
if(!LOOP_ON)
{
outfile.open("OUTPUT//OUTPUT_IonSpeciesInfo_36Ar_inloop_loopOFF.txt");//没开循环，准备输出文件
for(int i=0;i<NSpecies ;i++)
{
    ionspecies[i].PrintInfo(outfile);
}
outfile.close();
// svae to folder 

outfile.open(FILEPATH+"OUTPUT_IonSpeciesInfo_36Ar_inloop_loopOFF.txt");
for(int i=0;i<NSpecies ;i++)
{
    ionspecies[i].PrintInfo(outfile);
}
outfile.close();
}
//________________________________________________________________________________________________________
//cout<<" output ionspecies info in loop check ----- return at 5707"<<endl; return;



//========================= each ionspecies statistics =================================
 //20240709 再整理 在这里先进行Fill, 后面 canvas draw
// initialization

for(int i=0;i<NSpecies  ;i++)//分类别统计。初始化
{
    ionspecies[i].Tmin = 999999;ionspecies[i].Tmax = -999999;
    ionspecies[i].Cmin = 999999;ionspecies[i].Cmax = -999999;
    ionspecies[i].vmin = 999999;ionspecies[i].vmax = -999999;
    ionspecies[i].Bpmin = 999999;ionspecies[i].Bpmax = -999999;
    ionspecies[i].gtmin = 999999;ionspecies[i].gtmax = -999999;
}
for(int j=0;j<ions_n  ;j++)//统计数据的上下界
{
    int id=ions[j].Species;
    if(ionspecies[id].Tmin>ions[j].T){ionspecies[id].Tmin=ions[j].T;}if(ionspecies[id].Tmax<ions[j].T){ionspecies[id].Tmax=ions[j].T;}
    if(ionspecies[id].Cmin>ions[j].C){ionspecies[id].Cmin=ions[j].C;}if(ionspecies[id].Cmax<ions[j].C){ionspecies[id].Cmax=ions[j].C;}
    if(ionspecies[id].vmin>ions[j].v){ionspecies[id].vmin=ions[j].v;}if(ionspecies[id].vmax<ions[j].v){ionspecies[id].vmax=ions[j].v;}
    if(ionspecies[id].Bpmin>ions[j].Bp){ionspecies[id].Bpmin=ions[j].Bp;}if(ionspecies[id].Bpmax<ions[j].Bp){ionspecies[id].Bpmax=ions[j].Bp;}
    if(ionspecies[id].gtmin>ions[j].gammat){ionspecies[id].gtmin=ions[j].gammat;}if(ionspecies[id].gtmax<ions[j].gammat){ionspecies[id].gtmax=ions[j].gammat;}
}
/*debug
for(int i=0;i<NSpecies  ;i++)
{
    cout<<" debug! "<<ionspecies[i].Aname
    <<" "<<ionspecies[i].Tmin<<" "<<ionspecies[i].Tmax<<" "<<ionspecies[i].Cmin<<" "<<ionspecies[i].Cmax
    <<" "<<ionspecies[i].vmin<<" "<<ionspecies[i].vmax<<" "<<ionspecies[i].Bpmin<<" "<<ionspecies[i].Bpmax<<" "<<ionspecies[i].gtmin<<" "<<ionspecies[i].gtmax<<endl;
}*/


//Fill data
for(int j=0;j<ions_n  ;j++)//统计 γt-C   T
{
    int id = ions[j].Species;
    ionspecies[id].h2_gtC->Fill(ions[j].C,ions[j].gammat);
    ionspecies[id].h_C->Fill(ions[j].C);
    ionspecies[id].gr_gt_C->SetPoint(ionspecies[id].gr_gt_C->GetN() , ions[j].C, ions[j].gammat);
    ionspecies[id].h_T->Fill(ions[j].T);
}
// Draw will be in !LOOP_ON canvas
//______________________________ each ionspecies statistics


//20240702
double ISS_ave_gt=0;
for(int i=0;i<NSpecies  ;i++)//计算总的γt平均值
{
    //cout<<ionspecies[i].Aname<<"20240702---- ave gt = "<<ionspecies[i].Avegammat<<endl;
    ISS_ave_gt+=ionspecies[i].Avegammat;
}
ISS_ave_gt /= NSpecies;
cout<<"20240702---- ISS ALL ave gt = "<<ISS_ave_gt<<endl;
for(int i=0;i<NSpecies  ;i++)
{
    //if(ionspecies[i].Aname=="18Ne")ionspecies[i].GT_SHIFT_own = ionspecies[i].Avegammat - ISS_ave_gt   ;
    //if(ionspecies[i].Aname=="17Ne")ionspecies[i].GT_SHIFT_own = ionspecies[i].Avegammat - ISS_ave_gt +0.012  ;
    //if(ionspecies[i].Aname=="11C") ionspecies[i].GT_SHIFT_own = ionspecies[i].Avegammat - ISS_ave_gt +0.013;
    //if(ionspecies[i].Aname=="13N") ionspecies[i].GT_SHIFT_own = 0.002;
    
    //cout<<"20240702---- "<<ionspecies[i].Aname<<" GT_SHIFT_own = "<<ionspecies[i].GT_SHIFT_own<<endl;
}

//20240625 ------------- IsRef-2 -------------------
cout<<endl<<"\033[33m ===========  IsRef: final: ===========   \033[0m"<<endl;
outfile_logs<<"===========  IsRef: final: ==========="<<endl;

if(THIS_EXP=="2017_58Ni"){IsRef_MIN_N = 100;}
if(THIS_EXP=="2021_36Ar_SET2"){IsRef_MIN_N = 50;}
if(THIS_EXP=="2021_36Ar_SET3"){IsRef_MIN_N = 40;}

cout<<" RefIons should have counts > IsRef_MIN_N = "<<IsRef_MIN_N<<endl;
int IsRef_FINAL_N = 0;  
for(int i=0;i<NSpecies  ;i++)//参考核设置
{
    if(ionspecies[i].IsRef)//选择参考核
    {
        if(ionspecies[i].N>=IsRef_MIN_N)//初始粒子数足够
        {
            cout<<"   "<<ionspecies[i].Aname;
            outfile_logs<<" | "<<ionspecies[i].Aname<<" "<<ionspecies[i].A<<" "<<ionspecies[i].Z ;
            IsRef_FINAL_N++;
        }
        else{ionspecies[i].IsRef=0;}//粒子数不够多，不用
    }   
}
cout<<endl<<" choose ref species numer IsRef_FINAL_N = "<<IsRef_FINAL_N<<endl; 
outfile_logs<<endl<<" choose ref species numer IsRef_FINAL_N = "<<IsRef_FINAL_N<<endl; 


/*

//=============================== 20230531 Do_dA0T ================================
//if(Draw_each_T_C_ON)Draw_each_T_C(ionspecies);

double k_L, k_C, kgt, k_dLdC, k_dtddC,ktd = 0;
double chi_best = 0;
TGraph* gr_gtC_from_dA0T = new TGraph();
AxisFormat(gr_gtC_from_dA0T,"","C[m]","#gamma_{t}",kSpring);

TGraph* gr_LC_from_dA0T = new TGraph();
AxisFormat(gr_LC_from_dA0T,"","C[m]","L[m]",kAzure);

TGraph* gr_ddtC_from_dA0T = new TGraph();
AxisFormat(gr_ddtC_from_dA0T,"","C[m]","#Delta t_{d} [ns]",kAzure);

TGraph* gr_k1C_from_dA0T = new TGraph();
AxisFormat(gr_k1C_from_dA0T,"","C[m]","k1: dL/dC",1);

TGraph* gr_k2C_from_dA0T = new TGraph();
AxisFormat(gr_k2C_from_dA0T,"","C[m]","k2 d#Delta t_{d}/dC",1);

TGraph* gr_chi_from_dA0T = new TGraph();
AxisFormat(gr_chi_from_dA0T,"","C[m]","#chi_{n}",1);
TGraph* gr_CC_from_dA0T = new TGraph();
AxisFormat(gr_CC_from_dA0T,"","previous C[m] "," new C[m]",1);
//============================================ Do_dA0_T_ON ========================================
if(Do_dA0_T_ON)
{
    cout<<" -----------------Do_dA0_T_ON--------------------"<<endl;

int Do_dA0_T_times = 120; 

ii=0;
for(int i=0;i<Do_dA0_T_times;i++)
{
    double C_step=0.005;
    double dC_range = 0.05;
    double C_center = 128.808-0.3+C_step*i;  //step=0.01 dC = 0.005

    if(C_center<=128.95&&C_center>=128.55)
    {
//############################################
    Do_dA0_T(ionspecies, k_L, k_C, kgt, k_dLdC, k_dtddC,ktd,chi_best,C_center,dC_range , i);
//##########################################
    cout<<"k_L= "<<k_L<<" ktd= "<<ktd<<" k_C = "<<k_C<<" kgt= "<<kgt<<" k_dLdC= "<<k_dLdC<<" k_dtddC= "<<k_dtddC<<endl;
    gr_gtC_from_dA0T->SetPoint(ii,k_C,kgt );
    gr_LC_from_dA0T->SetPoint(ii,k_C,k_L );
    gr_ddtC_from_dA0T->SetPoint(ii,k_C,ktd );
    gr_k1C_from_dA0T->SetPoint(ii,k_C,k_dLdC );
    gr_k2C_from_dA0T->SetPoint(ii,k_C,k_dtddC );
    gr_chi_from_dA0T->SetPoint(ii,k_C,chi_best);
    gr_CC_from_dA0T->SetPoint(ii,C_center,k_C);
    ii++;
    }
}

if( 0 )
{///canvas draw
TCanvas* c_gr_gtC_from_dA0T = new TCanvas("c_gr_gtC_from_dA0T","c_gr_gtC_from_dA0T",1200,600);
gr_gtC_from_dA0T->Draw("ap");
c_gr_gtC_from_dA0T->Print(FILEPATH_DodA0T+"gr_gtC_from_dA0T.png");
TCanvas* c_gr_LC_from_dA0T = new TCanvas("c_gr_LC_from_dA0T","c_gr_LC_from_dA0T",1200,600);
gr_LC_from_dA0T->Draw("ap");
c_gr_LC_from_dA0T->Print(FILEPATH_DodA0T+"gr_LC_from_dA0T.png");
TCanvas* c_gr_ddtC_from_dA0T = new TCanvas("c_gr_ddtC_from_dA0T","c_gr_ddtC_from_dA0T",1200,600);
gr_ddtC_from_dA0T->Draw("ap");
c_gr_ddtC_from_dA0T->Print(FILEPATH_DodA0T+"gr_ddtC_from_dA0T.png");
TCanvas* c_gr_k1C_from_dA0T = new TCanvas("c_gr_k1C_from_dA0T","c_gr_k1C_from_dA0T",1200,600);
gr_k1C_from_dA0T->Draw("ap");
c_gr_k1C_from_dA0T->Print(FILEPATH_DodA0T+"gr_k1C_from_dA0T.png");
TCanvas* c_gr_k2C_from_dA0T = new TCanvas("c_gr_k2C_from_dA0T","c_gr_k2C_from_dA0T",1200,600);
gr_k2C_from_dA0T->Draw("ap");
c_gr_k2C_from_dA0T->Print(FILEPATH_DodA0T+"gr_k2C_from_dA0T.png");
TCanvas* c_gr_chi_from_dA0T = new TCanvas("c_gr_chi_from_dA0T","c_gr_chi_from_dA0T",1200,600);
gr_chi_from_dA0T->Draw("ap");
c_gr_chi_from_dA0T->Print(FILEPATH_DodA0T+"gr_chi_from_dA0T.png");
TCanvas* c_gr_CC_from_dA0T = new TCanvas("c_gr_CC_from_dA0T","c_gr_CC_from_dA0T",1200,600);
gr_CC_from_dA0T->Draw("ap");
c_gr_CC_from_dA0T->Print(FILEPATH_DodA0T+"gr_CC_from_dA0T.png");
}///canvas draw

TGraph_to_outfile( FILEPATH_DodA0T+"gr_gtC_from_dA0T_save.txt", gr_gtC_from_dA0T );
TGraph_to_outfile( FILEPATH_DodA0T+"gr_LC_from_dA0T_save.txt",  gr_LC_from_dA0T  );
TGraph_to_outfile( FILEPATH_DodA0T+"gr_ddtC_from_dA0T_save.txt",  gr_ddtC_from_dA0T  );
TGraph_to_outfile( FILEPATH_DodA0T+"gr_k1C_from_dA0T_save.txt",  gr_k1C_from_dA0T  );
TGraph_to_outfile( FILEPATH_DodA0T+"gr_k2C_from_dA0T_save.txt",  gr_k2C_from_dA0T  );
TGraph_to_outfile( FILEPATH_DodA0T+"gr_CC_from_dA0T_save.txt",  gr_CC_from_dA0T  );

}
//___________________________________  Do_dA0_T_ON end _______________________________________

//========================= 2023 0620 recalculate v C gt after Do_dA0_T ====================================
if(Recalculate_after_dA0T_ON)
{
cout<<"========================= recalculate v C gt after Do_dA0_T ===================================="<<endl;

    Recalculate_after_dA0T();

    C_ave_all /= ions_n;
    cout<<fixed<<setprecision(6)<<"\033[35m=== Cmin="<<Cmin<<"======== Cmax= "<<Cmax<<"====C_ave_all : "<<C_ave_all<<"\033[0m"<<endl;

}
/////////__________________________________ 
*/

if(Draw_each_T_C_ON)Draw_each_T_C(ionspecies);


/*
//1017 do fit after Record
if(Do_Fit_BpC_ON)
{
for(int i=0;i<NSpecies  ;i++)
{
    if(ionspecies[i].N>10)
    {
        ionspecies[i].Fit_BpC(0);  //0 BpC 1 lnBpC
        ionspecies[i].have_Fit_BpC = true;
    }
}
}
*/

//======================================================================================================
TH1F* h1 = new TH1F("h1","C_distribution",subregion_n,Cmin,Cmax); 
//!! warning!! 由于TH1F的bin是左闭右开，因此以Cmax作为右上界时，C=Cmax的那个离子的C是会超出最后一个Bin的！！
//解决方法： 发现超出时，把它归到最后一个bin-- 只有Cmax那一个离子
AxisFormat(h1,"C_distribution","C [m]","Counts");
//============================== 一些关于 gtC 的分析处理过程 ==============================
if(Draw_gtC_all_ON)
{
    Draw_gtC_all(0,0);
    Draw_gtC_all(1,6);
}
//====  20230224 ion gtC =================
for(int i=0;i<ions_n  ;i++)
{
    grerr_gtC_all->SetPoint(grerr_gtC_all_n++,ions[i].C,ions[i].gammat);//初步的γt(C)散点，带误差
}
GetGtC_Curve(grerr_gtC_all,grerr_gtC_ave, subregion_n,h1);//grerr_gtC_ave  每个 C 区间的 γt平均值和标准偏差

/// selection
for(int i=0;i<ions_n  ;i++)
{
    if(ions[i].gammat>1&&ions[i].gammat<3)//筛选条件
    {
        grerr_gtC_s->SetPoint(grerr_gtC_s_n++, ions[i].C,ions[i].gammat);//初步筛选后的 γt(C) 散点
    }
}
//==========  0929 gtC- Tz --0930 gtC for each ionspecies
for(int i=0;i<ions_n  ;i++)
{
    //if(ions[i].gammat >1.3&&ions[i].gammat <1.4)
    if(ions[i].gammat_err<=gtC_ERR_upper_bound)//γt误差足够小
    {
        for(int j=0;j<TzN  ;j++)
        {
            if( ionspecies[ions[i].Species].Tz==Tz_all_ions[j]  )//分类统计Tz
            {
                gtC_Tz[j][ii]+=ions[i].gammat;
                gtC2_Tz[j][ ii ] += ions[i].gammat*ions[i].gammat;        
                C_Division_n_Tz[j][ ii]++;
                count_Tz_n[j]++;
            }
        }
       
    }
}
//------------  0929 gtC- Tz --0930 gtC for each ionspecies



//----------------------           接下来进行构建 gtC 曲线    按照 steps 进行 -----------------------------------------------
//============================================ gtC STEP ===================================================
//+++++++++++++++++++++++ gtC step0 :original input all ions gt +++++++++++++++++++++++++++++++
if(!LOOP_ON&&SHOW_gtC_step0)//分步绘制
{
    //------------gtC step0 build:
    for(int i=0;i<ions_n;i++)
    {
        gr_gt_gterr_step0->SetPoint(gr_gt_gterr_step0_n++, ions[i].gammat, ions[i].gammat_err);
        gr_gtC_step0->SetPoint(gr_gtC_step0_n++, ions[i].C, ions[i].gammat);
        h1_gt_all_step0->Fill(ions[i].gammat);
        h1_gterr_all_step0->Fill(ions[i].gammat_err);
        h2_gtC_all_step0->Fill(ions[i].C,ions[i].gammat);
        h2_gt_gterr_step0->Fill(ions[i].gammat,ions[i].gammat_err);
        gtC_step0_n++;
    }
    //------------gtC step0 show:
    //Draw_one_TGraph(grerr_gtC_with_err , "c_grerr_gtC_with_err");
    //Draw_one_TGraph(grerr_gtC_chosen_v2 , "c_grerr_gtC_chosen_v2_before_smooth");
    
    //Draw_one_TGraph_yaxis_range(gr_gt_gterr_step0 , "c_gr_gt_gterr_step0",0,10);
    //Draw_one_TGraph_range(gr_gtC_step0, "c_gr_gtC_step0",2, -1,-1, 1.2,1.4 );
    //Draw_one_histogram(h1_gt_all_step0,"c_h1_gt_all_step0");
    //Draw_one_histogram(h1_gterr_all_step0,"c_h1_gterr_all_step0");

    gStyle->SetOptStat("ne");
    Draw_one_2Dhistogram(h2_gtC_all_step0, "c_h2_gtC_all_step0",1);
    Draw_one_2Dhistogram(h2_gt_gterr_step0, "c_h2_gt_gterr_step0");

    cout<<endl<<"+++++++ gtC step0 :original input all ions gt ++++++++  "<<endl;
    cout<<"    gtC_step0_n  = "<<gtC_step0_n<<endl;
}

//cout<<"!!!!!!!debug!!!!!! , return after gtC step0!!!"<<endl;return;//!!!!!!!!!!!!!!!!!!!!!!

//+++++++++++++++++++++++ gtC step1: gt_error filtering: filter out those ions whose gt_err are too large +++++++++++++++++++++++++++++++
if(!LOOP_ON&&SHOW_gtC_step1)
{
    //------------gtC step1 build:
    for(int i=0;i<ions_n;i++)
    {
        if(ions[i].gammat_err<gtC_ERR_upper_bound)
        {
            gr_gt_gterr_step1->SetPoint(gr_gt_gterr_step1_n++, ions[i].gammat, ions[i].gammat_err);
            gr_gtC_step1->SetPoint(gr_gtC_step1_n++, ions[i].C, ions[i].gammat);
            h1_gt_all_step1->Fill(ions[i].gammat);
            h1_gterr_all_step1->Fill(ions[i].gammat_err);
            h2_gtC_all_step1->Fill(ions[i].C,ions[i].gammat);
            h2_gt_gterr_step1->Fill(ions[i].gammat,ions[i].gammat_err);
            gtC_step1_n++;
        }
    }
    //------------gtC step1 show:
    //Draw_one_TGraph_yaxis_range(gr_gt_gterr_step1 , "c_gr_gt_gterr_step1",0,gtC_ERR_upper_bound);
    //Draw_one_TGraph_range(gr_gtC_step1, "c_gr_gtC_step1",2, -1,-1, 1.2,1.4 );
    //Draw_one_histogram(h1_gt_all_step1,"c_h1_gt_all_step1");
    //Draw_one_histogram(h1_gterr_all_step1,"c_h1_gterr_all_step1");
    Draw_one_2Dhistogram(h2_gtC_all_step1, "c_h2_gtC_all_step1",1);
    Draw_one_2Dhistogram(h2_gt_gterr_step1, "c_h2_gt_gterr_step1");

    cout<<endl<<"+++++++ gtC step1 :selection by gt error   ++++++++  "<<endl;
    cout<<"  ----------- gtC_ERR_upper_bound = "<<gtC_ERR_upper_bound<<" ----------- "<<endl;
    cout<<"    gtC_step1_n  = "<<gtC_step1_n<<endl;
    if(SHOW_gtC_step0)cout<<"    gtC_step0_n - gtC_step1_n  = "<<gtC_step0_n-gtC_step1_n<<" : "<<double(gtC_step0_n-gtC_step1_n)/gtC_step0_n*100<<" % "<<endl;
}
//cout<<"!!!!!!!debug!!!!!! , return after gtC step1!!!"<<endl;return;//!!!!!!!!!!!!!!!!!!!!!!


//+++++++++ h1_step0 C_distribution+++++++++++
TH1F* h1_STEP0 = new TH1F("h1_STEP0","h1_STEP0",subregion_n,Cmin,Cmax);
AxisFormat(h1_STEP0,"h1_STEP0_ C distribution"," C [m]", " Counts");
for(int i=0;i<ions_n  ;i++){h1_STEP0->Fill(ions[i].C);}
if(SHOW_h1_STEP0)Draw_one_histogram(h1_STEP0,"c_h1_STEP0");
 //+++++++++ h1_step1 C_distribution+++++++++++
if(Cmin<C_filter_min&&C_filter_ON)
{
    cout<<endl <<"!!!Cmin < C_filter_min change Cmin:"<<Cmin<<" to C_filter_min= "<<C_filter_min<<endl;
    Cmin = C_filter_min;
}
if(Cmax>C_filter_max&&C_filter_ON)
{
    cout<<endl <<"!!!Cmax > C_filter_max change Cmax:"<<Cmax<<" to C_filter_max= "<<C_filter_max<<endl;
    Cmax = C_filter_max;
}
TH1F* h1_STEP1 = new TH1F("h1_STEP1","h1_STEP1",subregion_n,Cmin,Cmax);

AxisFormat(h1_STEP1,"h1_STEP1_ C distribution"," C [m]", " Counts");
for(int i=0;i<ions_n  ;i++){if(ions[i].gammat_err<=gtC_ERR_upper_bound) {h1_STEP1->Fill(ions[i].C); } }

cout<<"+++h1_STEP1: gammat_err<=gtC_ERR_upper_bound , ion number = "<<h1_STEP1->GetEntries()<<endl;
if(SHOW_h1_STEP1)Draw_one_histogram(h1_STEP1,"c_h1_STEP1");



//----------- FIND effective left and right boundary of C distribution: EL, ER
// those ions with C out of this region[EL,ER] obviously have extreme C value and will not be included in the gt-C curve C region.
// the gtC curve subregioins will be the [EL,ER] dividing into subregion_n subregions.
double h1_STEP1_EL,h1_STEP1_ER;
//FIND_TH1F_effective_L_R(h1_STEP1, subregion_n,10, 3, h1_STEP1_EL,h1_STEP1_ER);

Get_h1_percent_region(h1_STEP1,h1_STEP1_EL_PERCENT,  h1_STEP1_EL,h1_STEP1_ER);

cout<<"\033[33m  effective C region:   \033[0m"<<"(h1_STEP1_EL,h1_STEP1_ER) = "<<h1_STEP1_EL<<" , "<<h1_STEP1_ER<<endl; 
//------


//=========================================================================
// use h1_STEP1_EL,h1_STEP1_ER to determine the region that will be divided into subregions 用h1_E 的分bin来作为标尺划分轨道小区间subregion
TH1F* h1_E = new TH1F("h1_E","effective C_distribution",subregion_n,h1_STEP1_EL,h1_STEP1_ER); 
AxisFormat(h1_E,"effective C_distribution","C [m]","Counts");

double subregion_Cmin = h1_STEP1_EL;
double subregion_Cmax = h1_STEP1_ER;


//====================================================================
if(Do_gtC_largeZline_ON)Do_gtC_largeZline(ions, ions_n, h1);
//______________________________________________________________________

//---------------------- 开始按轨道分小区间处理： divide ions into C subregions ---------
//===================== determine ions.C_region ========================================
// the region of subregions is determined by h1_E 在h1_E 的横轴范围内利用分bin进行小区间subregions 的划分
// 注意！！ TH1 的分bin, N 个Bin 从 1开始到N， bin 左闭右开。 最后1个Bin number N 在 new a[N] 数组是越界的！
int count_out_of_h1E_bin_L=0;int count_out_of_h1E_bin_R=0; int count_last_of_h1E_bin_R=0;
for(int i=0;i<ions_n;i++)
{
    ions[i].C_region = h1_E->FindBin(ions[i].C) - 1;//找到每个粒子其周长C对应的那个区间
    // ions[i].C_region 从 -1 ~ subregion_n 会有越界！ 以后每次用都要判断！
    if(ions[i].C_region==-1){count_out_of_h1E_bin_L++;}
    if(ions[i].C_region==subregion_n){count_out_of_h1E_bin_R++;} //
    if(ions[i].C_region==subregion_n-1){count_last_of_h1E_bin_R++;}
    ///! ions[i].C_region==subregion_n+1 会使得 [ions[i].C_region] 数组越界subregion_n
}
cout<<"----count_out_of_h1E_bin_L: "<<count_out_of_h1E_bin_L<<endl;//超h1E下界范围的粒子计数
cout<<"----count_out_of_h1E_bin_R: "<<count_out_of_h1E_bin_R<<endl;//超h1E上界范围的粒子计数
cout<<"----count_last_of_h1E_bin_R: "<<count_last_of_h1E_bin_R<<endl;//刚好在TH1 h1E 的最后一个bin粒子计数

//============== C_x initialization ================
double C_x[subregion_n]={0};   //center of each bin, the center of each C subregion
for(int i=0;i<subregion_n  ;i++)
{
    C_x[i] = h1_E->GetBinCenter(i+1);
    //cout<<"Cx[ "<<i<<" ] :"<<C_x[i]<<endl;  //debug
}
/*
//---0930 调用Fill_ion_gt()构建成员变量 每个种类自己的avegtC 数组，为 构建 avegtC 曲线做准备
// 关联 成员函数 BuildAvegtCCurve 
//2024 1113 换成外部函数实现的方式------ BuildAvegtCCurve_each()
for(int i=0;i<ions_n  ;i++)
{
    ionspecies[ions[i].Species].Fill_ion_gt(ions[i],subregion_n ,gtC_ERR_upper_bound);
}
for(int i=0;i<NSpecies  ;i++)
{
    if(ionspecies[i].N>100)
    {
        ionspecies[i].BuildAvegtCCurve(subregion_n,h1_E);   // build ionspecies[i].grerr_avegtC ionspecies[i].gr_gtC_own
        Grerr_sigma_to_gr(ionspecies[i].grerr_avegtC,ionspecies[i].gr_gtC_own_u,1);
        Grerr_sigma_to_gr(ionspecies[i].grerr_avegtC,ionspecies[i].gr_gtC_own_d,1);
        
    }
}
*/


//++++++++++++++++++++++++++ gtC step1 build ++++++++++++++++++++++
gtC_step1_n = BuildAvegtSR_step1(avegtSR_step1,avegtSR_sigma_step1,C_n_SR_step1,h1_E);
cout<<" after BuildAvegtSR_step1 gtC_step1_n= "<<gtC_step1_n<<endl;

TGraph* gr_gtC_chosen_step1 = new TGraph(subregion_n,C_x,avegtSR_step1);
TGraphErrors* grerr_gtC_chosen_step1 = new TGraphErrors(subregion_n,C_x,avegtSR_step1,0,avegtSR_sigma_step1);
AxisFormat(gr_gtC_chosen_step1,"gr_gtC_chosen_step1"," C[m] "," ave #gamma_{t} in each C interval", kAzure+1);
AxisFormat(grerr_gtC_chosen_step1,"grerr_gtC_chosen_step1"," C[m] "," ave #gamma_{t} in each C interval", kAzure+1);

if(SHOW_grerr_gtC_chosen_steps_ON && !LOOP_ON)
{
    
    //----------------------- show grerr_gtC_chosen_step1 ----------------------------
    //常驻显示 SHOW_grerr_gtC_chosen_steps_ON
    //Draw_one_TGraph_range(gr_gtC_chosen_step1, "c_gr_gtC_chosen_step1",2, -1,-1, 1.2,1.4 );

    //grerr_gtC_chosen_step1->Print();  //debug!!
    
    //Draw_one_TGraph_range(grerr_gtC_chosen_step1, "c_grerr_gtC_chosen_step1",2, -1,-1, 1.2,1.4 );
}


//+++++++++++++++++++++++ gtC step2: species selection : choose certain species(heavy ions) +++++++++++++++++++++++++++++++
//#################!!!!!! set gr_gtC_chosen selection!!!!!! ###############################
for(int i=0;i<NSpecies  ;i++)
{
    if( (ionspecies[i].A>CONDITION_gt_A) &&(ionspecies[i].MassUnknown==0)&&(ionspecies[i].Isomer_n==0) )//选择哪些核用于绘制γt(C)曲线
    {
        ionspecies[i].gtC_CHOSEN = true;
    }
    //##ARTIFICIAL
    //if(ionspecies[i].A==24&&ionspecies[i].Z==13){ionspecies[i].gtC_CHOSEN = 0;} //24Al              //特殊筛选，把某些核去掉或者加上
    //if(ionspecies[i].A==12&&ionspecies[i].Z==7){ionspecies[i].gtC_CHOSEN = 1;} //12N for SET3
    if(THIS_EXP=="2017_58Ni")
    {

    }
    if(THIS_EXP=="2021_36Ar_SET3")
    {
        if(ionspecies[i].A==14&&ionspecies[i].Z==8){ionspecies[i].gtC_CHOSEN = 0;}//14O
        if(ionspecies[i].A==10&&ionspecies[i].Z==6){ionspecies[i].gtC_CHOSEN = 0;} //10C for SET3
        if(ionspecies[i].A==20&&ionspecies[i].Z==12){ionspecies[i].gtC_CHOSEN = 0;} //20Mg for SET3
    }
    if (THIS_EXP == "2021_36Ar_SET1")
    {
        //未确认，先试试
        //if (ionspecies[i].A == 12 && ionspecies[i].Z == 7) { ionspecies[i].gtC_CHOSEN = true; }//12N可以用，再加上
        //if (ionspecies[i].A == 13 && ionspecies[i].Z == 7) { ionspecies[i].gtC_CHOSEN = false; }//13N gammat(C)曲线过低 
        //if (ionspecies[i].A == 11 && ionspecies[i].Z == 6) { ionspecies[i].gtC_CHOSEN = false; }//11C gammat(C)曲线过低
    }

}

//=================== gtC_STEP-2 species selection for gtC-curve ==========================
cout<<endl<<"=== \033[33m gtC_STEP-2 species selection for gtC-curve: \033[0m "<<endl;
outfile_logs<<endl<<"===  gtC_STEP-2 species selection for gtC-curve: "<<endl;
for(int i=0;i<NSpecies;i++)
{
    if(ionspecies[i].gtC_CHOSEN)//依次输出检查
    {
        cout<<"   "<<ionspecies[i].Aname;
        outfile_logs<<" | "<<ionspecies[i].Aname<<" "<<ionspecies[i].A<<" "<<ionspecies[i].Z;
    }
}
outfile_logs<<endl<<" C_DIVISION_CHOSEN_MIN =  "<<C_DIVISION_CHOSEN_MIN<<endl;
cout<<endl<<"----------"<<endl;
//------------------------- start building gtC_chosen -------------------------
// reset gtC_chosen[]
for(int ii=0;ii<subregion_n  ;ii++)
{
    gtC_chosen[ii]=0;gtC2_chosen[ii]=0;C_Division_n_chosen[ii]=0;
}
ii=0;
gtC_step2_n=0;
for(int i=0;i<ions_n;i++)
{
    ii=ions[i].C_region;//找到每个粒子其周长C对应的那个区间 
    if(ii<=0||ii>=subregion_n) continue;
    // 只使用 ii =1 ~ subregion_n -1 
    //#################!!!!!! use gr_gtC_chosen selection!!!!!! ###############################
    if( ionspecies[ions[i].Species].gtC_CHOSEN == true )
    {   
        if(ions[i].C>=(C_filter_min+0.0)&&ions[i].C<=C_filter_max-0.0 && ions[i].gammat_err<=gtC_ERR_upper_bound )//对用来绘制 γt(C) 核进行限制
        {
        //###########################   step2 selection here   ##########################################
            gtC_chosen[ii]+=ions[i].gammat;                                                    //γt平均值（等权重）
            gtC2_chosen[ii] += ions[i].gammat*ions[i].gammat;                                //γt平方（用作算方差）
            gtC_chosen_v2[ii]+=ions[i].gammat*(1.0/ions[i].gammat_err/ions[i].gammat_err);     //γt平均值（加权平均）
            gtC_chosen_v2_err[ii] +=(1.0/ions[i].gammat_err/ions[i].gammat_err);               //γt权重因子和（用来算方差）
            C_Division_n_chosen[ii] ++;                                                       //此周长区间的参考核数

            //------------gtC step2 build:
            if(!LOOP_ON&&SHOW_gtC_step2)   
            {
            gr_gt_gterr_step2->SetPoint(gr_gt_gterr_step2_n++, ions[i].gammat, ions[i].gammat_err);
            gr_gtC_step2->SetPoint(gr_gtC_step2_n++, ions[i].C, ions[i].gammat);
            h1_gt_all_step2->Fill(ions[i].gammat);
            h1_gterr_all_step2->Fill(ions[i].gammat_err);
            h2_gtC_all_step2->Fill(ions[i].C,ions[i].gammat);
            h2_gt_gterr_step2->Fill(ions[i].gammat,ions[i].gammat_err);
            gtC_step2_n++;
            }
            
            // step2: C_division_subregion counts in h1_E
            h1_E->Fill(ions[i].C);
        //____________________________#   step2 selection here   #_______________________________________________
        }
    }
    //gtC_chosen_v2[ii]+=ions[i].gammat*(1.0/ions[i].gammat_err/ions[i].gammat_err);
    //gtC_chosen_v2_err[ii] +=(1.0/ions[i].gammat_err/ions[i].gammat_err);
    
    gr_gtC_tmp->SetPoint(gr_gtC_tmp_n++,ions[i].C, ions[i].gammat );      // 可变内容
    // 20230707 add gt err
    grerr_gtC_with_err->SetPoint(grerr_gtC_with_err_n ,ions[i].C, ions[i].gammat  ); 
    grerr_gtC_with_err->SetPointError(grerr_gtC_with_err_n ,0, ions[i].gammat_err  );
    grerr_gtC_with_err_n++;
    
}

//========================   complete chosen gtC   =====================================

for(int i=0;i<subregion_n  ;i++)//初步计算γt（C）数组
{
    double C_center_i=h1_E->GetBinCenter(i+1);
    if(C_Division_n_chosen[i]>3)
    {
        gtC_chosen[i]/=C_Division_n_chosen[i];                                //计算中心值和标准差
        gtC2_chosen[i]/=C_Division_n_chosen[i];
        gtC_chosen_err[i]=sqrt(gtC2_chosen[i]-gtC_chosen[i]*gtC_chosen[i]);
        gtC_chosen_err[i] /= sqrt(C_Division_n_chosen[i] -1);
        //---------- gt err v2 20230707
        gtC_chosen_v2[i] /= gtC_chosen_v2_err[i];
        gtC_chosen_v2_err[i] = sqrt(1.0/gtC_chosen_v2_err[i] ) ;
        //cout<<"debug: "<<gtC_chosen_v2[i]<<" "<<gtC_chosen_v2_err[i]<<endl;
    }   
}

//cout<<" debug here 6237 "<<endl;

TGraph* gr_gtC_chosen_step2 = new TGraph();
TGraphErrors*  grerr_gtC_chosen_step2 = new TGraphErrors();
AxisFormat(gr_gtC_chosen_step2,"gr_gtC_chosen_step2"," C[m] "," ave #gamma_{t} in each subregion", kAzure+2);
AxisFormat(grerr_gtC_chosen_step2,"grerr_gtC_chosen_step2"," C[m] "," ave #gamma_{t} in each subregion", kAzure+2);
//----------------------- show grerr_gtC_chosen_step2 ----------------------------
//常驻显示 SHOW_grerr_gtC_chosen_steps_ON
if(SHOW_grerr_gtC_chosen_steps_ON && !LOOP_ON)
{
    for(int i=0;i<subregion_n  ;i++)
    {
        double C_center_i=h1_E->GetBinCenter(i+1);
        if(C_Division_n_chosen[i]>3)
        {
            //++++++++++++++++ build curve step2 +++++++++++++++++++
            gr_gtC_chosen_step2->SetPoint(gr_gtC_chosen_step2->GetN(),C_center_i, gtC_chosen[i] );
            grerr_gtC_chosen_step2->SetPoint(grerr_gtC_chosen_step2->GetN(),C_center_i, gtC_chosen[i] );
            grerr_gtC_chosen_step2->SetPointError(grerr_gtC_chosen_step2->GetN()-1,0, gtC_chosen_err[i] );
        }   
    }
    //Draw_one_TGraph_range(gr_gtC_chosen_step2, "c_gr_gtC_chosen_step2",2, -1,-1, 1.2,1.4 );
    //grerr_gtC_chosen_step2->Print();
    
    //Draw_one_TGraph_range(grerr_gtC_chosen_step2, "c_grerr_gtC_chosen_step2",2, -1,-1, 1.2,1.4 );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas* c_gtC_STEP1_and_2 = new TCanvas("c_gtC_STEP1_and_2","c_gtC_STEP1_and_2",1000,500);
    TMultiGraph* mg_gtC_chosen_step1and2 = new TMultiGraph();
    mg_gtC_chosen_step1and2->Add(grerr_gtC_chosen_step1);
    mg_gtC_chosen_step1and2->Add(grerr_gtC_chosen_step2);
    mg_gtC_chosen_step1and2->Draw("ap");
    c_gtC_STEP1_and_2->BuildLegend();
    //c_gtC_STEP1_and_2->Update();
    AxisFormat(mg_gtC_chosen_step1and2,"grerr_gtC_chosen step1 and step2"," C[m] ", " ave #gamma_{t} in each C interval ");
    mg_gtC_chosen_step1and2->GetYaxis()->SetRangeUser(1.2,1.4);


    //----------------------- gr_subregion_gt_err_step2 --------------------------
    // to see gammat error in each subregion:
    int gr_subregion_gt_err_step2_n=0;
    TGraph* gr_subregion_gt_err_step2 = new TGraph();
    AxisFormat(gr_subregion_gt_err_step2,"gr_subregion_gt_err_step2"," C[m] ", "subregion:  #sigma (#gamma_{t}) ");
    for(int i=0;i<subregion_n  ;i++)
    {
        if(C_Division_n_chosen[i]>C_DIVISION_CHOSEN_MIN)
        {
            gr_subregion_gt_err_step2->SetPoint(gr_subregion_gt_err_step2_n++,h1_E->GetBinCenter(i+1),gtC_chosen_err[i]) ;
        }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    TCanvas* c_gr_subregion_gt_err_step2 = new TCanvas("c_gr_subregion_gt_err_step2","c_gr_subregion_gt_err_step2",1000,500);
    gr_subregion_gt_err_step2->Draw("apl");

}

//----------------------- gr_subregion_gt_err_step2 ---------------------------

//+++++++++++++ ------------gtC step2 show:  ++++++++++++++++++++++++++++++++++++++++++++
if(!LOOP_ON&&SHOW_gtC_step2)
{
    Draw_one_TGraph_yaxis_range(gr_gt_gterr_step2 , "c_gr_gt_gterr_step2",0,gtC_ERR_upper_bound);
    Draw_one_TGraph_range(gr_gtC_step2, "c_gr_gtC_step2",2, -1,-1, 1.2,1.4 );

    double h1_gt_all_step2_90L,h1_gt_all_step2_90R;
    Get_h1_percent_region(h1_gt_all_step2,90, h1_gt_all_step2_90L,h1_gt_all_step2_90R);
    //Draw_one_histogram_range(h1_gt_all_step2,"c_h1_gt_all_step2",1,h1_gt_all_step2_90L,h1_gt_all_step2_90R,-1,-1);

    Draw_one_histogram(h1_gterr_all_step2,"c_h1_gterr_all_step2");
    Draw_one_2Dhistogram(h2_gtC_all_step2, "c_h2_gtC_all_step2",1);
    Draw_one_2Dhistogram(h2_gt_gterr_step2, "c_h2_gt_gterr_step2");
    cout<<endl<<"+++++++ gtC step2 :original input all ions gt ++++++++  "<<endl;
    cout<<"    gtC_step2_n  = "<<gtC_step2_n<<endl;
    if(SHOW_gtC_step0)cout<<"    gtC_step0_n - gtC_step2_n  = "<<gtC_step0_n-gtC_step2_n<<" : "<<double(gtC_step0_n-gtC_step2_n)/gtC_step0_n*100<<" % "<<endl;
    
}
if(SHOW_h1_E_step2)Draw_one_histogram(h1_E,"c_h1_E");

//cout<<" return at 6714 "<<endl;return;

//+++++++++++++++++++++++ gtC step3: n_sigma squeeze: filter out ions that is n_sigma away from central value of this subregion +++++++++++++++++++++++++++++++

//======================= squeeze_k ===========================
// Z-score 把偏离太大的点筛掉。重新计算γt(C)数组
// refresh gtC_chosen 刷新 gtC 数组
TGraphErrors* grerr_gtC_squeeze[squeeze_k]; //记录每一次刷新后的γt(C)曲线数组，刷新squeeze_k次
//--int grerr_gtC_squeeze_n[squeeze_k];
TGraph* gr_gtC_squeeze_zone[squeeze_k];  // 显示当前筛选的带状区域内留下的散点

if(squeeze_k>=1)
{
    for(int k=1;k<=squeeze_k;k++)
    {   
        gr_gtC_squeeze_zone[k-1] = new TGraph();
        AxisFormat(gr_gtC_squeeze_zone[k-1],strtmp.Format(" points zone :_%d",k-1),"C[m]","#gamma_{t}",kAzure-7);
        gr_gtC_squeeze_zone[k-1]->SetMarkerSize(1.0);
        gr_gtC_squeeze_zone[k-1]->SetMarkerStyle(24);  //空心圈
    
        //--grerr_gtC_squeeze_n[k-1]=0;
        grerr_gtC_squeeze[k-1] = new TGraphErrors();
        AxisFormat(grerr_gtC_squeeze[k-1],strtmp.Format("grerr_gtC_squeeze[%d]",k-1),"C[m]","#gamma_{t}",kGreen-3);
        cout<<" squeeze No. "<<k<<endl;
        squeeze_gtC(gtC_chosen,gtC_chosen_err , C_Division_n_chosen, n_sigma, h1_E,gr_gtC_squeeze_zone[k-1]);
        //调用squeeze_gtC 刷新 gtC_chosen[] ,gtC_chosen_err[]
        for(int i=0;i<subregion_n ;i++)
        {
            if(C_Division_n_chosen[i]>3)
            {
                //gr_gtC_chosen->SetPoint(grerr_n,h1_E->GetBinCenter(i+1),gtC_chosen[i]);
                grerr_gtC_squeeze[k-1]->SetPoint(grerr_gtC_squeeze[k-1]->GetN(),h1_E->GetBinCenter(i+1),gtC_chosen[i]);
                grerr_gtC_squeeze[k-1]->SetPointError(grerr_gtC_squeeze[k-1]->GetN()-1,0,gtC_chosen_err[i]);
                //--grerr_gtC_squeeze_n[k-1]++;
            }    
        }
        //grerr_gtC_squeeze[k-1]->Print();
        TGraphErrors_to_outfile(FILEPATH_chosen+strtmp.Format("grerr_gtC_squeeze[%d].txt",k-1),grerr_gtC_squeeze[k-1]);
        /*
        outfile_grerr_txt.open(FILEPATH_chosen+strtmp.Format("grerr_gtC_squeeze[%d].txt",k-1) );
        for(int i=0;i<grerr_gtC_squeeze[k-1]->GetN();i++)
        {
            grerr_gtC_squeeze[k-1]->GetPoint(i,xtmp,ytmp);
            yerrtmp=grerr_gtC_squeeze[k-1]->GetErrorY(i);
            outfile_grerr_txt<<xtmp<<" "<<ytmp<<" "<<yerrtmp<<endl;
        }
        cout<<"-------TGraphErrors_to_outfile : save "<<grerr_gtC_squeeze[k-1]->GetTitle()<<" to "<<FILEPATH_chosen+strtmp.Format("grerr_gtC_squeeze[%d].txt",k-1)<<endl;
        outfile_grerr_txt.close();
         */
        cout<<"||ions left:"<<gr_gtC_squeeze_zone[k-1]->GetN()<<endl;
        
    }
    if(!LOOP_ON&& Show_c_squeeze_tmp_ON)
    {
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TCanvas* c_squeeze_tmp = new TCanvas("c_squeeze_tmp","c_squeeze_tmp",1400,700);//绘制变化过程
        c_squeeze_tmp->Divide(2,2);
        
        c_squeeze_tmp->cd(1);
        gr_gtC_squeeze_zone[0]->Draw("ap");
        grerr_gtC_squeeze[0]->Draw("plsame");
        gr_gtC_squeeze_zone[0]->GetYaxis()->SetRangeUser(1.22,1.38);
        //c_squeeze_tmp->Print(FILEPATH_chosen+strtmp.Format("squeeze_loop%d_k_0.png",scan_loop_i ));
    
        if(squeeze_k>=2)
        {  
            c_squeeze_tmp->cd(2);
            gr_gtC_squeeze_zone[squeeze_k/2]->Draw("ap");
            grerr_gtC_squeeze[squeeze_k/2]->Draw("plsame");
            gr_gtC_squeeze_zone[squeeze_k/2]->GetYaxis()->SetRangeUser(1.22,1.38);
            
            //c_squeeze_tmp->Print(FILEPATH_chosen+strtmp.Format("squeeze_loop%d_k_%d.png",scan_loop_i,squeeze_k ));
                  
            c_squeeze_tmp->cd(3);
            gr_gtC_squeeze_zone[squeeze_k-1-1]->Draw("ap");
            grerr_gtC_squeeze[squeeze_k-1-1]->Draw("plsame");
            gr_gtC_squeeze_zone[squeeze_k-1-1]->GetYaxis()->SetRangeUser(1.22,1.38);
            //c_squeeze_tmp->Print(FILEPATH_chosen+strtmp.Format("squeeze_loop%d_k_%d.png",scan_loop_i,squeeze_k ));
               
            c_squeeze_tmp->cd(4);
            gr_gtC_squeeze_zone[squeeze_k-1]->Draw("ap");
            grerr_gtC_squeeze[squeeze_k-1]->Draw("plsame");
            gr_gtC_squeeze_zone[squeeze_k-1]->GetYaxis()->SetRangeUser(1.22,1.38);
            c_squeeze_tmp->Print(FILEPATH_chosen+strtmp.Format("squeeze_loop%d_k_%d.png",scan_loop_i,squeeze_k ));
        }   
    }
    if(LOOP_ON)
    {
        for(int k=1;k<=squeeze_k;k++)
        {
            delete gr_gtC_squeeze_zone[k-1] ;
            delete grerr_gtC_squeeze[k-1];
        }
    }
    

} //if (squeeze_k>=1)
//__________________________ squeeze _____________________________________ 

cout<<" termination debug after squeeze"<<endl;

//###################### set gr_gtC_chosen (GT-SHIFT) #############################
grerr_n=0;
if(Show_C_Division_n_chosen_ON)cout<<" now Set gr/grerr gtC chosen -----------"<<endl;
for(int i=0;i<subregion_n ;i++)
{   
    double Ci=h1_E->GetBinCenter(i+1); // 注意 第一个区间 gtC[0] 对应于 Bin(1)
    if(C_Division_n_chosen[i]>C_DIVISION_CHOSEN_MIN)//这个区间的可用于绘制γt(C)粒子数足够多
    {
        if(GT_SHIFT_Control_ON)
        {
            if(Ci>=GT_SHIFT_Left&&Ci<=GT_SHIFT_Right)//GT_SHIFT一个微小的偏移量，用于手动微调γt(C)曲线，一般为0
            {
                gr_gtC_chosen->SetPoint(grerr_n,Ci,gtC_chosen[i] + GT_SHIFT);
                grerr_gtC_chosen->SetPoint(grerr_n,Ci,gtC_chosen[i]+GT_SHIFT);
                grerr_gtC_chosen->SetPointError(grerr_n,0,gtC_chosen_err[i]);
            }
            else
            {
                gr_gtC_chosen->SetPoint(grerr_n,Ci,gtC_chosen[i] + 0);  //用TGraph保存筛选好的γt(C)曲线，实际计算使用的是gr_gtC_chosen
                grerr_gtC_chosen->SetPoint(grerr_n,Ci,gtC_chosen[i]+0);
                grerr_gtC_chosen->SetPointError(grerr_n,0,gtC_chosen_err[i]);
            }
        }
        else
        {
            gr_gtC_chosen->SetPoint(grerr_n,Ci,gtC_chosen[i] );          //用TGraph保存筛选好的γt(C)曲线，实际计算使用的是gr_gtC_chosen
            grerr_gtC_chosen->SetPoint(grerr_n,Ci,gtC_chosen[i]);
            grerr_gtC_chosen->SetPointError(grerr_n,0,gtC_chosen_err[i]);
        }

        //------ 20230707 v2
        grerr_gtC_chosen_v2->SetPoint(grerr_n,Ci,gtC_chosen_v2[i]+GT_SHIFT); //加权γt(C)曲线
        grerr_gtC_chosen_v2->SetPointError(grerr_n,0,gtC_chosen_v2_err[i]);
        grerr_n++;
        if(Show_C_Division_n_chosen_ON)cout<<" bin: "<<i+1<<" | C subregion center :"<<Ci<<" points:  "<<C_Division_n_chosen[i]<<endl;
    }    
    outfile_C_Division_n<<" bin: "<<i+1<<" C_subregion_center: "<<Ci<<" points: "<<C_Division_n_chosen[i]<<endl;
}
if(GT_SHIFT_Control_ON)
{
    cout<<" ---------- GT_SHIFT_Control_ON : GT_SHIFT = "<<GT_SHIFT<<" region: "<<GT_SHIFT_Left<<" ~~ "<<GT_SHIFT_Right<<endl;
    outfile_logs<<" ---------- GT_SHIFT_Control_ON : GT_SHIFT = "<<GT_SHIFT<<" region: "<<GT_SHIFT_Left<<" ~~ "<<GT_SHIFT_Right<<endl;
}
//======================= read in gtC =====================
cout<<" termination debug after readin gtC"<<endl;

TGraph* gr_gtC_readin = new TGraph();
AxisFormat(gr_gtC_readin,"gr_gtC_readin","C [m]","#gamma_{t}",3);
if(READIN_gtC_chosen_ON)
{
    cout<<" \033[33m  !!!!!! ---- READIN_gtC_chosen_ON ------ !!!!  \033[0m "<<endl;
    outfile_logs<<"          !!!!!! ---- READIN_gtC_chosen_ON ------ !!!!         "<<endl;
    TGraph_from_infile("INPUT//gtC_xym_20240625.txt",gr_gtC_readin);
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_SAVE_chosen_gtC = new TCanvas("c_SAVE_chosen_gtC","c_SAVE_chosen_gtC",1000,500);
if(READIN_gtC_chosen_ON)
{
    
    cout<<" \033[33m  !!!!!! ---- READIN_gtC_chosen_ON ------ !!!!  \033[0m "<<endl;
    outfile_logs<<"          !!!!!! ---- READIN_gtC_chosen_ON ------ !!!!         "<<endl;
    TGraph_from_infile("INPUT//gtC_xym_20240625.txt",gr_gtC_readin);
    gr_gtC_readin->DrawClone("Ap");
    c_SAVE_chosen_gtC->Print(FILEPATH_chosen+strtmp.Format("gtC_READIN_%d.png",scan_loop_i));
}

else 
{
    //这里的 gtC 还没有 smooth
    //gr_gtC_chosen->DrawClone("ap");
    //c_SAVE_chosen_gtC->Print(FILEPATH_chosen+strtmp.Format("gtC_%d.png",scan_loop_i));
}
delete c_SAVE_chosen_gtC;
*/


//20240709
//TGraphErrors_to_outfile(FILEPATH_chosen+strtmp.Format("gtC_data_%d.txt",scan_loop_i) , grerr_gtC_chosen);
/*
outfile_grerr_txt.open(FILEPATH_chosen+strtmp.Format("gtC_data_%d.txt",scan_loop_i));
for(int i=0;i<grerr_gtC_chosen->GetN();i++)
{
    grerr_gtC_chosen->GetPoint(i,xtmp,ytmp);
    yerrtmp=grerr_gtC_chosen->GetErrorY(i);
    outfile_grerr_txt<<xtmp<<" "<<ytmp<<" "<<yerrtmp<<endl;
}
cout<<"-------TGraphErrors_to_outfile : save "<<grerr_gtC_chosen->GetTitle()<<" to "<<FILEPATH_chosen+strtmp.Format("gtC_data_%d.txt",scan_loop_i)<<endl;
outfile_grerr_txt.close();
*/
//////////////////////////// _________ chosen gtC _______________________________


// ====================== gtC_ with time =======================
// 按照注入分别生成γt(C),在函数里面看一下γt(C)是否随时间变动， 内部画图
// 得到的 各个 grerr_avegtC_inj 在后面计算质量时 use_gtC_type ==3 的时候使用 
TGraphErrors* grerr_avegtC_inj_1 = new TGraphErrors();
TGraphErrors* grerr_avegtC_inj_2 = new TGraphErrors();
TGraphErrors* grerr_avegtC_inj_3 = new TGraphErrors();
TGraphErrors* grerr_avegtC_inj_4 = new TGraphErrors();

if(Do_gtC_with_time_ON)
{    
    if(THIS_EXP=="2017_58Ni")
    {
        Do_gtC_with_time(h1_E,grerr_avegtC_inj_1,grerr_avegtC_inj_2,grerr_avegtC_inj_3,grerr_avegtC_inj_4);
        cout<<" inj1: "<<gtC_inj_min_1<<" ~ "<<gtC_inj_max_1<<" inj2: "<<gtC_inj_min_2<<" ~ "<<gtC_inj_max_2
        <<" inj3: "<<gtC_inj_min_3<<" ~ "<<gtC_inj_max_3<<" inj4: "<<gtC_inj_min_4<<" ~ "<<gtC_inj_max_4<<endl;
    TGraphErrors_to_outfile(FILEPATH_chosen+"gtC_data_inj_1.txt" , grerr_avegtC_inj_1);
    TGraphErrors_to_outfile(FILEPATH_chosen+"gtC_data_inj_2.txt" , grerr_avegtC_inj_2);
    TGraphErrors_to_outfile(FILEPATH_chosen+"gtC_data_inj_3.txt" , grerr_avegtC_inj_3);
    TGraphErrors_to_outfile(FILEPATH_chosen+"gtC_data_inj_4.txt" , grerr_avegtC_inj_4);
    }
    
}


//if(Do_gtC_with_time_ON){cout<<"return after Do_gtC_with_time" <<endl;return;}
//_______________________ gtC with time ___________________________
cout<<" debug !! after gtC with time"<<endl;
//cout<<"return "<<endl;return;


/*
//=========== cubic spline interpolation
TSpline3* sp3 = new TSpline3( "sp3", gr_gtC_chosen);
TSpline3* sp3_smooth = new TSpline3( "sp3_smooth", gr_gtC_chosen_smooth);
*/
//+++++++++++++++++++++++++++++++++ gtC step4: smooth ++++++++++++++++++++++++++++++++++++++++++++++++
if(Do_gtC_smooth_ON)
{
cout<<" =============== Do_gtC_smooth_ON ==================================="<<endl;
    TGraph* gr_gtC_chosen_smooth = new TGraph();
    AxisFormat(gr_gtC_chosen_smooth,"gtC_smooth","C [m]","#gamma_{t}",kOrange);
    TGraphErrors* grerr_gtC_chosen_smooth = new TGraphErrors();
    AxisFormat(grerr_gtC_chosen_smooth,"gtC_smooth_err","C [m]","#gamma_{t}",kOrange);
    Smooth_DW(gr_gtC_chosen,gr_gtC_chosen_smooth,set_smooth_k, smooth_opt);              //平滑γt(C)曲线
    Smooth_DW_err(grerr_gtC_chosen,grerr_gtC_chosen_smooth,set_smooth_k, smooth_opt);    //平滑γt(C)曲线误差
    cout<<" termination debug ! after Smooth_DW "<<endl;
    if(Show_c_smooth_ON&&!LOOP_ON)
    {
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TCanvas* c_smooth = new TCanvas("c_smooth","c_smooth",1000,500);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        grerr_gtC_chosen->DrawClone("apl");
        grerr_gtC_chosen_smooth->DrawClone("plsame");
        c_smooth->Print(FILEPATH_chosen+strtmp.Format("smooth_%d.png",scan_loop_i));
        
    }
    

    //gr_gtC_chosen_smooth->Print();
    //grerr_gtC_chosen_smooth->Print();
    //gr_gtC_chosen_smooth->SetPoint(4,128.586,1.31);
    //grerr_gtC_chosen_smooth->SetPoint(4,128.586,1.31);
    /////##########          
    delete gr_gtC_chosen;//γt(C)曲线，计算用
    delete grerr_gtC_chosen;
    gr_gtC_chosen = gr_gtC_chosen_smooth;//直接替换为平滑后的γt(C)曲线
    grerr_gtC_chosen = grerr_gtC_chosen_smooth;
    AxisFormat(grerr_gtC_chosen,"#gamma_{t}-C curve with error bars","C [m]","#gamma_{t}",kRed);

}

TGraph* gr_gtC_chosen_u = new TGraph();
TGraph* gr_gtC_chosen_d = new TGraph();
AxisFormat(gr_gtC_chosen_u,"gr_gtC_chosen+#sigma","C [m]","#gamma_{t}",kAzure+1);
AxisFormat(gr_gtC_chosen_d,"gr_gtC_chosen-#sigma","C [m]","#gamma_{t}",kAzure+3);
Grerr_sigma_to_gr(grerr_gtC_chosen,gr_gtC_chosen_u,1);
Grerr_sigma_to_gr(grerr_gtC_chosen,gr_gtC_chosen_d,2);


//========================================= grerr_gtC_chosen  c u d=========================================================
//显示 最终的 gtC 以及上下界 保存到文件
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas* c_gr_gtC_chosen_ud = new TCanvas("c_gr_gtC_chosen_ud","c_gr_gtC_chosen_ud",1400,700);
grerr_gtC_chosen->Draw("apl");
gr_gtC_chosen_u->Draw("plsame");
gr_gtC_chosen_d->Draw("plsame");

if(THIS_EXP=="2021_36Ar_SET3")grerr_gtC_chosen->GetYaxis()->SetRangeUser(1.26,1.36);  //
if(THIS_EXP=="2021_36Ar_SET1")grerr_gtC_chosen->GetYaxis()->SetRangeUser(1.26,1.38);  //
if(THIS_EXP=="2017_58Ni")     grerr_gtC_chosen->GetYaxis()->SetRangeUser(1.33,1.38);
c_gr_gtC_chosen_ud->BuildLegend();
c_gr_gtC_chosen_ud->Print(FILEPATH_chosen+strtmp.Format("gr_gtC_chosen_ud_%d.png",scan_loop_i));
//____________________________________________________________________________________________

//cout<<" debug!! grerr_gtC_chosen->Print()"<<endl;
//grerr_gtC_chosen->Print();

//  有效 gtC 的轨道区间
gtC_chosen_CMIN = grerr_gtC_chosen->GetX()[0];
gtC_chosen_CMAX = grerr_gtC_chosen->GetX()[grerr_gtC_chosen->GetN()-1];

//  subregion_Cmax subregion_Cmin 是划分 subregion的左右端点

cout<<" ======== grerr_gtC_chosen first point: gtC_chosen_CMIN ==== "<<gtC_chosen_CMIN
    <<" ======== grerr_gtC_chosen last point : gtC_chosen_CMAX ==== "<<gtC_chosen_CMAX<<endl;

TGraph_to_outfile(FILEPATH_chosen+strtmp.Format("grr_gtC_u_data_%d.txt",scan_loop_i) , gr_gtC_chosen_u);
TGraph_to_outfile(FILEPATH_chosen+strtmp.Format("grr_gtC_d_data_%d.txt",scan_loop_i) , gr_gtC_chosen_d);
TGraphErrors_to_outfile(FILEPATH_chosen+strtmp.Format("grr_gtC_c_data_%d.txt",scan_loop_i) , grerr_gtC_chosen);

cout<<"\033[33m =========== gtC construction completed ==============  \033[0m"<<endl;


///===== own gtC
BuildAvegtCCurve_each(h1_E, gtC_ERR_upper_bound );

cout<<" ======== own gtC for each ionspecies is built ======"<<endl;


////新增部分 J.H.Lv
//========================================== from grr_gtC to ISF 从得到的gtC 曲线 计算 ISF =======================================
//ISF(Integral Step Factor) 积分步长因子  即 (1 + γt * γt * ΔC/ C)数组，在使用gtC时对所有离子是公共通用的，提前算好，加快运算速度

//double n_gtC_chosen_C_1 = (gtC_chosen_CMAX - gtC_chosen_CMIN) / ((subregion_Cmax - subregion_Cmin) / subregion_n) + 1 + 0.5; //有效γt(C)区间个数 +0.5四舍五入
//n_gtC_chosen_C = (int)n_gtC_chosen_C_1;//有效γt(C)区间个数

n_gtC_chosen_C = grerr_gtC_chosen->GetN();
double cal_dC = (gtC_chosen_CMAX - gtC_chosen_CMIN) / (n_gtC_chosen_C - 1.0) / dsubregion_n;

Calculate_GtC_ISF(gr_gtC_chosen, ISF);//得到(1 + γt * γt * ΔC/ C)数组，即 ISF 步长积分因子，加快运算速度
Calculate_GtC_ISF(gr_gtC_chosen_u, ISF_u);//上限
Calculate_GtC_ISF(gr_gtC_chosen_d, ISF_d);//下限
TH1D* h_locate = new TH1D("h_locate", "h_locate", n_gtC_chosen_C * dsubregion_n, 
                            gtC_chosen_CMIN - 0.5* ((subregion_Cmax - subregion_Cmin) / subregion_n), 
                            gtC_chosen_CMAX + 0.5* ((subregion_Cmax - subregion_Cmin) / subregion_n));
cout<<" h_locate: subregion_Cmin="<<subregion_Cmin<<" subregion_Cmax= "<<subregion_Cmax<<endl;

//这个直方图用于确定  各核周长所在区间
//这个直方图奇怪的范围是为了和ISF区间对齐,使得
//h_locate->GetBinCenter(h_locate->FindBin(C))       =        ISF[h_locate->FindBin(C)] 
//新增部分结束

//------------- 尝试 对于 特点种类离子使用单独的 ISF
int tmp_10C=-1;
int tmp_11C=-1;
if(THIS_EXP=="2021_36Ar_SET3")
{   
    
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ionspecies[i].Aname=="10C"){tmp_10C = ionspecies[i].Species;} 
        if(ionspecies[i].Aname=="11C"){tmp_11C = ionspecies[i].Species;} 
    }
    if(tmp_10C<0||tmp_11C<0){cout<<"error!! tmp_10C<0||tmp_11C<0"<<endl;}
    /*
    Calculate_GtC_ISF(ionspecies[tmp_10C].gr_gtC_own, ISF_1);
    Calculate_GtC_ISF(ionspecies[tmp_10C].gr_gtC_own_u, ISF_u_1);//上限
    Calculate_GtC_ISF(ionspecies[tmp_10C].gr_gtC_own_d, ISF_d_1);//下限
    Calculate_GtC_ISF(ionspecies[tmp_11C].gr_gtC_own, ISF_2);
    Calculate_GtC_ISF(ionspecies[tmp_11C].gr_gtC_own_u, ISF_u_2);//上限
    Calculate_GtC_ISF(ionspecies[tmp_11C].gr_gtC_own_d, ISF_d_2);//下限
    */
}
//----------------------------------------------------------------
/* shift ISF
TGraph* gr_gtC_chosen_shift_1 = new TGraph();
TGraph* gr_gtC_chosen_shift_1_u = new TGraph();
TGraph* gr_gtC_chosen_shift_1_d = new TGraph();
TGraph* gr_gtC_chosen_shift_2 = new TGraph();
TGraph* gr_gtC_chosen_shift_2_u = new TGraph();
TGraph* gr_gtC_chosen_shift_2_d = new TGraph();
TGraph_shift(gr_gtC_chosen, gr_gtC_chosen_shift_1,0.002);
TGraph_shift(gr_gtC_chosen_u, gr_gtC_chosen_shift_1_u,0.002);
TGraph_shift(gr_gtC_chosen_d, gr_gtC_chosen_shift_1_d,0.002);
TGraph_shift(gr_gtC_chosen, gr_gtC_chosen_shift_2,-0.001);
TGraph_shift(gr_gtC_chosen_u, gr_gtC_chosen_shift_2_u,-0.001);
TGraph_shift(gr_gtC_chosen_d, gr_gtC_chosen_shift_2_d,-0.001);

Calculate_GtC_ISF(gr_gtC_chosen_shift_1, ISF_1);
Calculate_GtC_ISF(gr_gtC_chosen_shift_1_u, ISF_u_1);//上限
Calculate_GtC_ISF(gr_gtC_chosen_shift_1_d, ISF_d_1);//下限
Calculate_GtC_ISF(gr_gtC_chosen_shift_2, ISF_2);
Calculate_GtC_ISF(gr_gtC_chosen_shift_2_u, ISF_u_2);//上限
Calculate_GtC_ISF(gr_gtC_chosen_shift_2_d, ISF_d_2);//下限
*/
//_____________________________________


for(int i=0;i<NSpecies  ;i++)
{
    ionspecies[i].Create_gr_gtC_shifted_own(gr_gtC_chosen,gr_gtC_chosen_u,gr_gtC_chosen_d);
}

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas* c_gtC_shifted_own = new TCanvas("c_gtC_shifted_own","c_gtC_shifted_own",1000,500);
TMultiGraph* mg_gtC_shifted_own = new TMultiGraph();
for(int i=0;i<NSpecies  ;i++)
{
    mg_gtC_shifted_own->Add(ionspecies[i].gr_gtC_shifted_own,"pl");
}
mg_gtC_shifted_own->Draw("apl");
*/

//----------- gtC_inj- ISF
if(THIS_EXP=="2017_58Ni"&&Do_gtC_with_time_ON)
{
    Calculate_GtC_ISF(grerr_avegtC_inj_1, ISF_inj1);//
    TH1D* h_locate_inj1; 
    From_grr_gtC_to_h_locate(grerr_avegtC_inj_1, h_locate_inj1,"h_locate_inj1",subregion_Cmin,subregion_Cmax);
    Calculate_GtC_ISF(grerr_avegtC_inj_2, ISF_inj2);//
    TH1D* h_locate_inj2; 
    From_grr_gtC_to_h_locate(grerr_avegtC_inj_2, h_locate_inj2,"h_locate_inj2",subregion_Cmin,subregion_Cmax);
    Calculate_GtC_ISF(grerr_avegtC_inj_3, ISF_inj3);//
    TH1D* h_locate_inj3; 
    From_grr_gtC_to_h_locate(grerr_avegtC_inj_3, h_locate_inj3,"h_locate_inj3",subregion_Cmin,subregion_Cmax);
    Calculate_GtC_ISF(grerr_avegtC_inj_4, ISF_inj4);//
    TH1D* h_locate_inj4; 
    From_grr_gtC_to_h_locate(grerr_avegtC_inj_4, h_locate_inj4,"h_locate_inj4",subregion_Cmin,subregion_Cmax);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    
if(INJ_SKIP_ON)Do_Set_INJ_SKIP();
   //===============================// calculate Mass for unknown ions //========================================
                                             //计算未知核质量 
TH1F* h1_refs = new TH1F("h1_refs","h1_refs",31,-0.5,30.5);
AxisFormat(h1_refs," "," reference ions used ", " count ");
TH1F* h_dC_ref_target = new TH1F("h_dC_ref_target","h_dC_ref_target",500,0,0.5);
AxisFormat(h_dC_ref_target, " #DeltaC of ref and this ion ", " #DeltaC [m] ", " count ");

int gr_err_dm_all_n=0;
TGraphErrors* gr_err_dm_all = new TGraphErrors();
AxisFormat(gr_err_dm_all,strtmp.Format(" error vs dm MASS_VER_%d",MASS_VER)," dm [keV] ", " error [keV]");
gr_err_dm_all->SetMarkerSize(1.2);
gr_err_dm_all->SetMarkerStyle(24);
TH2F* h2_err_dm_all = new TH2F("h2_err_dm_all","h2_err_dm_all", 100,-500,500,  200,0,2000);
AxisFormat(h2_err_dm_all,strtmp.Format(" h2 error vs dm MASS_VER_%d",MASS_VER)," dm [keV] ", " error [keV]");

TH2F* h2_massv1_C_all = new TH2F("h2_massv1_C_all","h2_massv1_C_all", 250,128.5,19.0,  200,-1000,1000);

calculated_unknown_ions_n = 0;
calculated_unknown_ions_n_v_others = 0;

ref_n=0;                        //how many ref ions will be used
j0 =0;                          //the starting subfix of ref ions
j_n=0;                          //the subfix of ref ions that put into new array
Bp_err2_1=0.0;                //sigma_Bpi square calculated using one ref, produced by AME_err

isomer_this_injection=0;
isomer_appeared=1;                      //isomer species appeared

grerr_n=0;

if(OUTPUT_singleion)
{
filename1=FILEPATH+"MassResult_ions_all"+strtmp.Format("%.3f_%.3f.txt",L,ddT);
outfile_single_ion.open( filename1);
outfile_single_ion<<"L= "<<L<<" ddT= "<<ddT<<" scan_k_ddtC= "<<scan_k_ddtC<<endl;
}

//----- 20230706 add v2 err gtC
if(OUTFILE_v2err_gtC_ON)outfile_v2error_gtC.open(FILEPATH+"v2error_gtC.txt");
if(INJ_sample_ON&&Do_PD_convergence_test) outfile_PD.open(FILEPATH+"PD_convergence_test.txt");
if(outfile_ERR_ANA_ON&&MASS_VER>2)
{
    outfile_ERR_ANA.open(FILEPATH+strtmp.Format("ERR_ANA_v%d_refweighted_%d.txt",MASS_VER, scan_loop_i));
}


int count_no_ref = 0;
int count_outof_gtC_ISF =0;
int count_outof_filter_chosenC = 0;
int count_outof_Cfilter_12N,count_outof_Cfilter_18Ne = 0;

int count_use_inj_1 =0;
int count_use_inj_2 =0;
int count_use_inj_3 =0;
int count_use_inj_4 =0;

cout<<endl<<"-------------------UNKNOWN ----------"<<endl
<<"use_gtC_type= "<<use_gtC_type<<endl
<<"---------"<<" MASS VERSION = "<< MASS_VER <<" : use_weighted_ON = "<<USE_WEIGHTED<<endl;
if(use_gtC_type==1)cout<<"fit on: "<<gtC_fit_ON<<" poln = "<<gtC_fit_pol<<endl;
if(C_filter_ON)cout<<"C_filter_ON on: "<<fixed<<setprecision(3)<<C_filter_min<<" ~ "<<C_filter_max<<endl;
cout<<" gt range :"<<gtC_chosen_gtMIN<<" ~ "<<gtC_chosen_gtMAX<<endl<<" Z min : >="<<choose_largeZ_min<<" Z max : <="<<choose_largeZ_max<<endl;
if(INJ_sample_ON)cout<<" \033[33m - - - - notice: INJ_sample_ON - - - - \033[0m "<<endl;

outfile_logs<<endl<<"-------------------UNKNOWN -------------------" <<endl
<<"use_gtC_type= "<<use_gtC_type<<endl
<<"---------"<<" MASS VERSION = "<< MASS_VER <<" : use_weighted_ON = "<<USE_WEIGHTED<<endl;
if(use_gtC_type==1)outfile_logs<<"fit on: "<<gtC_fit_ON<<" poln = "<<gtC_fit_pol<<endl;
if(C_filter_ON)outfile_logs<<"C_filter_ON on: "<<fixed<<setprecision(3)<<C_filter_min<<" ~ "<<C_filter_max<<endl;
outfile_logs<<" gt range :"<<gtC_chosen_gtMIN<<" ~ "<<gtC_chosen_gtMAX<<endl<<" Z min : >="<<choose_largeZ_min<<" Z max : <="<<choose_largeZ_max<<endl;
if(INJ_sample_ON)outfile_logs<<" \033[33m - - - - notice: INJ_sample_ON - - - - \033[0m "<<endl;

int DEALING_LINE = 0;
int SKIP_IONS_UNKNOWN_COUNT =0;
if(dA0T_ions_filter_ON){cout<<"----------- dA0T_ions_filter_ON -------------"<<endl;}
/*
//新增部分，计算 T v Bρ C误差平方以及协方差，最终计算每个核的影响因子    cal计算值  Sigma2方差 Cov协方差

    //思路：从 Cr实际值 -> Cr中心值 ----------------------> Ct中心值 -> Ct实际值
//误差影响：            dFr           γt(C)曲线的误差               dFt
//这里默认 Cr < Ct，反过来也差不多
double cal_T = 0, sigma2_T = 0, cal_v = 0, sigma2_v = 0, cov_T_v = 0, cal_gamma = 0;
double cal_Brou = 0, sigma2_Brou = 0, cal_C = 0, sigma2_C = 0, cov_Brou_C = 0;
double cov_C_v = 0;
*/
double dF_Cr_L = 0, dF_Cr_R = 0, dF_Ct = 0, dF_Ct_L = 0, dF_Ct_R = 0, dF_Cr = 0, dCr = 0, dCt = 0;//由于C误差导致的误差 （左 右 作为目标核 作为参考核）
double dFr = 0, dFt = 0;

for (int i = 0; i < ions_n; i++)
{
    /**/
    //cal_T = ions[i].A1 - 0.5 * ions[i].dA1;//ps
    //sigma2_T = ions[i].A1err + 0.25 * ions[i].dA1err - ions[i].cov16;//ps^2
    //cal_v = L / (ions[i].dA0 + ddT * 1000);//m/ps
    //sigma2_v = L * L / (pow((ions[i].dA0 + ddT * 1000), 4)) * ions[i].dA0err;//(m/ps)^2
    cal_gamma = 1 / sqrt(1 - pow((ions[i].v), 2) / pow(V_c, 2));
    //cov_T_v = L / (pow((ions[i].dA0 + ddT * 1000), 2)) * (-ions[i].cov15 + 0.5 * ions[i].cov56);//ps * m/ps

    //cal_Brou = (ionspecies[ions[i].Species].Mass / ions[i].Z) * cal_v * cal_gamma / V_c / V_c / 1000;//Tm   //M: keV  cal_v:m/ps  V_c:m/ns    
    //sigma2_Brou = pow((cal_Brou / ionspecies[ions[i].Species].Mass), 2) * pow(ionspecies[ions[i].Species].AME_err, 2) + cal_Brou * cal_Brou * pow(cal_gamma, 4) / cal_v / cal_v * sigma2_v;//Tm^2 这里考虑了质量的误差
    //cal_C = cal_v * cal_T;//m
    //sigma2_C = cal_v * cal_v * sigma2_T + cal_T * cal_T * sigma2_v + 2 * cal_C * cov_T_v;//m^2
    //cov_Brou_C = cal_Brou * cal_gamma * cal_gamma * cal_T / cal_v * sigma2_v + cal_Brou * cal_gamma * cal_gamma * cov_T_v;//Tm * m
    //cov_C_v = cal_T * sigma2_v + cal_v * cov_T_v;//m*m/ps

    //ions[i].v = cal_v * 1000.0;// [m / ns]
    //ions[i].Bp = cal_Brou;//[Tm]
    //ions[i].C = cal_C;//[m]
    //ions[i].T = cal_T/1000.0; //[ns]

    //ions[i].v_err2 = sigma2_v;
    //ions[i].Bp_err2 = sigma2_Brou;
    //ions[i].C_err2 = sigma2_C;
    //ions[i].cov_Brou_C = cov_Brou_C;
    //ions[i].cov_C_v = cov_C_v;


    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << ions[i].Bp << " " << cal_Brou  << endl;
    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << ions[i].v << " " << cal_v <<" "<< cal_gamma << endl;
    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << sqrt(ions[i].v_err2)<< " " << ions[i].v_err << endl;

    dF_Cr_L = Calculate_Unknown_Bp_dF(1, ions[i].C - sqrt(ions[i].C_err2), ions[i].C, gr_gtC_chosen, ISF, h_locate) - 1;// 无单位   这里是作为参考核的考虑，实际值到中心值
    dF_Cr_R = Calculate_Unknown_Bp_dF(1, ions[i].C + sqrt(ions[i].C_err2), ions[i].C, gr_gtC_chosen, ISF, h_locate) - 1;// 无单位

    if (abs(dF_Cr_L) < 0.9 && abs(dF_Cr_R) < 0.9)//都正常
    {
        if (abs(dF_Cr_L) > abs(dF_Cr_R)) { dF_Cr = dF_Cr_L; dCr = -sqrt(ions[i].C_err2); }// 无单位    取绝对值较大的那个作为 由于C误差导致的误差
        else { dF_Cr = dF_Cr_R; dCr = sqrt(ions[i].C_err2); }
    }
    else//有的核在边界附近，加上误差就出界了，导致计算结果为0
    {
        if (abs(dF_Cr_L) < 0.9) { dF_Cr = dF_Cr_L; dCr = -sqrt(ions[i].C_err2); }// 无单位    取有效的那个作为 由于C误差导致的误差
        else if (abs(dF_Cr_R) < 0.9) { dF_Cr = dF_Cr_R; dCr = sqrt(ions[i].C_err2); }
        else { dF_Cr = 0; dCr = sqrt(ions[i].C_err2); }//两个都为0，说明C中心值就已经出界了，后面不会以这个核为参考核，所以这里怎么处理都行
    }


    dF_Ct_L = 1.0 / (dF_Cr_L + 1.0) - 1.0; //这里是作为目标核时，中心值到实际值    dF_Cr与dF_Ct 绝对值差别  在千分之一以内?
    dF_Ct_R = 1.0 / (dF_Cr_R + 1.0) - 1.0; //这里是作为目标核时，中心值到实际值    dF_Cr与dF_Ct 绝对值差别  在千分之一以内?

    if (abs(dF_Ct_L) < 0.9 && abs(dF_Ct_R) < 0.9)//都正常
    {
        if (abs(dF_Ct_L) > abs(dF_Ct_R)) { dF_Ct = dF_Ct_L; dCt = -sqrt(ions[i].C_err2); }// 无单位    取绝对值较大的那个作为 由于C误差导致的误差
        else { dF_Ct = dF_Ct_R; dCt = sqrt(ions[i].C_err2); }
    }
    else//有的核在边界附近，加上误差就出界了，导致计算结果为0
    {
        if (abs(dF_Ct_L) < 0.9) { dF_Ct = dF_Ct_L; dCt = -sqrt(ions[i].C_err2); }// 无单位    取有效的那个作为 由于C误差导致的误差
        else if (abs(dF_Ct_R) < 0.9) { dF_Ct = dF_Ct_R; dCt = sqrt(ions[i].C_err2); }
        else { dF_Ct = 0; dCt = sqrt(ions[i].C_err2); }//两个都为0，说明C中心值就已经出界了，后面不会以这个核为参考核，所以这里怎么处理都行
    }



    //Bρ误差导致的误差约为 C 误差导致的误差大致相等，      Bρ（v也是）与C的相关系数接近于1，影响相互抵消，使得最终的dF因子很小       dF_Cr与dF_Ct
    ions[i].dFr = sqrt(ions[i].Bp_err2 / ions[i].Bp / ions[i].Bp + dF_Cr * dF_Cr + 2 * dF_Cr / dCr / ions[i].Bp * ions[i].cov_Brou_C);
    ions[i].dFt = sqrt(pow(cal_gamma, 4) / (ions[i].v/1000.0) / (ions[i].v/1000.0) * ions[i].v_err2 + dF_Ct * dF_Ct - 2 * dF_Ct / dCt * pow(cal_gamma, 2) / (ions[i].v/1000.0) * ions[i].cov_C_v);

    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << cov_Brou_C / sqrt(ions[i].C_err2) / sqrt(ions[i].Bp_err2) << " " <<sqrt(ions[i].Bp_err2) / cal_Brou <<" "<< sqrt(ions[i].C_err2) / cal_C << endl;


    //===============测试
    //testa = Calculate_Unknown_Bp_dF(1, 128.8, cal_C, gr_gtC_chosen, h_locate);
    //testb = Calculate_Unknown_Bp_dF(1, cal_C, 129.00, gr_gtC_chosen, h_locate);
    //testc = Calculate_Unknown_Bp_dF(1, 128.8, 129.00, gr_gtC_chosen, h_locate);
    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << testa << " " << testb <<" " << testc <<" " << testa * testb <<" "<< testc - testa * testb << endl;
    //在计算精度内， 任一段计算结果，严格等于内部 各分段结果相乘

    //if (ionspecies[ions[i].Species].IsRef == 1)
    //{
    //    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << ions[i].C << " " << cal_C << endl;
    //    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << ions[i].Bp << " " << cal_Brou << endl;
    //    //cout << cov_Brou_C <<"  " << cov_Brou_C / sqrt(ions[i].C_err2) / sqrt(ions[i].Bp_err2) << " " << cov_C_v <<"  "<< cov_C_v/ sqrt(ions[i].C_err2)/sqrt(sigma2_v) << endl;
    //    //cout << fixed << setiosflags(ios::left) << setw(15) << setprecision(14) << ions[i].dFr << " " << dF_Cr << " " << ions[i].dFt <<" " << dF_Ct << " " << sqrt(ions[i].Bp_err2) / cal_Brou << " " << sqrt(ions[i].Bp_err2) / cal_Brou / dF_Cr << " " << dF_Cr / dC << " " << dF_Ct / dC << endl;
    //    //cout << endl;
    //}
}
//double test[10];
//cout << test_C_min << " " << test_C_max << endl;
//TCanvas* c_sigma_C_test = new TCanvas("C_test", "C_test", 1200, 600);
//h_test->GetXaxis()->CenterTitle(1);
//h_test->GetXaxis()->SetTitle("sigma_C(m)");
//h_test->GetYaxis()->CenterTitle(1);
//h_test->GetYaxis()->SetTitle("Counts");
//h_test->Draw();
//新增部分结束



//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//正式开始计算
for(int i=0;i<ions_unknown_n  ;i++)//每个粒子，计算质量。
{
    if(i%1000==0){cout<<"now dealing with No. \033[4m\033[32m"<<i<<"\033[0m"<<endl<<"\033[1F"; }  //cmd 带颜色显示进度
    //if(i<800){    if(i==799){cout<<endl<<endl<<"!!!!!!  skip i<800 "<<endl<<endl;}   continue;}
    if(INJ_sample_ON)
    {   
        if( !(ions[i].inject_number%500==0) ) {continue;}
    }
    
    
    //###########################
    ions_unknown=ions[i];                             //运算符重载，此处是将该粒子的信息完全复制过去
    //############################
    int ts=ions_unknown.Species;  // for brevity, ts :this species 频繁使用ionspecies[ts]     此核种类的ID

    double this_iont_mass_v1=0;      //等权最终这个未知核的质量
    double delta_this_iont_mass_v1=0;//等权最终这个未知核的质量 - AME
    double this_iont_mass_VE=0;      //无论什么MASS_VER>2，最终这个未知核的质量
    double this_iont_mass_VE_err=1.14514;  //最终这个未知核的质量误差
    double delta_this_iont_mass_VE=0;//最终这个未知核的质量 - AME
    double this_iont_mvq_VE=0;          //最终这个未知核的m/q with err
    double this_iont_mvq_VE_err=0.114514;
    double this_iont_mvq_v1=0;


    //-----------------------skip this ion
    //if(ions[i].A==9&&ions[i].Z==6){continue;}//9C
    if(THIS_EXP=="2021_36Ar_SET3")
    {
        if(MASS_VER==5)if(ions[i].A==14&&ions[i].Z==8){continue;}//14O            //这两个核未区分
        if(MASS_VER==5)if(ions[i].A==21&&ions[i].Z==12){continue;}//21Mg
    }
    
    if(C_filter_ON)//开启轨道上下限制，超出范围的粒子舍弃
    {
        if( ions[i].C>C_filter_max||ions[i].C<C_filter_min)                    //目标核超出轨道限制，跳过
        {
            count_outof_filter_chosenC++;
            continue;
        }
        if(THIS_EXP=="2021_36Ar_SET3")
        {
        if(ions[i].A==12&&ions[i].Z==7){if( ions[i].C<128.68){count_outof_Cfilter_12N++;continue;} }  //某些特定核 划定了特殊的区间
        if(ions[i].A==18&&ions[i].Z==10){if( ions[i].C>128.97){count_outof_Cfilter_18Ne++;continue;} }
        }
    }
    //20230612 do_dA0T_FLAG
    if(dA0T_ions_filter_ON){if(ions[i].Do_dA0_T_flag==0){continue;}}

    //20230620 time_divide 900000 for part2
    if(time_divide_filter_ON){if(ions[i].time>900000){continue;}}

    if(choose_inj_option==1){if(ions_unknown.inject_number>7000)continue;}//分段（未开启）
    if(choose_inj_option==2){if(ions_unknown.inject_number<=7000){continue;}if(ions_unknown.inject_number>14000){continue;}}
    if(choose_inj_option==3){if(ions_unknown.inject_number<=14000){continue;}if(ions_unknown.inject_number>21000){continue;}}
    if(choose_inj_option==4){if(ions_unknown.inject_number<=21000){continue;} }
    //if(ions_unknown.inject_number<=14000)continue;
    
    //-------- 20240703
    //if(INJ_SKIP_ON) { if(INJ_SKIP[ions_unknown.inject_number] ){continue;} }
    if(INJ_SKIP_ON) //跳过某些注入
    { 
        INJ_SKIP_THIS = 0;
        for(int i=0;i<SKIP_n  ;i++)
        {
            if(SKIP_FILENAME[i]==ions_unknown.inject_filename){INJ_SKIP_THIS=true; break;}//这次注入有问题
        }
        if(INJ_SKIP_THIS)//跳过这次注入
        {
            //cout<<" ------------- skip this inj_file: "<<ions_unknown.inject_filename<<endl;
            SKIP_IONS_UNKNOWN_COUNT++;
            continue;
        } 
    }
    if (ISF[h_locate->FindBin(ions[i].C)] < 0.001)//超出范围的目标核粒子，ISF为0
    {
        //cout<<" debug !! "<<i<<endl;
        count_outof_gtC_ISF++;
        continue;//超出范围的粒子，Bp为0

    }

    //--------------------- do not skip this
    ions_t=ions_unknown;  // 整体复制

    ions_unknown.M_AME=ionspecies[ts].Mass;               //这种核的AME 裸核质量值
    ions_unknown.M_AME_err=ionspecies[ts].AME_err;        //AME的误差
    if(outputswitch2)  
    {
    cout<<endl;cout<<"-------- ion_unknown[ "<<i<<" ]---------------" <<endl;
    ions_unknown.PrintInfo();
    }

    j0=Injection_subfix[ions_unknown.inject_number];      //Injection_subfix 这一次注入的第一个粒子的 序号

    //initialization
    ions_this_injection =  Injection[ions_unknown.inject_number];  // 得到这次注入一共多少个离子
    
    if(ions_this_injection>=ref_n_MAX)//单次注入离子太多了
    {   
        cout<<"at inject "<<ions_unknown.inject_number<<"------ions_this_injection= "<<ions_this_injection<<" >"<<ref_n_MAX <<"!!!!-----"<<endl;
        ERROR_FLAG[0]=true;
        break;
    }
    
    if(ions_this_injection==1)//就一个核，算不了
    {
        if(outputswitch2)cout<<" there is no ref ion to use , skip-------"<<endl;
                 //do not use this ion_unknown
        count_no_ref++;
        continue;
    }

    //1021 C filtering

    if( (ions_unknown.C>C_filter_max||ions_unknown.C<C_filter_min)&&C_filter_ON )//目标核的轨道超出正常范围
    {
        continue;
    }
    if(THIS_EXP=="2021_36Ar_SET3")
    {
        // 单独对一个离子 卡轨道  涉嫌违规操作 放弃
        //if(ionspecies[ts].Aname=="10C")
        //if( (ions_unknown.C>128.85||ions_unknown.C<128.65)&&C_filter_ON )//目标核的轨道超出正常范围
        //{continue;}

        //单独挑选 gt 
        if(ionspecies[ts].Aname=="10C")
        {
            if(ions_unknown.gammat>1.366||ions_unknown.gammat<1.31 )
            //if(ions_unknown.gammat>1.365)
            //if(ions_unknown.gammat_err>gtC_ERR_upper_bound )
            {continue;}
        }
    }
    
    gammat_ave=0.0;
    Bp_ave_ref=0.0;
    //using ref ions
    j_n=0;
    m_tmp=0;
    ions_unknown.BpSum=0.0;ions_unknown.BpSum2=0.0;
    ions_unknown.Bp_err2 =0.0;

    for(int j=0;j< ions_this_injection  ;j++)//针对此核ions_unknown  ，使用这次注入的每个参考核进行一次计算
    {//==============loop : j  ions in this injection=================
        //ref_n >= 2
        // find aveBp using all the ref ions in this injection except itself

        //----------------------------------------------skip
        //20230612 do_dA0T_FLAG
        if(dA0T_ions_filter_ON)
            { if(ions[j0+j].Do_dA0_T_flag==0){continue;}}

        if(ions_unknown.ion_number == ions[j0+j].ion_number)//不用自己做参考核
        {   
            if(outputswitch2)cout<<"--itself---at "<<j+1 <<endl;
            continue;                                               
        }            //do not use itself as reference
        
        
        
        // choose ref ion ---  MassUnknown
        if( ionspecies[ions[j0+j].Species].MassUnknown)//这种核质量未知
        {
            if(outputswitch2)cout<<"-----------this ion"<<ions[j0+j].A<<ions[j0+j].name<<" MassUnknown, -----"<<endl;
            continue;
        }
        // choose ref ion --- IsRef
        if(ionspecies[ions[j0+j].Species].IsRef==0)//不是参考核
        {
            if(outputswitch2)cout<<"-----------this ion"<<ions[j0+j].A<<ions[j0+j].name<<" IsRef=false, -----"<<endl;
            count_use_massknown_but_notREF ++;
            //outfile_logs<<" not Ref debug 20240625"<<ions[j0+j].A<<ions[j0+j].name<<" when at ios_t"<<ions_unknown.A<<ions_unknown.name<<endl;
            continue;
        }
        // C_filter
        if(C_filter_ON)
        {if( ions[j0+j].C>C_filter_max||ions[j0+j].C<C_filter_min)//轨道限制，偏太远的不要
        {continue;}}

        if (ISF[h_locate->FindBin(ions[j0 + j].C)] < 0.001)//超出范围的参考核粒子，ISF为0
        {
            continue;//超出范围的参考核粒子，Bp为0
        }

        if(THIS_EXP=="2021_36Ar_SET3")
        {
            // 单独对一个离子 卡轨道  涉嫌违规操作 放弃
            //if(ionspecies[ions[j0+j].Species].Aname=="10C")
            //if( (ions[j0+j].C>128.85||ions[j0+j].C<128.65)&&C_filter_ON )//目标核的轨道超出正常范围
            //{continue;}
    
            //单独挑选 gt 
            if(ionspecies[ions[j0+j].Species].Aname=="10C")
            {
                if(ions[j0+j].gammat>1.366||ions[j0+j].gammat<1.31 )
                //if(ions[j0+j].gammat>1.365)
                //if(ions[j0+j].gammat_err>gtC_ERR_upper_bound )
                {continue;}
            }
        }
        
        //--------------------------------------------------do not skip

        ////////////////////////////////////////////////////////////for each ref ion    剩下的核都正常，可以用作参考核，开始计算
        // 历史遗留，由于是最核心最简单的等权质量计算v1，保留整个计算过程未封装
        //use tmp value, short name
        ions_r[j_n]=ions[j0+j];            //计算一个目标核，需要这一次注入的所有参考核，现在构建一个 参考核 数组
        j_n++;

        Ct=ions_unknown.C;
        Ci=ions[j0+j].C;
        Ct_err=sqrt(ions_unknown.C_err2);
        Ci_err=sqrt(ions[j0+j].C_err2);
        Bpi=ions[j0+j].Bp;
        vi_err=ions[j0+j].v_err;   
        vi=ions[j0+j].v;          
        ref_m=  ionspecies[ions[j0+j].Species].Mass;       //!!
        ref_m_err=  ionspecies[ions[j0+j].Species].AME_err;   


        Bp_ave_ref+=ions[ j0+j ].Bp;//本次注入其他参考核的平均Bρ
        
        //################# get one Bp value  for unknown ion   ###################################################gt
        if(use_gtC_type==1)//使用γt(C)曲线
        {
            if(!READIN_gtC_chosen_ON)// 使用计算出的 γt(C)曲线 ！！！！！最主要的分支（默认算法，等权重）
            {   
                //################## ### USE GTC HERE #############################################
                //之前的等步数算法-- 应当替换成 ISF
                //ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct ,  gr_gtC_chosen, gammat_ave,100,ions_unknown.inject_number);
            
                //20240702 尝试每种核用自己单独的gt-C 曲线
                //ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct , ionspecies[ts].gr_gtC_shifted_own , gammat_ave,100,ions_unknown.inject_number);
                //================ 2024L.JH 使用 ISF ==================
                ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, gr_gtC_chosen, ISF, h_locate);   //ISF!!!!!（默认算法，等权重）

                //58Ni红色ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, grerr_avegtC_inj_1, ISF_inj1, h_locate_inj1);
                //58Ni蓝色ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, grerr_avegtC_inj_2, ISF_inj2, h_locate_inj2);
                
                /*   一种离子单独使用特定的 gtC curve  ---- 尚需验证合理性
                if(THIS_EXP=="2021_36Ar_SET3"&&TEST_USE_OWN_gtC_ON)
                {
                    if(ionspecies[ts].Aname=="10C"||ionspecies[ts].Aname=="13O")
                    {
                        //ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, ionspecies[ts].gr_gtC_own, ISF_1, h_locate);
                        ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, gr_gtC_chosen_shift_1, ISF_1, h_locate);
                    }
                    if(ionspecies[ts].Aname=="11C"||ionspecies[ts].Aname=="13N")
                    {
                        //ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, ionspecies[ts].gr_gtC_own, ISF_2, h_locate);
                        ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, gr_gtC_chosen_shift_2, ISF_2, h_locate);
                    }
                }
                */
            }
            else //使用读取的  γt(C)曲线
            {
                ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct ,  gr_gtC_readin, gammat_ave,100,ions_unknown.inject_number);
            }
        }
        
        else if(use_gtC_type==2)ions_unknown.Bp = Calculate_Unknown_Bp_gt0(Bpi,Ci, Ct ,100,gt0);
        else if(use_gtC_type==3)
        {
            int t_inj = ions_unknown.inject_number;
            if     (t_inj<=gtC_inj_max_1&&t_inj>=gtC_inj_min_1){ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct ,grerr_avegtC_inj_1, gammat_ave,100,ions_unknown.inject_number);  count_use_inj_1++;}
            else if(t_inj<=gtC_inj_max_2&&t_inj>=gtC_inj_min_2){ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct ,grerr_avegtC_inj_2, gammat_ave,100,ions_unknown.inject_number);count_use_inj_2++;}  
            //else if(t_inj<=gtC_inj_max_3&&t_inj>=gtC_inj_min_3){ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct ,grerr_avegtC_inj_3, gammat_ave,100,ions_unknown.inject_number);count_use_inj_3++;}  
            //else if(t_inj<=gtC_inj_max_4&&t_inj>=gtC_inj_min_4){ions_unknown.Bp = Calculate_Unknown_Bp_1(Bpi,Ci, Ct ,grerr_avegtC_inj_4, gammat_ave,100,ions_unknown.inject_number);count_use_inj_4++;}  
            else    {ions_unknown.Bp = Calculate_Unknown_Bp_dF(Bpi, Ci, Ct, gr_gtC_chosen, ISF, h_locate);}   //（默认算法，等权重）
        }
        else if(use_gtC_type==4)
        {
            ions_unknown.Bp = Calculate_Unknown_Bp_4(Bpi,Ci, Ct ,  fitfun_gr_gtC_all, gammat_ave,100,ions_unknown.inject_number);
        }    
        //ions_unknown.Bp = Calculate_Unknown_Bp(Bpi,Ci, Ct ,  gr_gtC_chosen_smooth, gammat_ave,100,ions_unknown.inject_number);  
        //ions_unknown.Bp = Calculate_Unknown_Bp(Bpi,Ci, Ct ,  sp3_smooth , gammat_ave,100,ions_unknown.inject_number);

        //20230223-check one mass result ==========///debug
        //ions_unknown.Calculate_M_err();
//cout<<"DEBUG "<<i<<" "<<ions_unknown.A<<ions_unknown.name<<" "<<ions[ j0+j ].name<<" Bp= "<<ions_unknown.Bp<<" one mass : "<< ions_unknown.M_cal-ions_unknown.M_AME <<endl;

        //############################################## 得到了 ions_unknown.Bp #########################################################################
        ///cout<<"ref "<<j+1<<" gt ave = "<<gammat_ave<<endl; 

        ions_unknown.BpSum += ions_unknown.Bp;                        //ions_unknown.Bp 使用此参考核得到的结果，根据一次注入的所有参考核取平均值
        ions_unknown.BpSum2 += ions_unknown.Bp*ions_unknown.Bp;       //计算误差用
        
    if(outputswitch2)cout<<"ref ion : "<<j0+j<<"  "<<ions[ j0+j ].name<<" has Bp = "<<ions[ j0+j ].Bp<<fixed<<setprecision(10)<<" this time Bp = "<<ions_unknown.Bp<<endl;
        
        ions_unknown.Calculate_M_err();       //计算 质量M_cal 和误差（误差部分被注释掉了，没有计算）
        //############################################## 得到了 ions_unknown.M_cal , .Mvq ############################################################
        if(c_dmdC_ON)
        {
            gr_dmdC->SetPoint(gr_dmdC->GetN(),ions_unknown.C-ions[ j0+j ].C , ions_unknown.M_cal-ions_unknown.M_AME);
            h2_dmdC->Fill(ions_unknown.C-ions[ j0+j ].C , ions_unknown.M_cal-ions_unknown.M_AME);
        }
        h_dC_ref_target->Fill(abs(ions_unknown.C-ions[ j0+j ].C) );
        //20230419 err ana
        if(ionspecies[ts].MassUnknown==0 && !LOOP_ON)  //MassUnknown==0 此核是已知核
        {
            grerr_ERRYX->SetPoint(grerr_ERRYX_n++,  
                                pow( 1.0/(1 - pow(ions_unknown.v/V_c,2))-pow(gr_gtC_chosen->Eval(ions_unknown.C) ,2) , 2)     //(γ^2-γt^2)^2 此核
                                + pow( 1.0/(1 - pow(ions[ j0+j ].v/V_c,2))-pow(gr_gtC_chosen->Eval(ions[ j0+j ].C) ,2) , 2) , // + (γ^2-γt^2)^2 本次计算用的参考核
                                 pow( (ions_unknown.M_cal-ions_unknown.M_AME)/ions_unknown.M_AME, 1 )                         // (M计算值 - M_AME)/M_AME
                                );
            h2_ERRYX->Fill(pow( 1.0/(1 - pow(ions_unknown.v/V_c,2))-pow(gr_gtC_chosen->Eval(ions_unknown.C) ,2) , 2)
                                + pow( 1.0/(1 - pow(ions[ j0+j ].v/V_c,2))-pow(gr_gtC_chosen->Eval(ions[ j0+j ].C) ,2) , 2) , 
                                 pow( (ions_unknown.M_cal-ions_unknown.M_AME)/ions_unknown.M_AME, 1 ));
        }

        //20230420 gtC_v2
        if(DO_gtC_v2_ON&& ionspecies[ts].MassUnknown==0 && !LOOP_ON)
        {
            grerr_gtC_v2->SetPoint(grerr_gtC_v2_n,(ions_unknown.C+ions[ j0+j ].C)/2 , 
                                  sqrt( (log(ions_unknown.Bp_true/ions[j0+j].Bp_true))/(log(ions_unknown.C/ions[j0+j].C)) ));
            grerr_gtC_v2_n++;
            h2_gtC_v2->Fill((ions_unknown.C+ions[ j0+j ].C)/2 , 
                                  sqrt( (log(ions_unknown.Bp_true/ions[j0+j].Bp_true))/(log(ions_unknown.C/ions[j0+j].C)) ));
        }
        
    }//__________loop : j  ions in this injection___________________  针对此核ions_unknown  ，使用这次注入的每个参考核进行一次计算

    if(outputswitch2)  cout<<"-----------"<<endl;
    ref_n = j_n;   //obtain the ref ion numbers j_n,ref_n start from 1
    if(ref_n<=0)     //参考核全都有问题，都去掉了，没有可用的参考核 
    {
        if(outputswitch2)cout<<"no ion can be used in this injection"<<endl;
        count_no_ref++;
        continue;     //no ref to use for this unknown ion, skip to next unknown ion
    }

    h1_refs->Fill(ref_n);
    //################### get aveBp ,unknown ion #############################                                        
    ions_unknown.Bp=ions_unknown.BpSum /ref_n;    //set Bp to be Bp_ave    得到目标核的Bρ（各个参考核给出Bρ的 等权重 平均值）
    Bp_ave_ref/=ref_n;//参考核的平均Bρ
    if(outputswitch2)
    {
    cout<<"###################"
    <<fixed<<setprecision(10)<<"Bp_ave= "<<ions_unknown.Bp<<"  Bp_true= "<<ions_unknown.Bp_true<<endl;
    cout<<"delta Bp="<< (ions_unknown.Bp - ions_unknown.Bp_true)<<endl;
    }     
    //以上是历史遗留的 未经整合的 MASS_VER = 1 先得到 平均Bp 的过程       
//================================= 开始计算质量 MASS_VER ========================================================================

    //##################!! get unknown mass here!!!!!!!! #########################
    //if(MASS_VER==1) always 常驻
    {
        ions_unknown.Calculate_M_err();//calculate M (using Bp_ave)and M_err from AME_err
        this_iont_mass_v1 = ions_unknown.M_cal;
        delta_this_iont_mass_v1 = this_iont_mass_v1- ions_unknown.M_AME;  //ions_unknown.M_AME = ionspecies[i].Mass
        this_iont_mvq_v1 = ions_unknown.Mvq;
        ///cout<<ions_unknown.M_cal<<endl;
    }
    if(MASS_VER==2)
    {
       /////fit histogram:  this will be done after all the results from one species are filled into h1
       ///// based on single mass of MASS v1, no single mass value here for v3   
    }
    if(MASS_VER>=3)
    {// 有单个目标核的质量误差棒
        if(READIN_MERR_ON)
        {
            //read in err 外部读入每个核的质量误差，不需要计算
            this_iont_mass_VE = ReadIn_m_v2[ions_unknown.ID];
            this_iont_mass_VE_err = ReadIn_merr[ions_unknown.ID]; 
            /*
            outfile_errors_check<<ions_unknown.ID<<" "<<ions_unknown.ion_number<<" "<<ions_unknown.inject_number
            <<" "<<ions_unknown.A<<" "<<ions_unknown.Z<<" "<<this_iont_mass_VEerr <<endl;
            */
        }
        else
        {
            // 各种误差版本
            if(MASS_VER==3)//v3算法，各个参考核给出的误差不一致，采用加权平均
            {
                //#############!! get unknown mass  v3 here!!!!!!!!#########################################
                Calculate_Mass_with_err(ions_unknown , ions_r,ref_n,gr_gtC_chosen,gr_gtC_chosen_u,gr_gtC_chosen_d, this_iont_mass_VE,this_iont_mass_VE_err);
                //-------- this_iont_mass_VE is set --------
                
            }
            if(MASS_VER==4)//regard each calculation using one ref as independent 20231023
            {
                Calculate_Mass_with_err_v4(ions_unknown , ions_r,ref_n,gr_gtC_chosen,gr_gtC_chosen_u,gr_gtC_chosen_d, this_iont_mass_VE,this_iont_mass_VE_err);
                // v4 this function will directly set ions_unknown.M_cal_VE and this_iont_mass_VE
                
            }
            if(MASS_VER ==5)
            {   
                Calculate_iont_Mass_useDY(ions_unknown , ions_r,ref_n,gr_gtC_chosen,gr_gtC_chosen_u,gr_gtC_chosen_d, this_iont_mass_VE,this_iont_mass_VE_err);
                //debug!! if(i%100==0)cout<<" debug MASS_VER=5: dm = "<<delta_this_iont_mass_VE<<" +- "<<this_iont_mass_VE_err<<endl;
            }
            if (MASS_VER == 6)
            {
                Calculate_Mass_with_err_dF_all(ions_unknown, ions_r, ref_n, gr_gtC_chosen, gr_gtC_chosen_u, gr_gtC_chosen_d, h_locate, this_iont_mass_VE, this_iont_mass_VE_err);//##mass_dF 
            }
            
            // 得到 this_iont_mass_VE,this_iont_mass_VE_err
            // 整合 所有版本结果到 _VE 统一名称变量
            delta_this_iont_mass_VE   = this_iont_mass_VE - ions_unknown.M_AME;
            ions_unknown.M_cal_VE     = this_iont_mass_VE;
            ions_unknown.M_cal_err_VE = this_iont_mass_VE_err;
            ions_unknown.Mvq_VE       = this_iont_mass_VE/(ions_unknown.Z*u);
            ions_unknown.Mvq_err_VE   = this_iont_mass_VE_err/(ions_unknown.Z*u);
            this_iont_mvq_VE    = ions_unknown.Mvq_VE;
            this_iont_mvq_VE_err= ions_unknown.Mvq_err_VE;

        }
    }
    //else{ cout<<"error!! MASS_VER  is wrong  !!"<<endl; cout<<"return at MASS_VER error! "<<endl;return;}
    //Calculate_iont_Mass_v1(ions_unknown ,ions_r , ref_n,gr_gtC_chosen , 0 , 0);
    
//cout<<" [[Bp= "<<ions_unknown.Bp<<" ]]] "<<endl;
//cout<<" [[[[[[ "<<ions_unknown.M_cal<<" ]]] "<< delta_this_iont_mass_v1 << " ---- " <<( m_tmptmp-ions_unknown.M_AME)<<" +- "<<m_tmptmperr<< endl;  ///debug
//cout<<" debug!!  delta_this_iont_mass_VE= "<< delta_this_iont_mass_VE<<" +- "<< this_iont_mass_VE_err <<endl;  ///debug
    

    //____________________________________   this_iont_mass is set   __________________________________________________________________________



    if(!LOOP_ON)//记录数据，非循环情况下做一些统计分析
    {
        h_Mvq->Fill(this_iont_mvq_v1);
        h_Mvq_AME->Fill(ionspecies[ts].Mvq_AME);

        if(c_dmISO_ON)
        {
            if(ionspecies[ts].Aname=="14O")continue;
            if(ionspecies[ts].Isomer_n>0)continue;
            gr_dmISO->SetPoint(gr_dmISO->GetN(), pow(ions_unknown.Calculate_only_return_gamma(),2) - pow( gr_gtC_chosen->Eval(ions_unknown.C),2 ) , 
                ions_unknown.Mvq- ionspecies[ts].Mvq_AME);
            gr_TerrISO->SetPoint(gr_TerrISO->GetN(), pow(ions_unknown.Calculate_only_return_gamma(),2) - pow( gr_gtC_chosen->Eval(ions_unknown.C),2 ) , 
                ions_unknown.T_err);
            gr_verrISO->SetPoint(gr_verrISO->GetN(), pow(ions_unknown.Calculate_only_return_gamma(),2) - pow( gr_gtC_chosen->Eval(ions_unknown.C),2 ) , 
                ions_unknown.v_err);
            h2_dmvqISO->Fill(pow(ions_unknown.Calculate_only_return_gamma(),2) - pow( gr_gtC_chosen->Eval(ions_unknown.C),2 ) , 
                ions_unknown.Mvq- ionspecies[ts].Mvq_AME);
        }

        if(MASS_VER>=3)
        {
            if(Do_iont_error_out_ON&&MASS_VER>=3&&!READIN_MERR_ON)
            {
                outfile_errors<<ions_unknown.ID<<" "<<ions_unknown.ion_number<<" "<<ions_unknown.inject_number
                <<" "<<ions_unknown.A<<" "<<ions_unknown.Z<<" "<<this_iont_mass_VE<<" "<<this_iont_mass_VE_err<<endl;
            }
            if(c_dmdC_ON)
            {
                grerr_dmC_VE_all->SetPoint(grerr_dmC_VE_all_n,ions_unknown.C , delta_this_iont_mass_VE);
                grerr_dmC_VE_all->SetPointError(grerr_dmC_VE_all_n, 0 , ions_unknown.M_cal_err_VE);
                grerr_dmC_VE_all_n++;
            }
            
            gr_dm_dA0err_all->SetPoint(gr_dm_dA0err_all_n++,ions_unknown.T_err, delta_this_iont_mass_VE);
    
            gr_err_dm_all->SetPoint(gr_err_dm_all_n++,delta_this_iont_mass_VE,this_iont_mass_VE_err);
            h2_err_dm_all->Fill(delta_this_iont_mass_VE,this_iont_mass_VE_err);
        }
    }
    
//if(i%200==0)cout<<"DEBUG delta mass: "<<i<<" "<<abs(delta_this_iont_mass_v1)<<endl;
    

    if(abs(delta_this_iont_mass_v1)>500)//计算值与理论值差别过大，文件输出记录
    {
        if(ionspecies[ts].Isomer_n==0)outfile_large_dm<<delta_this_iont_mass_v1<<" "<<ions_unknown.A<<ions_unknown.name<<" inj: "<<ions_unknown.inject_number<<" "<<ions_unknown.inject_filename <<endl;
        else outfile_large_dm<<"   ##isomer "<<delta_this_iont_mass_v1<<" "<<ions_unknown.A<<ions_unknown.name<<" inj: "<<ions_unknown.inject_number<<" "<<ions_unknown.inject_filename <<endl;
    }

    ions_unknown.M_cal_err=0;
    
   

    
    //if(abs(ions_unknown.M_cal_VE-ions_unknown.M_AME)>500)
    //{outfile_large_dm<<"v2: "<<ions_unknown.M_cal_VE-ions_unknown.M_AME<<" "<<ions_unknown.A<<ions_unknown.name<<" inj: "<<ions_unknown.inject_number<< endl;}

//cout<<"DEBUG: "<<delta_this_iont_mass_v1<<endl;

    // if( abs(delta_this_iont_mass_v1) > 500 )
    // outfile1<<"####large deviation!! "<<ions_unknown.name<<" "<<delta_this_iont_mass_v1<<" "<<ions_unknown.inject_number<<" "<<ions_unknown.ion_number<<endl;

    //==============================================  ===============================================================================
    
    //if(abs(delta_this_iont_mass_v1)<1000)
    //!!!!!!!!!!!!!! // DM filter
    if(  DM_FILTER )
    {
        //if(ionspecies[ts].Aname=="10C"){if(delta_this_iont_mass_v1>350 ){DM_FILTER_n++; continue;}} 
        if(ionspecies[ts].Aname=="20Mg"){if(delta_this_iont_mass_v1<-400 ){DM_FILTER_n++; continue;}} 
        if(ionspecies[ts].Aname=="12N"){if(delta_this_iont_mass_v1<-300 ){DM_FILTER_n++; continue;}} 
        //if(ionspecies[ts].Aname=="17Ne"){if(delta_this_iont_mass_v1<-1000 ){DM_FILTER_n++; continue;}} 
    }// DM filter




    Ave_deltaM_all += delta_this_iont_mass_v1;             //所有核 统计计算值与AME的差别，此处为等权重结果
    calculated_unknown_ions_n++;                          // 记录计算的质量结果数  
    Injection_m[ions[i].inject_number]++;                 // 记录有质量结果的离子所在注入
    
    ionspecies[ts].HasResult=true;
    ionspecies[ts].N_unknown++;

    { // MASS_VER==1 常驻 always calculate massver1 : unweighted average, no error bar
        
        ionspecies[ts].vector_dm_v1.push_back(delta_this_iont_mass_v1);
        ionspecies[ts].Mass_cal +=ions_unknown.M_cal;
        ionspecies[ts].stdDeviation_V1 +=ions_unknown.M_cal*ions_unknown.M_cal;
        ionspecies[ts].Mvq_cal+=this_iont_mvq_v1;
        ions[i].mvq_v1 = ions_unknown.Mvq;
    }

    if(MASS_VER>=3)//以下为其他算法的结果
    {
        //if (ions_unknown.M_cal_VE > 0.1)
        {
        Ave_deltaM_all_v_others += ions_unknown.M_cal_VE-ions_unknown.M_AME;
        calculated_unknown_ions_n_v_others++;
        ionspecies[ts].Mass_cal_VE     +=1/(ions_unknown.M_cal_err_VE*ions_unknown.M_cal_err_VE)*ions_unknown.M_cal_VE;
        ionspecies[ts].Mass_cal_err_VE_ave +=1/(ions_unknown.M_cal_err_VE*ions_unknown.M_cal_err_VE);
        ionspecies[ts].stdDeviation_VE +=ions_unknown.M_cal_VE*ions_unknown.M_cal_VE;
        
        ions[i].m_VE = ions_unknown.M_cal_VE;
        ions[i].m_VE_err = ions_unknown.M_cal_err_VE;
        ionspecies[ts].h_iont_mass_err->Fill(ions_unknown.M_cal_err_VE);
        }
    }
    
    
    //20221016 gr_mvqC 每种 质荷比结果-C散点图
    ionspecies[ts].gr_mvqC->SetPoint(ionspecies[ts].gr_mvqC->GetN(),ions_unknown.C,this_iont_mvq_v1);
    //20240902h_mvq
    ionspecies[ts].h_mvq->Fill(this_iont_mvq_v1);
    if(MASS_VER>=3)
    {
        ionspecies[ts].grerr_mvqC_VE->SetPoint(ionspecies[ts].grerr_mvqC_VE->GetN(), ions_unknown.C,this_iont_mvq_VE);
        ionspecies[ts].grerr_mvqC_VE->SetPointError(ionspecies[ts].grerr_mvqC_VE->GetN()-1, 0,    this_iont_mvq_VE_err);
    }
    
    //20230814 gr_dmC
    ionspecies[ts].gr_dmC->SetPoint(ionspecies[ts].gr_dmC->GetN(),           ions_unknown.C, delta_this_iont_mass_v1);
    ionspecies[ts].grerr_dmC_VE->SetPoint(ionspecies[ts].grerr_dmC_VE->GetN(), ions_unknown.C, delta_this_iont_mass_VE  );
    ionspecies[ts].grerr_dmC_VE->SetPointError(ionspecies[ts].grerr_dmC_VE->GetN()-1,  0,   this_iont_mass_VE_err);
    
    //20240902 h_dm for showingg and also to be used in v2
    // fill h_dm after F-D bin width is determined
    //ionspecies[ts].h_dm->Fill(delta_this_iont_mass_v1);
    if(DoShow_mvqC_each_h2_ON)ionspecies[ts].h2_mvqC->Fill(ions_unknown.C,ions_unknown.Mvq);

    
    
    ///ions_unknown.PrintMassCalInfo();
    if(OUTPUT_singleion)
    {
    outfile_single_ion<<ions_unknown.A<<" "<<ions_unknown.Z<<" "<<ions_unknown.name<<" "<<ions_unknown.inject_number<<" "<<ions_unknown.ion_number
    <<" "<<ions_unknown.T<<"  "<<ions_unknown.v<<"  "<<ions_unknown.C<<" "<<ions_unknown.Bp
    <<" delta m= "<<fixed<<setprecision(10)<<delta_this_iont_mass_v1<<" +- "<<ions_unknown.M_cal_err
    <<" m= "<<ions_unknown.M_cal<<" mvq= "<<ions_unknown.Mvq<<endl;
    }
    //20240710
    if(outfile_each_mass_data_ON&&!LOOP_ON)
    {
        ionspecies[ts].outfile_mass<<ions_unknown.A<<" "<<ions_unknown.Z<<" "<<ions_unknown.inject_number<<" "<<ions_unknown.ion_number<<" "<<ions_unknown.inject_filename<<" "
        <<ions_unknown.C<<" "<<delta_this_iont_mass_v1<<" "<<ions_unknown.M_cal_err<<" "
        <<fixed<<setprecision(10)<<ions_unknown.Mvq<<" "<<ionspecies[ts].Mvq_AME<<" ";
        if(MASS_VER>=3)
        {
            ionspecies[ts].outfile_mass<<ions_unknown.M_cal_VE-ions_unknown.M_AME<<" "<<ions_unknown.M_cal_err_VE<<" "
            <<fixed<<setprecision(10)<<ions_unknown.Mvq_VE<<" ";
        }

        ionspecies[ts].outfile_mass<<endl;

    }
}  //end of : for(int i=0;i<ions_unknown_n  ;i++)               每个目标离子都计算完毕了
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-




if(OUTPUT_singleion)outfile_single_ion.close();

//---- 20230706 add v2 err gtC
if(OUTFILE_v2err_gtC_ON)outfile_v2error_gtC.close();
if(INJ_sample_ON&&Do_PD_convergence_test)outfile_PD.close();
if(outfile_ERR_ANA_ON&&MASS_VER>1)outfile_ERR_ANA.close();
//outfile1.close();



//==================================== COUT =============================
cout<<endl<<"----------- calculated ions : \033[31m"<<calculated_unknown_ions_n<<" / "<<ions_unknown_n
    <<"  ("<<double(calculated_unknown_ions_n)/double(ions_unknown_n)*100<<" %  )\033[0m-----------------"<<endl;

Ave_deltaM_all/=calculated_unknown_ions_n;
if(MASS_VER>1)Ave_deltaM_all_v_others/=calculated_unknown_ions_n_v_others;

cout<<"------------Ave_deltaM_all = \033[32m"<<Ave_deltaM_all<<"\033[0m-----------"<<endl;
if(MASS_VER>1)cout<<"------------Ave_deltaM_all_v"<<MASS_VER<<" = \033[32m"<<Ave_deltaM_all_v_others<<"\033[0m-----------"<<endl;

cout<<"------------count no ref = "<<count_no_ref<<"-----------"
<<"---count_outof_gtC_ISF  = "<<count_outof_gtC_ISF<<"-----------"<<endl;
if(C_filter_ON)
{
    cout<<"------!!C - filter - ON!!----- "<<C_filter_min<<" ~ "<<C_filter_max<<endl;
    cout<<"---count_outof_filter_C  = "<<count_outof_filter_chosenC<<"-----------"<<endl;  
    if(THIS_EXP=="2021_36Ar_SET3")
    {
        cout<<"in 36Ar-set3 skip---12N  = "<<count_outof_Cfilter_12N<<"-----------"<<endl;  
        cout<<"in 36Ar-set3 skip---18Ne  = "<<count_outof_Cfilter_18Ne<<"-----------"<<endl;       
    }  
    
}

cout<<"--20240626debug! : count_use_massknown_but_notREF: "<<count_use_massknown_but_notREF<<endl;
if(INJ_SKIP_ON)cout<<"--SKIP_IONS_UNKNOWN_COUNT by inject_filename: "<<SKIP_IONS_UNKNOWN_COUNT<<endl;
//==================================== logs =============================

outfile_logs<<endl<<"calculated ions : "<<calculated_unknown_ions_n<<" / "<<ions_unknown_n
    <<"  ("<<double(calculated_unknown_ions_n)/double(ions_unknown_n)*100<<" %  )-----------------"<<endl;

outfile_logs<<"---Ave_deltaM_all = "<<Ave_deltaM_all<<endl;
outfile_logs<<"------------Ave_deltaM_all_v"<<MASS_VER<<" = "<<Ave_deltaM_all_v_others<<"-----------"<<endl;

outfile_logs<<"---count no ref = "<<count_no_ref
<<"---count_outof_gtC_ISF  = "<<count_outof_gtC_ISF<<endl;
if(C_filter_ON)
{
    if(scan_loop_i<2)outfile_logs<<"C - filter - ON "<<C_filter_min<<" ~ "<<C_filter_max<<endl;
    outfile_logs<<"---count_outof_filter_C  = "<<count_outof_filter_chosenC<<endl;    
}
if(INJ_SKIP_ON)
{
    outfile_logs<<" INJ_SKIP_INFILENAME = "<<INJ_SKIP_INFILENAME<<endl;
    outfile_logs<<"--SKIP_IONS_UNKNOWN_COUNT by inject_filename : "<<SKIP_IONS_UNKNOWN_COUNT<<endl;
}
//20230707
if(!LOOP_ON &&(MASS_VER==3))
{
cout<<" $$$$3 check ----"<<endl
<<" count_gtC_masserr_condition_1= "<<count_gtC_masserr_condition_1<<" count_gtC_masserr_condition_2= "<<count_gtC_masserr_condition_2<<endl;
}

if(DM_FILTER)cout<<" DM_FILTER_ON: N="<<DM_FILTER_n<<endl;
if(TEST_USE_OWN_gtC_ON)cout<<" ===== TEST_USE_OWN_gtC_ON ====="<<endl;
//-------------------------------   IONSPECIES    ----------------------------


Ave_deltaM_all=0;

//initialization
if(Do_mass_err_scatter_ON&&MASS_VER>=3){outfile_sca_ave.open(FILEPATH+strtmp.Format("SCA_AVE_%d.txt",MASS_VER) );}


for(int i = 0; i<NSpecies; i++)
{ 
    if(ionspecies[i].HasResult)//各种核素依次进行计算
  
    {   
//------------ nuclear mass ---------------------------        
        ionspecies[i].Mass_cal /= ionspecies[i].N_unknown;
        ionspecies[i].Mass_cal_VE /= ionspecies[i].Mass_cal_err_VE_ave;
//------------ atomic mass excess---------------------------        

                       ionspecies[i].MassExcess_cal = ionspecies[i].GetMassExcess_cal(ionspecies[i].Mass_cal);
    
        if(MASS_VER>=3)ionspecies[i].MassExcess_cal_VE = ionspecies[i].GetMassExcess_cal(ionspecies[i].Mass_cal_VE);
//------------ mvq ----------------------------------------
        ionspecies[i].Mvq_cal /= ionspecies[i].N_unknown;
//------------ mass error ----------------------------------------
        
        ionspecies[i].stdDeviation_V1 = sqrt ( ionspecies[i].stdDeviation_V1/ionspecies[i].N_unknown - ionspecies[i].Mass_cal*ionspecies[i].Mass_cal );
        //ionspecies[i].Mass_cal_err =ionspecies[i].stdDeviation/sqrt(ionspecies[i].N_unknown-1);
        
        ionspecies[i].Mass_cal_err_VE_ave = sqrt(1/(ionspecies[i].Mass_cal_err_VE_ave)); 

        //====比较scattering error/average error选择大的作为Mass_cal_err_VE
        if(Do_mass_err_scatter_ON&&MASS_VER>=3)
        {
        
            for(int j = 0; j<ions_n; j++)
            {
                if(ions[j].Species!=i){continue;}
                if(ions[j].m_VE>0&&ions[j].m_VE_err>0){ionspecies[i].Mass_cal_err_VE_sca += pow(  (ions[j].m_VE-ionspecies[i].Mass_cal_VE)/ions[j].m_VE_err ,2);}            
            }
            if (ionspecies[i].N_unknown>1)
            {
                ionspecies[i].Mass_cal_err_VE_sca = sqrt(ionspecies[i].Mass_cal_err_VE_sca)/sqrt(ionspecies[i].N_unknown - 1)*ionspecies[i].Mass_cal_err_VE_ave;
            }
            else
            {
                ionspecies[i].Mass_cal_err_VE_sca = 0;
            }
        }
        if(ionspecies[i].Mass_cal_err_VE_sca>=ionspecies[i].Mass_cal_err_VE_ave)
        { ionspecies[i].Mass_cal_err_VE = ionspecies[i].Mass_cal_err_VE_sca;}
        else                                                                    
        { ionspecies[i].Mass_cal_err_VE = ionspecies[i].Mass_cal_err_VE_ave;}

        outfile_sca_ave<<ionspecies[i].Aname<<" "<<ionspecies[i].A<<" "<<ionspecies[i].Z<<" "<<ionspecies[i].N_unknown<<" "
        <<fixed<<setprecision(1)
        <<ionspecies[i].Mass_cal_err_VE_sca<<" "<<ionspecies[i].Mass_cal_err_VE_ave<<" "<<ionspecies[i].Mass_cal_err_VE<<" "
        <<(ionspecies[i].Mass_cal_err_VE_sca>=ionspecies[i].Mass_cal_err_VE_ave)<<endl;
        //______比较scattering error/average error_________ 

//------------deltaMass: mass - AME  ----------------------------------------
        ionspecies[i].deltaMass = ionspecies[i].Mass_cal - ionspecies[i].Mass ; //总的裸核质量差
        ionspecies[i].deltaMass_VE = ionspecies[i].Mass_cal_VE - ionspecies[i].Mass ;        
        if(ionspecies[i].N_unknown>1)
        {
            ionspecies[i].Mass_cal_err = ionspecies[i].stdDeviation_V1/sqrt(ionspecies[i].N_unknown - 1);
            ionspecies[i].deltaMass_err = sqrt( pow(ionspecies[i].Mass_cal_err,2)+pow(ionspecies[i].AME_err,2)  );
        }
        else 
        {
            ionspecies[i].deltaMass_err = 0; 
            ionspecies[i].Mass_cal_err = 0;
        }
        ionspecies[i].deltaMass_err_VE = sqrt( pow(ionspecies[i].Mass_cal_err_VE,2)+pow(ionspecies[i].AME_err,2)  );
        ionspecies[i].Mvq_cal_err = ionspecies[i].Mass_cal_err/ionspecies[i].Z/u;

        /*
        if(OUTPUT_MASSRESULT)
        {
        <<ionspecies[i].A<<" "<<ionspecies[i].name<<" "<<ionspecies[i].Z<<" "<<" | "<<ionspecies[i].MassUnknown<<" | "
        <<fixed<<setprecision(3)<<ionspecies[i].AveT
        <<" | nuc_mass= "<<ionspecies[i].Mass_cal
        <<" | AME20_ME= "<<ionspecies[i].AME<<" +- "<<ionspecies[i].AME_err
        <<fixed<<setprecision(1)
        <<" | ME= "<< ionspecies[i].MassExcess_cal<<" +- "<<ionspecies[i].Mass_cal_err 
        <<" | dm= "<<ionspecies[i].deltaMass<<" +- "<<ionspecies[i].deltaMass_err;
        
        if(MASS_VER>=3)outfile<<fixed<<setprecision(1)
        <<" | ME_VE= "<< ionspecies[i].MassExcess_cal_VE <<" +- "<<ionspecies[i].Mass_cal_err_VE
        <<" | dm_VE= "<<ionspecies[i].deltaMass_VE<<" +- "<<ionspecies[i].deltaMass_err_VE;
        else if(MASS_VER==1)outfile<<fixed<<setprecision(1)
        <<" | no_ME_VE= "<< 0 <<" +- "<<0
        <<" | no_dm_VE= "<<0<<" +- "<<0;
        else if(MASS_VER==2)outfile<<fixed<<setprecision(1)
        <<" | ME_v2= "<< ionspecies[i].MassExcess_cal_v2 <<" +- "<<ionspecies[i].Mass_cal_err_v2
        <<" | dm_v2= "<<ionspecies[i].deltaMass_v2<<" +- "<<ionspecies[i].deltaMass_err_v2;

        else outfile
        <<" | ME_v?= "<< 0 <<" +- "<<0
        <<" | dm_v?= "<<0<<" +- "<<0;
        
        outfile<<" | ions: "<<ionspecies[i].N_unknown<<endl;
        }
        */
     
    }
}
//------
if(Do_mass_err_scatter_ON&&MASS_VER>=3){outfile_sca_ave.close();}

//=============================== 生成 h_dm =============================
//h_dm 填入的是 v1 mass, 不是v2 也要有直方图， v2进行拟合
//为了 定义h_dm 首先需要知道如何分 bin, 以下采用 F-D 方法， 需要先算出 IQR

//粗略划分 14O和 21Mg 为了 计算F-D bin width
vector<double> vector_dm_v1_14Otmp; // @EXP@ THIS_EXP=="2021_36Ar_SET3"
vector<double> vector_dm_v1_21Mgtmp; 
double SET3_14O21Mg_dm_DIVIDE = -400;//##ARTIFICIAL
if(THIS_EXP=="2021_36Ar_SET3")
{
    //int N_14O_tmp=0;
    for(auto &j:ionspecies[ZN_ID[8][6]].vector_dm_v1)
    {
        if(j>SET3_14O21Mg_dm_DIVIDE){vector_dm_v1_14Otmp.push_back(j);}
        else{vector_dm_v1_21Mgtmp.push_back(j);}
    }
}

if(!LOOP_ON)
{
    // 输出文件
    outfile_VER2_FIT_info.open(FILEPATH+"VER2_FIT_info.txt");
    outfile_VER2_FIT_info<<" ---- IQR -----"<<endl;
    for(int i=0;i<NSpecies  ;i++)
    {  
        // 14O 未区分的情况下是双峰，如果用 vector_dm_v1 会按照双峰去算IQR 偏大
        if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O")
        {outfile_VER2_FIT_info<<ionspecies[i].Aname<<" "<<vector_dm_v1_14Otmp.size()<<" "<<Get_IQR(vector_dm_v1_14Otmp)<<endl;}
        else if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="21Mg")
        {outfile_VER2_FIT_info<<ionspecies[i].Aname<<" "<<vector_dm_v1_21Mgtmp.size()<<" "<<Get_IQR(vector_dm_v1_21Mgtmp)<<endl;} //
        
        else
        {outfile_VER2_FIT_info<<ionspecies[i].Aname<<" "<<ionspecies[i].N_unknown<<" "<<Get_IQR(ionspecies[i].vector_dm_v1)<<endl;}
    }      
}
if(!LOOP_ON)outfile_VER2_FIT_info<<" ---- F_D_binIwidth -----"<<endl;
//------------------ use Freedmann-Diaconis rule to get bin width:
for(int i = 0; i<NSpecies; i++)
{
    ionspecies[i].F_D_binIwidth = 2*Get_IQR(ionspecies[i].vector_dm_v1)/ pow( double(ionspecies[i].N_unknown) , 1.0/3.0) ;
    // 如果.vector_dm_v1 里面什么都没有 则 Get_IQR 直接返回 0 
    if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O")       //14O 单独处理 @EXP@
    {   
        ionspecies[i].F_D_binIwidth = 2*Get_IQR(vector_dm_v1_14Otmp)/ pow( double(vector_dm_v1_14Otmp.size()) , 1.0/3.0) ;
    }
    if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="21Mg") //21Mg
    {   
        ionspecies[i].F_D_binIwidth = 2*Get_IQR(vector_dm_v1_21Mgtmp)/ pow( double(vector_dm_v1_21Mgtmp.size()) , 1.0/3.0) ;
    }
    //cout<<ionspecies[i].Aname<<" h_dm: F_D_binIwidth = " <<ionspecies[i].F_D_binIwidth<<endl;
    if(!LOOP_ON)outfile_VER2_FIT_info<<ionspecies[i].Aname<<" "<<ionspecies[i].N_unknown<<" h_dm: F_D_binIwidth = " <<ionspecies[i].F_D_binIwidth<<endl;
    
    double FD_tmp = (ionspecies[i].F_D_binIwidth>0)?ionspecies[i].F_D_binIwidth:10;
    double mvq_FD_tmp = FD_tmp/(ionspecies[i].Z*u);
    //保证 FD_tmp 不是0

    //------------------ 定义 h_dm
    
    ionspecies[i].h_dm = new TH1F(ionspecies[i].Aname+" h_dm",ionspecies[i].Aname+" h_dm ", int(2000/FD_tmp),-1000,1000);   
    AxisFormat(ionspecies[i].h_dm,ionspecies[i].Aname+" h_dm ","#Deltam=M_{exp}-M_{AME}", "Counts");

    // Fill h_dm
    for(auto &j:ionspecies[i].vector_dm_v1){ionspecies[i].h_dm->Fill(j);}
         
}

//cout<<" return at 7210"<<endl;return;


//===============================  MASS VER 2  ============================

//---------- 2021 36Ar SET-3 14O/21Mg h_mvq_FD_14O21Mg fit 双高斯拟合
// 注意！！ 为了避免SET3的臃肿处理，可以输出两个峰的所有数据，用另外一个 Do_2Gaus 程序单独在外部做双高斯拟合
if(MASS_VER==2&&THIS_EXP == "2021_36Ar_SET3")
{
    int i_14O = ZN_ID[8][6];
    int i_21Mg = ZN_ID[12][9];
    // 定义 h_mvq_FD_14O21Mg 只用于需要双高斯拟合的 14O /21Mg @EXP@
    
    double FD_tmp_14O = (ionspecies[i_14O].F_D_binIwidth>0)?ionspecies[i_14O].F_D_binIwidth:10;
    double FD_tmp_21Mg = (ionspecies[i_21Mg].F_D_binIwidth>0)?ionspecies[i_21Mg].F_D_binIwidth:10;
    double mvq_FD_tmp = ( FD_tmp_14O/(8*u)+FD_tmp_21Mg/(12*u) )/2.0;
    cout<<" 2 gaus fit of 14O 21Mg: h_mvq_FD_14O21Mg bin width = "<<fixed<<setprecision(8)<<mvq_FD_tmp<<endl;
    
    TH1F* h_mvq_FD_14O21Mg = new TH1F(" h_mvq_FD_14O21Mg"," h_mvq_FD_14O21Mg", 
            int(2500/(8*u) /mvq_FD_tmp),ionspecies[i_14O].Mvq_AME-1500/(8*u),ionspecies[i_14O].Mvq_AME+1000/(8*u));
    AxisFormat(h_mvq_FD_14O21Mg, "h_mvq_FD_14O21Mg"," #frac{m}{q} [u/e]", "Counts");
    // Fill h_mvq_FD
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].A==14&&ions[j].Z==8){if(ions[j].mvq_v1>0){h_mvq_FD_14O21Mg->Fill(ions[j].mvq_v1);}}
    }
    

    double fit_mu_14O=0;double fit_mu_14O_err=0;
    double fit_sigma_14O=0;double fit_sigma_14O_err=0;
    double fit_mu_21Mg=0;double fit_mu_21Mg_err=0;
    double fit_sigma_21Mg=0;double fit_sigma_21Mg_err=0;
    double dm_fit_14O=0; double m_fit_14O_err=0;
    double dm_fit_21Mg=0; double m_fit_21Mg_err=0;

    //----------------------------- 2 GAUS FIT --------------------------------------------------    
    TF1* fitfun_2gaus_mvq_14O21Mg = new TF1(ionspecies[i_14O].Aname+"_fit_2gaus","gaus(0)+gaus(3)",1.74,1.76);
    fitfun_2gaus_mvq_14O21Mg->SetParameters(400, ionspecies[i_14O].Mvq_AME, 100/8/u,   200, ionspecies[i_21Mg].Mvq_AME, 100/12/u);//A1, mu1, sigma1, A2, mu2, sigma2,
    h_mvq_FD_14O21Mg->Fit(fitfun_2gaus_mvq_14O21Mg);
    
    fit_mu_14O = fitfun_2gaus_mvq_14O21Mg->GetParameter(1);
    fit_sigma_14O = fitfun_2gaus_mvq_14O21Mg->GetParameter(2);
    fit_mu_21Mg = fitfun_2gaus_mvq_14O21Mg->GetParameter(4);
    fit_sigma_21Mg = fitfun_2gaus_mvq_14O21Mg->GetParameter(5);

    fit_mu_14O_err = fitfun_2gaus_mvq_14O21Mg->GetParErrors()[1];
    fit_sigma_14O_err = fitfun_2gaus_mvq_14O21Mg->GetParErrors()[2];
    fit_mu_21Mg_err = fitfun_2gaus_mvq_14O21Mg->GetParErrors()[4];
    fit_sigma_21Mg_err = fitfun_2gaus_mvq_14O21Mg->GetParErrors()[5];

    cout<<" ----- h_mvq_FD_14O21Mg fit paras: ----- "<<endl;
    cout<<fixed<<setprecision(8)
    <<" 14O: mu= "<<fit_mu_14O<<" +- "<<fit_mu_14O_err<<" sigma= "<<fit_sigma_14O<<" +- "<<fit_sigma_14O_err<<endl
    <<" 21Mg: mu= "<<fit_mu_21Mg<<" +- "<<fit_mu_21Mg_err<<" sigma= "<<fit_sigma_21Mg<<" +- "<<fit_sigma_21Mg_err<<endl;

    if(!LOOP_ON) 
    {
        outfile_VER2_FIT_info<<" ----- h_mvq_FD_14O21Mg fit paras: ----- "<<endl;
        outfile_VER2_FIT_info<<" 14O: mu= "<<fit_mu_14O<<" +- "<<fit_mu_14O_err<<" sigma= "<<fit_sigma_14O<<" +- "<<fit_sigma_14O_err<<endl;
        outfile_VER2_FIT_info<<" 21Mg: mu= "<<fit_mu_21Mg<<" +- "<<fit_mu_21Mg_err<<" sigma= "<<fit_sigma_21Mg<<" +- "<<fit_sigma_21Mg_err<<endl;
    }

    dm_fit_14O = fit_mu_14O*8*u - ionspecies[i_14O].Mass;
    dm_fit_21Mg = fit_mu_21Mg*12*u - ionspecies[i_21Mg].Mass;
    m_fit_14O_err = fit_mu_14O_err*8*u;
    m_fit_21Mg_err = fit_mu_21Mg_err*12*u;

    if(!LOOP_ON) 
    {
    TCanvas * c_h_14O21Mg_2gaus = new TCanvas("c_h_14O21Mg_2gaus","c_h_14O21Mg_2gaus", 1000,500);
    h_mvq_FD_14O21Mg->DrawClone();
    }
    cout<<" ----- h_mvq_FD_14O21Mg fit results: "<<endl;
    cout<<" dm(14O) = "<<fixed<<setprecision(1)<< dm_fit_14O<<" +- "<<m_fit_14O_err<<" dm(21Mg) = "<<dm_fit_21Mg<<" +- "<<m_fit_21Mg_err<<endl;
    if(!LOOP_ON)
    {outfile_VER2_FIT_info<<endl<<" ----- h_mvq_FD_14O21Mg fit results: dm(14O) = "<<dm_fit_14O<<" +- "<<m_fit_14O_err<<" dm(21Mg) = "<<dm_fit_21Mg<<" +- "<<m_fit_21Mg_err<<endl;}

    ionspecies[i_14O ].deltaMass_v2 = dm_fit_14O;
    ionspecies[i_14O ].Mass_cal_err_v2 = m_fit_14O_err;
    ionspecies[i_14O ].deltaMass_err_v2 = sqrt( pow(ionspecies[i_14O].Mass_cal_err_v2,2)+pow(ionspecies[i_14O].AME_err,2)  );

    ionspecies[i_21Mg].deltaMass_v2 = dm_fit_21Mg;
    ionspecies[i_21Mg].Mass_cal_err_v2 = m_fit_21Mg_err;
    ionspecies[i_21Mg].deltaMass_err_v2 = sqrt( pow(ionspecies[i_21Mg].Mass_cal_err_v2,2)+pow(ionspecies[i_21Mg].AME_err,2)  );
    
    delete h_mvq_FD_14O21Mg;
    delete fitfun_2gaus_mvq_14O21Mg;
    //======================================= 2024 Re_IonIdentify  SET3==========================================
    vector<double> vector_re_id_14O_mvq;
    vector<double> vector_re_id_21Mg_mvq;
    if (Re_IonIdentify_on&&THIS_EXP=="2021_36Ar_SET3")   //针对 SET3 14O 和 21Mg 分开后重新制作新的input 文件，下次处理时，14O 和 21Mg 已鉴别
    {
        ionspecies[i_14O].ClearRecord();
        ionspecies[i_21Mg].ClearRecord();
        outfile.open("OUTPUT//ReIdentify_20230824_174150-del-2inj.txt"); //重新识别结果
        for (int i = 0; i < ions_n; i++)//在这里对没分开的核进行二次识别
        {
            if (ionspecies[ions[i].Species].Aname == "14O")
            {//14O 和 21Mg   在 m/q 双峰中心值3倍sigma以内进行划分
                //把原来认为是 14O 但是落在21Mg峰3sigma 范围内的， 重新识别为 21Mg。 没有得到质量的，以及中间区域的，舍弃
                if (abs(ions[i].mvq_v1 - fit_mu_14O)<3*fit_sigma_14O)
                {
                    ions[i].Species = ionspecies[i_14O].Species;
                    ions[i].A = 14;
                    ions[i].Z = 8;
                    ions[i].name = "O";
                    vector_re_id_14O_mvq.push_back(ions[i].mvq_v1);
                    ionspecies[i_14O].Record(ions[i]);   //对 T C v Bρ γt 总粒子数进行统计
                }
                else if (abs(ions[i].mvq_v1 - fit_mu_21Mg)<3*fit_sigma_21Mg)
                {
                    ions[i].Species = ionspecies[i_21Mg].Species;
                    ions[i].A = 21;
                    ions[i].Z = 12;
                    ions[i].name = "Mg";
                    vector_re_id_21Mg_mvq.push_back(ions[i].mvq_v1);
                    ionspecies[i_21Mg].Record(ions[i]);   //对 T C v Bρ γt 总粒子数进行统计
                }
                else
                {
                    continue; //不输出这个落在双峰范围外的
                }
            }
            //生成一个重新鉴别后的输入文件      ReIdentify_20230824_174150-del-2inj.txt
            outfile << fixed << setiosflags(ios::left)
                << setw(4) << ions[i].name << " "
                << setw(4) << ions[i].A << " "
                << setw(4) << ions[i].Z << " "
                << setw(13) << setprecision(9) << ions[i].T << " "
                << setw(13) << setprecision(9) << ions[i].T_err << " "
                << setw(6) << ions[i].ion_number << " "
                << setw(6) << ions[i].inject_number << " "
                << setw(18) << setprecision(9) << ions[i].A1 << " "
                << setw(18) << setprecision(9) << ions[i].A1err << " "
                << setw(18) << setprecision(12) << ions[i].A2 << " "
                << setw(18) << setprecision(12) << ions[i].A2err << " "
                << setw(18) << setprecision(12) << ions[i].A3 << " "
                << setw(18) << setprecision(12) << ions[i].A4 << " "
                << setw(18) << setprecision(9) << ions[i].dA0 / 1000.0 << " "
                << setw(18) << setprecision(9) << ions[i].dA0err << " "
                << setw(18) << setprecision(9) << ions[i].dA1 << " "
                << setw(18) << setprecision(9) << ions[i].dA1err << " "
                << setw(18) << setprecision(9) << ions[i].cov12 << " "
                << setw(18) << setprecision(9) << ions[i].cov15 << " "
                << setw(18) << setprecision(9) << ions[i].cov16 << " "
                << setw(18) << setprecision(9) << ions[i].cov25 << " "
                << setw(18) << setprecision(9) << ions[i].cov26 << " "
                << setw(18) << setprecision(9) << ions[i].cov56 << " "
                << setw(30) << ions[i].inject_filename << endl;
        }
        outfile.close();
    

        //------- 显示重新鉴别后的双峰 直方图 ------------
    
        double FD_tmp_re_id_14O21Mg = 2*Get_IQR(vector_re_id_14O_mvq)/ pow( double(vector_re_id_14O_mvq.size()) , 1.0/3.0) ;
        TH1F* h_mvq_re_id_14O21Mg = new TH1F("h_mvq_re_id_14O21Mg","h_mvq_re_id_14O21Mg",  
            int(2500/(8*u)/FD_tmp_re_id_14O21Mg), ionspecies[i_14O].Mvq_AME-1500/(8*u), ionspecies[i_14O].Mvq_AME+1000/(8*u)  );
        AxisFormat( h_mvq_re_id_14O21Mg,"h_mvq_re_id_14O21Mg"," #frac{m}{q} [u/e]", "Counts");
        TH1F* h_mvq_re_id_14O = new TH1F("h_mvq_re_id_14O","h_mvq_re_id_14O",  
            int(2000/(8*u)/FD_tmp_re_id_14O21Mg), ionspecies[i_14O].Mvq_AME-1000/(8*u), ionspecies[i_14O].Mvq_AME+1000/(8*u)  );
        AxisFormat( h_mvq_re_id_14O,"h_mvq_re_id_14O"," #frac{m}{q} [u/e]", "Counts");
        TH1F* h_mvq_re_id_21Mg = new TH1F("h_mvq_re_id_21Mg","h_mvq_re_id_21Mg",  
            int(2000/(12*u)/FD_tmp_re_id_14O21Mg), ionspecies[i_21Mg].Mvq_AME-1000/(12*u), ionspecies[i_21Mg].Mvq_AME+1000/(12*u)  );
        AxisFormat( h_mvq_re_id_21Mg,"h_mvq_re_id_21Mg"," #frac{m}{q} [u/e]", "Counts");

        ionspecies[i_14O].h_dm->Reset();
        ionspecies[i_21Mg].h_dm->Reset();
        for(auto&j :vector_re_id_14O_mvq )
        {
            h_mvq_re_id_14O21Mg->Fill(j);
            ionspecies[i_14O].h_dm->Fill(j*8*u-ionspecies[i_14O].Mass);
        }
        for(auto&j :vector_re_id_21Mg_mvq )
        {
            h_mvq_re_id_14O21Mg->Fill(j);
            ionspecies[i_21Mg].h_dm->Fill(j*12*u-ionspecies[i_21Mg].Mass);
        }

        if(!LOOP_ON)
        {
            TCanvas * c_h_re_id_14O21Mg = new TCanvas("c_h_re_id_14O21Mg","c_h_re_id_14O21Mg", 1000,500);
            h_mvq_re_id_14O21Mg->DrawClone();
            delete h_mvq_re_id_14O21Mg;
            TCanvas * c_h_re_id_14O = new TCanvas("c_h_re_id_14O","c_h_re_id_14O", 1000,500);
            //h_mvq_re_id_14O->Draw();
            ionspecies[i_14O].h_dm->Draw();
            TCanvas * c_h_re_id_21Mg = new TCanvas("c_h_re_id_21Mg","c_h_re_id_21Mg", 1000,500);
            //h_mvq_re_id_21Mg->Draw();
            ionspecies[i_21Mg].h_dm->Draw();
        }
        if(LOOP_ON)
        {
            delete h_mvq_re_id_14O21Mg;
            delete h_mvq_re_id_14O;
            delete h_mvq_re_id_21Mg;
        }


        //----------- 重新统计 14O/21Mg --------------

        ionspecies[i_14O].CalculateAve(); //determine average T C v Bp gammat 计算平均值
        ionspecies[i_21Mg].CalculateAve();
        
        for(int i=0;i<ions_n ;i++)
        {
            if(ions[i].A==14&&ions[i].Z==8) {ionspecies[i_14O].SigmaPreAdd(ions[i]);ionspecies[i_14O].SkewnessPreAdd(ions[i]);}
            if(ions[i].A==21&&ions[i].Z==12){ionspecies[i_21Mg].SigmaPreAdd(ions[i]);ionspecies[i_21Mg].SkewnessPreAdd(ions[i]);}    
        }
    
        ionspecies[i_14O].CalculateSigma();ionspecies[i_14O].CalculateSkewness();
        ionspecies[i_21Mg].CalculateSigma();ionspecies[i_21Mg].CalculateSkewness();
    
        ionspecies[i_14O].N_unknown = vector_re_id_14O_mvq.size();
        ionspecies[i_14O].Mass_cal = GetVectorMean(vector_re_id_14O_mvq)*8*u;
        ionspecies[i_14O].stdDeviation_V1 = GetVectorStdDev(vector_re_id_14O_mvq,GetVectorMean(vector_re_id_14O_mvq))*8*u;
        ionspecies[i_14O].Mass_cal_err = ionspecies[i_14O].stdDeviation_V1/(sqrt(ionspecies[i_14O].N_unknown));
        ionspecies[i_14O].deltaMass = ionspecies[i_14O].Mass_cal - ionspecies[i_14O].Mass;
        ionspecies[i_14O].deltaMass_err = sqrt( pow(ionspecies[i_14O].Mass_cal_err,2)+pow(ionspecies[i_14O].AME_err,2));
        
        ionspecies[i_21Mg].N_unknown = vector_re_id_21Mg_mvq.size();
        ionspecies[i_21Mg].Mass_cal = GetVectorMean(vector_re_id_21Mg_mvq)*12*u;
        ionspecies[i_21Mg].stdDeviation_V1 = GetVectorStdDev(vector_re_id_21Mg_mvq,GetVectorMean(vector_re_id_21Mg_mvq))*12*u;
        ionspecies[i_21Mg].Mass_cal_err = ionspecies[i_21Mg].stdDeviation_V1/(sqrt(ionspecies[i_21Mg].N_unknown));
        ionspecies[i_21Mg].deltaMass = ionspecies[i_21Mg].Mass_cal - ionspecies[i_21Mg].Mass;
        ionspecies[i_21Mg].deltaMass_err = sqrt( pow(ionspecies[i_14O].Mass_cal_err,2)+pow(ionspecies[i_14O].AME_err,2));

        ionspecies[i_21Mg].HasResult = true;


        cout<<fixed<<setprecision(4)<<" ---- after reidentification ----"
        <<" 14O N_unknown = "<< ionspecies[i_14O].N_unknown<<" AveT= "<<ionspecies[i_14O].AveT
        <<" 21Mg N_unknown = "<< ionspecies[i_21Mg].N_unknown<<" AveT= "<<ionspecies[i_21Mg].AveT<<fixed<<setprecision(1)<<endl;
    }
    //_______________________________________ 2024 Re_IonIdentify ___________________________________________
}




//========== mass ver 2 for 58Ni, deal with 3 nuclides that has isomer: 24Al , 44V, 52Co
// 注意！！ 为了避免臃肿处理，可以输出两个峰的所有数据，用另外一个 Do_2Gaus 程序单独在外部做双高斯拟合
if(MASS_VER==2&&THIS_EXP == "2017_58Ni")
{
    int i_24Al = ZN_ID[13][11];
    int i_44V = ZN_ID[23][21];
    int i_52Co = ZN_ID[27][25];
    int i_24Al_m=0;int i_44V_m=0;int i_52Co_m=0;
    for (int i = 0; i <NSpecies_Isomer; i++)
    {
       if(Isomer_ionspecis[i].Aname=="24Al_m"){i_24Al_m = i;}
       if(Isomer_ionspecis[i].Aname=="44V_m"){i_44V_m = i;}
       if(Isomer_ionspecis[i].Aname=="52Co_m"){i_52Co_m = i;}
    }
    cout<<"  debug ---- i_24Al_m = "<<i_24Al_m<<" i_44V_m= "<<i_44V_m<<" i_52Co_m = "<<i_52Co_m<<endl;
    //根据 M.Zhang 论文 使用同样的分BIN:
    // mvq: [u/e] 24Al 3.3e-6, 44V 1.4e-6, 52Co 1.4e-6 
    double h_mvq_24Al_isomer_MIN = 1.84556;
    double h_mvq_24Al_isomer_MAX = 1.84568;
    double h_mvq_24Al_isomer_bin = 3.0*0.000001;
    double h_mvq_44V_isomer_MIN = 1.91136;
    double h_mvq_44V_isomer_MAX = 1.91142;
    double h_mvq_44V_isomer_bin = 1.4*0.000001;
    double h_mvq_52Co_isomer_MIN = 1.92398;
    double h_mvq_52Co_isomer_MAX = 1.92406;
    double h_mvq_52Co_isomer_bin = 1.4*0.000001;
    
    TH1F* h_mvq_24Al_isomer = new TH1F(" h_mvq_24Al_isomer"," h_mvq_24Al_isomer", 
            int( (h_mvq_24Al_isomer_MAX -h_mvq_24Al_isomer_MIN)/h_mvq_24Al_isomer_bin ), h_mvq_24Al_isomer_MIN, h_mvq_24Al_isomer_MAX);
    AxisFormat(h_mvq_24Al_isomer, "h_mvq_24Al_isomer"," #frac{m}{q} [u/e]", "Counts");
    
    TH1F* h_mvq_44V_isomer = new TH1F(" h_mvq_44V_isomer"," h_mvq_44V_isomer", 
            int( (h_mvq_44V_isomer_MAX -h_mvq_44V_isomer_MIN)/h_mvq_44V_isomer_bin ), h_mvq_44V_isomer_MIN, h_mvq_44V_isomer_MAX);
    AxisFormat(h_mvq_44V_isomer, "h_mvq_44V_isomer"," #frac{m}{q} [u/e]", "Counts");
    
    TH1F* h_mvq_52Co_isomer = new TH1F(" h_mvq_52Co_isomer"," h_mvq_52Co_isomer", 
            int( (h_mvq_52Co_isomer_MAX -h_mvq_52Co_isomer_MIN)/h_mvq_52Co_isomer_bin ), h_mvq_52Co_isomer_MIN, h_mvq_52Co_isomer_MAX);
    AxisFormat(h_mvq_52Co_isomer, "h_mvq_52Co_isomer"," #frac{m}{q} [u/e]", "Counts");
    // Fill h_mvq
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].A==24&&ions[j].Z==13){if(ions[j].mvq_v1>0){h_mvq_24Al_isomer->Fill(ions[j].mvq_v1);}}
        if(ions[j].A==44&&ions[j].Z==23){if(ions[j].mvq_v1>0){h_mvq_44V_isomer->Fill(ions[j].mvq_v1);}}
        if(ions[j].A==52&&ions[j].Z==27){if(ions[j].mvq_v1>0){h_mvq_52Co_isomer->Fill(ions[j].mvq_v1);}}
    }
    
    double fit_mu_24Al=0;double fit_mu_24Al_err=0;
    double fit_sigma_24Al=0;double fit_sigma_24Al_err=0;
    double fit_mu_24Al_m=0;double fit_mu_24Al_m_err=0;
    double fit_sigma_24Al_m=0;double fit_sigma_24Al_m_err=0;
    double dm_fit_24Al=0; double m_fit_24Al_err=0;
    double dm_fit_24Al_m=0; double m_fit_24Al_m_err=0;

    double fit_mu_44V=0;double fit_mu_44V_err=0;
    double fit_sigma_44V=0;double fit_sigma_44V_err=0;
    double fit_mu_44V_m=0;double fit_mu_44V_m_err=0;
    double fit_sigma_44V_m=0;double fit_sigma_44V_m_err=0;
    double dm_fit_44V=0; double m_fit_44V_err=0;
    double dm_fit_44V_m=0; double m_fit_44V_m_err=0;

    double fit_mu_52Co=0;double fit_mu_52Co_err=0;
    double fit_sigma_52Co=0;double fit_sigma_52Co_err=0;
    double fit_mu_52Co_m=0;double fit_mu_52Co_m_err=0;
    double fit_sigma_52Co_m=0;double fit_sigma_52Co_m_err=0;
    double dm_fit_52Co=0; double m_fit_52Co_err=0;
    double dm_fit_52Co_m=0; double m_fit_52Co_m_err=0;

    Do_2gaus(h_mvq_24Al_isomer, 1.8, 1.9, 60, ionspecies[i_24Al].Mvq_AME, 10.0/13/u,   40, Isomer_ionspecis[i_24Al_m].Mvq_AME, 10.0/13/u,
    fit_mu_24Al,fit_mu_24Al_err,fit_mu_24Al_m,fit_mu_24Al_m_err,
    fit_sigma_24Al,fit_sigma_24Al_err,fit_sigma_24Al_m,fit_sigma_24Al_m_err,
    !LOOP_ON);

    Do_2gaus(h_mvq_44V_isomer, 1.8, 1.9, 40, ionspecies[i_44V].Mvq_AME, 10.0/23/u,   40, Isomer_ionspecis[i_44V_m].Mvq_AME, 10.0/23/u,
    fit_mu_44V,fit_mu_44V_err,fit_mu_44V_m,fit_mu_44V_m_err,
    fit_sigma_44V,fit_sigma_44V_err,fit_sigma_44V_m,fit_sigma_44V_m_err,
    !LOOP_ON);

    Do_2gaus(h_mvq_52Co_isomer, 1.8, 1.9, 80, ionspecies[i_52Co].Mvq_AME, 10.0/27/u,   40, Isomer_ionspecis[i_52Co_m].Mvq_AME, 10.0/27/u,
    fit_mu_52Co,fit_mu_52Co_err,fit_mu_52Co_m,fit_mu_52Co_m_err,
    fit_sigma_52Co,fit_sigma_52Co_err,fit_sigma_52Co_m,fit_sigma_52Co_m_err,
    !LOOP_ON);


    if(!LOOP_ON) 
    {
        outfile_VER2_FIT_info<<" ----- h_mvq_24Al_isomer fit paras: ----- "<<endl;
        outfile_VER2_FIT_info<<" 24Al: "<<" mu= "<<fit_mu_24Al<<" +- "<<fit_mu_24Al_err<<" sigma= "<<fit_sigma_24Al<<" +- "<<fit_sigma_24Al_err<<endl;
        outfile_VER2_FIT_info<<" 24Al_m: "<<" mu= "<<fit_mu_24Al_m<<" +- "<<fit_mu_24Al_m_err<<" sigma= "<<fit_sigma_24Al_m<<" +- "<<fit_sigma_24Al_m_err<<endl;
        outfile_VER2_FIT_info<<" ----- h_mvq_44V_isomer fit paras: ----- "<<endl;
        outfile_VER2_FIT_info<<" 44V: "<<" mu= "<<fit_mu_44V<<" +- "<<fit_mu_44V_err<<" sigma= "<<fit_sigma_44V<<" +- "<<fit_sigma_44V_err<<endl;
        outfile_VER2_FIT_info<<" 44V_m: "<<" mu= "<<fit_mu_44V_m<<" +- "<<fit_mu_44V_m_err<<" sigma= "<<fit_sigma_44V_m<<" +- "<<fit_sigma_44V_m_err<<endl;
        outfile_VER2_FIT_info<<" ----- h_mvq_52Co_isomer fit paras: ----- "<<endl;
        outfile_VER2_FIT_info<<" 52Co: "<<" mu= "<<fit_mu_52Co<<" +- "<<fit_mu_52Co_err<<" sigma= "<<fit_sigma_52Co<<" +- "<<fit_sigma_52Co_err<<endl;
        outfile_VER2_FIT_info<<" 52Co_m: "<<" mu= "<<fit_mu_52Co_m<<" +- "<<fit_mu_52Co_m_err<<" sigma= "<<fit_sigma_52Co_m<<" +- "<<fit_sigma_52Co_m_err<<endl;
    
    }

    dm_fit_24Al = fit_mu_24Al*ionspecies[i_24Al].Z*u - ionspecies[i_24Al].Mass;
    dm_fit_24Al_m = fit_mu_24Al_m*Isomer_ionspecis[i_24Al_m].Z*u - Isomer_ionspecis[i_24Al_m].Mass;
    m_fit_24Al_err = fit_mu_24Al_err*ionspecies[i_24Al].Z*u;
    m_fit_24Al_m_err = fit_mu_24Al_m_err*Isomer_ionspecis[i_24Al_m].Z*u;

    dm_fit_44V = fit_mu_44V*ionspecies[i_44V].Z*u - ionspecies[i_44V].Mass;
    dm_fit_44V_m = fit_mu_44V_m*Isomer_ionspecis[i_44V_m].Z*u - Isomer_ionspecis[i_44V_m].Mass;
    m_fit_44V_err = fit_mu_44V_err*ionspecies[i_44V].Z*u;
    m_fit_44V_m_err = fit_mu_44V_m_err*Isomer_ionspecis[i_44V_m].Z*u;
    
    dm_fit_52Co = fit_mu_52Co*ionspecies[i_52Co].Z*u - ionspecies[i_52Co].Mass;
    dm_fit_52Co_m = fit_mu_52Co_m*Isomer_ionspecis[i_52Co_m].Z*u - Isomer_ionspecis[i_52Co_m].Mass;
    m_fit_52Co_err = fit_mu_52Co_err*ionspecies[i_52Co].Z*u;
    m_fit_52Co_m_err = fit_mu_52Co_m_err*Isomer_ionspecis[i_52Co_m].Z*u;
    

    cout<<" ----- h_mvq_24Al_isomer fit results: "<<endl;
    cout<<" dm(24Al) = "<<fixed<<setprecision(1)<< dm_fit_24Al<<" +- "<<m_fit_24Al_err<<" dm(24Al_m) = "<<dm_fit_24Al_m<<" +- "<<m_fit_24Al_m_err<<endl;
    
    cout<<" ----- h_mvq_44V_isomer fit results: "<<endl;
    cout<<" dm(44V) = "<<fixed<<setprecision(1)<< dm_fit_44V<<" +- "<<m_fit_44V_err<<" dm(44V_m) = "<<dm_fit_44V_m<<" +- "<<m_fit_44V_m_err<<endl;
    
    cout<<" ----- h_mvq_52Co_isomer fit results: "<<endl;
    cout<<" dm(52Co) = "<<fixed<<setprecision(1)<< dm_fit_52Co<<" +- "<<m_fit_52Co_err<<" dm(52Co_m) = "<<dm_fit_52Co_m<<" +- "<<m_fit_52Co_m_err<<endl;
    
    if(!LOOP_ON)
    {
        outfile_VER2_FIT_info<<endl<<" ----- h_mvq_24Al_isomer fit results: dm(24Al) = "<<dm_fit_24Al<<" +- "<<m_fit_24Al_err<<" dm(24Al_m) = "<<dm_fit_24Al_m<<" +- "<<m_fit_24Al_m_err<<endl;
        outfile_VER2_FIT_info<<endl<<" ----- h_mvq_44V_isomer fit results: dm(44V) = "<<dm_fit_44V<<" +- "<<m_fit_44V_err<<" dm(44V_m) = "<<dm_fit_44V_m<<" +- "<<m_fit_44V_m_err<<endl;
        outfile_VER2_FIT_info<<endl<<" ----- h_mvq_52Co_isomer fit results: dm(52Co) = "<<dm_fit_52Co<<" +- "<<m_fit_52Co_err<<" dm(52Co_m) = "<<dm_fit_52Co_m<<" +- "<<m_fit_52Co_m_err<<endl;
    }


    /*
    ionspecies[i_24Al].deltaMass_v2 = dm_fit_24Al;
    ionspecies[i_24Al].Mass_cal_v2  = fit_mu_24Al*ionspecies[i_24Al].Z*u;
    ionspecies[i_24Al].MassExcess_cal_v2 = ionspecies[i_24Al].GetMassExcess_cal(ionspecies[i_24Al].Mass_cal_v2);
    ionspecies[i_24Al].Mass_cal_err_v2 = m_fit_24Al_err;
    ionspecies[i_24Al].deltaMass_err_v2 = sqrt( pow(ionspecies[i_24Al].Mass_cal_err_v2,2)+pow(ionspecies[i_24Al].AME_err,2)  );
    */
    
    /*
    Isomer_ionspecis[i_24Al_m].deltaMass_v2 = dm_fit_24Al_m;
    Isomer_ionspecis[i_24Al_m].Mass_cal_err_v2 = m_fit_24Al_m_err;
    Isomer_ionspecis[i_24Al_m].deltaMass_err_v2 = sqrt( pow(Isomer_ionspecis[i_24Al_m].Mass_cal_err_v2,2)+pow(Isomer_ionspecis[i_24Al_m].AME_err,2)  );
    
    ionspecies[i_44V].deltaMass_v2 = dm_fit_44V;
    ionspecies[i_44V].Mass_cal_err_v2 = m_fit_44V_err;
    ionspecies[i_44V].deltaMass_err_v2 = sqrt( pow(ionspecies[i_44V].Mass_cal_err_v2,2)+pow(ionspecies[i_44V].AME_err,2)  );
    Isomer_ionspecis[i_44V_m].deltaMass_v2 = dm_fit_44V_m;
    Isomer_ionspecis[i_44V_m].Mass_cal_err_v2 = m_fit_44V_m_err;
    Isomer_ionspecis[i_44V_m].deltaMass_err_v2 = sqrt( pow(Isomer_ionspecis[i_44V_m].Mass_cal_err_v2,2)+pow(Isomer_ionspecis[i_44V_m].AME_err,2)  );
    
    ionspecies[i_52Co].deltaMass_v2 = dm_fit_52Co;
    ionspecies[i_52Co].Mass_cal_err_v2 = m_fit_52Co_err;
    ionspecies[i_52Co].deltaMass_err_v2 = sqrt( pow(ionspecies[i_52Co].Mass_cal_err_v2,2)+pow(ionspecies[i_52Co].AME_err,2)  );
    Isomer_ionspecis[i_52Co_m].deltaMass_v2 = dm_fit_52Co_m;
    Isomer_ionspecis[i_52Co_m].Mass_cal_err_v2 = m_fit_52Co_m_err;
    Isomer_ionspecis[i_52Co_m].deltaMass_err_v2 = sqrt( pow(Isomer_ionspecis[i_52Co_m].Mass_cal_err_v2,2)+pow(Isomer_ionspecis[i_52Co_m].AME_err,2)  );
    */
    ionspecies[i_24Al].DoSetMass_v2(fit_mu_24Al,dm_fit_24Al,m_fit_24Al_err);
    Isomer_ionspecis[i_24Al_m].DoSetMass_v2(fit_mu_24Al_m,dm_fit_24Al_m,m_fit_24Al_m_err);
    ionspecies[i_44V].DoSetMass_v2(fit_mu_44V,dm_fit_44V,m_fit_44V_err);
    Isomer_ionspecis[i_44V_m].DoSetMass_v2(fit_mu_44V_m,dm_fit_44V_m,m_fit_44V_m_err);
    ionspecies[i_52Co].DoSetMass_v2(fit_mu_52Co,dm_fit_52Co,m_fit_52Co_err);
    Isomer_ionspecis[i_52Co_m].DoSetMass_v2(fit_mu_52Co_m,dm_fit_52Co_m,m_fit_52Co_m_err);
    
}





//=======================  mass ver 2 histogram fit ===============================
int h_m_Gauss_fit_min_count =400;   //##ARTIFICIAL

if(THIS_EXP=="2021_36Ar_SET3")h_m_Gauss_fit_min_count =150;
else if(THIS_EXP=="2017_58Ni")h_m_Gauss_fit_min_count =200;
else                          h_m_Gauss_fit_min_count =400;

if(MASS_VER==2)
{   
    if(mass_v2_show_each_fit_ON){cout<<" ================== mass_v2_show_each_fit_ON ==================="<<endl;}
    TCanvas* c_each_h_dm_fit[MAX_IONSPECIES];
    //只要数目足够 都先认为可以单gaus
    for(int i = 0; i<NSpecies; i++)
    {
        if(ionspecies[i].N_unknown>h_m_Gauss_fit_min_count){ionspecies[i].h_m_Gauss_fit_opt = 1;}
        //指定不能进行单高斯拟合的：
        if(ionspecies[i].Aname=="24Al")ionspecies[i].h_m_Gauss_fit_opt = 2;
        if(ionspecies[i].Aname=="44V" )ionspecies[i].h_m_Gauss_fit_opt = 2;
        if(ionspecies[i].Aname=="52Co")ionspecies[i].h_m_Gauss_fit_opt = 2;
    }
    


    for(int i = 0; i<NSpecies; i++)
    {
        if(mass_v2_show_each_fit_ON&&!LOOP_ON)c_each_h_dm_fit[i] = new TCanvas("c_h_dm_fit_"+ionspecies[i].Aname,"c_h_dm_fit_"+ionspecies[i].Aname, 1000,500);
        
        if(ionspecies[i].h_m_Gauss_fit_opt==1)
        {
            ionspecies[i].fitfun_gaus_dm = new TF1(ionspecies[i].Aname+"_fit_gaus","gaus(0)",-1000,1000);
            ionspecies[i].fitfun_gaus_dm->SetParameters(ionspecies[i].N_unknown, ionspecies[i].deltaMass, ionspecies[i].Mass_cal_err); // gaus(A,MIU,SIGMA)
            ionspecies[i].h_dm->Fit(ionspecies[i].fitfun_gaus_dm,"Q");
        }
        
        //###########################
        if(mass_v2_show_each_fit_ON&&!LOOP_ON)ionspecies[i].h_dm->DrawClone();
        //##############################
        //MASS VER 2 results initialization
        if(ionspecies[i].h_m_Gauss_fit_opt==0)
        {
            ionspecies[i].deltaMass_v2 = 0;
            ionspecies[i].Mass_cal_err_v2 = 0;
            ionspecies[i].deltaMass_err_v2 = sqrt( pow(ionspecies[i].Mass_cal_err_v2,2)+pow(ionspecies[i].AME_err,2)  );
        }
        if(ionspecies[i].h_m_Gauss_fit_opt==1)
        {
            ionspecies[i].deltaMass_v2 = ionspecies[i].fitfun_gaus_dm->GetParameter(1);
            ionspecies[i].Mass_cal_err_v2 = ionspecies[i].fitfun_gaus_dm->GetParErrors()[1];
            ionspecies[i].deltaMass_err_v2 = sqrt( pow(ionspecies[i].Mass_cal_err_v2,2)+pow(ionspecies[i].AME_err,2)  );
            ionspecies[i].fitfun_gaus_dm_A = ionspecies[i].fitfun_gaus_dm->GetParameter(0);    ionspecies[i].fitfun_gaus_dm_A_err = ionspecies[i].fitfun_gaus_dm->GetParErrors()[0];
            ionspecies[i].fitfun_gaus_dm_mu = ionspecies[i].fitfun_gaus_dm->GetParameter(1);   ionspecies[i].fitfun_gaus_dm_mu_err = ionspecies[i].fitfun_gaus_dm->GetParErrors()[1];
            ionspecies[i].fitfun_gaus_dm_sigma = ionspecies[i].fitfun_gaus_dm->GetParameter(2);ionspecies[i].fitfun_gaus_dm_sigma_err = ionspecies[i].fitfun_gaus_dm->GetParErrors()[2];
            ionspecies[i].Mass_cal_v2 = ionspecies[i].Mass + ionspecies[i].deltaMass_v2;
            ionspecies[i].MassExcess_cal_v2 = ionspecies[i].GetMassExcess_cal(ionspecies[i].Mass_cal_v2);
        }
        
    }
    if(LOOP_ON)
    {
        for(int i = 0; i<NSpecies; i++)
        {
            if(ionspecies[i].h_m_Gauss_fit_opt==1)
            {
                delete ionspecies[i].fitfun_gaus_dm;
            }
        }
    }

    cout<<endl<<"---------------------  MASS VER == 2 fit paras:"<<endl;
    if(!LOOP_ON)outfile_VER2_FIT_info<<"---------------------  MASS VER == 2 fit paras:"<<endl;
    for(int i = 0; i<NSpecies; i++)
    {
        ionspecies[i].Show_h_dm_fit_gaus_paras();
        if(!LOOP_ON)ionspecies[i].Show_h_dm_fit_gaus_paras(outfile_VER2_FIT_info);
    }

    cout<<endl<<"---------------------  MASS VER == 2 fit result:"<<endl;
    if(!LOOP_ON)outfile_VER2_FIT_info<<"---------------------  MASS VER == 2 fit result:"<<endl;
    for(int i = 0; i<NSpecies; i++)
    {
        cout<<ionspecies[i].Aname<<" "<<ionspecies[i].N_unknown<<" "<<ionspecies[i].h_m_Gauss_fit_opt<<" : "
        <<ionspecies[i].deltaMass_v2<<" +- "<<ionspecies[i].Mass_cal_err_v2<<endl;
        if(!LOOP_ON)outfile_VER2_FIT_info<<ionspecies[i].Aname<<" "<<ionspecies[i].N_unknown<<" "<<ionspecies[i].h_m_Gauss_fit_opt<<" : "
        <<ionspecies[i].deltaMass_v2<<" +- "<<ionspecies[i].Mass_cal_err_v2<<endl;
    }



    //对于 v2 中不能做拟合的，计数太少： opt=0 使用 v1 结果

    for(int i = 0; i<NSpecies; i++)
    {
        if(ionspecies[i].h_m_Gauss_fit_opt==0)
        {
            ionspecies[i].deltaMass_v2 = ionspecies[i].deltaMass;
            ionspecies[i].Mass_cal_err_v2 = ionspecies[i].Mass_cal_err;
            ionspecies[i].deltaMass_err_v2 = ionspecies[i].deltaMass_err ; 
            ionspecies[i].Mass_cal_v2 = ionspecies[i].Mass_cal;
            ionspecies[i].MassExcess_cal_v2 = ionspecies[i].MassExcess_cal;
        }
    }
    cout<<endl<<"---------------------  MASS VER == 2 final  results:"<<endl;
    if(!LOOP_ON)outfile_VER2_FIT_info<<"---------------------  MASS VER == 2 final  results:"<<endl;
    for(int i = 0; i<NSpecies; i++)
    {
        cout<<ionspecies[i].Aname<<" "<<ionspecies[i].N_unknown<<" "<<ionspecies[i].h_m_Gauss_fit_opt<<" : "
            <<ionspecies[i].deltaMass_v2<<" +- "<<ionspecies[i].Mass_cal_err_v2<<endl;
        if(!LOOP_ON)outfile_VER2_FIT_info<<ionspecies[i].Aname<<" "<<ionspecies[i].N_unknown<<" "<<ionspecies[i].h_m_Gauss_fit_opt<<" : "
            <<ionspecies[i].deltaMass_v2<<" +- "<<ionspecies[i].Mass_cal_err_v2<<endl;
    }

if(!LOOP_ON)outfile_VER2_FIT_info.close();

}








//以上做完了MASS VER2 再统一输出
//======================================= 输出质量结果文件 ===========================
if(OUTPUT_MASSRESULT)//输出质量结果
{
    filename2=FILEPATH_m_s+"m_s_v2024_"+scan_loop_i+".txt";
    outfile.open( filename2);
    outfile<<" MASS_VER= "<<MASS_VER<<endl;
    outfile<<fixed<<setprecision(4)<<"L= "<<L<<" ddT= "<<ddT<<endl;

    for(int i = 0; i<NSpecies; i++)
    {

    outfile<<ionspecies[i].A<<" "<<ionspecies[i].name<<" "<<ionspecies[i].Z<<" "<<" | "<<ionspecies[i].MassUnknown<<" | "
    <<fixed<<setprecision(3)<<ionspecies[i].AveT
    <<" | nuc_mass= "<<ionspecies[i].Mass_cal
    <<" | AME20_ME= "<<ionspecies[i].AME<<" +- "<<ionspecies[i].AME_err
    <<fixed<<setprecision(4)
    <<" | ME= "<< ionspecies[i].MassExcess_cal<<" +- "<<ionspecies[i].Mass_cal_err 
    <<" | dm= "<<ionspecies[i].deltaMass<<" +- "<<ionspecies[i].deltaMass_err;
         
    if(MASS_VER>=3)outfile<<fixed<<setprecision(4)
    <<" | ME_VE= "<< ionspecies[i].MassExcess_cal_VE <<" +- "<<ionspecies[i].Mass_cal_err_VE
    <<" | dm_VE= "<<ionspecies[i].deltaMass_VE<<" +- "<<ionspecies[i].deltaMass_err_VE;
    else if(MASS_VER==1)outfile<<fixed<<setprecision(4)
    <<" | no_ME_VE= "<< 0 <<" +- "<<0
    <<" | no_dm_VE= "<<0<<" +- "<<0;
    else if(MASS_VER==2)outfile<<fixed<<setprecision(4)
    <<" | ME_v2= "<< ionspecies[i].MassExcess_cal_v2 <<" +- "<<ionspecies[i].Mass_cal_err_v2
    <<" | dm_v2= "<<ionspecies[i].deltaMass_v2<<" +- "<<ionspecies[i].deltaMass_err_v2;

    else outfile
    <<" | ME_v?= "<< 0 <<" +- "<<0
    <<" | dm_v?= "<<0<<" +- "<<0;
             
    outfile<<" | ions: "<<ionspecies[i].N_unknown<<endl;
    
    }
}

if(OUTPUT_MASSRESULT)outfile.close();


if(MASS_VER==2)
{
    cout<<endl<<"---- for mass ver2 cout Isomer ---"<<endl;
    for(int i=0;i<NSpecies_Isomer  ;i++)
    {
        cout<<Isomer_ionspecis[i].Aname<<" "<<fixed<<setprecision(1)
        <<" | AME20_ME= "<<Isomer_ionspecis[i].AME<<" +- "<<Isomer_ionspecis[i].AME_err
        <<" | ME_v2= "<< Isomer_ionspecis[i].MassExcess_cal_v2 <<" +- "<<Isomer_ionspecis[i].Mass_cal_err_v2
        <<" | dm_v2= "<<Isomer_ionspecis[i].deltaMass_v2<<" +- "<<Isomer_ionspecis[i].deltaMass_err_v2<<endl;
    }
    
    
}
    
//cout<<" return at 7306"<<endl;return;
///////////////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//       所有的主要计算过程已经结束 
//       下面开始各种画图 
//       根据是否有循环扫描 LOOP_ON 分成两大部分 
//       如果没有循环扫描 则打开一系列 canvas


// =====================================   Show Mass Result ========================================
for(int i=0;i<Show_Mass_Result_opt_n;i++){Show_Mass_Result_opt[i]=0;}




                                                 if(LOOP_ON)
                                                { // loop on 循环打开 需要重置很多New 节省内存


gr2d_LddtC->SetPoint(scan_loop_i-1,L,ddT,C_intersection);


if(MASS_VER==1)Show_Mass_Result_opt[0]=1;
if(MASS_VER==2)Show_Mass_Result_opt[1]=1;
if(MASS_VER>=3)Show_Mass_Result_opt[2]=1;
if(MASS_VER>=3&&Show_MASS_VER1_ON)Show_Mass_Result_opt[0]=1;
if(Only_draw_mass_VE){ Show_Mass_Result_opt[0]=0;Show_Mass_Result_opt[1]=0;  Show_Mass_Result_opt[2]=1;  }
double Save_YAxis_UP, Save_YAxis_LOW=0;
Save_YAxis_LOW = -40;
Save_YAxis_UP = 40;
if(L_ddt_correlation_ON){Save_YAxis_LOW = -40;Save_YAxis_UP = 40;}
if(Show_Mass_Result_ON)Show_Mass_Result(ionspecies,NSpecies,0,ThisParaInfo,true,true,Save_YAxis_LOW,Save_YAxis_UP,false,0, Show_Mass_Result_opt,Show_Mass_Result_opt_n,false,false,true);
else cout<<" --------------Show_Mass_Result OFF----------"<<endl;

gr2d_Lddt_Xn_v1->SetPoint(scan_loop_i-1,L,ddT,L_ddT_Xn2[scan_loop_i]);  //L_ddT_Xn2 就是从 scan_loop_i =1 开始有数据

//
//for (int i = 0; i < NSpecies; i++)
//{
//    if (ionspecies[i].h_each_ref_cal_mass_err != NULL) delete ionspecies[i].h_each_ref_cal_mass_err;
//    if (ionspecies[i].h_iont_mass_err != NULL)delete ionspecies[i].h_iont_mass_err;
//    if (ionspecies[i].h_iont_chi_n != NULL)delete ionspecies[i].h_iont_chi_n;
//    if (ionspecies[i].h2_gtC != NULL)delete ionspecies[i].h2_gtC;
//    if (ionspecies[i].h_C != NULL)delete ionspecies[i].h_C;
//    if (ionspecies[i].h_T != NULL)delete ionspecies[i].h_T;
//
//}

//====================================  reset   =====================================================
for(int i =0;i<NSpecies; i++)
{

    ionspecies[i].Mass_cal=0;ionspecies[i].Mass_cal_err=0;ionspecies[i].N_unknown=0;
    ionspecies[i].stdDeviation_V1=0;
    //v2
    ionspecies[i].Mass_cal_v2=0;ionspecies[i].Mass_cal_err_v2=0;
    ionspecies[i].stdDeviation_VE=0;
    ionspecies[i].ClearRecord();

    
    

    for(int j=ionspecies[i].gr_BpC->GetN();j>=0  ;j--){ionspecies[i].gr_BpC->RemovePoint(j);}
    for(int j=ionspecies[i].gr_lnBpC->GetN();j>=0  ;j--){ionspecies[i].gr_lnBpC->RemovePoint(j);}
    for(int j=ionspecies[i].gr_mvqC->GetN();j>=0  ;j--){ionspecies[i].gr_mvqC->RemovePoint(j);}
    for(int j=ionspecies[i].gr_dmC->GetN();j>=0  ;j--){ionspecies[i].gr_dmC->RemovePoint(j);}
    for(int j=ionspecies[i].grerr_mvqC_VE->GetN();j>=0  ;j--){ionspecies[i].grerr_mvqC_VE->RemovePoint(j);}
    for(int j=ionspecies[i].grerr_dmC_VE->GetN();j>=0  ;j--){ionspecies[i].grerr_dmC_VE->RemovePoint(j);}
    for(int j=ionspecies[i].grerr_avegtC->GetN();j>=0  ;j--){ionspecies[i].grerr_avegtC->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gtC_own->GetN();j>=0  ;j--){ionspecies[i].gr_gtC_own->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gtC_own_u->GetN();j>=0  ;j--){ionspecies[i].gr_gtC_own_u->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gtC_own_d->GetN();j>=0  ;j--){ionspecies[i].gr_gtC_own_d->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gtC_shifted_own->GetN();j>=0  ;j--){ionspecies[i].gr_gtC_shifted_own->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gtC_shifted_own_u->GetN();j>=0  ;j--){ionspecies[i].gr_gtC_shifted_own_u->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gtC_shifted_own_d->GetN();j>=0  ;j--){ionspecies[i].gr_gtC_shifted_own_d->RemovePoint(j);}
    for(int j=ionspecies[i].gr_gt_C->GetN();j>=0  ;j--){ionspecies[i].gr_gt_C->RemovePoint(j);}

    //h1->Reset
    ionspecies[i].h2_gtC->Reset();
    ionspecies[i].h_C->Reset();
    ionspecies[i].h_T->Reset();
    ionspecies[i].h_mvq->Reset();
    
    ionspecies[i].h_each_ref_cal_mass_err->Reset();
    ionspecies[i].h_iont_mass_err->Reset(); 
    ionspecies[i].h_iont_chi_n->Reset();

            
} 
for(int i=0;i<NSpecies_IAS ;i++)
{
    //IAS_ionspecis[i].ClearRecord();
    IAS_ionspecis[i].Mass_cal = 0;
    IAS_ionspecis[i].Mass_cal_err = 0 ;
    
    IAS_ionspecis[i].N_unknown=0;   
}
for(int i=0;i<NSpecies_Isomer ;i++)
{
    //Isomer_ionspecis[i].ClearRecord();
    Isomer_ionspecis[i].Mass_cal = 0;
    Isomer_ionspecis[i].Mass_cal_err = 0 ;
    
    Isomer_ionspecis[i].N_unknown=0;   
}
//------------ TGraph remove point
for(int i=gr1->GetN();i>=0;i--){gr1->RemovePoint(i);} //由于每次删除一个点后，会自动缩进，因此必须从后往前删除

//循环需要清掉
//---- 注意gr_gtC_chosen 和  grerr_gtC_chosen在之前平滑的时候已经delete 掉了，之后和  gr_gtC_chosen_smooth  grerr_gtC_chosen_smooth 指向同一片空间
//因此 gr_gtC_chosen  和 gr_gtC_chosen_smooth 只需要delete 掉其中一个就行，不能再次delete
//grerr_gtC_chosen 和 grerr_gtC_chosen_smooth  只需要delete 掉其中一个就行，不能再次delete

delete gr_gtC_readin;
delete gr_gtC_chosen; //gr_gtC_chosen  = gr_gtC_chosen_smooth 
delete grerr_gtC_chosen; //grerr_gtC_chosen = grerr_gtC_chosen_smooth

delete gr_gtC_chosen_u;
delete gr_gtC_chosen_d;

delete grerr_gtC_chosen_v2;
delete grerr_gtC_with_err;
delete grerr_gtC_all;
delete grerr_gtC_s;
delete grerr_gtC_ave;
delete grerr_gtC_v2;
delete grerr_avegtC_inj_1;
delete grerr_avegtC_inj_2;
delete grerr_avegtC_inj_3;
delete grerr_avegtC_inj_4;
delete gr_err_dm_all;

delete gr_gtC_chosen_step1;
delete grerr_gtC_chosen_step1;
delete gr_gtC_chosen_step2;
delete grerr_gtC_chosen_step2;


delete h2_gtC_v2;
delete h1_STEP0;
delete h1_STEP1;
delete h1;
delete h1_E;
delete h_locate;
delete h2_err_dm_all;
delete h2_massv1_C_all;

for(int i =0;i<NSpecies; i++){delete ionspecies[i].h_dm;}  // 在循环中每次new 循环结束 delete

delete fitfun_gtC ;
delete fitfun_gtC_s ;
delete fitfun_gtC0_s ;

//delete sp3;
//delete sp3_smooth ;


//_________________________________________  reset   ___________________________________________ 



                                                }//loop on








                                                  if(!LOOP_ON)
                                                  { // DRAW ON  保证没有循环的情况下才画图，只运行一遍，所以不用考虑重复使用
//================================================        DRAW            ====================================================================================
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);

// =====================================   Show Mass Result ========================================
if(Show_Mass_Result_ON)
{
    gr_errsys_chin = new TGraph(); //在Showmass函数里面使用，记录卡方随着添加不同系统误差的变化 global, for the usage in Show_Mass_Result
    gr_errsys_chin_n = 0;
    AxisFormat(gr_errsys_chin,"","#sigma_{sys} [keV]","#chi_{n}");
    gr_errsys_chin->SetMarkerColor(kAzure+2);
    //old Show_Mass_Result(ionspecies,NSpecies,ThisParaInfo,true,true,-200,200,false,0);
    Show_Mass_Result_opt[0]=1;
    // [0] 等权 只有已知核 刻度
    Show_Mass_Result(ionspecies,NSpecies,0,ThisParaInfo,true,true,-40,40,false,0,Show_Mass_Result_opt,Show_Mass_Result_opt_n,false,false,true);
    // [1] 等权 所有核 save another graph showing all ions
    Show_Mass_Result(ionspecies,NSpecies,1,ThisParaInfo,true,true,-200,200,false,0,Show_Mass_Result_opt,Show_Mass_Result_opt_n,false,true,true);


    Show_Mass_Result_opt[0]=0;
    if(Show_MASS_VER1_ON||MASS_VER==1)Show_Mass_Result_opt[0]=1; 
    if(Only_draw_mass_VE){ Show_Mass_Result_opt[0]=0;Show_Mass_Result_opt[1]=0;  Show_Mass_Result_opt[2]=1;  }  
    if(MASS_VER==2)Show_Mass_Result_opt[1]=1;
    if(MASS_VER>=3)Show_Mass_Result_opt[2]=1;
    // [2] 误差加权 只有已知核 刻度
    Show_Mass_Result(ionspecies,NSpecies,2,ThisParaInfo,true,true,-60,60,false,0,Show_Mass_Result_opt,Show_Mass_Result_opt_n,false,false,true);
    
    // [3] 误差加权 所有核 save another graph showing all ions
    Show_Mass_Result(ionspecies,NSpecies,3,ThisParaInfo,true,true,-200,200,false,0,Show_Mass_Result_opt,Show_Mass_Result_opt_n,false,true,true);

    // ERR_SYS_ON
    //for(int i=0;i<20;i++)
    //{
        //Show_Mass_Result(ionspecies,NSpecies,i,ThisParaInfo,false,true,-200,200,true,i*0.5,Show_Mass_Result_opt,Show_Mass_Result_opt_n);
    //}
    //Show_Mass_Result(ionspecies,NSpecies,21,ThisParaInfo,true,true,-200,200,true,6,Show_Mass_Result_opt,Show_Mass_Result_opt_n);
    //Draw_one_TGraph(gr_errsys_chin,"c_gr_errsys_chin");

}
else cout<<" --------------Show_Mass_Result OFF----------"<<endl;
//------------------------------------------
//Draw_one_TGraph(gr_dm_dA0err_all,"c_gr_dm_dA0err_all");

//========================================= 1、 c0=========================================================
if(c0_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c0 = new TCanvas("c0","c0 title",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gr_n=0;
for(int i=0;i<NSpecies  ;i++)
{
    grerr0->SetPoint(grerr0_n, ionspecies[i].AveT*1000 , ionspecies[i].SigmaT*1000);
    grerr0->SetPointError(grerr0_n, 0 , ionspecies[i].SigmaT*1000/(sqrt(2*ionspecies[i].N -2) ) );
    grerr0_n++;
    //cout<<ionspecies[i].AveC<<"&&&&&&&&7"<<ionspecies[i].Avevgammat<<endl;   
    //cout<<ionspecies[i].Aname<<" "<<ionspecies[i].AveT*1000<<" "<<ionspecies[i].SigmaT*1000<<" +- "<<ionspecies[i].SigmaT*1000/(sqrt(2*ionspecies[i].N -2) )
    //<<" "<<ionspecies[i].N<<endl;
}
AxisFormat(grerr0,"SigmaT vs AveT for all ionspecies","AveT [ps]","SigmaT [ps]",2);
grerr0->Draw("ap");

lat_n->SetTextColor(1);
lat_n->SetTextFont(43);
lat_n->SetTextSize(20);
for(int i=0;i<NSpecies  ;i++)
{
    
    //if(ionspecies[i].N<200)
    {
        lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].N);
        lat_n->DrawLatex(ionspecies[i].AveT*1000,ionspecies[i].SigmaT*1000,lat_text);
    }   
}
}//____________________________________________c0________________________________________

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_C_distribution = new TCanvas("c_C_distribution","c_C_distribution",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h1->DrawClone();
*/


if(c_h1_CDistribution_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_h1_CDistribution = new TCanvas("c_h1_CDistribution","c_h1_CDistribution",1000,500);
c_h1_CDistribution->cd(1);
for(int i=0;i<ions_n  ;i++)
{        
    h1->Fill(ions[i].C);

}
AxisFormat(h1,"C distribution ddt0","C [m]","Counts",1);

h1->DrawClone("");
c_h1_CDistribution->BuildLegend();
double Bpmax=-1;
double Bpmin=99;
TCanvas *c_h_Bp = new TCanvas("c_h_Bp","c_h_Bp",1000,500);
TH1F* h_Bp = new TH1F("h_Bp","h_Bp",100,4.8,4.85);
AxisFormat(h_Bp,"Bp distribution ","Bp [Tm]","Counts",1);
for(int i=0;i<ions_n  ;i++)
{        
    if(ionspecies[ions[i].Species].MassUnknown==0)
    {
        h_Bp->Fill(ions[i].Bp);
        if(ions[i].Bp>Bpmax)Bpmax = ions[i].Bp;
        if(ions[i].Bp<Bpmin)Bpmin = ions[i].Bp;
    }
}
h_Bp->Draw();
cout<<" Bpmax ="<<Bpmax<<" Bpmin= "<<Bpmin<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
}


//=============================  c1  =========================================================
if(c1_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c1_Bp_C = new TCanvas("c1_Bp_C","c1_Bp_C",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//0928------------
TGraph* gr_parchange = new TGraph();
int gr_parchange_n=0;
AxisFormat(gr_parchange,"","ions ","k_lnBpC_all");
TF1* fit_tmp = new TF1("fit_tmp","pol1",0,10);
TGraph* gr1_part = new TGraph();
int gr1_part_n=0;
AxisFormat(gr1_part,"","lnC ","lnBp",kAzure+7);
//0928---------------

//TH2F* hh_BpC = new TH2F("hh_BpC", "all ions Bp_C ",100,128.65,129.10,500, 5.46, 5.48 );
gr_n=0;
for(int i=0;i<ions_n  ;i++)
{
    //hh_BpC->Fill();
    gr1->SetPoint(i,(ions[i].C),( ions[i].Bp) );
    //if(ions[i].A==22&&ions[i].Z==14)cout<<i<<" is Si 22 "<<endl;
    //0928----
    if(i<10000)gr1_part->SetPoint(i,(ions[i].C),( ions[i].Bp) );
    
    //if(i%100==0){ gr1->Fit(fit_tmp,"q");gr_parchange->SetPoint(gr_parchange_n++,i,fit_tmp->GetParameter(0)); }
    //0928-----
    //gr1->SetPoint(i,(ions[i].Z),( ions[i].C) );
}
//hh_BpC->Draw("colz");
AxisFormat(gr1,"Bp_C distribution of all ions"," C [m]"," Bp[T]",1);
gr1->SetMarkerColor(kRed+1);
//AxisFormat(gr1,"Bp_C distribution of all ions"," lnC [C:m]","lnBp [Bp:Tm]",1);
//AxisFormat(gr1,"Bp_C distribution of all ions"," Z","C [m]",1);
gr1->SetMarkerSize(1);
gr1->Draw("ap");
gr1_part->SetMarkerSize(1);
//gr1_part->Draw("samep"); //0928

/* 0928
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c1_ss = new TCanvas("c1_ss","c1_ss",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gr_parchange->Draw("ap");
*/
}//____________________________  c1  ______________________________________



if(c1_Bp_C_inj_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c1_Bp_C_inj = new TCanvas("c1_Bp_C_inj","c1_Bp_C_inj",1600,800);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int gr2d_BpCInj_n=0;
TGraph2D *gr2d_BpCInj = new TGraph2D();
for(int i=0;i<ions_n  ;i++)
{
    //if(ions[i].inject_number<4000&&ions[i].inject_number>3900)
    {
        gr2d_BpCInj->SetPoint(gr2d_BpCInj_n++,ions[i].C,ions[i].Bp,ions[i].inject_number);
    }
}
gr2d_BpCInj->Draw("pcol");
}//________if(c1_Bp_C_inj_ON)

//1013================================c1_Bp_C_ionspecies===============================================
if(c1_Bp_C_ionspecies_ON)
{

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c1_Bp_C_ionspecies = new TCanvas("c1_Bp_C_ionspecies","c1_Bp_C_ionspecies",1600,800);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ionspecies[0].gr_BpC->Draw("ap");
for(int i=0;i<NSpecies  ;i++)
{
    ionspecies[i].gr_BpC->Draw("samep");
}
auto legend_BpC_all = new TLegend(0.80,0.10,0.95,0.35); //downright
    legend_BpC_all->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"C");
for(int i=0;i<NSpecies  ;i++)
{
    //legend_BpC_all->AddEntry(ionspecies[i].gr_lnBpC,ionspecies[i].Aname,"p" );
    legend_BpC_all->AddEntry(ionspecies[i].gr_BpC,ionspecies[i].Aname,"p" );
}
legend_BpC_all->Draw("same");


}
//______________________________c1_Bp_C_ionspecies_______________________________________________________________

// ==================================== c2 =======================================
if(c2_ON)
{
cout<<" ================ c2 ON ======================="<<endl;

TCanvas *c8 = new TCanvas("c2","c2",2000,1000);
if(SHOW_gtC_step0)h2_gtC_all_step0->Draw("colz");
//if(SHOW_gtC_step2)h2_gtC_all_step2->Draw("colz");
grerr_gtC_chosen->Draw("plsame");
//c8->Update();
//h2_gtC_all_step0->GetXaxis()->SetRangeUser(128.69, 129.01);

}//____________________________________ c2 ________________________________________


//============================= c_h2_gtC_all_ON =======================================
if(c_h2_gtC_all_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_h2_gammat_C = new TCanvas("c_h2_gammat_C","c_h2_gammat_C",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TH2F* h2_gammat_C = new TH2F ("h2_gammat_C","h2_gammat_C", 250,128.5,129.0, 100,1.3,1.4);
AxisFormat(h2_gammat_C,""," C [m] ","#gamma_{t} v2");
for(int i=0;i< ions_n ;i++)
{
    h2_gammat_C->Fill(ions[i].C,ions[i].gammat);
}
h2_gammat_C->Draw("colz");

}//_____________________________ c_gtC_all ____________________________________________

/* 
======================== two ways of calculating gt ==========================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_h_dgt = new TCanvas("c_h_dgt","c_h_dgt",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1F* h_dgt = new TH1F("h_dgt","h_dgt",2000,-1,1);
AxisFormat(h_dgt,"","gammatv2 - gammat","Counts");
for(int i=0;i<ions_n  ;i++)
{
    h_dgt->Fill(ions[i].gammat - ions[i].gammat_v2);
    
}
h_dgt->Draw();
_________________________ two ways of calculating gt ____________________________
*/

if(c_gtC_chosen_ON)
{//============================= c_gtC_chosen_ON =======================================
    cout<<"============================ c_gtC_chosen_ON ========================"<<endl;
/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_grerr_gtC_all = new TCanvas("c_grerr_gtC_all","c_grerr_gtC_all",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grerr_gtC_all->Draw("ap");
grerr_gtC_ave->Draw("psame");
//cout<<"grerr_gtC_ave->Print(); "<<endl;
//grerr_gtC_ave->Print();
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_grerr_gtC_chosen = new TCanvas("c_grerr_gtC_chosen","c_grerr_gtC_chosen",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grerr_gtC_chosen->Draw("apl");
c_grerr_gtC_chosen->Print(FILEPATH_chosen+"grerr_gtC_chosen.png");


lat_n->SetTextSize(20);
for(int i=0;i<subregion_n  ;i++)
{
    if(C_Division_n_chosen[i]<1.2*C_DIVISION_CHOSEN_MIN)
    {
        lat_text = lat_text.Format("  %d",C_Division_n_chosen[i]);
        //lat_n->DrawLatex(h1_E->GetBinCenter(i+1),gtC_chosen[i],lat_text);
    }
}
//TGraphErrors_to_outfile(FILEPATH+"grerr_gtC_chosen-data.txt",grerr_gtC_chosen);

outfile_grerr_txt.open(FILEPATH+"grerr_gtC_chosen-data.txt");
for(int i=0;i<grerr_gtC_chosen->GetN();i++)
{
    grerr_gtC_chosen->GetPoint(i,xtmp,ytmp);
    yerrtmp=grerr_gtC_chosen->GetErrorY(i);
    outfile_grerr_txt<<xtmp<<" "<<ytmp<<" "<<yerrtmp<<endl;
}
cout<<"-------TGraphErrors_to_outfile : save "<<grerr_gtC_chosen->GetTitle()<<" to "<<FILEPATH+"grerr_gtC_chosen-data.txt"<<endl;
outfile_grerr_txt.close();


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//TCanvas *c_gr_gtC_chosen = new TCanvas("c_gr_gtC_chosen","c_gr_gtC_chosen",1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//gr_gtC_chosen->Draw("apl");


//gr_gtC_chosen_smooth->Draw("psame");
/*
sp3->SetLineColor(kAzure-7);
sp3->SetLineWidth(3);
//sp3->Draw("Csame");
sp3_smooth->SetLineColor(kOrange-3);
sp3_smooth->SetLineWidth(3);
//sp3_smooth->Draw("Csame");
*/
//c_gr_gtC_chosen->BuildLegend();


}//____________________________________c_gtC_chosen_ON___________________________________________________


//===================================== c_gtC_v2 ==================================================
if(c_gtC_v2_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas* c_gtC_v2 = new TCanvas("c_gtC_v2","c_gtC_v2",1800,1000);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AxisFormat(grerr_gtC_v2,"","C[m]","#gamma_{t}");
//grerr_gtC_v2->Draw("Ap");

TGraphErrors* grerr_gtC_v2_make = new TGraphErrors();
GetGtC_Curve(grerr_gtC_v2,grerr_gtC_v2_make,subregion_n,h1);
AxisFormat(grerr_gtC_v2_make,"","C[m]","#gamma_{t}");
grerr_gtC_v2_make->Draw("ap");
TGraph_to_outfile("OUTPUT//grerr_gtC_v2_make.txt",grerr_gtC_v2_make);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas* c_gtC_v2_2 = new TCanvas("c_gtC_v2_2","c_gtC_v2_2",1800,1000);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AxisFormat(h2_gtC_v2,"","C[m]","#gamma_{t}");
//h2_gtC_v2->Draw("colz");
}
//_____________________________  c_gtC_v2  ________________________________





//============================= c_gtC_Tz_ON =======================================
//------------0929
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(c_gtC_Tz_ON)
{

TCanvas *c_0929 = new TCanvas("c_0929","gt-C Tz "  ,1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TGraphErrors* grerr_avegtC_Tz[TzN];
int grerr_avegtC_Tz_n[TzN];
for(int i=0;i<TzN  ;i++){grerr_avegtC_Tz_n[i]=0;}  // initialization
for(int j=0;j<TzN  ;j++)
{   
    grerr_avegtC_Tz[j] = new TGraphErrors();
    AxisFormat(grerr_avegtC_Tz[j],"","C","#gamma_{t}",2+j);
}
for(int i=0;i<TzN  ;i++)
{
    
    for(int j=0;j<subregion_n  ;j++)
    {
        if(C_Division_n_Tz[i][j]>3)
        {
            gtC_Tz[i][j] /= C_Division_n_Tz[i][j];
            gtC2_Tz[i][j] /= C_Division_n_Tz[i][j];
            sigma_gtC_Tz[i][j] = sqrt (gtC2_Tz[i][j]-gtC_Tz[i][j]*gtC_Tz[i][j]) ;
            grerr_avegtC_Tz[i]->SetPoint( grerr_avegtC_Tz_n[i],h1->GetBinCenter(j+1), gtC_Tz[i][j] );
            grerr_avegtC_Tz[i]->SetPointError(grerr_avegtC_Tz_n[i], 0, sigma_gtC_Tz[i][j]/sqrt( C_Division_n_Tz[i][j]-1) ) ;
            grerr_avegtC_Tz_n[i]++;
        }        
    }
    //cout<<"grerr_avegtC_Tz_n["<<i<<"] "<<grerr_avegtC_Tz_n[i]<<endl;
}


for(int i=0;i<TzN  ;i++)
{
    if(i==0)grerr_avegtC_Tz[i]->Draw("ap");
    else grerr_avegtC_Tz[i]->Draw("samep");
    //grerr_avegtC_Tz[i]->Print();
   
}
grerr_avegtC_Tz[0]->GetYaxis()->SetRangeUser(1.3,1.4);

auto legend0929 = new TLegend(0.80,0.10,0.95,0.25); //downright
    legend0929->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"C");
for(int i=0;i<TzN  ;i++)
{
    legend0929->AddEntry(grerr_avegtC_Tz[i], strtmp.Format("  Tz = %.2f   %d",Tz_all_ions[i], count_Tz_n[i]),"P");
}
legend0929->Draw("same");

    
//------------------------0929
}//__________________________________________ c_gtC_Tz_ON _______________________________________



if(c_avegtC_ionspecies_ON)
{//============================= c_avegtC_ionspecies =======================================
//------------0930
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_avegtC_ionspecies = new TCanvas("c_avegtC_ionspecies","gt-C for all ionspecies "  ,1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//c_avegtC_ionspecies->cd(1);
TMultiGraph* mg_avegtC_ionspecies = new TMultiGraph();
AxisFormat(mg_avegtC_ionspecies,"avegtC_ionspecies","C [m]","ave #gamma_{t}");
for(int i=0;i<NSpecies  ;i++)
{
    if(ionspecies[i].Has_own_avegtC)
    {
        if(THIS_EXP == "2017_58Ni")
        {
            if(ionspecies[i].Aname=="20Na"){continue;}
            if(ionspecies[i].Aname=="14O"){continue;}
            if(ionspecies[i].Aname=="22Mg"){continue;}
            if(ionspecies[i].Aname=="18Ne"){continue;}
            if(ionspecies[i].Aname=="13N"){continue;}
        }
        
        //mg_avegtC_ionspecies->Add( ionspecies[i].grerr_avegtC,"pl");
        mg_avegtC_ionspecies->Add( ionspecies[i].gr_gtC_own,"pl");
        /*
        mg_avegtC_ionspecies->Add( ionspecies[i].gr_gtC_own,"pl");
        mg_avegtC_ionspecies->Add( ionspecies[i].gr_gtC_own_u,"pl");
        mg_avegtC_ionspecies->Add( ionspecies[i].gr_gtC_own_d,"pl");
        */
    
    }
}

mg_avegtC_ionspecies->Draw("ap ");
mg_avegtC_ionspecies->GetXaxis()->SetRangeUser(128.6,128.9);
mg_avegtC_ionspecies->GetYaxis()->SetRangeUser(gtC_chosen_gtMIN,gtC_chosen_gtMAX);
//c_avegtC_ionspecies->BuildLegend();
auto legend_gtC_species = new TLegend(0.80,0.20,0.95,0.9); //downright
legend_gtC_species->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"");
for(int i=0;i<NSpecies  ;i++)
{
    if(ionspecies[i].Has_own_avegtC)
    {
        //legend_gtC_species->AddEntry(ionspecies[i].grerr_avegtC, ionspecies[i].Aname+strtmp.Format(" %d", ionspecies[i].count_avegtC_n),"P");
        legend_gtC_species->AddEntry(ionspecies[i].gr_gtC_own, ionspecies[i].Aname+strtmp.Format(" %d", ionspecies[i].count_avegtC_n),"P");
    }
}
legend_gtC_species->Draw("same");

if(Do_gtC_divide_Z)
{
    
    TGraphErrors* grerr_avegtC_l = new TGraphErrors();// 0~low Z
    TGraphErrors* grerr_avegtC_h = new TGraphErrors();// > high Z
    AxisFormat(grerr_avegtC_l,"","C[m]","ave #gamma_{t}",2);
    AxisFormat(grerr_avegtC_h,"","C[m]","ave #gamma_{t}",4);
    
    int tmp_count_lZ=0;int tmp_count_hZ=0;
    int division_Z = CONDITION_gt_Z;
    GetGtC_Curve_Z_divide( grerr_avegtC_l,  h1_E,  0, 99999,   gtC_chosen_gtMIN,gtC_chosen_gtMAX,  0,division_Z,   tmp_count_lZ);
    GetGtC_Curve_Z_divide( grerr_avegtC_h,  h1_E,  0, 99999,   gtC_chosen_gtMIN,gtC_chosen_gtMAX,  division_Z+1,999, tmp_count_hZ);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c_avegtC_ionspecies_Z = new TCanvas("c_avegtC_ionspecies_Z","avegtC light&heavy ions by Z"  ,1000,500);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grerr_avegtC_l->GetXaxis()->SetRangeUser(128.5,128.9);
    grerr_avegtC_l->GetYaxis()->SetRangeUser(gtC_chosen_gtMIN,gtC_chosen_gtMAX);
    grerr_avegtC_l->Draw("ap");
    grerr_avegtC_h->Draw("psame");
    
    auto legend_avegtC_spe_Z = new TLegend(0.80,0.10,0.95,0.25); //downright
    legend_avegtC_spe_Z->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"C");
    
    legend_avegtC_spe_Z->AddEntry(grerr_avegtC_l, strtmp.Format("Z<= %d  %d", division_Z , tmp_count_lZ),"P");
    legend_avegtC_spe_Z->AddEntry(grerr_avegtC_h, strtmp.Format("Z>= %d  %d", division_Z+1, tmp_count_hZ),"P");
    
    legend_avegtC_spe_Z->Draw("same");
}
///////////////////////////////////////
if(Do_gtC_divide_A)
{
    TGraphErrors* grerr_avegtC_l_A = new TGraphErrors(); // 0~low A
    TGraphErrors* grerr_avegtC_h_A = new TGraphErrors(); // > high A
    AxisFormat(grerr_avegtC_l_A,"","C[m]","ave #gamma_{t}",2);
    AxisFormat(grerr_avegtC_h_A,"","C[m]","ave #gamma_{t}",4);
    
    int tmp_count_lA=0;int tmp_count_hA=0;
    int division_A = CONDITION_gt_A;
    GetGtC_Curve_A_divide( grerr_avegtC_l_A,  h1_E,  0, 99999,   gtC_chosen_gtMIN,gtC_chosen_gtMAX,  0,division_A,   tmp_count_lA);
    GetGtC_Curve_A_divide( grerr_avegtC_h_A,  h1_E,  0, 99999,   gtC_chosen_gtMIN,gtC_chosen_gtMAX,  division_A+1,999, tmp_count_hA);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c_avegtC_ionspecies_A = new TCanvas("c_avegtC_ionspecies_A","avegtC light&heavy ions by A"  ,1000,500);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grerr_avegtC_l_A->GetXaxis()->SetRangeUser(128.6,128.9);
    grerr_avegtC_l_A->GetYaxis()->SetRangeUser(gtC_chosen_gtMIN,gtC_chosen_gtMAX);
    grerr_avegtC_l_A->Draw("ap");
    grerr_avegtC_h_A->Draw("psame");
    
    auto legend_avegtC_spe_A = new TLegend(0.80,0.10,0.95,0.25); //downright
    legend_avegtC_spe_A->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"C");
    
    legend_avegtC_spe_A->AddEntry(grerr_avegtC_l_A, strtmp.Format("A<= %d  %d", division_A , tmp_count_lA),"P");
    legend_avegtC_spe_A->AddEntry(grerr_avegtC_h_A, strtmp.Format("A>= %d  %d", division_A+1, tmp_count_hA),"P");
    
    legend_avegtC_spe_A->Draw("same");
}





}//_______________________________c_avegtC_ionspecies_______________________________________________


//================================= each ionspecies =======================================
if(c_each_ionspecies_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(show_each_h_C_gtC_ON)
{

//2023 所有种类的 canvas divide上下 ：上面是gammat-C 二维直方图 下面是 C分布直方图
TCanvas*c_gtC_ionspecies[30];
for(int i=0;i<NSpecies  ;i++)
{
    c_gtC_ionspecies[i] = new TCanvas("c_gtC_ionspecies_"+ionspecies[i].Aname,"c_gtC_ionspecies"+ionspecies[i].Aname  ,800,1400);
    c_gtC_ionspecies[i]->Divide(1,3,0.01,0.01);

    c_gtC_ionspecies[i]->cd(1);
    ionspecies[i].gr_gt_C->GetXaxis()->SetRangeUser(128.5,129.0);
    ionspecies[i].gr_gt_C->GetYaxis()->SetRangeUser(1.2,1.4);
    ionspecies[i].gr_gt_C->Draw("ap");
    

    c_gtC_ionspecies[i]->cd(2);
    ionspecies[i].h2_gtC->Draw("colz");
    
    c_gtC_ionspecies[i]->cd(3);
    ionspecies[i].h_C->Draw();

    if(savefile_each_h_C_ON)c_gtC_ionspecies[i]->Print(ionspecies[i].folder_path+"h_C_and_h2_gtC.png");

}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// T distribution
if(show_each_h_T_ON)
{
TCanvas*c_ionspecies_h_T[MAX_IONSPECIES];
for(int i=0;i<NSpecies  ;i++)
{
    c_ionspecies_h_T[i] = new TCanvas("c_ionspecies_h_T_"+ionspecies[i].Aname,"c_ionspecies_h_T"+ionspecies[i].Aname  ,1000,500);
    ionspecies[i].h_T->Draw();
    if(savefile_each_h_T_ON)c_ionspecies_h_T[i]->Print(ionspecies[i].folder_path+"h_T.png");
}
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(show_each_h_mvq_ON)
{
//20230811 h_mvq
TCanvas* c_h_mvq_ionspecies[MAX_IONSPECIES];
for(int i=0;i<NSpecies  ;i++)
{   
    if(outfile_each_h_mvq_ON)TGraph_to_outfile(FILEPATH+ionspecies[i].Aname+"_gr_mvqC_out.txt",ionspecies[i].gr_mvqC);
    if(ionspecies[i].h_mvq->GetEntries()<100)
    {
        cout<<"-- skip "<<ionspecies[i].Aname
        <<" h_mvq->GetEntries()= "<<ionspecies[i].h_mvq->GetEntries()<<" <100 --"<<endl;
        continue;
    }
    c_h_mvq_ionspecies[i] = new TCanvas("c_h_mvq_ionspecies_"+ionspecies[i].Aname,"c_h_mvq_ionspecies"+ionspecies[i].Aname  ,1000,500);
    
    //ionspecies[i].h_mvq = new TH1F(ionspecies[i].Aname+" h_mvq",ionspecies[i].Aname+" h_mvq ", 50000,1.5,2.0);
    //AxisFormat(ionspecies[i].h_mvq,ionspecies[i].Aname+" h_mvq ","#frac{m}{q} [u/e]", "Counts");
    //TGraph_to_TH1F(ionspecies[i].gr_mvqC,ionspecies[i].h_mvq);    
    ionspecies[i].h_mvq->DrawClone();
}

}//if(show_each_h_mvq_ON)


/*
//============= 20230419 ==========
int grerr_Cfilter_n =0;
TGraphErrors* grerr_Cfilter = new TGraphErrors();
AxisFormat(grerr_Cfilter,"","ave C [m]","ave #gamma_{t} ");

TGraphErrors* grerr_Cfilter_2 = new TGraphErrors();
AxisFormat(grerr_Cfilter_2,"","ave v [m/ns]","ave #gamma_{t} ");

TGraphErrors* grerr_Cfilter_3 = new TGraphErrors();
AxisFormat(grerr_Cfilter_3,"","ave T [ns]","ave #gamma_{t} ");

for(int j=0;j<ions_n  ;j++)
{
    if( (ions[j].C>128.7&&ions[j].C<128.9)&&(ions[j].gammat>0&&ions[j].gammat<20) )
    {
        ionspecies[ions[j].Species].Cfilter_n++;
        ionspecies[ions[j].Species].Cfilter_aveC += ions[j].C;
        ionspecies[ions[j].Species].Cfilter_avegt += ions[j].gammat;
        ionspecies[ions[j].Species].Cfilter_aveCstd += ions[j].C*ions[j].C;
        ionspecies[ions[j].Species].Cfilter_avegtstd += ions[j].gammat*ions[j].gammat;
        ionspecies[ions[j].Species].Cfilter_avev += ions[j].v;
        ionspecies[ions[j].Species].Cfilter_aveT += ions[j].T;
        ionspecies[ions[j].Species].Cfilter_avevstd += ions[j].v*ions[j].v;
        ionspecies[ions[j].Species].Cfilter_aveTstd += ions[j].T*ions[j].T;
    }
    
}
for(int i=0;i<NSpecies  ;i++)
{
    ionspecies[i].Cfilter_aveC /= ionspecies[i].Cfilter_n;
    ionspecies[i].Cfilter_avegt /= ionspecies[i].Cfilter_n;
    ionspecies[i].Cfilter_aveCstd = sqrt(ionspecies[i].Cfilter_aveCstd/ionspecies[i].Cfilter_n - ionspecies[i].Cfilter_aveC*ionspecies[i].Cfilter_aveC);
    ionspecies[i].Cfilter_avegtstd = sqrt(ionspecies[i].Cfilter_avegtstd/ionspecies[i].Cfilter_n - ionspecies[i].Cfilter_avegt*ionspecies[i].Cfilter_avegt);
    ionspecies[i].Cfilter_avev /= ionspecies[i].Cfilter_n;
    ionspecies[i].Cfilter_aveT /= ionspecies[i].Cfilter_n;
    ionspecies[i].Cfilter_avevstd = sqrt(ionspecies[i].Cfilter_avevstd/ionspecies[i].Cfilter_n - ionspecies[i].Cfilter_avev*ionspecies[i].Cfilter_avev);
    ionspecies[i].Cfilter_aveTstd = sqrt(ionspecies[i].Cfilter_aveTstd/ionspecies[i].Cfilter_n - ionspecies[i].Cfilter_aveT*ionspecies[i].Cfilter_aveT);
    if(ionspecies[i].Cfilter_n>100)
    {
        grerr_Cfilter->SetPoint(grerr_Cfilter_n,ionspecies[i].Cfilter_aveC,ionspecies[i].Cfilter_avegt);
        grerr_Cfilter_2->SetPoint(grerr_Cfilter_n,ionspecies[i].Cfilter_avev,ionspecies[i].Cfilter_avegt);
        grerr_Cfilter_3->SetPoint(grerr_Cfilter_n,ionspecies[i].Cfilter_aveT,ionspecies[i].Cfilter_avegt);
        //grerr_Cfilter->SetPointError(grerr_Cfilter_n,0,0);
        //grerr_Cfilter_2->SetPointError(grerr_Cfilter_n,0,0);
        //grerr_Cfilter_3->SetPointError(grerr_Cfilter_n,0,0);
        grerr_Cfilter->SetPointError(grerr_Cfilter_n, ionspecies[i].Cfilter_aveCstd,ionspecies[i].Cfilter_avegtstd);
        grerr_Cfilter_2->SetPointError(grerr_Cfilter_n, ionspecies[i].Cfilter_avevstd,ionspecies[i].Cfilter_avegtstd);
        grerr_Cfilter_3->SetPointError(grerr_Cfilter_n, ionspecies[i].Cfilter_aveTstd,ionspecies[i].Cfilter_avegtstd);
        grerr_Cfilter_n++;
    }    
}
lat_n->SetTextColor(2);
lat_n->SetTextFont(43);
lat_n->SetTextSize(40);

TCanvas*c_0419=new TCanvas("c_0419","c_0419" ,1000,500);
grerr_Cfilter->Draw("ap");
for(int i=0;i<NSpecies  ;i++)
{
    lat_n->SetTextColor(A_color[ionspecies[i].A]);
    lat_text = ionspecies[i].Aname+ strtmp.Format( ": %d", ionspecies[i].Cfilter_n);
    if(ionspecies[i].Cfilter_n>100)lat_n->DrawLatex(ionspecies[i].Cfilter_aveC,ionspecies[i].Cfilter_avegt, lat_text);
}

TCanvas*c_0419_2=new TCanvas("c_0419_2","c_0419_2" ,1000,500);
grerr_Cfilter_2->Draw("ap");
for(int i=0;i<NSpecies  ;i++)
{
    lat_n->SetTextColor(A_color[ionspecies[i].A]);
    lat_text = ionspecies[i].Aname+ strtmp.Format( ": %d", ionspecies[i].Cfilter_n);
    if(ionspecies[i].Cfilter_n>100)lat_n->DrawLatex(ionspecies[i].Cfilter_avev,ionspecies[i].Cfilter_avegt, lat_text);
}

TCanvas*c_0419_3=new TCanvas("c_0419_3","c_0419_3" ,1000,500);
grerr_Cfilter_3->Draw("ap");
for(int i=0;i<NSpecies  ;i++)
{
    lat_n->SetTextColor(A_color[ionspecies[i].A]);
    lat_text = ionspecies[i].Aname+ strtmp.Format( ": %d", ionspecies[i].Cfilter_n);
    if(ionspecies[i].Cfilter_n>100)lat_n->DrawLatex(ionspecies[i].Cfilter_aveT,ionspecies[i].Cfilter_avegt, lat_text);
}

//_________________ 20230419 ________________
*/

}//---------------------------------- each ionspecies  ---------------------------------------------


if(c3_h_INJ_ON)
{
    TH1F* h_ions_inj = new TH1F ("h_ions_inj","h_ions_inj", 31, -0.5,30.5); // 120 ns  step 50 ps
    AxisFormat(h_ions_inj,""," ions in one injection ","count");
    h_ions_inj->SetLineWidth(2);

    TH1F* h_ions_inj_m = new TH1F ("h_ions_inj_m","h_ions_inj_m", 31, -0.5,30.5); //
    AxisFormat(h_ions_inj_m,""," mass results in one injection ","count",kRed);
    h_ions_inj_m->SetLineWidth(2);
    int tmp_sum =0 ;
    int ions_more_than_1=0;
    int ions_more_than_2=0;
    int ions_more_than_3=0;
    for(int i=1;i<=total_injection  ;i++)
    {
        tmp_sum+=Injection[i];
        h_ions_inj->Fill(Injection[i]);
        h_ions_inj_m->Fill(Injection_m[i]);
        if(Injection[i]>1)ions_more_than_1++;
        if(Injection[i]>2)ions_more_than_2++;
        if(Injection[i]>3)ions_more_than_3++;
    }
    cout<<"    sum of ions from each injection = "<< tmp_sum<< " ions_n = "<<ions_n<<endl;
    cout<<" injections with ions > 1: "<<ions_more_than_1<<" "<<double(ions_more_than_1)/total_injection*100.0<<" % "<<endl;
    cout<<" injections with ions > 2: "<<ions_more_than_2<<" "<<double(ions_more_than_2)/total_injection*100.0<<" % "<<endl;
    cout<<" injections with ions > 3: "<<ions_more_than_3<<" "<<double(ions_more_than_3)/total_injection*100.0<<" % "<<endl;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas*c3_h_ions_inj=new TCanvas("c3_h_ions_inj","c3_h_ions_inj" ,1000,500);
    h_ions_inj->GetYaxis()->SetTitleOffset(1.0);
    h_ions_inj->GetXaxis()->SetNdivisions(120);
    h_ions_inj->GetXaxis()->SetRangeUser(-0.5,13.5);
    if(THIS_EXP=="2017_58Ni")h_ions_inj->GetXaxis()->SetRangeUser(-0.5,22.5);
    if(THIS_EXP=="2021_36Ar_SET3")h_ions_inj->GetXaxis()->SetRangeUser(-0.5,28.5);
    h_ions_inj->Draw("HIST TEXT0");
    //Tfile_save_hist->cd(); // 切换到 Tfile 但是会让之后的tcanvas全部空白！所以不使用
    Tfile_save_hist = new TFile(FILEPATH+"save_hist.root", "UPDATE");
    h_ions_inj->Write(THIS_EXP+"_h_ions_inj"); // 写入到 Tfile
    Tfile_save_hist->Close(); // Close 大写C
    h_ions_inj->SaveAs(FILEPATH+THIS_EXP+"_h_ions_inj.txt");
    
    //h_ions_inj_remain->Draw("same");!!

    //---- save file
    if(Do_outfile_inj_ions)
    {
        ofstream outfile_inj_ions;
        outfile_inj_ions.open(FILEPATH+"inj_ions_data.txt");
        for(int i=1;i<=total_injection  ;i++)
        {
            outfile_inj_ions<<" injection_number: "<<i<<" ions_in_this_injection: "<<Injection[i]<<endl;
        }
        outfile_inj_ions.close();
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas*c3_h_ions_inj_m=new TCanvas("c3_h_ions_inj_m","c3_h_ions_inj_m" ,1000,500);
    h_ions_inj_m->GetYaxis()->SetTitleOffset(1.0);
    h_ions_inj_m->GetXaxis()->SetNdivisions(120);
    h_ions_inj_m->GetXaxis()->SetRangeUser(-0.5,13.5);
    if(THIS_EXP=="2017_58Ni")h_ions_inj_m->GetXaxis()->SetRangeUser(-0.5,22.5);
    if(THIS_EXP=="2021_36Ar_SET3")h_ions_inj_m->GetXaxis()->SetRangeUser(-0.5,28.5);
    h_ions_inj_m->Draw("HIST TEXT0");

    Tfile_save_hist = new TFile(FILEPATH+"save_hist.root", "UPDATE");
    h_ions_inj_m->Write(THIS_EXP+"_h_ions_inj_m"); // 写入到 Tfile
    Tfile_save_hist->Close(); // Close 大写C
    h_ions_inj_m->SaveAs(FILEPATH+THIS_EXP+"_h_ions_inj_m.txt");

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas*c3_h_refs=new TCanvas("c3_h_refs","c3_h_refs" ,1000,500);
    h1_refs->GetXaxis()->SetNdivisions(120);
    h1_refs->GetXaxis()->SetRangeUser(-0.5,13.5);
    if(THIS_EXP=="2017_58Ni")h1_refs->GetXaxis()->SetRangeUser(-0.5,22.5);
    if(THIS_EXP=="2021_36Ar_SET3")h1_refs->GetXaxis()->SetRangeUser(-0.5,22.5);
    h1_refs->Draw("HIST TEXT0");
    
    Tfile_save_hist = new TFile(FILEPATH+"save_hist.root", "UPDATE");
    h1_refs->Write(THIS_EXP+"_h1_refs"); // 写入到 Tfile
    Tfile_save_hist->Close(); // Close 大写C
    h1_refs->SaveAs(FILEPATH+THIS_EXP+"_h1_refs.txt");

    //输出到文件看某些注入情况
    /*
    ofstream outfile_inj_has_one_ion;
    outfile_inj_has_one_ion.open(FILEPATH+"inj_ions_analysis.txt");
    for(int i=0;i<ions_n  ;i++)
    {
        if(Injection[ions[i].inject_number] ==2&&Injection_m[ions[i].inject_number] ==1 )
        {
            outfile_inj_has_one_ion<<ionspecies[ions[i].Species].Aname<<" "<<ions[i].inject_filename<<endl;
        }
    }
    outfile_inj_has_one_ion.close();
    */
}



// =====================  c4 磁刚度 周期 和 injection time 的关系 ======================
if(c4_Do_analysis_with_time_ON)
{

// 校准 连续的injection 和 时间
TGraph* gr_INJ_TIME = new TGraph();
AxisFormat(gr_INJ_TIME,"INJ_TIME"," time [s] ", " injection number");
for(int j=0;j<ions_n  ;j++)
{
    gr_INJ_TIME->SetPoint(gr_INJ_TIME->GetN(),ions[j].time,ions[j].inject_number );
}
TCanvas* c_INJ_TIME = new TCanvas("c_INJ_TIME","c_INJ_TIME",1800,600);
gr_INJ_TIME->SetMarkerSize(1);
/*
gr_INJ_TIME->GetXaxis()->SetTimeDisplay(1);
gr_INJ_TIME->GetXaxis()->SetNdivisions(515);
if(THIS_EXP=="2017_58Ni")gr_INJ_TIME->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2017-12-27 12:00:00");  // 58Ni
*/
gr_INJ_TIME->Draw("Ap");


//=========== 所有种类 T - time   ========== 
if(Generate_each_T_time)
{

    for(int i=0;i<NSpecies  ;i++)
    {
        ionspecies[i].grerr_T_t_n = 0;
        ionspecies[i].grerr_T_t  = new TGraphErrors();
        AxisFormat(ionspecies[i].grerr_T_t,ionspecies[i].Aname+" T_time"," time [s] ", "T[ns]");
        ionspecies[i].grerr_T_t->SetMarkerSize(1.0);
        ionspecies[i].grerr_T_t->SetLineWidth(1.0);
        ionspecies[i].grerr_T_t->GetYaxis()->SetTitleOffset(0.8);
        ionspecies[i].grerr_T_t->GetXaxis()->SetTimeDisplay(1);
        ionspecies[i].grerr_T_t->GetXaxis()->SetNdivisions(515);
        ionspecies[i].grerr_T_t->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2021-10-23 16:00:01");  // part-2 beginning
    }
    for(int j=0;j<ions_n  ;j++)
    {
        ionspecies[ions[j].Species].grerr_T_t->SetPoint(ionspecies[ions[j].Species].grerr_T_t_n , ions[j].time,ions[j].T );
        ionspecies[ions[j].Species].grerr_T_t->SetPointError(ionspecies[ions[j].Species].grerr_T_t_n , 0, ions[j].T_err);
        ionspecies[ions[j].Species].grerr_T_t_n++;
    }
    
    TCanvas*c_T_t_ionspecies[100];
    for(int i=0;i<NSpecies  ;i++)
    {
        c_T_t_ionspecies[i] = new TCanvas("c_T_t_ionspecies_"+ionspecies[i].Aname,"c_T_t_ionspecies"+ionspecies[i].Aname  ,2200,400);
        ionspecies[i].grerr_T_t->Draw("ap");
        //c_T_t_ionspecies[i]->Print(ionspecies[i].folder_path+"dA0_T.png");
        c_T_t_ionspecies[i]->Print(FILEPATH+ ionspecies[i].Aname+"_T_t.png");
    }

}//__________________________所有种类 T - time_________________________


//=========== 所有种类 mvq - time   ========== 
if(Generate_each_mvq_time)
{
    for(int i=0;i<NSpecies  ;i++)
    {
        ionspecies[i].grr_mvq_time  = new TGraphErrors();
        AxisFormat(ionspecies[i].grr_mvq_time,ionspecies[i].Aname+" T_time"," time [s] ", "mvq (u/e)");
        ionspecies[i].grr_mvq_time->SetMarkerSize(1.0);
        ionspecies[i].grr_mvq_time->SetLineWidth(1.0);
        ionspecies[i].grr_mvq_time->GetYaxis()->SetTitleOffset(0.8);
        ionspecies[i].grr_mvq_time->GetXaxis()->SetTimeDisplay(1);
        ionspecies[i].grr_mvq_time->GetXaxis()->SetLabelSize(0.02);
        ionspecies[i].grr_mvq_time->GetXaxis()->SetNdivisions(515);
        //if(THIS_EXP=="2021_36Ar_SET2")ionspecies[i].grr_mvq_time->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2021-10-24 00:00:01");  // set2  
        if(THIS_EXP=="2017_58Ni")ionspecies[i].grr_mvq_time->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2017-12-27 12:00:00");  // 58Ni 
    }
    for(int j=0;j<ions_n  ;j++)
    {
        if(ions[j].mvq_v1<1)continue;
        ionspecies[ions[j].Species].grr_mvq_time->SetPoint(ionspecies[ions[j].Species].grr_mvq_time->GetN() , ions[j].time,ions[j].mvq_v1 );
        //ionspecies[ions[j].Species].grr_mvq_time->SetPointError(ionspecies[ions[j].Species].grr_mvq_time->GetN()-1, 0, 0);
    }
    
    TCanvas*c_mvq_time_each[MAX_IONSPECIES];
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ionspecies[i].MassUnknown||ionspecies[i].N<100){continue;}
        c_mvq_time_each[i] = new TCanvas("c_mvq_time_each_"+ionspecies[i].Aname,"c_mvq_time_each"+ionspecies[i].Aname  ,2200,400);
        ionspecies[i].grr_mvq_time->Draw("ap");
        //c_mvq_time_each[i]->Print(ionspecies[i].folder_path+"dA0_T.png");
        c_mvq_time_each[i]->Print(FILEPATH+ ionspecies[i].Aname+"_mvq_time.png");
    }
}
//__________________________所有种类 mvq - time_________________________

//=========== 所有种类 Bp - time   ========== 
if(Generate_each_Bp_time)
{
    // 只有已知核能一开始算出 磁刚度
    for(int i=0;i<NSpecies  ;i++)
    {
        
        ionspecies[i].grr_Bp_time  = new TGraphErrors();
        AxisFormat(ionspecies[i].grr_Bp_time,ionspecies[i].Aname+" B#rho varying with time"," time [s] ", "B#rho [Tm]");
        //AxisFormat(ionspecies[i].grr_Bp_inj,ionspecies[i].Aname+" B#rho varying with injection"," injection number ", "B#rho [Tm]");
        ionspecies[i].grr_Bp_time->SetMarkerSize(1.0);
        ionspecies[i].grr_Bp_time->SetLineWidth(1.0);
        ionspecies[i].grr_Bp_time->GetYaxis()->SetTitleOffset(0.8);
        ionspecies[i].grr_Bp_time->GetXaxis()->SetTimeDisplay(1);
        ionspecies[i].grr_Bp_time->GetXaxis()->SetNdivisions(515);
        //if(THIS_EXP=="2021_36Ar_SET2")ionspecies[i].grr_Bp_time->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2021-10-23 16:00:01");  // 2021 36Ar 22Si  
        //if(THIS_EXP=="2021_36Ar_SET3")ionspecies[i].grr_Bp_time->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2021-11-06 00:00:01");  // 2021 36Ar 26P 
        if(THIS_EXP=="2017_58Ni")ionspecies[i].grr_Bp_time->GetXaxis()->SetTimeFormat("#frac{%m%d}{%H:%M:%S} %F2017-12-27 12:00:00");  // 58Ni 
    }
    for(int j=0;j<ions_n  ;j++)
    {
        if(ionspecies[ions[j].Species].IsRef)
        {
            ionspecies[ions[j].Species].grr_Bp_time->SetPoint(ionspecies[ions[j].Species].grr_Bp_time->GetN() , ions[j].time,ions[j].Bp );
            //ionspecies[ions[j].Species].grr_Bp_time->SetPointError(ionspecies[ions[j].Species].grr_Bp_time->GetN()-1 , 0, 0);
        }
        
    }
    
    TCanvas*c_Bp_t_ionspecies[100];
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ionspecies[i].IsRef&&ionspecies[i].N>3000)
        {
            c_Bp_t_ionspecies[i] = new TCanvas("c_Bp_t_ionspecies_"+ionspecies[i].Aname,"c_Bp_t_ionspecies"+ionspecies[i].Aname  ,2200,400);
            ionspecies[i].grr_Bp_time->Draw("ap");
            c_Bp_t_ionspecies[i]->Print(FILEPATH+ ionspecies[i].Aname+"_Bp_t.png");
        }
        
    }

}//__________________________所有种类 Bp - time_________________________

//=========== 所有种类 Bp - inj   ========== 
if(Generate_each_Bp_inj)
{
    // 只有已知核能一开始算出 磁刚度
    for(int i=0;i<NSpecies  ;i++)
    {
        
        ionspecies[i].grr_Bp_inj  = new TGraphErrors();
        AxisFormat(ionspecies[i].grr_Bp_inj,ionspecies[i].Aname+" B#rho varying with injection"," injection number ", "B#rho [Tm]");
        ionspecies[i].grr_Bp_inj->SetMarkerSize(1.0);
        ionspecies[i].grr_Bp_inj->SetLineWidth(1.0);
        ionspecies[i].grr_Bp_inj->GetYaxis()->SetTitleOffset(0.8);
          
    }
    for(int j=0;j<ions_n  ;j++)
    {
        if(ionspecies[ions[j].Species].IsRef)
        {
            ionspecies[ions[j].Species].grr_Bp_inj->SetPoint(ionspecies[ions[j].Species].grr_Bp_inj->GetN() , ions[j].inject_number,ions[j].Bp );
            //ionspecies[ions[j].Species].grr_Bp_inj->SetPointError(ionspecies[ions[j].Species].grr_Bp_inj->GetN()-1 , 0, 0);
        }
        
    }
    
    TCanvas*c_Bp_inj_ionspecies[100];
    for(int i=0;i<NSpecies  ;i++)
    {
        if(ionspecies[i].IsRef&&ionspecies[i].N>3000)
        {
            c_Bp_inj_ionspecies[i] = new TCanvas("c_Bp_inj_ionspecies_"+ionspecies[i].Aname,"c_Bp_inj_ionspecies"+ionspecies[i].Aname  ,2200,400);
            ionspecies[i].grr_Bp_inj->Draw("ap");
            c_Bp_inj_ionspecies[i]->Print(FILEPATH+ ionspecies[i].Aname+"_Bp_t.png");
        }
    }

}//__________________________所有种类 Bp - inj_________________________


}//__________________________ c4 磁刚度 周期 和 injection time 的关系 _________________________



if(c5_Do_each_T_fix_ON)
{
//===================================================================================================================
//=================================== C_T for each ionspecies and T_fix 20230908  ===================================
    cout<<endl<<" ========================= Do_each_T_fix_ON ===================== "<<endl;
    //Build_and_Draw_each_T_C(ionspecies);
    
    //Draw_each_T_histogram(ionspecies);
    
    int h_Tfix_bins_n = 250;  // 1000/bins = bin width
    int h_m0_bins_n  = 200;   // 4000/200 = 20 keV per bin 
    for(int i=0;i<NSpecies  ;i++)
    {
        ionspecies[i].h_Tfix0 = new TH1F(ionspecies[i].Aname+" h_Tfix0 ",ionspecies[i].Aname+" h_Tfix0 ",h_Tfix_bins_n,ionspecies[i].AveT-0.5,ionspecies[i].AveT+0.5);
        AxisFormat(ionspecies[i].h_Tfix0,ionspecies[i].Aname+" h_Tfix0 ","T[ns]","count",1);
        ionspecies[i].h_Tfix1 = new TH1F(ionspecies[i].Aname+" h_Tfix1 ",ionspecies[i].Aname+" h_Tfix1 ",h_Tfix_bins_n,ionspecies[i].AveT-0.5,ionspecies[i].AveT+0.5);
        AxisFormat(ionspecies[i].h_Tfix1,ionspecies[i].Aname+" h_Tfix1 ","T[ns]","count",kRed);
        ionspecies[i].h_Tfix2 = new TH1F(ionspecies[i].Aname+" h_Tfix2 ",ionspecies[i].Aname+" h_Tfix2 ",h_Tfix_bins_n,ionspecies[i].AveT-0.5,ionspecies[i].AveT+0.5);
        AxisFormat(ionspecies[i].h_Tfix2,ionspecies[i].Aname+" h_Tfix2 ","T[ns]","count",kRed);
        ionspecies[i].h_Tfix3 = new TH1F(ionspecies[i].Aname+" h_Tfix3 ",ionspecies[i].Aname+" h_Tfix3 ",h_Tfix_bins_n,ionspecies[i].AveT-0.5,ionspecies[i].AveT+0.5);
        AxisFormat(ionspecies[i].h_Tfix3,ionspecies[i].Aname+" h_Tfix3 ","T[ns]","count",kRed);

        ionspecies[i].h_m0_v1 = new TH1F(ionspecies[i].Aname+" h_m0_v1 ",ionspecies[i].Aname+" h_m0_v1 ",h_m0_bins_n,ionspecies[i].Mass_cal-2000,ionspecies[i].Mass_cal+2000);
        AxisFormat(ionspecies[i].h_m0_v1,ionspecies[i].Aname+" h_m0_v1 ","mass[keV]","count",1);
        ionspecies[i].h_m0_v2 = new TH1F(ionspecies[i].Aname+" h_m0_v2 ",ionspecies[i].Aname+" h_m0_v2 ",h_m0_bins_n,ionspecies[i].Mass_cal-2000,ionspecies[i].Mass_cal+2000);
        AxisFormat(ionspecies[i].h_m0_v2,ionspecies[i].Aname+" h_m0_v2 ","mass[keV]","count",kAzure);
        
        //v3 的算法中心值依赖于Bpfix Cfix,如果没给对，那么算出来的分布可能无法落在h_m0_v3 的范围内所以改用直接算std的方式，不完全依赖直方图
        ionspecies[i].h_m0_v3 = new TH1F(ionspecies[i].Aname+" h_m0_v3 ",ionspecies[i].Aname+" h_m0_v3 ",h_m0_bins_n,ionspecies[i].Mass_cal-2000,ionspecies[i].Mass_cal+2000);
        AxisFormat(ionspecies[i].h_m0_v3,ionspecies[i].Aname+" h_m0_v3 ","mass[keV]","count",kTeal);
        ionspecies[i].h_m1 = new TH1F(ionspecies[i].Aname+" h_m1 ",ionspecies[i].Aname+" h_m1 ",h_m0_bins_n,ionspecies[i].Mass_cal-2000,ionspecies[i].Mass_cal+2000);
        AxisFormat(ionspecies[i].h_m1,ionspecies[i].Aname+" h_m1 ","mass[keV]","count",kRed);
        
    }
    
    

    Build_gr_TC0(ionspecies);  // 必须有这一步创建 h_Tfix0 为了后续使用
    //Build_gr_TC1(ionspecies);  // gr_TC1 and h_Tfix1
    //Build_gr_TC2(ionspecies,gr_gtC_chosen);  // gr_TC2 and h_Tfix2
    
    Build_gr_TC3(ionspecies);
    //cout<<"debug"<<Get_Tfix2(ions[8888],128.808,gr_gtC_chosen);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas*c_TC0_ISS[MAX_IONSPECIES];
    for(int i=0;i<NSpecies  ;i++)
    {
        //if(ionspecies[i].Aname=="31S"&&THIS_EXP=="2017_58Ni")
        if(ionspecies[i].Aname=="11C"&&THIS_EXP=="2021_36Ar_SET3")
        //if( (ionspecies[i].Aname=="22Si"||ionspecies[i].Aname=="7Be")&&THIS_EXP=="2021_36Ar_SET2")
        {
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        c_TC0_ISS[i] = new TCanvas("c_TC0_"+ionspecies[i].Aname,"c_TC0_"+ionspecies[i].Aname  ,1000,500);
        ionspecies[i].gr_TC0->DrawClone("ap");
        //ionspecies[i].gr_TC1->DrawClone("psame");
        //ionspecies[i].gr_TC2->DrawClone("psame");
        ionspecies[i].gr_TC3->DrawClone("psame");
        }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas*c_T_ionspecies[MAX_IONSPECIES];

    for(int i=0;i<NSpecies  ;i++)
    {
        ionspecies[i].sigma_Tfix0 = ionspecies[i].h_Tfix0->GetStdDev();
        ionspecies[i].sigma_Tfix1 = ionspecies[i].h_Tfix1->GetStdDev();
        ionspecies[i].sigma_Tfix2 = ionspecies[i].h_Tfix2->GetStdDev();
        ionspecies[i].sigma_Tfix3 = ionspecies[i].h_Tfix3->GetStdDev();
        ionspecies[i].sigma_Tfix3_v3 = GetVectorStdDev(ionspecies[i].vector_Tfix3_v3, GetVectorMean(ionspecies[i].vector_Tfix3_v3));   //固定Bp C
        if(ionspecies[i].gr_TC3_n>10) 
        {   
            /*
            ionspecies[i].fitfun_h_Tfix3 = new TF1("fitfun_h_Tfix3_"+ionspecies[i].Aname,"gaus",500,700);
            ionspecies[i].h_Tfix3->Fit(ionspecies[i].fitfun_h_Tfix3,"qn");
            ionspecies[i].sigma_Tfix3 = ionspecies[i].fitfun_h_Tfix3->GetParameter(2);
            delete ionspecies[i].fitfun_h_Tfix3;
            */
        }
        cout<<ionspecies[i].Aname<<" sigma_Tfix0[ps]: "<<1000*ionspecies[i].sigma_Tfix0<<endl
            //<<ionspecies[i].Aname<<" sigma_Tfix1[ps]:"<<1000*ionspecies[i].sigma_Tfix1<<endl
            //<<ionspecies[i].Aname<<" sigma_Tfix2[ps]: "<<1000*ionspecies[i].sigma_Tfix2<<endl
            <<ionspecies[i].Aname<<" sigma_Tfix3[ps]: "<<1000*ionspecies[i].sigma_Tfix3<<endl
            <<ionspecies[i].Aname<<" sigma_Tfix3_v3[ps]: "<<1000*ionspecies[i].sigma_Tfix3_v3<<endl
            <<endl;
    
        //if(ionspecies[i].Aname=="31S"&&THIS_EXP=="2017_58Ni")
        if(ionspecies[i].Aname=="11C"&&THIS_EXP=="2021_36Ar_SET3")
        //if( (ionspecies[i].Aname=="22Si"||ionspecies[i].Aname=="7Be")&&THIS_EXP=="2021_36Ar_SET2")
        {
            
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        c_T_ionspecies[i] = new TCanvas("c_T_"+ionspecies[i].Aname,"c_T_"+ionspecies[i].Aname  ,1000,500);
        
        ionspecies[i].h_Tfix3->DrawClone("");
        //ionspecies[i].h_Tfix2->DrawClone("");
        //ionspecies[i].h_Tfix1->DrawClone("");
        ionspecies[i].h_Tfix0->DrawClone("same");
        
        }    
    }
    

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c_Isochronous = new TCanvas("c_Isochronous","c_Isochronous",1000,500);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TGraphErrors* grerr_sigmaTfix0 = new TGraphErrors();int grerr_sigmaTfix0_n=0;
    TGraphErrors* grerr_sigmaTfix1 = new TGraphErrors();int grerr_sigmaTfix1_n=0;
    TGraphErrors* grerr_sigmaTfix2 = new TGraphErrors();int grerr_sigmaTfix2_n=0;
    TGraphErrors* grerr_sigmaTfix3 = new TGraphErrors();int grerr_sigmaTfix3_n=0;
    TGraphErrors* grerr_sigmaTfix3_v3 = new TGraphErrors();int grerr_sigmaTfix3_v3_n=0;
    TGraphErrors* grr_sigmaTfix3_fit = new TGraphErrors();
    AxisFormat(grerr_sigmaTfix0,"std of Tfix0 for all ionspecies","T [ns]","#sigma_{T} [ps]",1);
    AxisFormat(grerr_sigmaTfix1,"std of Tfix1 for all ionspecies","T [ns]","#sigma_{Tfix1} [ps]",kRed);
    AxisFormat(grerr_sigmaTfix2,"std of Tfix2 for all ionspecies","T [ns]","#sigma_{Tfix2} [ps]",kAzure);
    AxisFormat(grerr_sigmaTfix3,"std of Tfix3 for all ionspecies","T [ns]","#sigma_{Tfix3} [ps]",kAzure+7);
    AxisFormat(grerr_sigmaTfix3_v3,"std of Tfix3_v3 for all ionspecies","T [ns]","#sigma_{Tfix3} [ps]",kTeal);
    AxisFormat(grr_sigmaTfix3_fit,"grr fit by pol2","T [ns]","#sigma_{Tfix3} [ps]",1);
    int N_plotmin = 10;
    int N_fit_min = 100;
    cout<<"----- the min number to be shown in grerr_sigmaTfix is N_plotmin = "<<N_plotmin<<endl;
    cout<<"----- the min number to be included in the grerr for fitting by pol2 is N_fit_min = "<<N_fit_min<<endl;
    ofstream OUTFILE_Isochronous_fix0;
    OUTFILE_Isochronous_fix0.open(FILEPATH+"Isochronous_fix0.txt");
    ofstream OUTFILE_Isochronous_fix3;
    OUTFILE_Isochronous_fix3.open(FILEPATH+"Isochronous_fix3.txt");
    for(int i=0;i<NSpecies  ;i++)
    {
        if(THIS_EXP=="2017_58Ni"&&ionspecies[i].Aname=="14O"){continue;}  // 58Ni 14O 存在问题
        if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O"){continue;}  // 36Ar 14O 未分辨 存在问题
        if(ionspecies[i].gr_TC0_n>N_plotmin&&ionspecies[i].Isomer_n==0)
        {
            grerr_sigmaTfix0->SetPoint(grerr_sigmaTfix0_n, ionspecies[i].AveT , ionspecies[i].sigma_Tfix0*1000);
            grerr_sigmaTfix0->SetPointError(grerr_sigmaTfix0_n, 0 , ionspecies[i].sigma_Tfix0*1000/(sqrt(2*ionspecies[i].gr_TC0_n -2) ) );
            grerr_sigmaTfix0_n++;
            grerr_sigmaTfix1->SetPoint(grerr_sigmaTfix1_n, ionspecies[i].AveT , ionspecies[i].sigma_Tfix1*1000);
            grerr_sigmaTfix1->SetPointError(grerr_sigmaTfix1_n, 0 , ionspecies[i].sigma_Tfix1*1000/(sqrt(2*ionspecies[i].gr_TC1_n -2) ) );
            grerr_sigmaTfix1_n++;
            grerr_sigmaTfix2->SetPoint(grerr_sigmaTfix2_n, ionspecies[i].AveT , ionspecies[i].sigma_Tfix2*1000);
            grerr_sigmaTfix2->SetPointError(grerr_sigmaTfix2_n, 0 , ionspecies[i].sigma_Tfix2*1000/(sqrt(2*ionspecies[i].gr_TC2_n -2) ) );
            grerr_sigmaTfix2_n++;
            OUTFILE_Isochronous_fix0<<ionspecies[i].Aname<<" "<<ionspecies[i].AveT<<" "<<ionspecies[i].Mvq_AME<<" "
            <<ionspecies[i].sigma_Tfix0*1000<<" "<<ionspecies[i].sigma_Tfix0*1000/(sqrt(2*ionspecies[i].gr_TC0_n -2) )<<endl;
        }    
        if(ionspecies[i].gr_TC0_n>N_plotmin&&ionspecies[i].Isomer_n==0)
        {
            
            grerr_sigmaTfix3->SetPoint(grerr_sigmaTfix3_n, ionspecies[i].AveT , ionspecies[i].sigma_Tfix3*1000);
            grerr_sigmaTfix3->SetPointError(grerr_sigmaTfix3_n, 0 , ionspecies[i].sigma_Tfix3*1000/(sqrt(2*ionspecies[i].gr_TC3_n -2) ) );
            grerr_sigmaTfix3_n++;
            grerr_sigmaTfix3_v3->SetPoint(grerr_sigmaTfix3_v3_n, ionspecies[i].AveT , ionspecies[i].sigma_Tfix3_v3*1000);
            grerr_sigmaTfix3_v3->SetPointError(grerr_sigmaTfix3_v3_n, 0 , ionspecies[i].sigma_Tfix3_v3*1000/(sqrt(2*ionspecies[i].gr_TC3_n -2) ) );
            grerr_sigmaTfix3_v3_n++;
            if(ionspecies[i].gr_TC3_n>N_fit_min)
            {
                grr_sigmaTfix3_fit->SetPoint(grr_sigmaTfix3_fit->GetN(), ionspecies[i].AveT , ionspecies[i].sigma_Tfix3*1000);
                grr_sigmaTfix3_fit->SetPointError(grr_sigmaTfix3_fit->GetN()-1, 0 , ionspecies[i].sigma_Tfix3*1000/(sqrt(2*ionspecies[i].gr_TC3_n -2) ) );
            }

            OUTFILE_Isochronous_fix3<<ionspecies[i].Aname<<" "<<ionspecies[i].AveT<<" "<<ionspecies[i].Mvq_AME<<" "
            <<ionspecies[i].sigma_Tfix3*1000<<" "<<ionspecies[i].sigma_Tfix3*1000/(sqrt(2*ionspecies[i].gr_TC3_n -2) )<<" "
            <<ionspecies[i].sigma_Tfix3_v3*1000<<" "<<ionspecies[i].sigma_Tfix3_v3*1000/(sqrt(2*ionspecies[i].gr_TC3_n -2) )
            <<endl;
        }
    }
    OUTFILE_Isochronous_fix0.close();
    OUTFILE_Isochronous_fix3.close();
    bool Draw_TFIX3=0;
    grerr_sigmaTfix0->Draw("ap");Draw_TFIX3=0;
    //grerr_sigmaTfix1->Draw("psame");
    //grerr_sigmaTfix2->Draw("psame");
    grerr_sigmaTfix3->Draw("psame");   Draw_TFIX3=0;//显示名字
    grerr_sigmaTfix3_v3->Draw("psame");
    grerr_sigmaTfix0->GetYaxis()->SetRangeUser(0,100);
    
    lat_n->SetTextAngle(90);
    lat_n->SetTextFont(43);
    lat_n->SetTextSize(20);
    for(int i=0;i<NSpecies  ;i++)
    {
        if(THIS_EXP=="2017_58Ni"&&ionspecies[i].Aname=="14O"){continue;}  // 58Ni 14O 存在问题
        if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O"){continue;}
        if(ionspecies[i].gr_TC0_n >N_plotmin&&ionspecies[i].Isomer_n==0)
        {
            lat_n->SetTextColor(1);
            lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].gr_TC0_n);
            lat_n->DrawLatex(ionspecies[i].AveT,ionspecies[i].sigma_Tfix0*1000,lat_text);
            if(Draw_TFIX3==1)
            {
                lat_n->SetTextColor(kAzure+7);
                lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].gr_TC3_n);
                lat_n->DrawLatex(ionspecies[i].AveT,ionspecies[i].sigma_Tfix3*1000,lat_text);
            }
        }   
    }
    //___________________________________ C_T for each ionspecies and Tfix 20230908 ___________________________________

    //Fit Isochronous curve:
    TF1* fitful_pol2_ISC = new TF1("fitful_pol2_ISC","pol2",550,700);
    fitful_pol2_ISC->SetParameters(2435,-7.702, 0.0061);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c_ISC_fit = new TCanvas("c_ISC_fit","c_ISC_fit",1000,500);
    grr_sigmaTfix3_fit ->Draw("ap");
    grerr_sigmaTfix3->Draw("psame");
    grr_sigmaTfix3_fit ->Fit(fitful_pol2_ISC);

    cout<<" chi square = "<<fitful_pol2_ISC->GetChisquare()<<endl;
    // TLatex 标名字
    for(int i=0;i<NSpecies  ;i++)
    {
        if(THIS_EXP=="2017_58Ni"&&ionspecies[i].Aname=="14O"){continue;}  // 58Ni 14O 存在问题
        if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O"){continue;}
        if(ionspecies[i].gr_TC0_n>N_plotmin&&ionspecies[i].Isomer_n==0)
        {
            if(ionspecies[i].gr_TC3_n >N_fit_min){lat_n->SetTextColor(1);}
            else {lat_n->SetTextColor(kAzure+7);}
            lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].gr_TC3_n);
            lat_n->DrawLatex(ionspecies[i].AveT,ionspecies[i].sigma_Tfix3*1000,lat_text);
            
        }   
    }



    //=========== 2025 0210 
    if(Do_Isochronous_mass_ON)
    {
        Build_h_m0(ionspecies);
        
        TCanvas*c_h_m0_ionspecies[MAX_IONSPECIES];

        for(int i=0;i<NSpecies  ;i++)
        {
            ionspecies[i].stddev_m0_v1 = ionspecies[i].h_m0_v1->GetStdDev();
            //ionspecies[i].stddev_m0_v2 = ionspecies[i].h_m0_v2->GetStdDev();
            
            //ionspecies[i].stddev_m0_v3 = ionspecies[i].h_m0_v3->GetStdDev();//选定的Bpfix Cfix 不一定能使结果落在直方图范围
            ionspecies[i].stddev_m0_v3 = GetVectorStdDev(ionspecies[i].vector_m0_v3, GetVectorMean(ionspecies[i].vector_m0_v3));
            
            cout<<ionspecies[i].Aname<<" stddev_m0_v1: "<<ionspecies[i].stddev_m0_v1<<" stddev_m0_v1/Z: "<<ionspecies[i].stddev_m0_v1/ionspecies[i].Z<<endl
                <<ionspecies[i].Aname<<" stddev_m0_v2: "<<ionspecies[i].stddev_m0_v2<<" stddev_m0_v2/Z: "<<ionspecies[i].stddev_m0_v2/ionspecies[i].Z<<endl
                <<ionspecies[i].Aname<<" stddev_m0_v3: "<<ionspecies[i].stddev_m0_v3<<" stddev_m0_v3/Z: "<<ionspecies[i].stddev_m0_v3/ionspecies[i].Z<<endl
                <<ionspecies[i].Aname<<" stddev_m1: "<<ionspecies[i].stdDeviation_V1<<" stddev_m1/Z: "<<ionspecies[i].stdDeviation_V1/ionspecies[i].Z<<endl
                <<endl;
            cout<<" m0_N= "<<ionspecies[i].m0_n<<"  N_unknown= "<<ionspecies[i].N_unknown<<endl;
        
            if(ionspecies[i].Aname=="31S"&&THIS_EXP=="2017_58Ni")
            //if(ionspecies[i].Aname=="11C"&&THIS_EXP=="2021_36Ar_SET3")
            //if( (ionspecies[i].Aname=="22Si"||ionspecies[i].Aname=="7Be")&&THIS_EXP=="2021_36Ar_SET2")
            {
                
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            c_h_m0_ionspecies[i] = new TCanvas("c_h_m0_"+ionspecies[i].Aname,"c_h_m0_"+ionspecies[i].Aname  ,1000,500);
            
            ionspecies[i].h_m0_v1->DrawClone("");
            //ionspecies[i].h_m0_v2->DrawClone("same");
            ionspecies[i].h_m0_v3->DrawClone("same");
            ionspecies[i].h_m1->DrawClone("same");

            
            }    
        }
        //=================20250210
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TCanvas *c_Isochronous_m = new TCanvas("c_Isochronous_m","c_Isochronous_m",1000,500);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TGraphErrors* grr_stddev_m1 = new TGraphErrors();int grr_stddev_m1_n=0;
        TGraphErrors* grr_stddev_m0_v1 = new TGraphErrors();int grr_stddev_m0_v1_n=0;
        //TGraphErrors* grr_stddev_m0_v2 = new TGraphErrors();int grr_stddev_m0_v2_n=0;
        TGraphErrors* grr_stddev_m0_v3 = new TGraphErrors();int grr_stddev_m0_v3_n=0;
        TGraphErrors* grr_stddev_m1_fit = new TGraphErrors();
        AxisFormat(grr_stddev_m1,"standard deviation of mass for all ionspecies","T [ns]","#sigma_{m1} [keV]",kRed);
        AxisFormat(grr_stddev_m0_v1,"standard deviation of mass for all ionspecies","T [ns]","#sigma_{m0_v1} [keV]",kBlack);
        //AxisFormat(grr_stddev_m0_v2,"standard deviation of mass for all ionspecies","T [ns]","#sigma_{m0_v2} [keV]",kBlue);
        AxisFormat(grr_stddev_m0_v3,"standard deviation of mass for all ionspecies","T [ns]","#sigma_{m0_v3} [keV]",kTeal);
        grr_stddev_m0_v1->SetMarkerStyle(24);
        //grr_stddev_m0_v2->SetMarkerStyle(25);
        grr_stddev_m0_v3->SetMarkerStyle(26);
        AxisFormat(grr_stddev_m1_fit,"grr fit by pol2","T [ns]","#sigma_{Tfix3} [ps]",1);
        N_plotmin = 10;
        N_fit_min = 100;
        cout<<"----- the min number to be shown in grr_stddev_m1 is N_plotmin = "<<N_plotmin<<endl;
        cout<<"----- the min number to be included in the grerr for fitting by pol2 is N_fit_min = "<<N_fit_min<<endl;
        ofstream OUTFILE_Isochronous_stddev_m;
        OUTFILE_Isochronous_stddev_m.open(FILEPATH+"Isochronous_stddev_m.txt");
        
        for(int i=0;i<NSpecies  ;i++)
        {
            if(THIS_EXP=="2017_58Ni"&&ionspecies[i].Aname=="14O"){continue;}  // 58Ni 14O 存在问题
            if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O"){continue;}  // 36Ar 14O 未分辨 存在问题
            if(ionspecies[i].N_unknown>N_plotmin&&ionspecies[i].Isomer_n==0)
            {
                grr_stddev_m1->SetPoint(grr_stddev_m1_n, ionspecies[i].AveT , ionspecies[i].stdDeviation_V1/ionspecies[i].Z);
                grr_stddev_m1->SetPointError(grr_stddev_m1_n, 0 , ionspecies[i].stdDeviation_V1/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) ) );
                grr_stddev_m1_n++;

                grr_stddev_m0_v1->SetPoint(grr_stddev_m0_v1_n, ionspecies[i].AveT , ionspecies[i].stddev_m0_v1/ionspecies[i].Z);
                grr_stddev_m0_v1->SetPointError(grr_stddev_m0_v1_n, 0 , ionspecies[i].stddev_m0_v1/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) ) );
                grr_stddev_m0_v1_n++;
                //grr_stddev_m0_v2->SetPoint(grr_stddev_m0_v2_n, ionspecies[i].AveT , ionspecies[i].stddev_m0_v2/ionspecies[i].Z);
                //grr_stddev_m0_v2->SetPointError(grr_stddev_m0_v2_n, 0 , ionspecies[i].stddev_m0_v2/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) ) );
                //grr_stddev_m0_v2_n++;
                grr_stddev_m0_v3->SetPoint(grr_stddev_m0_v3_n, ionspecies[i].AveT , ionspecies[i].stddev_m0_v3/ionspecies[i].Z);
                grr_stddev_m0_v3->SetPointError(grr_stddev_m0_v3_n, 0 , ionspecies[i].stddev_m0_v3/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) ) );
                grr_stddev_m0_v3_n++;
                if(ionspecies[i].N_unknown>N_fit_min)
                {
                    grr_stddev_m1_fit->SetPoint(grr_stddev_m1_fit->GetN(), ionspecies[i].AveT , ionspecies[i].stdDeviation_V1/ionspecies[i].Z);
                    grr_stddev_m1_fit->SetPointError(grr_stddev_m1_fit->GetN()-1, 0 , ionspecies[i].stdDeviation_V1/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) ) );
                }
                OUTFILE_Isochronous_stddev_m<<ionspecies[i].Aname<<" "<<ionspecies[i].AveT<<" "<<ionspecies[i].Mvq_AME<<" "
                <<ionspecies[i].stddev_m0_v1/ionspecies[i].Z<<" "<<ionspecies[i].stddev_m0_v1/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) )<<" "
                <<ionspecies[i].stddev_m0_v3/ionspecies[i].Z<<" "<<ionspecies[i].stddev_m0_v3/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) )<<" "
                <<ionspecies[i].stdDeviation_V1/ionspecies[i].Z<<" "<<ionspecies[i].stdDeviation_V1/ionspecies[i].Z/(sqrt(2*ionspecies[i].m0_n-2) )<<" "
                <<ionspecies[i].m0_n<<" "
                <<endl;
            }    
            
        }
        OUTFILE_Isochronous_stddev_m.close();
        
        bool Draw_steddev_m1=0;
        
        grr_stddev_m0_v1->Draw("ap");Draw_steddev_m1=0;
        //grr_stddev_m0_v2->Draw("psame");
        grr_stddev_m0_v3->Draw("psame");
        grr_stddev_m1->Draw("psame");Draw_steddev_m1=1;
        grr_stddev_m0_v1->GetYaxis()->SetRangeUser(0,100);
        
        lat_n->SetTextAngle(90);
        lat_n->SetTextFont(43);
        lat_n->SetTextSize(20);
        /*
        for(int i=0;i<NSpecies  ;i++)
        {
            if(THIS_EXP=="2017_58Ni"&&ionspecies[i].Aname=="14O"){continue;}  // 58Ni 14O 存在问题
            if(THIS_EXP=="2021_36Ar_SET3"&&ionspecies[i].Aname=="14O"){continue;}
            if(ionspecies[i].m0_n>N_plotmin&&ionspecies[i].Isomer_n==0)
            {
                lat_n->SetTextColor(1);
                lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].m0_n);
                lat_n->DrawLatex(ionspecies[i].AveT, ionspecies[i].stddev_m0_v1/ionspecies[i].Z, lat_text);
                if(Draw_steddev_m1==1)
                {
                    lat_n->SetTextColor(kAzure+7);
                    lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].m0_n);
                    lat_n->DrawLatex(ionspecies[i].AveT, ionspecies[i].stdDeviation_V1/ionspecies[i].Z, lat_text);
                }
            }   
        }
        */
        


    } //if(Do_Isochronous_mass_ON)

}// end of Do Tfix
//______________________________________________________________________________________________________________________________


if(c6_CdC_ON)
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c6_CdC = new TCanvas("c6_CdC","c6_CdC"  ,1000,500);
    
    h_dC_ref_target->Draw();

}

//==================================================================================================================
//----------------------------
if(Show_mvqC_each_ON)
{
    cout<<" \033[32m ---------- Show_mvqC_each_ON ------------- \033[0m"<<endl;
    Show_mvqC_each(ionspecies);
}
if(Show_dmC_each_ON)
{
    cout<<" \033[32m ---------- Show_dmC_each_ON ------------- \033[0m"<<endl;
    Show_dmC_each(ionspecies);
}
if(DoShow_mvqC_each_h2_ON)Show_mvqC_each_h2(ionspecies);
if(Show_h_iont_merr_each_ON)
{
    cout<<" \033[32m ---------- Show_h_iont_merr_each_ON ------------- \033[0m"<<endl;
    Show_h_iont_merr_each(ionspecies);
}
if(Show_h_refcal_merr_each_ON)
{
    cout<<" \033[32m ---------- Show_h_refcal_merr_each_ON ------------- \033[0m"<<endl;
    Show_h_refcal_merr_each(ionspecies);
}
if(Show_h_iont_refs_chi_each_ON)
{
    cout<<" \033[32m ---------- Show_h_iont_refs_chi_each ------------- \033[0m"<<endl;
    Show_h_iont_refs_chi_each(ionspecies);
}
//-----------------------------

//===============================================================
if(c_avegt_ionspecies_ON)
{
//------------0930
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_avegt_ionspecies = new TCanvas("c_avegt_ionspecies","avegt - ionspecies "  ,1000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int grerr_avegt_A_n = 0;
TGraphErrors* grerr_avegt_A = new TGraphErrors();
AxisFormat(grerr_avegt_A,"grerr_avegt_A","Tz","ave #gamma_{t}");
for(int i=0;i<NSpecies  ;i++)
{
    if(ionspecies[i].N>10)
    {
        grerr_avegt_A->SetPoint(grerr_avegt_A_n,ionspecies[i].Tz,ionspecies[i].Avegammat);
        grerr_avegt_A->SetPointError(grerr_avegt_A_n, 0, ionspecies[i].Sigmagammat);
        grerr_avegt_A_n++;
    }
}
grerr_avegt_A->Draw("ap");
}
//___________________________________________________________________


if(c7_ON)
{//=====================================================

cout<<" ================ c7 ON ======================="<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c7 = new TCanvas("c7","c7 h_T_all",2000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double T_bin_width = 10; //[ps]
TH1F* h_T_all = new TH1F("h_T_all","h_T_all",(T_RANGE_MAX- T_RANGE_MIN)/T_bin_width*1000,T_RANGE_MIN,T_RANGE_MAX); 
AxisFormat(h_T_all,"h_T_all"," Revolution Time T (ns) ", strtmp.Format(" Counts (per %.1f ps) ",T_bin_width));
h_T_all->GetYaxis()->SetTitleOffset(0.4);
h_T_all->GetYaxis()->SetTickLength(0.015);

TH1F** h_T_all_Tz = new TH1F* [All_Tz_n];
// note for 2017_58Ni: Tz[0] = -2, Tz[1] = -1.5 ...
for(int j=0;j<All_Tz_n  ;j++)
{
    h_T_all_Tz[j] = new TH1F ("h_T_all_Tz="+All_Tz_str[j],"h_T_all_Tz="+All_Tz_str[j], (T_RANGE_MAX- T_RANGE_MIN)/T_bin_width*1000,T_RANGE_MIN,T_RANGE_MAX);
    AxisFormat(h_T_all_Tz[j],""," Revolution Time T (ns) ", strtmp.Format(" Counts (per %.1f ps) ",T_bin_width));
    h_T_all_Tz[j]->SetLineWidth(1);
    h_T_all_Tz[j]->SetLineColor(GetTzColor_v2(All_Tz[j]));
    
}

for(int i=0;i<ions_n  ;i++)
{
    if(ions[i].T>=T_RANGE_MIN&&ions[i].T<T_RANGE_MAX&&ionspecies[ions[i].Species].N>0)
    {   
        h_T_all->Fill(ions[i].T);
        for(int j=0;j<All_Tz_n  ;j++)
        {
            if(ionspecies[ions[i].Species].Tz == All_Tz[j])
            {
                h_T_all_Tz[j]->Fill(ions[i].T);
                break;
            }
        }
    }    
}

//##############################
//h_T_all->DrawClone();
h_T_all_Tz[All_Tz_n-1]->DrawClone();
h_T_all_Tz[All_Tz_n-1]->GetYaxis()->SetRangeUser(0.5,20000);
h_T_all_Tz[All_Tz_n-1]->GetXaxis()->SetRangeUser(610,644);
h_T_all_Tz[All_Tz_n-1]->GetYaxis()->SetTickLength(0.02);
h_T_all_Tz[All_Tz_n-1]->GetYaxis()->SetTitleOffset(0.5);
for(int j=All_Tz_n-2;j>=1  ;j--)  // 最后
{
    h_T_all_Tz[j]->DrawClone("same");
    //h_T_all_Tz[j]->GetYaxis()->SetRangeUser(0.5,7000);
}
////////////////////////////////////
//c7->SetGridy(1);
c7->SetLogy(1);
//在峰上标名字
lat_n->SetTextColor(1);
lat_n->SetTextFont(43);
lat_n->SetTextSize(35);
lat_n->SetTextAlign(12);
lat_n->SetTextAngle(90);
for(int i=0;i<NSpecies;i++)
{
    if(ionspecies[i].AveT>=T_RANGE_MIN&&ionspecies[i].AveT<T_RANGE_MAX&&ionspecies[i].N>0)
    {
        //lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].N);    
        lat_n->SetTextColor(GetTzColor_v2( ionspecies[i].Tz));
        lat_text = ionspecies[i].name_latex;  
        if(ionspecies[i].Tz>-2)lat_n->DrawLatex(ionspecies[i].AveT,ionspecies[i].N,lat_text);
    }    
}
//delete h_T_all;
}//____________________________________ c7 ________________________________________

if(c7_b_ON)
{//======================= mvq spectrum ==============================

cout<<" ================ c7_b mvq ON ======================="<<endl;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c7_b = new TCanvas("c7_b","c7_b h_mvq_all",2000,500);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double mvq_bin_width = 0.0001; //[ps]
TH1F* h_mvq_all = new TH1F("h_mvq_all","h_mvq_all",(mvq_RANGE_MAX- mvq_RANGE_MIN)/mvq_bin_width,mvq_RANGE_MIN,mvq_RANGE_MAX); 
AxisFormat(h_mvq_all,"h_mvq_all"," m/q  (u/e) ", " Counts (per 10^{-4} u/e) ");
h_mvq_all->GetYaxis()->SetTitleOffset(0.4);
h_mvq_all->GetYaxis()->SetTickLength(0.015);

TH1F** h_mvq_all_Tz = new TH1F* [All_Tz_n];
// note for 2017_58Ni: Tz[0] = -2, Tz[1] = -1.5 ...
for(int j=0;j<All_Tz_n  ;j++)
{
    h_mvq_all_Tz[j] = new TH1F ("h_mvq_all_Tz="+All_Tz_str[j],"h_mvq_all_Tz="+All_Tz_str[j], (mvq_RANGE_MAX- mvq_RANGE_MIN)/mvq_bin_width,mvq_RANGE_MIN,mvq_RANGE_MAX);
    AxisFormat(h_mvq_all_Tz[j],""," m/q  (u/e) ", " Counts (per 10^{-4} u/e) ");
    h_mvq_all_Tz[j]->SetLineWidth(1);
    h_mvq_all_Tz[j]->SetLineColor(GetTzColor_v2(All_Tz[j]));
    
}

for(int i=0;i<ions_n  ;i++)
{
    if(ions[i].mvq_v1>=mvq_RANGE_MIN&&ions[i].mvq_v1<mvq_RANGE_MAX&&ionspecies[ions[i].Species].N_unknown>0)
    {   
        h_mvq_all->Fill(ions[i].mvq_v1);
        for(int j=0;j<All_Tz_n  ;j++)
        {
            if(ionspecies[ions[i].Species].Tz == All_Tz[j])
            {
                h_mvq_all_Tz[j]->Fill(ions[i].mvq_v1);
                break;
            }
        }
    }    
}

//##############################
h_mvq_all->DrawClone();
h_mvq_all->GetYaxis()->SetRangeUser(0.5,7000);
h_mvq_all->GetYaxis()->SetTitleOffset(0.7);
/*
h_mvq_all_Tz[All_Tz_n-1]->DrawClone();
h_mvq_all_Tz[All_Tz_n-1]->GetYaxis()->SetRangeUser(0.5,10000);
//h_mvq_all_Tz[All_Tz_n-1]->GetXaxis()->SetRangeUser(610,644);
h_mvq_all_Tz[All_Tz_n-1]->GetYaxis()->SetTickLength(0.02);
h_mvq_all_Tz[All_Tz_n-1]->GetYaxis()->SetTitleOffset(0.5);
for(int j=All_Tz_n-2;j>=1  ;j--)  // 最后
{
    h_mvq_all_Tz[j]->DrawClone("same");
    //h_mvq_all_Tz[j]->GetYaxis()->SetRangeUser(0.5,7000);
}
*/
////////////////////////////////////
//c7_b->SetGridy(1);
c7_b->SetLogy(1);
//在峰上标名字
lat_n->SetTextColor(1);
lat_n->SetTextFont(43);
lat_n->SetTextSize(35);
lat_n->SetTextAlign(12);
lat_n->SetTextAngle(90);
for(int i=0;i<NSpecies;i++)
{
    if(ionspecies[i].Mvq_cal>=mvq_RANGE_MIN&&ionspecies[i].Mvq_cal<mvq_RANGE_MAX&&ionspecies[i].N_unknown>0)
    {
        //lat_text = ionspecies[i].name_latex+ strtmp.Format(" : %d ",ionspecies[i].N);    
        lat_text = ionspecies[i].name_latex;


        //lat_n->SetTextColor(GetTzColor_v2( ionspecies[i].Tz));
        if(ionspecies[i].MassUnknown)lat_n->SetTextColor(kRed);
        else
        {
            if(ionspecies[i].IsRef)lat_n->SetTextColor(1);
            else lat_n->SetTextColor(kBlue);
        }
          
        //if(ionspecies[i].Tz>-2)lat_n->DrawLatex(ionspecies[i].Mvq_cal,ionspecies[i].N_unknown,lat_text);

        lat_n->DrawLatex(ionspecies[i].Mvq_cal,ionspecies[i].N_unknown,lat_text);

    }    
}

}//____________________________________ c7_b ________________________________________

// ==================================== c8 =======================================
if(c8_ON)
{

TGraphErrors* grr_gt_mass = new TGraphErrors();
TH2F* h2_gt_mass_10C = new TH2F("h2_gt_mass_10C","h2_gt_mass_10C",100,1.32,1.4, 100,1.6687,1.6691);

TGraphErrors* grr_gt_mass_12N = new TGraphErrors();
AxisFormat(grr_gt_mass_12N,"grr_gt_mass_12N","#gamma_{t}"," mvq_v1");
grr_gt_mass_12N->SetMarkerSize(1.0);
grr_gt_mass_12N->SetMarkerStyle(24);
TH2F* h2_gt_mass_12N = new TH2F("h2_gt_mass_12N","h2_gt_mass_12N",100,1.31,1.37, 100,1.71634,1.71644);

TGraphErrors* grr_gt_mass_20Mg = new TGraphErrors();
AxisFormat(grr_gt_mass_20Mg,"grr_gt_mass_20Mg","#gamma_{t}"," mvq_v1");
grr_gt_mass_20Mg->SetMarkerSize(1.0);
grr_gt_mass_20Mg->SetMarkerStyle(24);
TH2F* h2_gt_mass_20Mg = new TH2F("h2_gt_mass_20Mg","h2_gt_mass_20Mg",100,1.31,1.37, 100,1.6676,1.6678);

if(THIS_EXP=="2021_36Ar_SET3")
{
    
    for(int j=0;j<ions_n  ;j++)
    {
        if( ionspecies[ ions[j].Species ].Aname =="10C" )
        {
            //if(ions[j].mvq_v1>1&&ions[j].gammat>1.3&&ions[j].gammat<1.4&&ions[j].gammat_err<0.05)
            if(ions[j].mvq_v1>1)
            {
                grr_gt_mass->SetPoint(grr_gt_mass->GetN(),ions[j].gammat, ions[j].mvq_v1);
                //grr_gt_mass->SetPointError(grr_gt_mass->GetN()-1,ions[j].gammat_err, 0);
                grr_gt_mass->SetPointError(grr_gt_mass->GetN()-1,0, 0);

                h2_gt_mass_10C->Fill(ions[j].gammat, ions[j].mvq_v1);
            }
            
        }
        if( ionspecies[ ions[j].Species ].Aname =="12N" )
        {
            if(ions[j].mvq_v1>1)
            {
                grr_gt_mass_12N->SetPoint(grr_gt_mass_12N->GetN(),ions[j].gammat, ions[j].mvq_v1);
                h2_gt_mass_12N->Fill(ions[j].gammat, ions[j].mvq_v1);
            }
        }
        if( ionspecies[ ions[j].Species ].Aname =="20Mg" )
        {
            if(ions[j].mvq_v1>1)
            {
                grr_gt_mass_20Mg->SetPoint(grr_gt_mass_20Mg->GetN(),ions[j].gammat, ions[j].mvq_v1);
                h2_gt_mass_20Mg->Fill(ions[j].gammat, ions[j].mvq_v1);
            }
        }
    }
    AxisFormat(grr_gt_mass,"grr_gt_mass-10C","#gamma_{t}"," mvq_v1");
    grr_gt_mass->SetMarkerSize(1.0);
    grr_gt_mass->SetMarkerStyle(24);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c8_gt_mass = new TCanvas("c8_gt_mass","c8_gt_mass",800,1200);
c8_gt_mass->Divide(1,2);
c8_gt_mass->cd(1);
grr_gt_mass->Draw("Ap");
for(int i=0;i<NSpecies;i++){ if(ionspecies[i].Aname=="10C"){ionspecies[i].f0_mvq_AME->Draw("same");}}

c8_gt_mass->cd(2);
h2_gt_mass_10C->Draw("colz");
for(int i=0;i<NSpecies;i++){ if(ionspecies[i].Aname=="10C"){ionspecies[i].f0_mvq_AME->Draw("same");}}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c8_gt_mass_12N = new TCanvas("c8_gt_mass_12N","c8_gt_mass_12N",800,1200);
c8_gt_mass_12N->Divide(1,2);
c8_gt_mass_12N->cd(1);
grr_gt_mass_12N->Draw("Ap");
for(int i=0;i<NSpecies;i++){ if(ionspecies[i].Aname=="12N"){ionspecies[i].f0_mvq_AME->Draw("same");}}

c8_gt_mass_12N->cd(2);
h2_gt_mass_12N->Draw("colz");
for(int i=0;i<NSpecies;i++){ if(ionspecies[i].Aname=="12N"){ionspecies[i].f0_mvq_AME->Draw("same");}}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c8_gt_mass_20Mg = new TCanvas("c8_gt_mass_20Mg","c8_gt_mass_20Mg",800,1200);
c8_gt_mass_20Mg->Divide(1,2);
c8_gt_mass_20Mg->cd(1);
grr_gt_mass_20Mg->Draw("Ap");
for(int i=0;i<NSpecies;i++){ if(ionspecies[i].Aname=="20Mg"){ionspecies[i].f0_mvq_AME->Draw("same");}}

c8_gt_mass_20Mg->cd(2);
h2_gt_mass_20Mg->Draw("colz");
for(int i=0;i<NSpecies;i++){ if(ionspecies[i].Aname=="20Mg"){ionspecies[i].f0_mvq_AME->Draw("same");}}


}
// ==================================== c8 =======================================



if(c9_h_Mvq_ON)
{//=====================================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c9_h_Mvq = new TCanvas("c9_h_Mvq","c9_h_Mvq Mvq",2000,1000);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool c9_show_count_ON = 0;

AxisFormat(h_Mvq,"Mvq Calculated","m/q","Counts",2);
h_Mvq->SetLineColor(kAzure);
h_Mvq->SetLineWidth(4);
h_Mvq->GetXaxis()->SetTitleOffset(1.0);
h_Mvq->GetXaxis()->SetTitleSize(0.05);
h_Mvq->GetXaxis()->SetLabelSize(0.04);
h_Mvq->GetYaxis()->SetTitleOffset(0.7);
h_Mvq->GetYaxis()->SetTitleSize(0.05);
h_Mvq->GetYaxis()->SetLabelSize(0.04);
c9_h_Mvq->SetLogy(1);

h_Mvq->Draw();
//在峰上标名字
lat_n->SetTextColor(1);
lat_n->SetTextFont(43);
lat_n->SetTextSize(35);
lat_n->SetTextAngle(90);   // 竖着写名字

cout<<" c9_h_Mvq ON---------------------------"<<endl;
for(int i = 0; i<NSpecies; i++)
{
    if(c9_show_count_ON)lat_text = ionspecies[i].name_latex + strtmp.Format(" : %d",ionspecies[i].N_unknown);
    else                lat_text = ionspecies[i].name_latex ;
    
    lat_n->DrawLatex(  ionspecies[i].Mvq_cal, 
                       (ionspecies[i].N_unknown>1000)?1000:ionspecies[i].N_unknown,
                        lat_text);
    cout<<ionspecies[i].A<<" "<<ionspecies[i].name<<" "
    <<fixed<<setprecision(8)
    <<" mvq_cal_v1(average) = "<<ionspecies[i].Mvq_cal<<" AME_mvq= "<<ionspecies[i].Mvq_AME<<endl;   
}
lat_n->SetTextAngle(0);    //reset

}//__________________________c9_h_Mvq______________________________


//=========================c10=========================================
if(c10_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c10 = new TCanvas("c10","c10 Fit BpC for all ionspecies",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(int i=0;i<NSpecies  ;i++)
{
    if(ionspecies[i].N>10)
    {
        ionspecies[i].Fit_BpC(0);  //0 BpC 1 lnBpC
    }
}
ionspecies[0].fitfun_pol1_BpC->DrawClone("");
for(int i=0;i<NSpecies  ;i++)
{
    ionspecies[i].fitfun_pol1_BpC->DrawClone("same");
}
auto legend_c10 = new TLegend(0.80,0.10,0.95,0.35); //downright
     legend_c10->SetHeader(strtmp.Format("L = %.3f,ddt = %.3f",L,ddT),"C");
for(int i=0;i<NSpecies  ;i++)
{
    legend_c10->AddEntry(ionspecies[i].fitfun_pol1_BpC,ionspecies[i].Aname,"l" );
}
legend_c10->Draw("same");

double C_tmp;
//Generate_BpC_dispersion("c_BpC_dispersion",ionspecies,NSpecies, 0,false,false,C_tmp);


}//_________________________c10____________________________________

//=========================c10=========================================
if(c11_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c11 = new TCanvas("c11","c11 k_lnBpC for all ionspecies",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int gr_k_lnBpC_n =0;
TGraph* gr_k_lnBpC = new TGraph();
AxisFormat(gr_k_lnBpC,"k_lnBpC for all ionspecies","AveT (ns)","k_lnBpC");
int gr_k_BpC_n =0;
TGraph* gr_k_BpC = new TGraph();
AxisFormat(gr_k_BpC,"k_BpC for all ionspecies","AveT (ns)","k_BpC");
for(int i=0;i<NSpecies;i++)
{
    if(ionspecies[i].N>30)
    {
       gr_k_BpC->SetPoint(gr_k_lnBpC_n++,ionspecies[i].AveT*0.001,ionspecies[i].Get_k_BpC());
       //gr_k_lnBpC->SetPoint(gr_k_lnBpC_n++,ionspecies[i].AveT*0.001,ionspecies[i].Get_k_lnBpC());
    }
    
}
gr_k_BpC->Draw("ap");
//gr_k_lnBpC->Draw("ap");

}//_________________________c11____________________________________

//===================================   mass deviation analysis =============================

if(c_gr_dmdC_ON)
{//=========================c_dmdc=========================================
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_dmdc = new TCanvas("c_dmdc","c_dmdc dmdC",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TF1* fit_dmdc = new TF1("fit_dmdc","pol1",0,100);
TLatex* lat_n_dmdc = new TLatex();
cout<<"dmdC = "<<gr_dmdC->GetN()<<endl;
AxisFormat(gr_dmdC,"dmdC for all mass calculations","#DeltaC = C_x - C_ref [m]","m_{cal}-m_{AME} [keV]");
gr_dmdC->SetMarkerSize(1.0);
gr_dmdC->SetMarkerStyle(24);
gr_dmdC->GetYaxis()->SetRangeUser(-3000,3000);
gr_dmdC->Fit(fit_dmdc);
gr_dmdC->Draw("ap");

lat_text = Info_fitfun_pol1(fit_dmdc);
lat_n_dmdc->SetTextFont(22);
lat_n_dmdc->SetTextSize(0.06);
lat_n_dmdc->SetTextColor(2);
lat_n_dmdc->DrawLatex(0,-2000,lat_text);

if(MASS_VER>=3)
{
    TCanvas *c_dmC_VE_all = new TCanvas("c_dmC_VE_all","c_dmC_VE_all",1200,600);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cout<<"grerr_dmC_VE_all_n = "<<grerr_dmC_VE_all_n<<endl;
    AxisFormat(grerr_dmC_VE_all,"C dm for each ion_t mass calculations","C [m]","m_{calv2}-m_{AME} [keV]",4);
    
    grerr_dmC_VE_all->SetMarkerSize(1.2);
    grerr_dmC_VE_all->GetYaxis()->SetRangeUser(-3000,3000);
    grerr_dmC_VE_all->Draw("ap");
}

}//_________________________c_dmdc____________________________________
if(c_h2_dmdC_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_h2_dmdc = new TCanvas("c_h2_dmdc","c_h2_dmdc",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
AxisFormat(h2_dmdC,"dmdC for all mass calculations","#DeltaC = C_x - C_ref [m]","m_{cal}-m_{AME} [keV]");
h2_dmdC->Draw("colz");
}

if(c_dmISO_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_dmISO = new TCanvas("c_dmISO","c_dmISO",1200,1200);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TF1* fit_dmISO = new TF1("fit_dmISO","pol1",0,100);
TLatex* lat_n_dmISO = new TLatex();
cout<<"dmISO = "<<gr_dmISO->GetN()<<endl;
// gr_dmISO --- 实际是mvq
AxisFormat(gr_dmISO,"dmISO for all mass calculations"," #gamma^{2} - #gamma_{t}^{2}(iont.C)","iont.mvq_{cal}-mvq_{AME} [keV]");
AxisFormat(h2_dmvqISO,"h2_dmvqISO for all mass calculations"," #gamma^{2} - #gamma_{t}^{2}(iont.C)","iont.mvq_{cal}-mvq_{AME} [keV]");
gr_dmISO->SetMarkerSize(1.0);
gr_dmISO->SetMarkerStyle(24);
//gr_dmISO->GetYaxis()->SetRangeUser(-3000,3000);
//gr_dmISO->Fit(fit_dmdc);

//gr_dmISO->Draw("ap");

c_dmISO->Divide(1,2);
c_dmISO->cd(1);
h2_dmvqISO->Draw("colz");

TH1F* h1_dmvqISO_StdDeV = new TH1F("h1_dmvqISO_StdDeV","h1_dmvqISO_StdDeV", 400,-0.1,0.2);  // 用于分区间 和 h2_dmvqISO 横轴保持一致
TGraphErrors* grr_dmvqStdDev_ISO = new TGraphErrors();
AxisFormat(grr_dmvqStdDev_ISO,"grr_dmvqStdDev_ISO","#gamma^{2} - #gamma_{t}^{2}(iont.C)","StdDeV of #Delta (iont.mvq)");
Get_subregion_StdDev_from_gr(gr_dmISO, h1_dmvqISO_StdDeV,400, grr_dmvqStdDev_ISO);
//grr_dmvqStdDev_ISO->Print();
c_dmISO->cd(2);
grr_dmvqStdDev_ISO->Draw("Ap");
//================= subregionn StdDev ==================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_errISO = new TCanvas("c_errISO","c_errISO",1200,1200);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AxisFormat(gr_TerrISO,"TerrISO for all mass calculations"," #gamma^{2} - #gamma_{t}^{2}(iont.C)","Terr");
AxisFormat(gr_verrISO,"verrISO for all mass calculations"," #gamma^{2} - #gamma_{t}^{2}(iont.C)","verr");
gr_TerrISO->SetMarkerSize(1.0);
gr_TerrISO->SetMarkerStyle(24);
gr_verrISO->SetMarkerSize(1.0);
gr_verrISO->SetMarkerStyle(24);
c_errISO->Divide(1,2);
c_errISO->cd(1);
gr_TerrISO->Draw("Ap");
c_errISO->cd(2);
gr_verrISO->Draw("Ap");
}

//==================================== 20230419 ERR ANA =======================================
if(c_ERRANA_on)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_ERRANA = new TCanvas("c_ERRANA","c_ERRANA",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AxisFormat(grerr_ERRYX," "," (#gamma_{i}^{2} - #gamma_{ti}^{2})^{2} + (#gamma_{j}^{2} - #gamma_{tj}^{2})^{2}  "," Y = (#Deltam/m)^{2}");
grerr_ERRYX->Draw("ap");
TF1* f_ERRYX = new TF1("f_ERRYX","pol1",0,10);
grerr_ERRYX->Fit(f_ERRYX);
lat_text = Info_fitfun_pol1(f_ERRYX);
lat_n->SetTextFont(22);
lat_n->SetTextSize(0.06);
lat_n->SetTextColor(2);
lat_n->DrawLatex(0.5,0.000001,lat_text);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *c_ERRANA_2 = new TCanvas("c_ERRANA_2","c_ERRANA_2",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h2_ERRYX->Draw("colz");
}
//_______________________________________ 20230419 ERR ANA ___________________________________________


//======================================= 2024 c_h_Mvq_2gaus ==========================================
if (c_h_Mvq_2gaus_on)
{
    TCanvas* c_Mvq_test = new TCanvas("c_Mvq_test", "c_Mvq_test", 1200, 600);
    AxisFormat(h_Mvq_14O21Mg_ReIdentify, "^{14}O/^{21}Mg Mvq Calculated", "#fraq{m}{q} [u/e]", "Counts", 2);//名称自定义
    /*
    h_Mvq_14O21Mg_ReIdentify->SetLineColor(kAzure);
    h_Mvq_14O21Mg_ReIdentify->SetLineWidth(4);
    AxisFormat(h_Mvq_14O21Mg_ReIdentify," ^{14}O/^{21}Mg Mvq reidentified ","#fraq{m}{q} [u/e]", "Counts",kAzure);
    h_Mvq_14O21Mg_ReIdentify->GetXaxis()->SetRangeUser(1.75038, 1.75058);
    h_Mvq_14O21Mg_ReIdentify->GetYaxis()->SetRangeUser(0, 400);
    h_Mvq_14O21Mg_ReIdentify->Draw();
    */
}
//_______________________________________ 2024 c_Mvq_test ___________________________________________

//======================================= 2024 ERR ANA ==========================================
if(Do_Draw_2024ERRANA_ON&&MASS_VER>=3)
{
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //TCanvas *c_gr_err_dm_all = new TCanvas("c_gr_err_dm_all","c_gr_err_dm_all",1200,600);
    //gr_err_dm_all->Draw("Ap");
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCanvas *c_h2_err_dm_all = new TCanvas("c_h2_err_dm_all","c_h2_err_dm_all",1400,1200);
    h2_err_dm_all->Draw("colz");
}

//_______________________________________ 2024 ERR ANA ___________________________________________

/////////////////////////////////___________________________   DRAW _____________________________________ /////////////////////////////////////////


                                                    } // DRAW ON







//___________________________________________________________________________________________________________________________________________//











if(Do_iont_error_out_ON&&MASS_VER>=3){outfile_errors.close();}
outfile_errors_check.close();
outfile_C_Division_n.close();


cout<<"  ------------ L-ddT-scan for loops one circle finish ------------------- "<<endl;
}////______________________________L-ddT-scan for loops end ________________________________________________

/////////////////////////////////////////////////////////////////////////
// 已经出了 L ddT 循环



outfile_large_dm.close();

if(LOOP_ON)cout<<"\033[32m||LOOP ON "<<" ||\033[0m"<<endl;
else        cout<<"\033[33m||LOOP OFF "<<" ||\033[0m"<<endl;
cout<<" \033[32m -------------------------now out of L-ddT-scan for loops  -------------------------------------- \033[0m"<<endl;

outfile_logs.close();

//Tfile_save_hist->Close(); // Close 大写C

if(cc1_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *cc1 = new TCanvas("cc1","cc1 gr2D L ddT C_intersection",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AxisFormat(gr2d_LddtC,"C_intersection","L[m]","ddT[ns]","C_intersection [m]");
gr2d_LddtC->Draw("pcol");

TGraph2D_to_outfile(FILEPATH+"gr2d_LddtC.txt",gr2d_LddtC); 

//Show_gr2d_LddtC("./INPUT/gr2d_LddtC.txt");
}
if(LOOP_ON&&cc_xn_v1_ON)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCanvas *cc_xn_v1 = new TCanvas("cc_xn_v1"," gr2D L ddT #chi_{n} v1",1200,600);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AxisFormat(gr2d_Lddt_Xn_v1,"#chi_{n} vs L ddT for mass v1","L[m]","ddT[ns]","#Chi_{n} v1 ");
gr2d_Lddt_Xn_v1->Draw("pcol");

TGraph2D_to_outfile(FILEPATH+"gr2d_Lddt_Xn_v1.txt",gr2d_Lddt_Xn_v1); 


}


for(int i=0;i<5;i++){if(ERROR_FLAG[i]){cout<<endl<<endl<<"\033[7m\033[31m  ERROR!!!! FLAG:  \033[0m "<<i<<endl;  }  }

//==========================

outfile_Lddt_Xn2.close();
/*delete L_ddT_Xn2;
delete L_ddT_Xn2_v2;
delete L_ddT_Xn2_sys;
delete L_ddT_Xn2_v2_sys;
*/
cout<<"\033[32m _______________________________FINISH__________________________________ \033[0m"<<endl;

}
//end of void main 
//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
//结束
//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
//L ddt:
//2021_36Ar_SET3 INFILE: INPUT//2021_36Ar_SET3//IdentifyResult_20230824_174150-del-2inj.txt |  18.057 0.0965
//2017_58Ni      58Ni_gtCalculator_outfile_20231121.txt   | 18.044 0.1372