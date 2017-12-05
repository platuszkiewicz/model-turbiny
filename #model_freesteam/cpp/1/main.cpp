#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip> // do setprecision

using namespace std;

    extern "C"
{
#include "freesteam/steam_ps.h"
#include "freesteam/steam_pT.h"
#include "freesteam/steam_ph.h"
#include "freesteam/region4.h"
#include "freesteam/steam_pv.h"
#include "freesteam/derivs.h"
#include "freesteam/solver2.h"
}

double PRZEL(double m, double p_omega,double tau, double p_alfa_0,      // alfa - wlot, omega - wylot, 
            double p_omega_0, double m_0)                               // 0 - warunki obliczeniowe       
{                                                                                                
       double p_alfa; // ciœnienie wylotowe
       p_alfa=pow(   pow((m/m_0),2)*tau*(pow(p_alfa_0,2)-pow(p_omega_0,2))+pow(p_omega,2)   ,0.5);
       return p_alfa;
}

double funkcja_KROK_16 (double m_g8[], double p_7[], double p_k[], double tau_7_p[],int i,double p_x) // zwraca deltê dla dowolnego p_7=p_x
{
   double eps_0=p_k[0]/p_7[0];
   double beta=0.546*0.6*0.7*0.8;
   double delta=m_g8[i]-m_g8[0]/pow(1-pow((eps_0-beta)/(1-beta),2),0.5)*p_x/p_7[0]*pow(1/tau_7_p[i],0.5)*pow(1-pow((p_k[i]/p_x-beta)/(1-beta),2),0.5);
   return delta;
}
double KROK_16(double m_g8[], double p_7[], double p_k[], double tau_7_p[],int i) // funkcja zwraca ciœnienie przed ostatni¹ grup¹ stopni
{                                                                                 // ciœnienie wyznaczane jest metod¹ bisekcji
       double pL=p_k[i]; // przeszukiwanie OD p_k Pa
       double pR=200e3;  //                 DO 200 kPa
       double pM;
       double eps=0.005; // po¿¹dana dok³adnoœæ
     
       while(fabs(pL-pR)>eps)
       {                     
             pM=(pL+pR)/2;
             if(funkcja_KROK_16(m_g8,p_7,p_k,tau_7_p,i,pL)*funkcja_KROK_16(m_g8,p_7,p_k,tau_7_p,i,pM) < 0){
                pR = pM;}
             else {
                  pL=pM;}
       }
       
       return pM;
}  

void ROZ(double m[],double p_0[],double T_0[],double p_r[],
        double m_z1[],double m_z2[],double m_z3[],double m_z4[],
        double p_zh[],
        double p_z1[],double p_z2[],double p_z3[],double p_z4[],
        double i_z[],double v_z[],
        int i,
        int zawor_czesc_otw[],
        int opcja [],
        double gx[])      
{
        
       /********** WARUNKI OBLICZENIOWE ***************************************************************************/
       m_z1[0]=0.4*m[0];  // przep³ywy w dyszy na podstawie Rys 6.8 str 97
       m_z2[0]=0.25*m[0];
       m_z3[0]=0.2*m[0];
       m_z4[0]=0.15*m[0];
                              
       p_z1[0]=12.3e6;        // cisnienie za zaworami (nale¿y za³o¿yæ) str 118
       p_z2[0]=12.3e6;
       p_z3[0]=12.3e6;
       p_z4[0]=12.3e6;
       
       double eps_1_0,eps_2_0,eps_3_0,eps_4_0; // epsilon(stosunek cieœnieñ) obliczeniowy 
       eps_1_0=p_r[0]/p_z1[0];
       eps_2_0=p_r[0]/p_z2[0];
       eps_3_0=p_r[0]/p_z3[0];
       eps_4_0=p_r[0]/p_z4[0];
       
       double beta_d=0.546; // tablica 2.2
       
       double E_1_0,E_2_0,E_3_0,E_4_0;
       
       if (eps_1_0>beta_d){
          E_1_0=pow(1-pow((eps_1_0-beta_d)/(1-beta_d),2),0.5);}
       else{
          E_1_0=1;}
              
       if (eps_2_0>beta_d){
          E_2_0=pow(1-pow((eps_2_0-beta_d)/(1-beta_d),2),0.5);}
       else{
          E_2_0=1;}        
          
       if (eps_3_0>beta_d){
          E_3_0=pow(1-pow((eps_3_0-beta_d)/(1-beta_d),2),0.5);}
       else{
          E_3_0=1;}
          
       if (eps_4_0>beta_d){
          E_4_0=pow(1-pow((eps_4_0-beta_d)/(1-beta_d),2),0.5);}
       else{
          E_4_0=1;}
                           
                 cout << "    E_1_0 " << E_1_0 << endl;
                 cout << "    E_2_0 " << E_2_0 << endl;
                 cout << "    E_3_0 " << E_3_0 << endl;
                 cout << "    E_4_0 " << E_4_0 << endl;
       
       double m_kr1_0,m_kr2_0,m_kr3_0,m_kr4_0; // przep³ywy krytyczne
       
       m_kr1_0=m_z1[0]/E_1_0; cout << "mkr1 " << m_kr1_0 << endl;
       m_kr2_0=m_z2[0]/E_2_0; cout << "mkr2 " << m_kr2_0 << endl;
       m_kr3_0=m_z3[0]/E_3_0; cout << "mkr3 " << m_kr3_0 << endl;
       m_kr4_0=m_z4[0]/E_4_0; cout << "mkr4 " << m_kr4_0 << endl << endl;
       
       /********** WARUNKI ZMIENIONE ******************************************************************************************************/      
       double a=0.96; // wspó³czynnik strat ciœnienia na zaworze szybkozamykaj¹cym i regulacyjnym, przy ca³kowitym otwarciu   
       
       p_zh[i]=a*p_0[i]; // ciœnienie za ca³kowicie otwartym zaworem
       p_zh[0]=a*p_0[0];
       double eps_h=p_r[i]/p_zh[i]; cout << " eps_h " << eps_h << endl; // stosunek ciœnienia za SR do ciœnienia przed SR przy ca³kowitym otwarciu zaworu
       
       
       SteamState S1=freesteam_set_pT(p_0[i],T_0[i]+273.15);                       // entalpia i objêtoœæ w³aœciwa za zaworami
       i_z[i]=freesteam_h(S1); cout << "i_z[" << i << "] " << i_z[i] << endl; 
       SteamState S2=freesteam_set_ph(p_zh[i],i_z[i]);
       v_z[i]=freesteam_v(S2); cout << "v_z[" << i << "] " << v_z[i] << endl;
      
       SteamState S3=freesteam_set_pT(p_0[0],T_0[0]+273.15);                       // OBLICZENIOWE entalpia i objêtoœæ w³aœciwa za zaworami
       i_z[0]=freesteam_h(S3); cout << "i_z[0] " << i_z[0] << endl; 
       SteamState S4=freesteam_set_ph(p_zh[0],i_z[0]);
       v_z[0]=freesteam_v(S4); cout << "v_z[0] " << v_z[0] << endl;
        
            /*  SteamState S2=freesteam_set_pv(13000,0.02);
              SteamState S = freesteam_set_pT(13e4,300+273.5);
              char der[]={'p','v','s'};
              double deriv = freesteam_deriv(S, der); cout << endl << " deriv " << deriv << endl << endl;
              double kappa = -freesteam_v(S)/13e6*deriv; cout << " kappa %%% " << kappa << endl << endl; */
                  
       double beta_dz=0.546; // tablica 2.2 , beta_d w warunkach zmienionych
       
       double E_h;
       
       if (eps_h>beta_dz){
          E_h=pow(1-pow((eps_h-beta_dz)/(1-beta_dz),2),0.5);}
       else{
          E_h=1;}
       
       double m_zh1,m_zh2,m_zh3,m_zh4;  // przep³ywy przez ca³kowicie otwarty zawór
       m_zh1=m_kr1_0*p_zh[i]/p_zh[0]*pow((p_zh[0]*v_z[0])/(p_zh[i]*v_z[i]),0.5)*E_h; cout << "m_zh1 " << m_zh1 << endl;
       m_zh2=m_kr2_0*p_zh[i]/p_zh[0]*pow((p_zh[0]*v_z[0])/(p_zh[i]*v_z[i]),0.5)*E_h; cout << "m_zh2 " << m_zh2 << endl;
       m_zh3=m_kr3_0*p_zh[i]/p_zh[0]*pow((p_zh[0]*v_z[0])/(p_zh[i]*v_z[i]),0.5)*E_h; cout << "m_zh3 " << m_zh3 << endl;
       m_zh4=m_kr4_0*p_zh[i]/p_zh[0]*pow((p_zh[0]*v_z[0])/(p_zh[i]*v_z[i]),0.5)*E_h; cout << "m_zh4 " << m_zh4 << endl;
       cout << "SUMA-- " << m_zh1+m_zh2+m_zh3+m_zh4 << endl << endl;
       
       zawor_czesc_otw[i]=0;           // badanie który zawór jest czêœciowo otwarty                                                                                     
        if (m[i]<m_zh1){
          zawor_czesc_otw[i]=1;}
       else{
            if(m[i]<(m_zh1+m_zh2)){
               zawor_czesc_otw[i]=2;}
            else{
                 if(m[i]<(m_zh1+m_zh2+m_zh3)){
                    zawor_czesc_otw[i]=3;}
                 else{ 
                    if(m[i]<(m_zh1+m_zh2+m_zh3+m_zh4)){   
                        zawor_czesc_otw[i]=4;}
                    else{
                        zawor_czesc_otw[i]=5;} // wszystkie zawory otwarte 
                     }
                }   
            }       
            
       cout << "Zawor czescowo otwarty to " << zawor_czesc_otw[i] << endl << endl;

      
       double m_zk,m_zhk; // przep³yw faktyczny i maksymalny w czêœciowo otwartym zaworze

       switch(zawor_czesc_otw[i])
       {
         case 0:
              cout << "blad!" << endl;
         case 1:
              m_z1[i]=m[i];
              m_z2[i]=0; 
              m_z3[i]=0; 
              m_z4[i]=0;  
              m_zk=m_z1[i];
              m_zhk=m_zh1;
              break; 
         case 2:
              m_z1[i]=m_zh1;                  //  p_z1[i]=p_zh;            
              m_z2[i]=m[i]-m_zh1;             //  p_z2[i]=p_zk;
              m_z3[i]=0;                      //  p_z3[i]=p_r[i];
              m_z4[i]=0;                      //  p_z4[i]=p_r[i];
              m_zk=m_z2[i];    
              m_zhk=m_zh2;                     
              break;
         case 3:
              m_z1[i]=m_zh1;                  //  p_z1[i]=p_zh;                   
              m_z2[i]=m_zh2;                  //  p_z2[i]=p_zh;           
              m_z3[i]=m[i]-m_zh1-m_zh2;       //  p_z3[i]=p_zk;
              m_z4[i]=0;                      //  p_z4[i]=p_r[i];
              m_zk=m_z3[i]; 
              m_zhk=m_zh3; 
              break;
         case 4: 
              m_z1[i]=m_zh1;                  //  p_z1[i]=p_zh;                      
              m_z2[i]=m_zh2;                  //  p_z2[i]=p_zh;                       
              m_z3[i]=m_zh3;                  //  p_z3[i]=p_zh;            
              m_z4[i]=m[i]-m_zh1-m_zh2-m_zh3; //  p_z4[i]=p_zk;
              m_zk=m_z4[i];
              m_zhk=m_zh4;                           
              break;
         case 5:
              m_z1[i]=m_zh1;                  //  p_z1[i]=p_zh;                      
              m_z2[i]=m_zh2;                  //  p_z2[i]=p_zh;                       
              m_z3[i]=m_zh3;                  //  p_z3[i]=p_zh;            
              m_z4[i]=m_zh4;                  //  p_z4[i]=p_zh;   
              m_zk=m_z4[i];
              m_zhk=m_zh4;     
              cout << "Wszystkie zawory otwarte!" << endl;          
         }
         
         //************** ciœienie za czêœciowo otwartym zaworem p_zk ************//
       double p_zk=0;
       double g=m_zk/m_zhk; cout << " g " << g << endl; gx[i]=g;
       double A,B;
       
       if(eps_h>beta_dz){
           A=(pow(p_r[i],2))*(1-beta_dz);
           B=pow(  pow(1-beta_dz,2) - pow(eps_h-beta_dz,2) , 0.5 )/(1-beta_dz);             
                      cout << "eps_h/(beta_dz*B) " << eps_h/(beta_dz*B) << endl;    
                      if(g<=(eps_h/(beta_dz*B))){
                           p_zk=p_r[i]*(-beta_dz+(1-beta_dz)*pow(1+(1-2*beta_dz)*pow(g*p_zh[i]/p_r[i],2),0.5))/(1-2*beta_dz);
                           cout << " p_r " << p_r[i] << endl << "beta_dz " << beta_dz << endl << "p_zh " << p_zh[i] << endl;
                           
                           
                           cout << "p_zk " << p_zk << endl << endl;
                           cout<<"#1"<<endl;opcja[i]=1;}
                      else{
                           p_zk=g*p_zh[i]*B;cout<<"#2 "<< endl << "g " << g << endl << "p_r " << p_r[i] << endl << "B " << B << endl << "p_zk " << p_zk << endl << endl ;opcja[i]=2;}                                                          
          }
       else{
           A=(pow(p_r[i],2))*(1-beta_dz);
           B=1;    
                      cout << "eps_h/(beta_dz*B) " << eps_h/(beta_dz*B) << endl;     
                      if(g<=eps_h/(beta_dz*B)){
                           p_zk=p_r[i]*(-beta_dz+(1-beta_dz)*pow(1+(1-2*beta_dz)*pow(g*p_zh[i]/p_r[i],2),0.5))/(1-2*beta_dz);
                           cout<<"#3"<<endl;opcja[i]=3;}
                      else{
                           p_zk=g*p_zh[i]*B;cout<<"#4 "<< endl << "g " << g << endl << "p_r " << p_r[i] << endl << "B " << B << endl << "p_zk " << p_zk << endl << endl;opcja[i]=4;}     
            }
          cout << " A " << A << endl;
          cout << " B " << B << endl;
          
           switch(zawor_czesc_otw[i])
       {
         case 0:
              cout << "blad!" << endl;
         case 1:
              p_z1[i]=p_zk;
              p_z2[i]=p_r[i];
              p_z3[i]=p_r[i];
              p_z4[i]=p_r[i];
              break; // zaleznosci wyznaczono dalej
         case 2:
              p_z1[i]=p_zh[i];            
              p_z2[i]=p_zk;
              p_z3[i]=p_r[i];
              p_z4[i]=p_r[i];                  
              break;
         case 3:
              p_z1[i]=p_zh[i];                   
              p_z2[i]=p_zh[i];           
              p_z3[i]=p_zk;
              p_z4[i]=p_r[i];
              break;
         case 4: 
              p_z1[i]=p_zh[i];                      
              p_z2[i]=p_zh[i];                       
              p_z3[i]=p_zh[i];            
              p_z4[i]=p_zk;                           
              break;
         case 5:
              p_z1[i]=p_zh[i];                      
              p_z2[i]=p_zh[i];                       
              p_z3[i]=p_zh[i];            
              p_z4[i]=p_zh[i];        
              cout << "Wszystkie zawory otwarte!" << endl;          
         }
       

         cout << "Zawor 1 m_z1 " << m_z1[i] << " // p_z1 " << p_z1[i] << endl;
         cout << "Zawor 2 m_z2 " << m_z2[i] << " // p_z2 " << p_z2[i] << endl;  
         cout << "Zawor 3 m_z3 " << m_z3[i] << " // p_z3 " << p_z3[i] << endl; 
         cout << "Zawor 4 m_z4 " << m_z4[i] << " // p_z4 " << p_z4[i] << endl;                
}

double SPADEK_IZENTROPOWY_i(double p1, double i1, double p2) // oblicza spadek izentropowy od punktu 1 do punktu 2
{
     SteamState S1=freesteam_set_ph(p1,i1); // punkt 1
     double s12=freesteam_s(S1);            // entropia
     SteamState S2=freesteam_set_ps(p2,s12);// punkt 2     
     double H=i1-freesteam_h(S2);
	 
	 return H;
}

double SPADEK_IZENTROPOWY_T(double p1, double T, double p2) // oblicza spadek izentropowy od punktu 1 do punktu 2, temp w Kelwinach
{
       cout << "   spadek izentropowy " << endl;
     SteamState S1=freesteam_set_pT(p1,T); cout << "   p1 " << p1 << " T1 " << T-273.15 <<endl;// punkt 1
     double i1=freesteam_h(S1); cout << "   i1 " << i1 << endl;
     double s12=freesteam_s(S1); cout << "    s12 " << s12 << endl;                    // entropia
     SteamState S2=freesteam_set_ps(p2,s12); cout << "   i2 " << freesteam_h(S2)<<endl;// punkt 2     
     double H=i1-freesteam_h(S2);
	 cout << "      H " << H << endl;
	 return H;
}

double v_pT(double p, double T)
{
       SteamState S1=freesteam_set_pT(p,T+273.15);
       return freesteam_v(S1);
}

double v_pi(double p, double i)
{
       SteamState S1=freesteam_set_ph(p,i);
       return freesteam_v(S1);
}

void REG(   double m_z1[], double m_z2[], double m_z3[], double m_z4[], 
            double p_z1[], double p_z2[], double p_z3[], double p_z4[], 
            double H_z1[], double H_z2[], double H_z3[], double H_z4[], 
            double i_z[], 
            double p_r[],
            double p_zh[],
            int i,
            int zawor_czesc_otw[],
            double P_ru[], double delta_P_k[], double delta_P_tw[], double P_r[], 
            double T_r[], double v_r[], double i_r[], ofstream& spr_reg)
{
      H_z1[i]=SPADEK_IZENTROPOWY_i(p_z1[i],i_z[i],p_r[i]); cout << "H_1 " << H_z1[i] << endl;
      H_z2[i]=SPADEK_IZENTROPOWY_i(p_z2[i],i_z[i],p_r[i]); cout << "H_2 " << H_z2[i] << endl;
      H_z3[i]=SPADEK_IZENTROPOWY_i(p_z3[i],i_z[i],p_r[i]); cout << "H_3 " << H_z3[i] << endl;
      H_z4[i]=SPADEK_IZENTROPOWY_i(p_z4[i],i_z[i],p_r[i]); cout << "H_4 " << H_z4[i] << endl;
      
      cout << " Mz_1 * Hz1 " << m_z1[i]*H_z1[i] << endl;
      cout << " Mz_2 * Hz2 " << m_z2[i]*H_z2[i] << endl;
      cout << " Mz_3 * Hz3 " << m_z3[i]*H_z3[i] << endl;
      cout << " Mz_4 * Hz4 " << m_z4[i]*H_z4[i] << endl;
      
      double H_h0=SPADEK_IZENTROPOWY_i(p_zh[i],i_z[i],p_r[i]);
      
      double n=3600;        // obr/min
      double D=1.500;       // œrednia œrednica stopnia regulacyjnego
      double u=M_PI*D*n/60; // prêdkoœæ obwodowa stopnia
      
      double x_h0=u/pow(2*H_h0,0.5); cout <<"x_h0 " << x_h0 << endl;
      double x_1=0,x_2=0,x_3=0,x_4=0;
      
      switch (zawor_czesc_otw[i])
      {
      case 1: x_1=x_h0*pow(H_h0/H_z1[i],0.5); 
      x_2=0; 
      x_3=0;
      x_4=0;
      break;
      
      case 2: x_1=x_h0*pow(H_h0/H_z1[i],0.5); 
      x_2=x_h0*pow(H_h0/H_z2[i],0.5); 
      x_3=0;
      x_4=0;
      break;
      
      case 3: x_1=x_h0*pow(H_h0/H_z1[i],0.5); 
      x_2=x_h0*pow(H_h0/H_z2[i],0.5); 
      x_3=x_h0*pow(H_h0/H_z3[i],0.5); 
      x_4=0;
      break;
      
      case 4: x_1=x_h0*pow(H_h0/H_z1[i],0.5); 
      x_2=x_h0*pow(H_h0/H_z2[i],0.5); 
      x_3=x_h0*pow(H_h0/H_z3[i],0.5); 
      x_4=x_h0*pow(H_h0/H_z4[i],0.5);
      break;
      
      case 5: x_1=x_h0*pow(H_h0/H_z1[i],0.5); 
      x_2=x_h0*pow(H_h0/H_z2[i],0.5); 
      x_3=x_h0*pow(H_h0/H_z3[i],0.5); 
      x_4=x_h0*pow(H_h0/H_z4[i],0.5);
      
      }

        cout << "Zawor czescowo otwarty to " << zawor_czesc_otw[i] << endl << endl;
        cout <<"x_1 " << x_1 << " x_1/x0 " << x_1/x_h0 << endl;
        cout <<"x_2 " << x_2 << " x_2/x0 " << x_2/x_h0 << endl;
        cout <<"x_3 " << x_3 << " x_3/x0 " << x_3/x_h0 << endl;
        cout <<"x_4 " << x_4 << " x_4/x0 " << x_4/x_h0 <<  endl;
        
      double eta_u_obl_1=0.75, eta_u_obl_2=0.75, eta_u_obl_3=0.75, eta_u_obl_4=0.75;
      double eta_u1=eta_u_obl_1*(-1*pow(x_1/x_h0,2)+2*(x_1/x_h0)); cout << "eta_u1  " << eta_u1 << endl;
      double eta_u2=eta_u_obl_2*(-1*pow(x_2/x_h0,2)+2*(x_2/x_h0)); cout << "eta_u2  " << eta_u2 << endl;
      double eta_u3=eta_u_obl_3*(-1*pow(x_3/x_h0,2)+2*(x_3/x_h0)); cout << "eta_u3  " << eta_u3 << endl;
      double eta_u4=eta_u_obl_4*(-1*pow(x_4/x_h0,2)+2*(x_4/x_h0)); cout << "eta_u4  " << eta_u4 << endl;
      
      P_ru[i]=m_z1[i]/3.6*H_z1[i]*eta_u1+m_z2[i]/3.6*H_z2[i]*eta_u2+m_z3[i]/3.6*H_z3[i]*eta_u3+m_z4[i]/3.6*H_z4[i]*eta_u4;
      
      double B=100; // szerokoœæ wieñca ³opatkowego 50
      double l=125; // wysokoœæ ³opatek wiruj¹cych 20
      double x_h=x_h0; // wskaŸnik prêdkoœci dla strumienia przep³ywaj¹cego przez zawór ca³kowicie otwarty
      double z=zawor_czesc_otw[i]; // liczba par koñców segmentów dyszowych
      double F=150000; // minimalny przekrój czynnych dysz 2500
      double E=0; // wspó³czynnik przyjmuj¹cy wartoœci jak ni¿ej
      if (p_r[i]/p_zh[i]<0.546){
               E=1;}
      else {
           E=1.3033-0.5555*p_r[i]/p_zh[i];}
      
      delta_P_k[i]=0.11*B*l*x_h*P_ru[i]*z/(E*F);
      
      i_r[i]=i_z[i]-(P_ru[i]-delta_P_k[i])/(m_z1[i]/3.6+m_z2[i]/3.6+m_z3[i]/3.6+m_z4[i]/3.6); // m_zj jest podawane w t/h!!!
      
      double eps_v=0.001; // dok³adnoœæ obliczeñ
      
      double vr_p=v_pi(p_r[i],i_r[i]); 
      double vr=0; double vr_poprzednie;
      cout << "vr_p start " << vr_p << endl;
      
      do 
      {
          double fi=0.755; // wspó³czynnik zwi¹zany z geometri¹ stopnia
          double e=zawor_czesc_otw[i]/4; // ³uk zasilania
          
          delta_P_tw[i]=(fi*pow(D-1,2)+0.344e-3*(1-e))*D*pow(l,1.5)*pow(u,3)/(1e6*vr_p);
          
          P_r[i]=P_ru[i]-delta_P_k[i]-delta_P_tw[i];
          
          i_r[i]=i_z[i]-P_r[i]/(m_z1[i]/3.6+m_z2[i]/3.6+m_z3[i]/3.6+m_z4[i]/3.6);
          
          
          vr=v_pi(p_r[i],i_r[i]); cout << "vr " << vr << endl;
          
          vr_poprzednie=vr_p; cout << " vr_poprzednie " << vr_poprzednie << endl;
          vr_p=vr;
      } while ((   (fabs(vr_poprzednie-vr))/vr)>(eps_v));
      
      cout << "P_ru " << P_ru[i] << endl << "----delta_P_k " << delta_P_k[i] << endl << "----delta_P_tw " << delta_P_tw[i] << endl << "P_r " << P_r[i] << endl;
      v_r[i]=vr;      
      SteamState S3=freesteam_set_ph(p_r[i],i_r[i]);
      cout << "p_r " << p_r[i] << endl;
      T_r[i]=freesteam_T(S3)-273.15;
      cout << "T_r " << T_r[i]<< endl;
      cout << "v_r " << v_r[i] << endl;
      cout << "i_r " << i_r[i] << endl;
      spr_reg.setf(ios::fixed);

      spr_reg << H_z1[i] << " " << H_z2[i] << " " << H_z3[i] << " " << H_z4[i] << " " << m_z1[i]+m_z2[i]+m_z3[i]+m_z4[i] << " " << eta_u1 << " " << eta_u2 << " " << eta_u3 << " " << eta_u4 << " " << P_r[i] << " " << P_ru[i] << " " << delta_P_k[i] << " " << delta_P_tw[i] << " " << m_z1[i]/3.6*H_z1[i]+m_z2[i]/3.6*H_z2[i]+m_z3[i]/3.6*H_z3[i]+m_z4[i]/3.6*H_z4[i]<< endl;      
}

void OSI(double m_g, double p_alfa, double T_alfa, double p_omega, double eta,               // wjeœcia
         double H_x[], double P_x[], double T_x[], double i_x[], double v_x[],               // wyjœcia
         int i)
{
      H_x[i]=SPADEK_IZENTROPOWY_T(p_alfa,T_alfa+273.15,p_omega); 
      P_x[i]=m_g/3.6*H_x[i]*eta; 
      SteamState S_alfa=freesteam_set_pT(p_alfa, T_alfa+273.15); 
      double i_alfa=freesteam_h(S_alfa); 
     
      double i_omega=i_alfa-P_x[i]/(m_g/3.6);
      i_x[i]=i_omega;
      SteamState S_omega=freesteam_set_ph(p_omega,i_omega);
      T_x[i]=freesteam_T(S_omega)-273.15;
      v_x[i]=freesteam_v(S_omega);  
}

int main(void)
{
    const int K=1157; // liczba punktów pracy (bez pkt obliczeniowego!) 
    const int N=K+1; // rozmiar tabel(K+pkt obliczeniowy)
    
    /*################ D E K L A R A C J E ##############*/
	
	double m_0[N],m[N],m_r[N]; // para œwie¿a, rozrz¹d
	double m_z1[N],m_z2[N],m_z3[N],m_z4[N],m_zh1[N],m_zh2[N],m_zh3[N],m_zh4[N];
	double m_g1[N],m_u1[N],m_g2[N],m_u2[N]; // st I oraz II
	double m_p[N]; // przegrzew wtórny
	double m_g3[N],m_u3[N],m_g4[N],m_u4[N],m_g5[N],m_u5[N],m_g6[N],m_u6[N],m_g7[N],m_u7[N],m_g8[N],m_k[N]; // st III-VII
	double m_wt[N]; // wtrysk do pary wtornej
	double m_przec[N]; // przeciek z czêœci WP
	
	int zawor_czesc_otw[N];
	
	double tau_r_p[N],tau_1_p[N],tau_p_p[N],tau_3_p[N],tau_4_p[N],tau_5_p[N],tau_N_p[N],tau_7_p[N]; // stosunki temperatur
	
	double T_0[N],T_r[N],T_1[N],T_2[N],T_p[N],T_3[N],T_4[N],T_5[N],T_6[N],T_N[N],T_7[N],T_k[N]; // temperatury
	
	double p_0[N],p_r[N],p_1[N],p_2[N],p_p[N],p_3[N],p_4[N],p_5[N],p_6[N],p_N[N],p_7[N],p_k[N]; // ciœnienia
	
    double p_zh[N],p_z1[N],p_z2[N],p_z3[N],p_z4[N],p_zh1[N],p_zh2[N],p_zh3[N],p_zh4[N];	// rozrz¹d
    double gx[N];
    double P_ru[N], delta_P_k[N], delta_P_tw[N], P_r[N];
    int opcja[N]; // opisuje sposób liczenia p_zk

    
    double i_0[N],i_r[N],i_1[N],i_2[N],i_p[N],i_3[N],i_4[N],i_5[N],i_6[N],i_N[N],i_7[N],i_k[N]; // entalpie
    double i_z[N]; // rozrz¹d
	
    double v_0[N],v_z[N],v_r[N],v_1[N],v_2[N],v_p[N],v_3[N],v_4[N],v_5[N],v_6[N],v_N[N],v_7[N],v_k[N]; // objêtoœæ w³aœciwa
    
    double H_z1[N],H_z2[N],H_z3[N],H_z4[N]; // osi¹gi st. reg.
    
    double P_1[N], P_2[N], P_3[N], P_4[N], P_5[N], P_6[N], P_7[N], P_8[N]; // OSI¥GI STOPNI
    double P[N], P_el[N], delta_P_el[N];
    double H_1[N], H_2[N], H_3[N], H_4[N], H_5[N], H_6[N], H_7[N], H_8[N]; 

    double eta_1_0, eta_2_0, eta_3_0, eta_4_0, eta_5_0, eta_6_0;
    eta_1_0=0.82; eta_2_0=0.82; eta_3_0=0.84; eta_4_0=0.845; eta_5_0=0.848; eta_6_0=0.854; // tylko grupy 1-6 !
    
    double eta_io_7; 
    eta_io_7=0.85;
    
    double eps_8[N];
    double eta_i8_s[N];
    
    double y_7[N], y_7_p[N], y_k[N], y_k_p[N]; // dotyczy grupy 7 i 8;
    
    double a=0.95; // straty ciœnienia w przelotni
    int iter_tau[N]; iter_tau[0]=0;// dla ka¿dego punktu obliczeniowego zlicza iloœæ iteracji koniecznych do dok³adnego wyznaczenia tau 
    /*################ ZAPIS DO PLIKU ZEWNÊTRZNEGO ###########################*/
   // ofstream st_reg; 
   //st_reg.open("st_reg.txt");
    
    ofstream spr_reg;
    //spr_reg.open("spr_reg.txt"); spr_reg << "H_z1 " << "H_z2 " << "H_z3 " << "H_z4 " << "m_0 " << "eta_u1 " << "eta_u2 " << "eta_u3 " << "eta_u4 " << "P_r " << "P_ru " << "P_k " << "P_tw " << "P_0 " << endl;
    
    //ofstream osiagi; 
    //osiagi.open("osiagi.txt");    
    
    //ofstream testing; 
    //testing.open("wyjscie_testing.txt");   
    	
	/*################ W A R U N K I      O B L I C Z E N I O W E ############*/
	P_el[0]=220e6;
	m_u1[0]=32; m_u2[0]=58; m_u3[0]=25; m_u4[0]=15; m_u5[0]=30; m_u6[0]=30; m_u7[0]=6; 
	m_przec[0]=10; m_wt[0]=28; // 10 - 28
	
                                      m_0[0]=m[0]=656;
	T_0[0]=538; p_0[0]=12.8e6;        m_r[0]=m_0[0];
	T_r[0]=477; p_r[0]=8.41e6;        m_g1[0]=m_0[0];
	T_1[0]=383;	p_1[0]=4.24e6;        m_g2[0]=m_0[0]-m_u1[0];
	T_2[0]=332; p_2[0]=2.82e6; 
	T_p[0]=535; p_p[0]=2.49e6;        m_g3[0]=m_p[0]=m_0[0]-m_u1[0]-m_u2[0]+m_wt[0]-m_przec[0];
	T_3[0]=410; p_3[0]=1.00e6;        m_g4[0]=m_0[0]-m_u1[0]-m_u2[0]-m_u3[0];              
	T_4[0]=323;	p_4[0]=0.49e6;        m_g5[0]=m_0[0]-m_u1[0]-m_u2[0]-m_u3[0]-m_u4[0];
	T_5[0]=249;	p_5[0]=249e3;         m_g6[0]=m_0[0]-m_u1[0]-m_u2[0]-m_u3[0]-m_u4[0]-m_u5[0];
	T_6[0]=190;	p_6[0]=137e3;  
	T_N[0]=190;	p_N[0]=136.5e3;       m_g7[0]=m_0[0]-m_u1[0]-m_u2[0]-m_u3[0]-m_u4[0]-m_u5[0]-m_u6[0];
	T_7[0]=73;  p_7[0]=32e3;          m_g8[0]=m_0[0]-m_u1[0]-m_u2[0]-m_u3[0]-m_u4[0]-m_u5[0]-m_u6[0]-m_u7[0];
	T_k[0]=30;	p_k[0]=4.75e3; 
	
	v_0[0]=v_pT(p_0[0],T_0[0]); 
	v_r[0]=v_pT(p_r[0],T_r[0]); 
	v_1[0]=v_pT(p_1[0],T_1[0]); 
	v_2[0]=v_pT(p_2[0],T_2[0]);
	v_p[0]=v_pT(p_p[0],T_p[0]);
	v_3[0]=v_pT(p_3[0],T_3[0]);
	v_4[0]=v_pT(p_4[0],T_4[0]);
	v_5[0]=v_pT(p_5[0],T_5[0]);
	v_6[0]=v_pT(p_6[0],T_6[0]);
	v_N[0]=v_pT(p_N[0],T_N[0]);
	v_7[0]=v_pT(p_7[0],T_7[0]);
	v_k[0]=v_pT(p_k[0],T_k[0]);
	
	/*################ W E J Œ C I A    M O D E L U ##########################*/
	ifstream wejscie;
	wejscie.open("wejscie_testing_n=1157.txt");
	for (int i=1; i<=K; i++) 	// indeks 0 tabeli odnosi siê do wielkoœci obliczeniowej
	{
        wejscie >> m_0[i]; 
        wejscie >> T_0[i]; 
        wejscie >> p_0[i]; 
        wejscie >> T_p[i]; 
        wejscie >> p_k[i]; 
        wejscie >> m_u1[i];
        wejscie >> m_u2[i];
        wejscie >> m_u3[i];
        wejscie >> m_u4[i];
        wejscie >> m_u5[i];
        wejscie >> m_u6[i];
        wejscie >> m_u7[i];     
        
        wejscie >> m_wt[i];//=0.0871*m_0[i]-23.5192;         
        wejscie >> m_przec[i];//=0.0140*m_0[i]+0.1044; 
   }

	/*################ G £ Ó W N Y    A L G O R Y T M  ##########################*/
	for (int i=1; i<=K; i++) // pêtla po kolejnych punktach pracy
	{
        cout << endl << "***************************** Punkt pracy numer " << i << " *****************************" << endl; iter_tau[i]=0;
        
        cout << "Obliczenia bilansu masy (1-12)" << endl; SteamState S0=freesteam_set_pT(p_0[i],T_0[i]+273.15); i_0[i]=freesteam_h(S0);
        m[i]=m_0[i];             cout << "m    " << m[i]    << " t/h" << endl; // kroki 1-12
        m_r[i]=m[i];             cout << "m_r  " << m_r[i]  << endl;
        m_g1[i]=m_r[i];          cout << "m_g1 " << m_g1[i] << endl;
        m_g2[i]=m_r[i]-m_u1[i];  cout << "m_g2 " << m_g2[i] << endl;
        m_p[i]=m_g2[i]-m_u2[i]+m_wt[i]-m_przec[i];  cout << "m_p  " << m_p[i]  << endl;
        m_g3[i]=m_p[i];          cout << "m_g3 " << m_g3[i] << endl;
        m_g4[i]=m_g3[i]-m_u3[i]; cout << "m_g4 " << m_g4[i] << endl;
        m_g5[i]=m_g4[i]-m_u4[i]; cout << "m_g5 " << m_g5[i] << endl;
        m_g6[i]=m_g5[i]-m_u5[i]; cout << "m_g6 " << m_g6[i] << endl;
        m_g7[i]=m_g6[i]-m_u6[i]; cout << "m_g7 " << m_g7[i] << endl;
        m_g8[i]=m_g7[i]-m_u7[i]; cout << "m_g8 " << m_g8[i] << endl << endl;
        m_k[i]=m_g8[i];
               
        tau_1_p[i]=1;             // krok 13
        tau_r_p[i]=tau_1_p[i];    // krok 13
        tau_p_p[i]=T_p[i]/T_p[0]; // krok 14
        tau_3_p[i]=tau_p_p[i];    // krok 15
        tau_4_p[i]=tau_p_p[i];
        tau_5_p[i]=tau_p_p[i];
        tau_N_p[i]=tau_p_p[i];
        tau_7_p[i]=tau_p_p[i];

double eps_tau=0.01;       
double tau_r, tau_1, tau_p, tau_3, tau_4, tau_5, tau_N, tau_7;  
tau_r=0, tau_1=0, tau_p=0, tau_3=0, tau_4=0, tau_5=0, tau_N=0, tau_7=0;
double tau_r_poprzednie, tau_1_poprzednie, tau_p_poprzednie, tau_3_poprzednie, tau_4_poprzednie, tau_5_poprzednie, tau_N_poprzednie, tau_7_poprzednie;         

do
{
        cout << "Obliczenia cisnien (16-25)" << endl;
        p_7[i]=KROK_16( m_g8,  p_7,  p_k,  tau_7_p, i);                cout << "p_7 " << p_7[i] << " Pa" << endl;    
        p_N[i]=PRZEL(m_g7[i],p_7[i],tau_7_p[i],p_N[0],p_7[0],m_g7[0]); cout << "p_N " << p_N[i] << endl;
        p_6[i]=p_N[i]/a;                                               cout << "p_6 " << p_6[i] << endl;                                                    
        p_5[i]=PRZEL(m_g6[i],p_6[i],tau_5_p[i],p_5[0],p_6[0],m_g6[0]); cout << "p_5 " << p_5[i] << endl;
        p_4[i]=PRZEL(m_g5[i],p_5[i],tau_4_p[i],p_4[0],p_5[0],m_g5[0]); cout << "p_4 " << p_4[i] << endl;
        p_3[i]=PRZEL(m_g4[i],p_4[i],tau_3_p[i],p_3[0],p_4[0],m_g4[0]); cout << "p_3 " << p_3[i] << endl;
        p_p[i]=PRZEL(m_g3[i],p_3[i],tau_p_p[i],p_2[0],p_3[0],m_g3[0]); cout << "p_p " << p_p[i] << endl;
        p_2[i]=p_p[i]+(p_2[0]-p_p[0])*pow(m_p[i]/m_p[0],2);            cout << "p_2 " << p_2[i] << endl;
        p_1[i]=PRZEL(m_g2[i],p_2[i],tau_1_p[i],p_1[0],p_2[0],m_g2[0]); cout << "p_1 " << p_1[i] << endl;
        p_r[i]=PRZEL(m_g1[i],p_1[i],tau_r_p[i],p_r[0],p_1[0],m_g1[0]); cout << "p_r " << p_r[i] << endl << endl;

        cout << "********************Obliczenia rozrzadu pary (26)" << endl;
        
        ROZ(m,p_0,T_0,p_r,
        m_z1,m_z2,m_z3,m_z4,
        p_zh,
        p_z1,p_z2,p_z3,p_z4,
        i_z,v_z,i,zawor_czesc_otw,opcja,gx);
        
        
        cout << "********************Obliczenia stopnia regulacyjnego (27)" << endl;
      
        REG( m_z1, m_z2, m_z3, m_z4, 
             p_z1, p_z2, p_z3, p_z4, 
             H_z1, H_z2, H_z3, H_z4, 
             i_z, 
             p_r,
             p_zh,
             i,
             zawor_czesc_otw,
             P_ru, delta_P_k, delta_P_tw, P_r, 
             T_r,  v_r,  i_r,spr_reg );                  /* spr_reg << m_z1[i] << " ";       
                                                        spr_reg << m_z2[i] << " ";   
                                                         spr_reg << m_z3[i] << " ";   
                                                         spr_reg << m_z4[i] << " " << endl;  */
                                                         
                                                         
    // cout << "m_g1 " << m_g1[i] << " p_alfa " << p_r[i] << " T alfa " << T_r[i] << " p omega " << p_1[i] << endl;
     
     OSI(m_g1[i], p_r[i], T_r[i], p_1[i], eta_1_0, // 1. grupa stopni
         H_1,  P_1,  T_1, i_1,  v_1,               
         i);
     cout << " P_1 " << P_1[i] << " T_1 " << T_1[i] << " i_1 " << i_1[i] << " v_1 " << v_1[i] << endl;  
     
     OSI(m_g2[i], p_1[i], T_1[i], p_2[i], eta_2_0, // 2. grupa stopni
         H_2,  P_2,  T_2, i_2,  v_2,               
         i);              
     cout << " P_2 " << P_2[i] << " T_2 " << T_2[i] << " i_2 " << i_2[i] << " v_2 " << v_2[i] << endl;
     
     SteamState P=freesteam_set_pT(p_p[i],T_p[i]+273.15); // parametry pary przegrzanej
     i_p[i]=freesteam_h(P); cout << i_p[i] << endl;
     v_p[i]=freesteam_v(P); cout << v_p[i] << endl;
     
     OSI(m_g3[i], p_p[i], T_p[i], p_3[i], eta_3_0, // 3. grupa stopni
         H_3,  P_3,  T_3, i_3,  v_3,               
         i);              
     cout << " P_3 " << P_3[i] << " T_3 " << T_3[i] << " i_3 " << i_3[i] << " v_3 " << v_3[i] << endl;     
     
     OSI(m_g4[i], p_3[i], T_3[i], p_4[i], eta_4_0, // 4. grupa stopni
         H_4,  P_4,  T_4, i_4,  v_4,               
         i);              
     cout << " P_4 " << P_4[i] << " T_4 " << T_4[i] << " i_4 " << i_4[i] << " v_4 " << v_4[i] << endl;     
     
     OSI(m_g5[i], p_4[i], T_4[i], p_5[i], eta_5_0, // 5. grupa stopni
         H_5,  P_5,  T_5, i_5,  v_5,               
         i);              
     cout << " P_5 " << P_5[i] << " T_5 " << T_5[i] << " i_5 " << i_5[i] << " v_5 " << v_5[i] << endl; 
     
     OSI(m_g6[i], p_5[i], T_5[i], p_6[i], eta_6_0, // 6. grupa stopni
         H_6,  P_6,  T_6, i_6,  v_6,               
         i);              
     cout << " P_6 " << P_6[i] << " T_6 " << T_6[i] << " i_6 " << i_6[i] << " v_6 " << v_6[i] << endl; 
     
     // ****************************** P R Z E L O T N I A **********************************************// 
     SteamState N=freesteam_set_ph(p_N[i],i_6[i]); // parametry pary w przelotni                         //   
     T_N[i]=freesteam_T(N)-273.15; cout << T_N[i] << endl;                                               //      
     v_N[i]=freesteam_v(N); cout << v_N[i] << endl;                                                      //
     // *************************************************************************************************//
     
     // ****************************** G R U P A    7 ***************************************************// 
     double eps_y_7=0.001;                                                                               //
     double y_7_poprzednie; // zmienna pomocnicza do okreœlenia warunku while                            //
     y_7_p[i]=1;                                                                                         //
     y_7[i]=0;                                                                                           //
     do                                                                                                  //
     {                                                                                                   //
         H_7[i]=SPADEK_IZENTROPOWY_i(p_N[i],i_6[i],p_7[i]);                                              //
         P_7[i]=m_g6[i]/3.6*H_7[i]*eta_io_7*y_7_p[i];                                                    //
         i_7[i]=i_6[i]-P_7[i]/(m_g7[i]/3.6);                                                             //
         SteamState S7=freesteam_set_ph(p_7[i],i_7[i]);                                                  //
         y_7[i]=freesteam_x(S7);                                                                         //
         y_7_poprzednie=y_7_p[i];                                                                        //
         y_7_p[i]=y_7[i];                                                                                //
     } while(fabs((y_7_poprzednie-y_7[i])/y_7[i])>eps_y_7);                                              //
     cout << "y_7[i] " << y_7[i] << endl;                                                                //
                                                                                                         //
     SteamState S7=freesteam_set_ph(p_7[i],i_7[i]);                                                      //
     v_7[i]=freesteam_v(S7);                                                                             //
     T_7[i]=freesteam_T(S7)-273.15;                                                                      //
     cout << "p_7[i] " << p_7[i] << endl;                                                                //
     cout << "i_7[i] " << i_7[i] << endl;                                                                //
     cout << "v_7[i] " << v_7[i] << endl;                                                                //
     cout << "T_7[i] " << T_7[i] << endl;   
     cout << " P_7 " << P_7[i] << " T_7 " << T_7[i] << " i_7 " << i_7[i] << " v_7 " << v_7[i] << endl;                                                             //
     // *************************************************************************************************// 
     
     // ****************************** W E R Y F I K A C J A   T A U ************************************//  
     
     
     cout << "p=100 000 (1 bar) T=300C v= " << v_pT(100000,50) << endl; 
     
     tau_r=(p_r[i]*v_r[i])/(p_r[0]*v_r[0]); cout << tau_r << endl;
     tau_1=(p_1[i]*v_1[i])/(p_1[0]*v_1[0]); cout << tau_1 << endl;
     tau_p=(p_p[i]*v_p[i])/(p_p[0]*v_p[0]); cout << tau_p << endl;
     tau_3=(p_3[i]*v_3[i])/(p_3[0]*v_3[0]); cout << tau_3 << endl;
     tau_4=(p_4[i]*v_4[i])/(p_4[0]*v_4[0]); cout << tau_4 << endl;
     tau_5=(p_5[i]*v_5[i])/(p_5[0]*v_5[0]); cout << tau_5 << endl;
     tau_N=(p_N[i]*v_N[i])/(p_N[0]*v_N[0]); cout << tau_N << endl;
     tau_7=(p_7[i]*v_7[i])/(p_7[0]*v_7[0]); cout << tau_7 << endl;
     
     tau_r_poprzednie=tau_r_p[i];
     tau_1_poprzednie=tau_1_p[i];
     tau_p_poprzednie=tau_p_p[i];
     tau_3_poprzednie=tau_3_p[i];
     tau_4_poprzednie=tau_4_p[i];
     tau_5_poprzednie=tau_5_p[i];
     tau_N_poprzednie=tau_N_p[i];
     tau_7_poprzednie=tau_7_p[i];
     
     tau_r_p[i]=tau_r;
     tau_1_p[i]=tau_1;
     tau_p_p[i]=tau_p;
     tau_3_p[i]=tau_3;
     tau_4_p[i]=tau_4;
     tau_5_p[i]=tau_5;
     tau_N_p[i]=tau_N;
     tau_7_p[i]=tau_7;
     
     cout << "tau_r_poprzednie " << tau_r_poprzednie << "tau_r " << tau_r << endl;
     cout << "tau_1_poprzednie " << tau_1_poprzednie << "tau_1 " << tau_1 << endl;
     cout << "tau_p_poprzednie " << tau_p_poprzednie << "tau_p " << tau_p << endl;
     cout << "tau_3_poprzednie " << tau_3_poprzednie << "tau_3 " << tau_3 << endl;
     cout << "tau_4_poprzednie " << tau_4_poprzednie << "tau_4 " << tau_4 << endl;
     cout << "tau_5_poprzednie " << tau_5_poprzednie << "tau_5 " << tau_5 << endl;
     cout << "tau_N_poprzednie " << tau_N_poprzednie << "tau_N " << tau_N << endl;
     cout << "tau_7_poprzednie " << tau_7_poprzednie << "tau_7 " << tau_7 << endl; iter_tau[i]++;

} while (fabs( tau_r_poprzednie-tau_r)/tau_r>eps_tau ||
         fabs( tau_1_poprzednie-tau_1)/tau_1>eps_tau ||
         fabs( tau_p_poprzednie-tau_p)/tau_p>eps_tau ||
         fabs( tau_3_poprzednie-tau_3)/tau_3>eps_tau || 
         fabs( tau_4_poprzednie-tau_4)/tau_4>eps_tau ||
         fabs( tau_5_poprzednie-tau_5)/tau_5>eps_tau ||
         fabs( tau_N_poprzednie-tau_N)/tau_N>eps_tau ||
         fabs( tau_7_poprzednie-tau_7)/tau_7>eps_tau ) ;  
     
     // ************* P O P R A W I O N E    T E M P E R A T U R Y *************************************************//
     T_r[i]=tau_r*T_r[0];
     T_1[i]=tau_1*T_1[0];
     //T_p[i]=tau_p*T_p[0];
     T_3[i]=tau_3*T_3[0];
     T_4[i]=tau_4*T_4[0];
     T_5[i]=tau_5*T_5[0];
     T_N[i]=tau_N*T_N[0];
     T_7[i]=tau_7*T_7[0];  
     
     // ************* G R U P A    8 ******************************************************************************//
     double eps_y_k=0.001;
     SteamState SK0=freesteam_set_ph(p_k[0],2430e3);  cout << "p_k " << p_k[0] << "i_k   " << i_k[0] << endl;                                              
     y_k[0]=freesteam_x(SK0);  cout << "y_k[0] " << y_k[0] << endl;
     
     y_k_p[i]=y_k[0];
     H_8[i]=SPADEK_IZENTROPOWY_i(p_7[i],i_7[i],p_k[i]);
     eps_8[i]=p_k[i]/p_7[i];
     eta_i8_s[i]=1.6954*pow(eps_8[i],2)-1.3486*eps_8[i]+0.9061;
     
     double y_k_poprzednie;
     do
     {
         P_8[i]=m_g8[i]/3.6*H_8[i]*eta_i8_s[i]*y_k_p[i];
         i_k[i]=i_7[i]-P_8[i]/(m_g8[i]/3.6); cout << "i_k[i] " << i_k[i] << endl;
         
         SteamState SK=freesteam_set_ph(p_k[i],i_k[i]);                                                           
         y_k[i]=freesteam_x(SK); cout << "y_k[i] " << y_k[i] << endl; 
         y_k_poprzednie=y_k_p[i]; cout << "y_k_poprzednie " << y_k_poprzednie << endl << endl;
         y_k_p[i]=y_k[i];
     }while (fabs((y_k_poprzednie-y_k[i])/y_k[i])>eps_y_k);
     
     SteamState SK=freesteam_set_ph(p_k[i],i_k[i]);                                                                
     v_k[i]=freesteam_v(SK); cout << "v_k[i] " << v_k[i] << endl;
     T_k[i]=freesteam_T(SK)-273.15; cout << "T_k[i] " << T_k[i] << endl;          
     cout << " P_8 " << P_8[i] << " T_k " << T_k[i] << " i_k " << i_k[i] << " v_k " << v_k[i] << "H_8 " << H_8[i] << endl;                                                                                          
     // ***********************************************************************************************************//
     
     // ***************************** O B L I C Z E N I A    M O C Y **********************************************// 
     double eps_P_el=0.01;
     
     double delta_P_m0=7.5e6; // straty mechaniczne
     P[i]=P_r[i]+P_1[i]+P_2[i]+P_3[i]+P_4[i]+P_5[i]+P_6[i]+P_7[i]+P_8[i]-delta_P_m0; cout << "P[i] " << P[i] << endl;
     
     double P_el_p=P[i];
     double P_el_poprzednie;
     
     do
     {
            delta_P_el[i]=1000*(1247+1163*pow(P_el_p/P_el[0],2)); cout << "delta_P_el " << delta_P_el[i] << endl;
            P_el[i]=P[i]-delta_P_el[i]; cout << "P_el " << P_el[i] << endl;
            
            P_el_poprzednie=P_el_p; cout << "P_el_poprzednie " << P_el_poprzednie << endl;
            P_el_p=P_el[i];
     }while( fabs( (P_el_poprzednie-P_el[i])/P_el[i])>eps_P_el);
     
     } // koniec pêtli po punktach pracy
     
     /*################ ZAPIS DO PLIKU ZEWNÊTRZNEGO ###########################*/
     /*
     st_reg << "m_0 " << "p_0 " << "p_zh " << "p_r " << "p_z1 " << "p_z2 " << "p_z3 " << "p_z4 "  
     << "m_z1 " << "m_z2 " << "m_z3 " << "m_z4 " << "gx " << "P_r " << endl;
     
     for (int j=K;j>0;j--)
     {
         st_reg.setf(ios::fixed);
         st_reg << m_0[j] << " ";
         st_reg << p_0[j] << " ";
         st_reg << p_zh[j] << " ";
         st_reg << p_r[j] << " ";
         st_reg << p_z1[j] << " ";
         st_reg << p_z2[j] << " ";
         st_reg << p_z3[j] << " ";
         st_reg << p_z4[j] << " ";
         st_reg << m_z1[j] << " ";
         st_reg << m_z2[j] << " ";
         st_reg << m_z3[j] << " ";
         st_reg << m_z4[j] << " ";
         st_reg << gx[j] << " ";
         st_reg << P_r[j]  << " ";
         st_reg << opcja[j] << endl;  
 
     }
     
    osiagi << setprecision(2) << fixed;
     
     osiagi << "p_r "  << "p_1 "  << "p_2 "  << "p_p "  << "p_3 "  << "p_4 "  << "p_5 "  << "p_6 "  << "p_7 "  << "p_k "  <<
               "T_r "  << "T_1 "  << "T_2 "  << "T_p "  << "T_3 "  << "T_4 "  << "T_5 "  << "T_6 "  << "T_7 "  << "T_k "  <<
               "i_0 " <<"i_1 "  << "i_2 "  << "i_p "  << "i_3 "  << "i_4 "  << "i_5 "  << "i_6 "  << "i_7 "  << "i_k "  <<
                          "m_g1 " << "m_g2 "            << "m_g3 " << "m_g4 " << "m_g5 " << "m_g6 " << "m_g7 " << "m_g8 " <<
                          "m_u1 " << "m_u2 "            << "m_u3 " << "m_u4 " << "m_u5 " << "m_u6 " << "m_u7 " << 
                          "m_wt " << "m_przec " <<
  "P_r0 " <<"P_ru " << "delta_P_k " << "delta_P_tw "   <<"P_r " << "P_1 "  << "P_2 "             << "P_3 "  << "P_4 "  << "P_5 "  << "P_6 "  << "P_7 "  << "P_8 "  << "P-deltaP_m " << "delta_P_el " << "P_el " << 
                "x_7 " << "x_k " << "eps " << "eta_8 " << "iter_tau " << " jzc" << endl;    
               
     for (int j=K;j>0;j--)
     {
         //osiagi.setf(ios::fixed);
         osiagi << p_r[j]/1e6   << " "; // MPa
         osiagi << p_1[j]/1e6   << " ";
         osiagi << p_2[j]/1e6   << " ";
         osiagi << p_p[j]/1e6   << " ";
         osiagi << p_3[j]/1e6   << " ";
         osiagi << p_4[j]/1e6   << " ";
         osiagi << p_5[j]/1e6   << " ";
         osiagi << p_6[j]/1e3   << " ";
         osiagi << p_7[j]/1e3   << " "; // kPa
         osiagi << p_k[j]/1e3   << " ";
         osiagi << T_r[j] << " ";
         osiagi << T_1[j] << " ";
         osiagi << T_2[j] << " ";
         osiagi << T_p[j] << " ";
         osiagi << T_3[j] << " ";
         osiagi << T_4[j] << " ";
         osiagi << T_5[j] << " ";
         osiagi << T_6[j] << " ";
         osiagi << T_7[j] << " ";
         osiagi << T_k[j] << " ";
         osiagi << i_0[j]/1e3 << " ";
         osiagi << i_1[j]/1e3 << " "; // kJ/kg
         osiagi << i_2[j]/1e3 << " ";
         osiagi << i_p[j]/1e3 << " ";
         osiagi << i_3[j]/1e3 << " ";
         osiagi << i_4[j]/1e3 << " ";
         osiagi << i_5[j]/1e3 << " ";
         osiagi << i_6[j]/1e3 << " ";
         osiagi << i_7[j]/1e3 << " ";
         osiagi << i_k[j]/1e3 << " ";
         osiagi << m_g1[j] << " ";
         osiagi << m_g2[j] << " ";
         osiagi << m_g3[j] << " ";
         osiagi << m_g4[j] << " ";
         osiagi << m_g5[j] << " ";
         osiagi << m_g6[j] << " ";
         osiagi << m_g7[j] << " ";
         osiagi << m_g8[j] << " ";
         osiagi << m_u1[j] << " ";
         osiagi << m_u2[j] << " ";
         osiagi << m_u3[j] << " ";
         osiagi << m_u4[j] << " ";
         osiagi << m_u5[j] << " ";
         osiagi << m_u6[j] << " ";
         osiagi << m_u7[j] << " ";
         osiagi << m_wt[j] << " ";
         osiagi << m_przec[j] << " ";
         osiagi << (m_z1[j]/3.6*H_z1[j] + m_z2[j]/3.6*H_z2[j] + m_z3[j]/3.6*H_z3[j] + m_z4[j]/3.6*H_z4[j])/1e6 << " "; // ca³a moc oddana do SR
         osiagi << P_ru[j]/1e6 << " "; // moc wew
         osiagi << delta_P_k[j]/1e6 << " "; // straty
         osiagi << delta_P_tw[j]/1e6 << " "; // straty
         osiagi << P_r[j]/1e6 << " "; // MW
         osiagi << P_1[j]/1e6  << " ";
         osiagi << P_2[j]/1e6  << " ";
         osiagi << P_3[j]/1e6  << " ";
         osiagi << P_4[j]/1e6  << " ";
         osiagi << P_5[j]/1e6  << " ";
         osiagi << P_6[j]/1e6  << " ";
         osiagi << P_7[j]/1e6  << " ";
         osiagi << P_8[j]/1e6  << " ";
         osiagi << P[j]/1e6  << " ";
         osiagi << delta_P_el[j]/1e6  << " ";
         osiagi << P_el[j]/1e6  << " ";
         osiagi << y_7[j] << " ";
         osiagi << y_k[j] << " ";
         osiagi << eps_8[j]  << " ";
         osiagi << eta_i8_s[j]*y_k[j]<< " "; 
         osiagi << iter_tau[j] << " "; 
         osiagi << (
         m_0[j]/3.6*(i_0[j]/1e3-(-0.0007*pow(m_0[j],2)+1.181*m_0[j]+612.49))
         +(m_g2[j]-m_u2[j])/3.6*(i_p[j]/1e3-i_2[j]/1e3)
         +m_wt[j]/3.6*(i_p[j]/1e3-(-0.0007*pow(m_0[j],2)+1.181*m_0[j]+612.49))
         +m_przec[j]/3.6*(i_2[j]/1e3-(-0.0007*pow(m_0[j],2)+1.181*m_0[j]+612.49))
         )*3600/(P_el[j]/1e3) << endl;

     }     
     
     
    
	st_reg.close();
	spr_reg.close();
	osiagi.close(); 
	
	
	
	*/
	
	
	
    ofstream testing;
	testing.open("wyjscie_testing.txt");
    for (int j=1;j<=K;j++)
     {
         //osiagi.setf(ios::fixed);
         testing << m_g1[j] << " ";
         testing << P_el[j]/1e6  << endl;
     } 
    //cin.get(); 
	return 0;



}



