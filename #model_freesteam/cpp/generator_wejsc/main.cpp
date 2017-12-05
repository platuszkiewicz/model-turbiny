/* Plik zwraca K wejœæ do modelu w postaci pliku .txt o K linijkach. Dane wystêpuj¹ w kolejnoœci                       
   m_0, T_0, p_0, T_p, p_k, m_u1, m_u2, m_u3, m_u4, m_u5, m_5, m_u6, m_u7                               
   Nale¿y ustawiæ: 
   -------------iloœæ punktów pracy (K)
   -------------wartoœæ pocz¹tkow¹ (np. m_0_start)
   -------------interwa³ (np. m_0_delta)    
   Interwa³ dodatni powoduje zwiêkszanie wielkoœci, a ujemny zmniejszanie                               */

#include <cstdlib>  
#include <iostream> 
#include <fstream>
#include <math.h>
 
using namespace std;

int main(int argc, char *argv[])
{
    ofstream zapis;
    zapis.open("D:\\##STUDIA\\+III sem mgr\\#Praca magisterska\\#model_freesteam\\cpp\\1\\wejscie.txt"); 
    
    int K=341; // iloœæ punktów pracy
    double m_0[K-1],T_0[K-1],p_0[K-1], T_p[K-1], p_k[K-1];
    double m_u1[K-1], m_u2[K-1], m_u3[K-1], m_u4[K-1], m_u5[K-1], m_u6[K-1], m_u7[K-1];
    double m_wt[K-1], m_przec[K-1];
    
    // ********************** TEN FRAGMENT SIÊ MODYFIKUJE ******** //
    double m_0_start=655, T_0_start=538, p_0_start=12.8e6, T_p_start=535, p_k_start=4.75e3;
    double m_u1_start=0.0616*m_0_start-7.6766, m_u2_start=0.0664*m_0_start+14.900, m_u3_start=0.0876*m_0_start-31.749,
     m_u4_start=0.0292*m_0_start-3.7444, m_u5_start=0.0524*m_0_start-4.4238, m_u6_start=0.0471*m_0_start-2.4041, m_u7_start=0.0173*m_0_start-5.7293;
    
    double m_0_delta=-1, T_0_delta=0, p_0_delta=0, T_p_delta=0, p_k_delta=0;
    
    // ********************** KONIEC MODYFIKACJI ***************** //
    
    m_0[0]=m_0_start; T_0[0]=T_0_start; p_0[0]=p_0_start; T_p[0]=T_p_start; p_k[0]=p_k_start;
    m_u1[0]=m_u1_start; m_u2[0]=m_u2_start; m_u3[0]=m_u3_start; m_u4[0]=m_u4_start; m_u5[0]=m_u5_start; m_u6[0]=m_u6_start; m_u7[0]=m_u7_start;
                        m_wt[0]=0.0871*m_0[0]-23.5192; 
       	                m_przec[0]=0.0140*m_0[0]+0.1044; 

    for (int i=1; i<=K-1; i++)  // okreœlenie parametrów wejœciowych
    {
        m_0[i]=m_0[i-1]+m_0_delta; 
        T_0[i]=T_0[i-1]+T_0_delta;
        p_0[i]=p_0[i-1]+p_0_delta;
        T_p[i]=T_p[i-1]+T_p_delta;
        p_k[i]=p_k[i-1]+p_k_delta;
        
       /* m_u1[i]=m_u1[i-1]+m_u1_delta; 
        m_u2[i]=m_u2[i-1]+m_u2_delta; 
        m_u3[i]=m_u3[i-1]+m_u3_delta; 
        m_u4[i]=m_u4[i-1]+m_u4_delta; 
        m_u5[i]=m_u5[i-1]+m_u5_delta; 
        m_u6[i]=m_u6[i-1]+m_u6_delta; 
        m_u7[i]=m_u7[i-1]+m_u7_delta; */
        
        m_u1[i]=0.059*m_0[i]-6.891;
        m_u2[i]=0.066*m_0[i]+15.016;
        if(m_0[i]<375){m_u3[i]=0;}else
        m_u3[i]=0.0876*m_0[i]-31.749; // zmodyfikowane (przy zmianie zmieniæ te¿ m_u3_start)
        m_u4[i]=0.026*m_0[i]-2.194;
        m_u5[i]=0.051*m_0[i]-3.159;
        m_u6[i]=0.048*m_0[i]-1.113;
        if(m_0[i]<320){m_u7[i]=0;}else
        m_u7[i]=0.018*m_0[i]-5.729; // zmodyfikowane (przy zmianie zmieniæ te¿ m_u7_start)
        m_wt[i]=0.0871*m_0[i]-23.5192; 
       	m_przec[i]=0.0140*m_0[i]+0.1044; 
    }

    zapis.setf(ios::fixed); // notacja zapisu
    zapis.precision(2);     // precyzja zapisu
   
    
    for (int i=0; i<K; i++) // zapis do pliku
    {  
        zapis << m_0[i] << " "; cout << m_0[i] << " ";
        zapis << T_0[i] << " "; cout << T_0[i] << " ";
        zapis << p_0[i] << " "; cout << p_0[i] << " ";
        zapis << T_p[i] << " "; cout << T_p[i] << " ";
        zapis << p_k[i] << " "; cout << p_k[i] << " ";
        zapis << m_u1[i] << " "; cout << m_u1[i] << " ";
        zapis << m_u2[i] << " "; cout << m_u2[i] << " ";
        zapis << m_u3[i] << " "; cout << m_u3[i] << " ";
        zapis << m_u4[i] << " "; cout << m_u4[i] << " ";
        zapis << m_u5[i] << " "; cout << m_u5[i] << " ";
        zapis << m_u6[i] << " "; cout << m_u6[i] << " ";
        zapis << m_u7[i] << " "; cout << m_u7[i] << " ";
        zapis << m_wt[i] << " "; cout << m_wt[i] << " ";
        zapis << m_przec[i] << endl; cout << m_przec[i] << endl;
    }
        
    system("PAUSE");
    return EXIT_SUCCESS;
}
