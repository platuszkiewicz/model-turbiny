#include <cstdlib>
#include <iostream>
#include<fstream>
#include<math.h>

using namespace std;

double d = 0;
double det(int n, double mat[6][6])
{
    int c, subi, i, j, subj;
    double submat[6][6];  
    if (n == 2) 
    {
        return( (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    }
    else
    {  
        for(c = 0; c < n; c++)
        {  
            subi = 0;  
            for(i = 1; i < n; i++)
            {  
                subj = 0;
                for(j = 0; j < n; j++)
                {    
                    if (j == c)
                    {
                        continue;
                    }
                    submat[subi][subj] = mat[i][j];
                    subj++;
                }
                subi++;
            }
        d = d + (pow(-1 ,c) * mat[0][c] * det(n - 1 ,submat));

        }
    }
    return d;
}

int main(int argc, char *argv[])
{
    ifstream wejscia;
	wejscia.open("wej.txt"); // wszystkich wejsc jest 3855
	int N=3;
	double m_0[N];
	double T_0[N];
	double p_0[N];
	double T_p[N];
	double p_k[N];
	
	cout << "N = " << N << endl;cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
	
	for (int i=0; i<=N-1; i++)
	{
        wejscia >> m_0[i];
        wejscia >> T_0[i];
        wejscia >> p_0[i];
        wejscia >> T_p[i]; 
        wejscia >> p_k[i];        
    }
   
    double Z[6][N]; // wype³nianie macierzy Z[l. wierszy][l.kolumn]
    for (int k=0; k<N; k++)
    {
                Z[0][k]=m_0[k];
                Z[1][k]=T_0[k];
                Z[2][k]=p_0[k];
                Z[3][k]=T_p[k];
                Z[4][k]=p_k[k];
                Z[5][k]=1;
    }
    
    cout << "MACIERZ Z " << endl;
    for (int w=0; w<6; w++)
    {
        for (int k=0; k<N; k++)
        { 
            cout << Z[w][k] << " " ;
        }
        cout << endl;
    }cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
    
    double Z_T[N][6]; // wype³nienie macierzy Z transponowanej
    for (int w=0; w<N; w++)
    {
                Z_T[w][0]=Z[0][w];
                Z_T[w][1]=Z[1][w];
                Z_T[w][2]=Z[2][w];
                Z_T[w][3]=Z[3][w];
                Z_T[w][4]=Z[4][w];
                Z_T[w][5]=Z[5][w];
    }
    
    cout << "MACIERZ Z_T " << endl;
    for (int w=0; w<N; w++)
    {
        for (int k=0; k<6; k++)
        { 
            cout << Z_T[w][k] << " " ;
        }
        cout << endl;
    }cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
    
    double ZZ_T[6][6]; // mno¿enie macierzy Z oraz Z_T
    for (int w=0; w<6; w++)
    {
        for (int k=0; k<6; k++)
        {
            ZZ_T[w][k]=0;
            for (int i=0; i<N; i++)
            { ZZ_T[w][k]=ZZ_T[w][k]+Z[w][i]*Z_T[i][k];}
        }
    }
    
     ofstream kalk;
	kalk.open("kalk.txt");
    cout << "MACIERZ Z*Z_T " << endl;
    for (int w=0; w<6; w++)
    {
        for (int k=0; k<6; k++)
        { 
            cout << ZZ_T[w][k] << " " ;
            kalk << ZZ_T[w][k] << " " ;
        }
        cout << endl;
        kalk << endl;
    }cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
    
    cout << det(3,ZZ_T);
    system("PAUSE");
    return EXIT_SUCCESS;
}
