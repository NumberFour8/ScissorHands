#include <iostream>
#include <stdio.h>
#include <string>

#include "kernel.h"

using namespace std;

int main()
{
	string cf,X,Y,Z,T,N;
	mpz_t zcf,zX,zY,zZ,zT,zN;

	cout << "N = ?" << endl;
	cin >> N;
	cout << "X = ? " << endl;
	cin >> X;
	cout << "Y = ?" << endl;
	cin >> Y;
	cout << "Z = ?" << endl;
	cin >> Z;
	cout << "T = ?" << endl;
	cin >> T;
	cout << "k = ?" << endl;
	cin >> cf;
	
	mpz_init_set_str(zN,N.c_str(),10);
	mpz_init_set_str(zcf,cf.c_str(),10);
	mpz_init_set_str(zX,X.c_str(),10);
	mpz_init_set_str(zY,Y.c_str(),10);
	mpz_init_set_str(zZ,Z.c_str(),10);
	mpz_init_set_str(zT,T.c_str(),10);
	
	NAF coeffNaf(2,zcf);
	mpz_clear(zcf);
	
	to_mont_repr(zX,zN);
	to_mont_repr(zY,zN);
	to_mont_repr(zZ,zN);
	to_mont_repr(zT,zN);
	
	h_ExtendedPoint pts;
	pts.fromMPZ(zX,zY,zZ,zT);

	h_Aux ax;
	mpz_to_biguint(ax.N,zN);
	// TODO: Přidat výpočet 3N a invN!!!!!

	cudaError_t cudaStatus = computeExtended(ax,&pts,coeffNaf);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
    
	mpz_clears(zX,zY,zZ,zT,zN);
	
	char c;
	cin >> c;
    return 0;
}
