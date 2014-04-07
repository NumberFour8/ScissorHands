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
	
	// Koeficient převést do NAF
	mpz_init_str(zN,N.c_str(),10);
	mpz_init_str(zcf,cf.c_str(),10);
	mpz_init_str(zX,X.c_str(),10);
	mpz_init_str(zY,Y.c_str(),10);
	mpz_init_str(zZ,Z.c_str(),10);
	mpz_init_str(zT,T.c_str(),10);
	
	NAF coeffNaf(2,zcf);
	mpz_clear(zcf);
	
	to_mont_repr(zX,zN);
	to_mont_repr(zY,zN);
	to_mont_repr(zZ,zN);
	to_mont_repr(zT,zN);
	
	h_ExtendedPoint pts;
	pts.fromMPZ(zX,zY,zZ,zT);

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
    
	mpz_clear(zX);
	mpz_clear(zY);
	mpz_clear(zZ);
	mpz_clear(zT);
	mpz_clear(zN);

	char c;
	cin >> c;
    return 0;
}
