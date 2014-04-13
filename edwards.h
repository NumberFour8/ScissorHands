#define curvesAdd edwardsAdd
#define curvesDbl edwardsDbl
#define curvesSub edwardsSub

__global__ void edwardsAdd(void* R, void *P, void *Q,void* aux)
{
	Aux *ax = (Aux*)aux;
	
	// Prom�nn� ve sd�len� pam�ti pro bod P
    __shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];

	__shared__ VOL digit_t x2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Pomocn� prom�nn� ve sd�len� pam�ti pro p�enos a t0,t1,t2
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// V�sek pro konkr�tn� k�ivku
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];

	VOL digit_t* c_x2 = x2[threadIdx.y];
	VOL digit_t* c_y2 = y2[threadIdx.y];
	VOL digit_t* c_z2 = z2[threadIdx.y];
	VOL digit_t* c_t2 = t2[threadIdx.y];
		
	// Pomocn� prom�nn� a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	VOL carry_t* _CARRY = carry[threadIdx.y];  // p�enos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocn� prom�nn� pro n�soben�
	
	const digit_t _N    = ax->N[threadIdx.x];	// x-t� cifra N
	const digit_t _3N   = ax->N3[threadIdx.x];  // x-t� cifra 3*N
	const digit_t _INVN = ax->invN;				// -N^(-1) mod W
	
	// Na��t�n� dat (4 sou�adnice po MAX_BYTES bajtech)
	const digit_t idx = 4*NB_DIGITS*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Te� m��eme p�e��st spr�vn� bod P
	VOL digit_t* Qd   = ((digit_t*)Q)+idx; // Te� m��eme p�e��st spr�vn� bod P


	// Nakop�rov�n� pracovn�ch dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x+0*NB_DIGITS); // prvn�ch 32 cifer pat�� k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+1*NB_DIGITS); // dal��ch 32 cifer pat�� k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+2*NB_DIGITS); // dal��ch 32 k sou�adnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+3*NB_DIGITS); // ... a posledn� k sou�adnici T

	c_x2[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvn�ch 32 cifer pat�� k X
	c_y2[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dal��ch 32 cifer pat�� k Y
	c_z2[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dal��ch 32 k sou�adnici Z
	c_t2[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a posledn� k sou�adnici T

	c_tcy[threadIdx.x] = 0;
	c_tt0[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 
	_AUX[threadIdx.x]  = 0; 

	// Edwards Extended (add-2008-hwcd-2), a = 1, independent of d,incomplete
	/////////////////////////////////////////	
	
	MUL_MOD(c_tt0,c_x1,c_x2);
	MUL_MOD(c_tt1,c_y1,c_y2);
	
	MUL_MOD(c_t2,c_t2,c_z1);
	MUL_MOD(c_t1,c_t1,c_z2);
	
	ADD_MOD(c_z1,c_x2,c_y2);
	SUB_MOD(c_z2,c_x1,c_y1);
	
	ADD_MOD(c_x2,c_t1,c_t2);
	ADD_MOD(c_y2,c_t1,c_t2);
	
	MUL_MOD(c_z2,c_z2,c_z1);
	ADD_MOD(c_z1,c_z2,c_tt1);
	
	SUB_MOD(c_t1,c_z1,c_tt0);
	ADD_MOD(c_z2,c_tt1,c_tt0);
	
	MUL_MOD(c_x1,c_x2,c_t1);
	MUL_MOD(c_y1,c_z2,c_y2);
	MUL_MOD(c_z1,c_t1,c_z2);
	MUL_MOD(c_t1,c_x2,c_y2);
	
	/////////////////////////////////////////
	
	VOL digit_t* Rd   = ((digit_t*)R)+idx; // Te� m��eme p�e��st spr�vn� bod R

	// Nakop�rov�n� pracovn�ch dat zp�tky
	*(Rd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvn�ch 32 cifer pat�� k X
	*(Rd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dal��ch 32 cifer pat�� k Y
	*(Rd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dal��ch 32 k sou�adnici Z
	*(Rd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a posledn� k sou�adnici T

	__syncthreads();
}

__global__ void edwardsSub(void* R, void *P, void *Q,void* aux)
{
	Aux *ax = (Aux*)aux;
	
	// Prom�nn� ve sd�len� pam�ti pro bod P
    __shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];

	__shared__ VOL digit_t x2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z2[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Pomocn� prom�nn� ve sd�len� pam�ti pro p�enos a t0,t1,t2
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// V�sek pro konkr�tn� k�ivku
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];

	VOL digit_t* c_x2 = x2[threadIdx.y];
	VOL digit_t* c_y2 = y2[threadIdx.y];
	VOL digit_t* c_z2 = z2[threadIdx.y];
	VOL digit_t* c_t2 = t2[threadIdx.y];
		
	// Pomocn� prom�nn� a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	VOL carry_t* _CARRY = carry[threadIdx.y];  // p�enos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocn� prom�nn� pro n�soben�
	
	const digit_t _N    = ax->N[threadIdx.x];	// x-t� cifra N
	const digit_t _3N   = ax->N3[threadIdx.x];  // x-t� cifra 3*N
	const digit_t _INVN = ax->invN;				// -N^(-1) mod W
	
	// Na��t�n� dat (4 sou�adnice po MAX_BYTES bajtech)
	const digit_t idx = 4*NB_DIGITS*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Te� m��eme p�e��st spr�vn� bod P
	VOL digit_t* Qd   = ((digit_t*)Q)+idx; // Te� m��eme p�e��st spr�vn� bod P


	// Nakop�rov�n� pracovn�ch dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x+0*NB_DIGITS); // prvn�ch 32 cifer pat�� k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+1*NB_DIGITS); // dal��ch 32 cifer pat�� k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+2*NB_DIGITS); // dal��ch 32 k sou�adnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+3*NB_DIGITS); // ... a posledn� k sou�adnici T

	c_x2[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS); // prvn�ch 32 cifer pat�� k X
	c_y2[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS); // dal��ch 32 cifer pat�� k Y
	c_z2[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS); // dal��ch 32 k sou�adnici Z
	c_t2[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); // ... a posledn� k sou�adnici T

	c_tcy[threadIdx.x] = 0;
	c_tt0[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 
	_AUX[threadIdx.x]  = 0; 

	// Edwards Extended (add-2008-hwcd-2), a = 1, independent of d,incomplete
	/////////////////////////////////////////	
	
	MUL_MOD(c_tt0,c_x1,c_x2);
	MUL_MOD(c_tt1,c_y1,c_y2);
	
	MUL_MOD(c_t2,c_z1,c_t2);
	MUL_MOD(c_t1,c_t1,c_z2);
	
	SUB_MOD(c_z1,c_y2,c_x2); 
	SUB_MOD(c_z2,c_x1,c_y1);
	
	SUB_MOD(c_x2,c_t1,c_t2); 
	ADD_MOD(c_y2,c_t1,c_t2); 
	
	MUL_MOD(c_z2,c_z2,c_z1);
	ADD_MOD(c_z1,c_z2,c_tt1);
	
	ADD_MOD(c_t1,c_z1,c_tt0); 
	SUB_MOD(c_z2,c_tt1,c_tt0);
	
	MUL_MOD(c_x1,c_x2,c_t1); 
	MUL_MOD(c_y1,c_z2,c_y2);
	MUL_MOD(c_z1,c_t1,c_z2);
	MUL_MOD(c_t1,c_x2,c_y2);	
	
	/////////////////////////////////////////
	
	VOL digit_t* Rd   = ((digit_t*)R)+idx; // Te� m��eme p�e��st spr�vn� bod R

	// Nakop�rov�n� pracovn�ch dat zp�tky
	*(Rd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvn�ch 32 cifer pat�� k X
	*(Rd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dal��ch 32 cifer pat�� k Y
	*(Rd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dal��ch 32 k sou�adnici Z
	*(Rd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a posledn� k sou�adnici T

	__syncthreads();
}

__global__ void edwardsDbl(void* R,void* P,void* aux)
{
	Aux *ax = (Aux*)aux;
	
	// Prom�nn� ve sd�len� pam�ti pro bod P
    __shared__ VOL digit_t x1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t y1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t z1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t t1[CURVES_PER_BLOCK][NB_DIGITS];
	
	// Pomocn� prom�nn� ve sd�len� pam�ti pro p�enos a t0,t1,t2
	__shared__ VOL carry_t carry[CURVES_PER_BLOCK][NB_DIGITS]; 
	__shared__ VOL digit_t temp0[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp1[CURVES_PER_BLOCK][NB_DIGITS];
	__shared__ VOL digit_t temp2[CURVES_PER_BLOCK][NB_DIGITS];
	
	// V�sek pro konkr�tn� k�ivku
	VOL digit_t* c_x1 = x1[threadIdx.y];
	VOL digit_t* c_y1 = y1[threadIdx.y];
	VOL digit_t* c_z1 = z1[threadIdx.y];
	VOL digit_t* c_t1 = t1[threadIdx.y];
		
	// Pomocn� prom�nn� a konstanty
	VOL digit_t* c_tt0  = temp0[threadIdx.y];   // t0
	VOL digit_t* c_tt1  = temp1[threadIdx.y];   // t1
	VOL carry_t* _CARRY = carry[threadIdx.y];  // p�enos
	VOL digit_t* _AUX   = temp2[threadIdx.y];  // pomocn� prom�nn� pro n�soben�
	
	const digit_t _N    = ax->N[threadIdx.x];	// x-t� cifra N
	const digit_t _3N   = ax->N3[threadIdx.x];  // x-t� cifra 3*N
	const digit_t _INVN = ax->invN;				// -N^(-1) mod W
	
	// Na��t�n� dat (4 sou�adnice po MAX_BYTES bajtech)
	const digit_t idx = 4*NB_DIGITS*(blockIdx.x*blockDim.y + threadIdx.y);
	VOL digit_t* Pd   = ((digit_t*)P)+idx; // Te� m��eme p�e��st spr�vn� bod P

	
	// Nakop�rov�n� pracovn�ch dat	(celkem 32*4 = 128 cifer)
	c_x1[threadIdx.x] = *(Pd+threadIdx.x+0*NB_DIGITS); // prvn�ch 32 cifer pat�� k X
	c_y1[threadIdx.x] = *(Pd+threadIdx.x+1*NB_DIGITS); // dal��ch 32 cifer pat�� k Y
	c_z1[threadIdx.x] = *(Pd+threadIdx.x+2*NB_DIGITS); // dal��ch 32 k sou�adnici Z
	c_t1[threadIdx.x] = *(Pd+threadIdx.x+3*NB_DIGITS); // ... a posledn� k sou�adnici T

	c_tcy[threadIdx.x] = 0;
	c_tt0[threadIdx.x] = 0;
	c_tt1[threadIdx.x] = 0; 
	_AUX[threadIdx.x]  = 0; 
 
	// Edwards Extended (dbl-2008-hwcd), a = 1, independent of d,incomplete
	/////////////////////////////////////////
	
	SQR_MOD(c_tt0,c_y1);
	SQR_MOD(c_tt1,c_x1);
	
	ADD_MOD(c_t1,c_x1,c_y1);
	SQR_MOD(c_t1,c_t1);
	
	SUB_MOD(c_y1,c_t1,c_tt1);
	SQR_MOD(c_z1,c_z1);
	
	DBE_MOD(c_t1,c_z1);
	SUB_MOD(c_x1,c_y1,c_tt0);
	
	ADD_MOD(c_z1,c_tt1,c_tt0);
	SUB_MOD(c_y1,c_tt1,c_tt0);
	
	SUB_MOD(c_tt0,c_z1,c_t1);
	MUL_MOD(c_t1,c_x1,c_y1);
	MUL_MOD(c_x1,c_x1,c_tt0);
	MUL_MOD(c_y1,c_y1,c_z1);
	MUL_MOD(c_z1,c_z1,c_tt0);
	
	////////////////////////////////////////

	VOL digit_t* Rd   = ((digit_t*)R)+idx; // Te� m��eme p�e��st spr�vn� bod R

	// Nakop�rov�n� pracovn�ch dat zp�tky
	*(Rd+threadIdx.x+0*NB_DIGITS) = c_x1[threadIdx.x];  // prvn�ch 32 cifer pat�� k X
	*(Rd+threadIdx.x+1*NB_DIGITS) = c_y1[threadIdx.x];  // dal��ch 32 cifer pat�� k Y
	*(Rd+threadIdx.x+2*NB_DIGITS) = c_z1[threadIdx.x];  // dal��ch 32 k sou�adnici Z
	*(Rd+threadIdx.x+3*NB_DIGITS) = c_t1[threadIdx.x];  // ... a posledn� k sou�adnici T

	__syncthreads();
}
