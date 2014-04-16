// Vycisteni docasnych promennych
#define CLEAR_TEMP() c_tcy[threadIdx.x] = 0;	\
					 c_tt0[threadIdx.x] = 0;	\
					 c_tt1[threadIdx.x] = 0;	\
					 _AUX[threadIdx.x]  = 0; 
 
// Nakopírování pracovních dat bodu Q (celkem 32*4 = 128 cifer) 
#define COPY_Q()	c_x2[threadIdx.x] = *(Qd+threadIdx.x+0*NB_DIGITS);  \
					c_y2[threadIdx.x] = *(Qd+threadIdx.x+1*NB_DIGITS);  \
					c_z2[threadIdx.x] = *(Qd+threadIdx.x+2*NB_DIGITS);  \
					c_t2[threadIdx.x] = *(Qd+threadIdx.x+3*NB_DIGITS); 
	

// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
#define curvesAdd() COPY_Q();					\
					CLEAR_TEMP();				\
					SUB_MOD(c_tt0,c_y1,c_x1);	\
					ADD_MOD(c_tt1,c_y2,c_x2);	\
					MUL_MOD(c_tt0,c_tt0,c_tt1);	\
					ADD_MOD(c_tt1,c_y1,c_x1);	\
					SUB_MOD(c_x1,c_y2,c_x2);	\
					MUL_MOD(c_tt1,c_tt1,c_x1);  \
					DBL_MOD(c_z2);				\
					DBL_MOD(c_t2);				\
					MUL_MOD(c_z1,c_z1,c_t2);	\
					MUL_MOD(c_z2,c_z2,c_t1);	\
					ADD_MOD(c_y2,c_z2,c_z1);	\
					SUB_MOD(c_x2,c_z2,c_z1);	\
					SUB_MOD(c_z2,c_tt1,c_tt0);	\
					ADD_MOD(c_t2,c_tt1,c_tt0);	\
					MUL_MOD(c_x1,c_y2,c_z2);	\
					MUL_MOD(c_y1,c_t2,c_x2);	\
					MUL_MOD(c_t1,c_y2,c_x2);	\
					MUL_MOD(c_z1,c_z2,c_t2);	
	
/////////////////////////////////////////
	
// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
#define curvesSub() COPY_Q();					\
					CLEAR_TEMP();				\
					SUB_MOD(c_tt0,c_y1,c_x1);	\
					SUB_MOD(c_tt1,c_y2,c_x2);	\
					MUL_MOD(c_tt0,c_tt0,c_tt1);	\
					ADD_MOD(c_tt1,c_y1,c_x1);	\
					ADD_MOD(c_x1,c_y2,c_x2);	\
					MUL_MOD(c_tt1,c_tt1,c_x1);	\
					DBL_MOD(c_z2);				\
					DBL_MOD(c_t2);				\
					MUL_MOD(c_z1,c_z1,c_t2);	\
					MUL_MOD(c_z2,c_z2,c_t1);	\
					SUB_MOD(c_y2,c_z2,c_z1);	\
					ADD_MOD(c_x2,c_z2,c_z1);	\
					SUB_MOD(c_z2,c_tt1,c_tt0);	\
					ADD_MOD(c_t2,c_tt1,c_tt0);	\
					MUL_MOD(c_x1,c_y2,c_z2);	\
					MUL_MOD(c_y1,c_t2,c_x2);	\
					MUL_MOD(c_t1,c_y2,c_x2);	\
					MUL_MOD(c_z1,c_z2,c_t2);
	
/////////////////////////////////////////
	
// Twisted Edwards Extended (dbl-2008-hwcd-4), a = -1, independent of d,incomplete	
#define curvesDbl()	CLEAR_TEMP()				\
					ADD_MOD(c_tt0,c_x1,c_y1);	\
					SQR_MOD(c_tt1,c_x1);		\
					SQR_MOD(c_x1,c_y1);			\
					SQR_MOD(c_y1,c_z1);			\
					ADD_MOD(c_t1,c_tt1,c_x1);	\
					SUB_MOD(c_z1,c_tt1,c_x1);	\
					SQR_MOD(c_tt1,c_tt0);		\
					DBL_MOD(c_y1);				\
					SUB_MOD(c_tt0,c_t1,c_tt1);	\
					ADD_MOD(c_tt1,c_y1,c_z1);	\
					MUL_MOD(c_x1,c_tt1,c_tt0);	\
					MUL_MOD(c_y1,c_t1,c_z1);	\
					MUL_MOD(c_t1,c_t1,c_tt0);	\
					MUL_MOD(c_z1,c_z1,c_tt1);

/////////////////////////////////////////