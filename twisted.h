// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
#define twistedAdd() COPY_Q();					\
					 CLEAR_TEMP();				\
					 SUB_MOD(c_tt0,c_y1,c_x1);	\
					 ADD_MOD(c_tt1,c_y2,c_x2);	\
					 MUL_MOD(c_tt0,c_tt0,c_tt1);\
					 ADD_MOD(c_tt1,c_y1,c_x1);	\
					 SUB_MOD(c_x1,c_y2,c_x2);	\
					 MUL_MOD(c_tt1,c_tt1,c_x1); \
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
////////////////////////////////////////////////////
	
// Twisted Edwards Extended (add-2008-hwcd-4), a = -1, independent of d,incomplete
#define twistedSub() COPY_Q();					\
					 CLEAR_TEMP();				\
					 SUB_MOD(c_tt0,c_y1,c_x1);	\
					 SUB_MOD(c_tt1,c_y2,c_x2);	\
					 MUL_MOD(c_tt0,c_tt0,c_tt1);\
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
////////////////////////////////////////////////////
	
// Twisted Edwards Extended (dbl-2008-hwcd-4), a = -1, independent of d,incomplete	
#define twistedDbl() CLEAR_TEMP();				\
					 ADD_MOD(c_tt0,c_x1,c_y1);	\
					 SQR_MOD(c_tt1,c_x1);		\
					 SQR_MOD(c_x1,c_y1);		\
					 SQR_MOD(c_y1,c_z1);		\
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
////////////////////////////////////////////////////
