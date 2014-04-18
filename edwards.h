// Edwards Extended (add-2008-hwcd-2), a = 1, independent of d,incomplete
#define edwardsAdd() COPY_Q();					\
					 CLEAR_TEMP();				\
					 MUL_MOD(c_tt0,c_x1,c_x2);	\
					 MUL_MOD(c_tt1,c_y1,c_y2);	\
					 MUL_MOD(c_t2,c_t2,c_z1);	\
					 MUL_MOD(c_t1,c_t1,c_z2);	\
					 ADD_MOD(c_z1,c_x2,c_y2);	\
					 SUB_MOD(c_z2,c_x1,c_y1);	\
					 ADD_MOD(c_x2,c_t1,c_t2);	\
					 SUB_MOD(c_y2,c_t1,c_t2);	\
					 MUL_MOD(c_z2,c_z2,c_z1);	\
					 ADD_MOD(c_z1,c_z2,c_tt1);	\
					 SUB_MOD(c_t1,c_z1,c_tt0);	\
					 ADD_MOD(c_z2,c_tt1,c_tt0);	\
					 MUL_MOD(c_x1,c_x2,c_t1);	\
					 MUL_MOD(c_y1,c_z2,c_y2);	\
					 MUL_MOD(c_z1,c_t1,c_z2);	\
					 MUL_MOD(c_t1,c_x2,c_y2);   

////////////////////////////////////////////////////

// Edwards Extended (add-2008-hwcd-2), a = 1, independent of d,incomplete
#define edwardsSub() COPY_Q();					\
					 CLEAR_TEMP();				\
					 MUL_MOD(c_tt0,c_x1,c_x2);	\
					 MUL_MOD(c_tt1,c_y1,c_y2);	\
					 MUL_MOD(c_t2,c_z1,c_t2);	\
					 MUL_MOD(c_t1,c_t1,c_z2);	\
					 SUB_MOD(c_z1,c_y2,c_x2); 	\
					 SUB_MOD(c_z2,c_x1,c_y1);	\
					 SUB_MOD(c_x2,c_t1,c_t2); 	\
					 ADD_MOD(c_y2,c_t1,c_t2); 	\
					 MUL_MOD(c_z2,c_z2,c_z1);	\
					 ADD_MOD(c_z1,c_z2,c_tt1);	\
					 ADD_MOD(c_t1,c_z1,c_tt0); 	\
					 SUB_MOD(c_z2,c_tt1,c_tt0);	\
					 MUL_MOD(c_x1,c_x2,c_t1); 	\
					 MUL_MOD(c_y1,c_z2,c_y2);	\
					 MUL_MOD(c_z1,c_t1,c_z2);	\
					 MUL_MOD(c_t1,c_x2,c_y2);	
////////////////////////////////////////////////////


// Edwards Extended (dbl-2008-hwcd), a = 1, independent of d,incomplete
#define edwardsDbl() CLEAR_TEMP();				\
					 SQR_MOD(c_tt0,c_y1);		\
					 SQR_MOD(c_tt1,c_x1);		\
					 ADD_MOD(c_t1,c_x1,c_y1);	\
					 SQR_MOD(c_t1,c_t1);		\
					 SUB_MOD(c_y1,c_t1,c_tt1);	\
					 SQR_MOD(c_z1,c_z1);		\
					 DBE_MOD(c_t1,c_z1);		\
					 SUB_MOD(c_x1,c_y1,c_tt0);	\
					 ADD_MOD(c_z1,c_tt1,c_tt0);	\
					 SUB_MOD(c_y1,c_tt1,c_tt0);	\
					 SUB_MOD(c_tt0,c_z1,c_t1);	\
					 MUL_MOD(c_t1,c_x1,c_y1);	\
					 MUL_MOD(c_x1,c_x1,c_tt0);	\
					 MUL_MOD(c_y1,c_y1,c_z1);	\
					 MUL_MOD(c_z1,c_z1,c_tt0);
////////////////////////////////////////////////////
