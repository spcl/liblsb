/*
 * liblsb_fortran.cpp
 *
 *  Created on: Dec 18, 2011
 *      Author: htor
 */

#include "liblsb_internal.hpp"

// TODO: should be auto-generated like in libnbc!
#define F77_FUNC_(name,NAME) name ## _

extern "C" { void F77_FUNC_(lsb_init,LSB_INIT)(const char *projname, int *length, int *autoprof_interval); }
void F77_FUNC_(lsb_init,LSB_INIT)(const char *projname, int *length, int *autoprof_interval) {

  char *c=(char*)calloc((*length+1)*sizeof(char), 1);
  memcpy(c, projname, *length);

  LSB_Init(c, *autoprof_interval);

  free(c);
}

extern "C" { void F77_FUNC_(lsb_rec,LSB_REC)(unsigned int *id); }
void F77_FUNC_(lsb_rec,LSB_REC)(unsigned int *id) {
   LSB_Rec(*id);
}

/*extern "C" { void F77_FUNC_(lsb_rec_ints,LSB_REC_INTS)(unsigned int *id, int *int1, int *int2); }
void F77_FUNC_(lsb_rec_ints,LSB_REC_INTS)(unsigned int *id, int *int1, int *int2) {
   LSB_Rec_ints(*id, *int1, *int2);
}

extern "C" { void F77_FUNC_(lsb_rec_intdbl,LSB_REC_INTDBL)(unsigned int *id, int *int1, double *dbl); }
void F77_FUNC_(lsb_rec_intdbl,LSB_REC_INTDBL)(unsigned int *id, int *int1, double *dbl) {
   LSB_Rec_intdbl(*id, *int1, *dbl);
}*/


extern "C" { void F77_FUNC_(lsb_reg_param_int,LSB_REG_PARAM_INT)(char *name, int *length, int *value); }
void F77_FUNC_(lsb_reg_param_int,LSB_REG_PARAM_INT)(char *name, int *length, int *value) {
  char *c=(char*)calloc((*length+1)*sizeof(char), 1);
  memcpy(c, name, *length);

  LSB_Reg_param("%s = %i", c, *value);

  free(c);

}

extern "C" { void F77_FUNC_(lsb_reg_id,LSB_REG_ID)(char *name, int *length, int *value); }
void F77_FUNC_(lsb_reg_id,LSB_REG_ID)(char *name, int *length, int *value) {
  char *c=(char*)calloc((*length+1)*sizeof(char), 1);
  memcpy(c, name, *length);

  LSB_Reg_id("%s = %i", c, *value);

  free(c);

}


extern "C" { void F77_FUNC_(lsb_set_rparam_int, LSB_SET_RPARAM_INT)(const char **index, int *value); }
void F77_FUNC_(lsb_set_rparam, LSB_SET_RPARAM_INT)(const char **index, int *value) {
  LSB_Set_Rparam_int(*index, *value);
}  

extern "C" { void F77_FUNC_(lsb_set_rparam_str, LSB_SET_RPARAM_STR)(const char **index, const char **value); }
void F77_FUNC_(lsb_set_rparam, LSB_SET_RPARAM_STR)(const char **index, const char **value) {
  LSB_Set_Rparam_string(*index, *value);
}  


extern "C" { void F77_FUNC_(lsb_flush,LSB_FLUSH)(void); }
void F77_FUNC_(lsb_flush,LSB_FLUSH)(void) { LSB_Flush(); }

extern "C" { void F77_FUNC_(lsb_finalize,LSB_FINALIZE)(void); }
void F77_FUNC_(lsb_finalize,LSB_FINALIZE)(void) { LSB_Finalize(); }

extern "C" { void F77_FUNC_(lsb_res,LSB_RES)(void); }
void F77_FUNC_(lsb_res,LSB_RES)(void) { LSB_Res(); }

extern "C" { void F77_FUNC_(lsb_rec_enable,LSB_REC_ENABLE)(void); }
void F77_FUNC_(lsb_rec_enable,LSB_REC_ENABLE)(void) { LSB_Rec_enable(); }

extern "C" { void F77_FUNC_(lsb_rec_disable,LSB_REC_DISABLE)(void); }
void F77_FUNC_(lsb_rec_disable,LSB_REC_DISABLE)(void) { LSB_Rec_disable(); }
