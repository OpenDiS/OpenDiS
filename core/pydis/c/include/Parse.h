#ifndef _PARSE_H
#define _PARSE_H
/****************************************************************************
 *
 *      Module:      Parse.h  
 *      Description: Contains definitions, structures and prototypes
 *                   for code used in parsing the parameter files.
 *
 ***************************************************************************/
#include "Home.h"

/*
 *      Define a set values that may be returned by
 *      functions used in parsing the user supplied
 *      values in the control file.
 */
#define TOKEN_ERR            -1
#define TOKEN_NULL            0
#define TOKEN_GENERIC         1
#define TOKEN_EQUAL           2
#define TOKEN_BEGIN_VAL_LIST  3
#define TOKEN_END_VAL_LIST    4

/*
 *      Define the variable types that may be associated with input
 *      parameters
 */
#define V_NULL    0
#define V_DBL     1
#define V_INT     2
#define V_STRING  3
#define V_COMMENT 4

#define VFLAG_NULL        0x00
#define VFLAG_ALIAS       0x01
#define VFLAG_INITIALIZED 0x02
#define VFLAG_DISABLED    0x04
#define VFLAG_SET_BY_USER 0x08


#ifdef __cplusplus
extern "C" void BindVar(ParamList_t *list, const char *name, void *addr,
                    int type, int cnt, int flags);
extern "C" int  GetNextToken(FILE *fp, char *token, int maxTokenSize);
extern "C" int  GetParamVals(FILE *fp, int valType, int valsExpected,
                    void *valList);
extern "C" int  LookupParam(ParamList_t *list, char *token);
extern "C" void WriteParam(ParamList_t *list, int index, FILE *fp);
#else
void BindVar(ParamList_t *list, const char *name, void *addr, int type,
         int cnt, int flags);
void DisableUnneededParams(Home_t *home);
int  GetNextToken(FILE *fp, char *token, int maxTokenSize);
int  GetParamVals(FILE *fp, int valType, int valsExpected,
        void *valList);
int  LookupParam(ParamList_t *list, char *token);
void MarkParamDisabled(ParamList_t *CPList, char *name);
void MarkParamEnabled(ParamList_t *CPList, char *name);
void WriteParam(ParamList_t *list, int index, FILE *fp);
#endif

#endif  /* _PARSE_H */
