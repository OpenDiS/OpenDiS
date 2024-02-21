/**************************************************************************
 *
 *      Module:       Parse.c
 *      Description:  Contains functions used to bind parameter names
 *                    to code variables, parsing files containing
 *                    parameter/value pairs, and managing lists of
 *                    bound parameters.
 *
 *      Public functions:
 *          BindVar()
 *          GetNextToken()
 *          GetParamVals()
 *          LookupParam()
 *          WriteParam()
 *
 *************************************************************************/
#include <stdarg.h>
#include <ctype.h>
#include "Home.h"
#include "Parse.h"


/*---------------------------------------------------------------------------
 *
 *      Function:       BindVar
 *      Description:    Register the identifier of a parameter along with
 *                      the address of its corresponding variable
 *                      and some ancillary info about the parameter.
 *                      Used mainly in the parsing of control files
 *                      containing parameter id/value pairs.
 *
 *      Arguments:
 *          list   pointer to the parameter list to which the
 *                 specified parameter name and info should be added
 *          name   identifier associated with the parameter.
 *                 Note: <name> must be a character string less
 *                 than MAX_STRING_LEN bytes long.
 *          addr   address of the variable associated with the
 *                 parameter
 *          type   specifies the type of variable that is associated
 *                 with the specified parameter
 *          cnt    Indicates the number of values of type <type>
 *                 associated with the parameter.
 *          flags  bit-field containing flags indicating additional
 *                 traits associated with the parameter.
 *
 *
 *      NOTE: For strings (type == V_STRING) if cnt == 1, the
 *            addr is a (char *) pointing to the string.  If
 *            cnt > 1, addr is actually a (char **) pointing to
 *            an array of length cnt of pointers each pointing to a
 *            character string previously allocated by the caller.
 *
 *-------------------------------------------------------------------------*/
void BindVar(ParamList_t *list, const char *name, void *addr, int type,
             int cnt, int flags)
{
        int       pCnt;
        VarData_t *vList;

        vList = list->varList;
        pCnt = list->paramCnt;
        list->paramCnt += 1;

        vList = (VarData_t *)realloc(vList, list->paramCnt *
                                     sizeof(VarData_t));
        
        memset(vList[pCnt].varName, 0, MAX_STRING_LEN);

        strncpy(vList[pCnt].varName, name, MAX_STRING_LEN-1);
        vList[pCnt].valList = addr;
        vList[pCnt].valType = type;
        vList[pCnt].valCnt  = cnt;
        vList[pCnt].flags   = flags;

        list->varList = vList;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       WriteParam
 *      Description:    Writes to the specified file descriptor, one
 *                      or more parameter identifiers and associated
 *                      values.
 *      Arguments:
 *          list    Pointer to the list of parameters
 *          index   index in <list> of the parameter and value to
 *                  be printed (if <index> is less than zero, all
 *                  parameters/values will be printed.
 *          fp      file descriptor to which to write the data.
 *
 *      Last Modified: 04/18/2008 gh  - bug fix to prevent aliases from
 *                                being written along with normal params.
 *
 *-------------------------------------------------------------------------*/
void WriteParam(ParamList_t *list, int index, FILE *fp)
{
        int  i, j, minIndex, maxIndex, type, count;
        char *lineFeed, *lBracket, *rBracket;
        char format[MAX_STRING_LEN];

/*
 *      We'll be looping over all parameters to be printed.  If
 *      a single parameter is selected, set loop indices so we
 *      only iterate once for the specific parameter, otherwise
 *      loop over all parameters.
 */
        if (index >= 0) {
            minIndex = index;
            maxIndex = index+1;
        } else {
            minIndex = 0;
            maxIndex = list->paramCnt;
        }

        for (i = minIndex; i < maxIndex; i++) {

/*
 *          If we're dumping the whole list of parameters, don't
 *          bother dumping out a parameter that is an alias for
 *          for another.
 */
            if (((list->varList[i].flags & VFLAG_ALIAS) != 0) && (index < 0)) {
                continue;
            }

/*
 *          Use of some parameters precludes the use of others, in which
 *          case, unused parameters may have been flagged as irrelevant
 *          for this particular simulation.  If the current parameter has
 *          has been flagged, there's no need to write it out.
 *
 *          However, just as a safety net, if value was explicitly supplied
 *          by the user in the restart file, go ahead and write the value
 *          out anyway.
 */
            if (((list->varList[i].flags & VFLAG_DISABLED) == VFLAG_DISABLED) &&
                ((list->varList[i].flags & VFLAG_SET_BY_USER) == 0)) {
                continue;
            }

            type = list->varList[i].valType;
            count = list->varList[i].valCnt;

            if (count > 1) {
                lineFeed = "\n";
                lBracket = "[";
                rBracket = "]";
            } else {
                lineFeed = "";
                lBracket = "";
                rBracket = "";
            }

/*
 *          If the parameter type is not a 'comment', start by
 *          writing the identifier/name.
 */
            if (type != V_COMMENT) {
                fprintf(fp, "%s = %s%s", list->varList[i].varName,
                        lBracket, lineFeed);
            }

            for (j = 0; j < count; j++) {
                switch(type) {
                    case V_INT:
                        fprintf(fp, "  %d%s",
                               ((int *)(list->varList[i].valList))[j],
                               lineFeed);
                        break;
                    case V_DBL:
                        fprintf(fp, "  %e%s",
                               ((double *)(list->varList[i].valList))[j],
                               lineFeed);
                        break;
                    case V_STRING:
                        if (list->varList[i].valCnt == 1) {
                            fprintf(fp, "  \"%s\"%s",
                                   (char *)list->varList[i].valList, lineFeed);
                        } else {
                            fprintf(fp, "  \"%s\"%s",
                                   ((char **)(list->varList[i].valList))[j],
                                   lineFeed);
                        }
                        break;
                    case V_COMMENT:
                        fprintf(fp, "#\n#  %s\n#",
                               (char *)list->varList[i].varName);
                        break;

                }  /* end switch */

            }  /* for (j = 0; j < count; ...) */

            fprintf(fp, "  %s\n", rBracket);

        }  /* for (i = minIndex; i < maxIndex; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       GetParamVals
 *      Description:    Given a variable type and the number of expected
 *                      associated values, reads the values from the
 *                      file stream and (potentially) returns them to the
 *                      caller.  
 *
 *                      Note: This function assumes the format of the
 *                      data to be read in will be one of:
 *                            "= <value>" 
 *                        or  "= [<value_list>]"
 *
 *      Arguments:
 *          fp           File descriptor from which to read.
 *          valType      Type of values to read (i.e. V_INT, V_DBL, etc).
 *          valsExpected Number of values we expect to read.
 *          valList      Address at which to store the values read.
 *                       Note: the caller is responsible for insuring
 *                       that the storage at <valList> is of sufficient
 *                       size to accept the values read in.  If we're
 *                       expected to read an array of strings, then
 *                       <valList> should be an array of pointers to
 *                       character strings, each character string having
 *                       been pre-allocated by the caller.  If this
 *                       pointer is NULL, the function will read
 *                       the values but discard them rather than attempt
 *                       to return them to the caller.
 *
 *      Returns:  0 on error
 *                1 if the proper numbers were read from <fp>
 *
 *-------------------------------------------------------------------------*/
int GetParamVals(FILE *fp, int valType, int valsExpected, void *valList)
{
        int  valCnt, tokenType;
        int  lBracketCnt, isValList, itemCnt;
        char token[256];

        isValList     = 0;
        lBracketCnt   = 0;

/*
 *      Always check for "=" before value or value_list.
 */
        tokenType = GetNextToken(fp, token, sizeof(token));

        if (tokenType != TOKEN_EQUAL) {
            printf("Error: Expected '=' in parameter assignment\n");
            return(0);
        }

/*
 *      Loop until we've parsed the next value/value_list, or hit
 *      an error...
 */
        for (itemCnt = 0, valCnt = 0; ; valCnt++, itemCnt++) {

            tokenType = GetNextToken(fp, token, sizeof(token));

/*
 *          If the first thing we find is a '[', expect a list of
 *          values contained within []'s.
 */
            if ((tokenType == TOKEN_BEGIN_VAL_LIST) && (itemCnt == 0)) {
                lBracketCnt = 1;
                isValList = 1;
                valCnt--;
                continue;
            }

/*
 *          Take action based on the type of token found
 */
            switch (tokenType) {

                case TOKEN_ERR:
                case TOKEN_NULL:
                    return(0);

/*
 *              Nested "[".  ignore it but bump up the nesting count
 */
                case TOKEN_BEGIN_VAL_LIST:
                    lBracketCnt++;
                    valCnt--;
                    break;

/*
 *              Decrement the nesting count for [] pairs.  If this
 *              brings the count negative, there's an error, otherwise
 *              return to the caller if we found the expected number
 *              of values or we were parsing an unknown number of values.
 */
                case TOKEN_END_VAL_LIST:
                    lBracketCnt--;
                    if (lBracketCnt < 0) return(0);
                    if (lBracketCnt == 0) {
                        if ((valsExpected == 0) ||
                            (valCnt == valsExpected)) {
                            return(1);
                        }
                        if (valCnt < valsExpected) return(0);
                    }
                    valCnt--;
                    break;

/*
 *              An "=" is just skipped.
 */
                case TOKEN_EQUAL:
                    valCnt--;
                    continue;
                    break;

/*
 *              Anything else is treated as a value.
 */
                default:
/*
 *                  If the caller did not provide a location in which
 *                  to return the results, just continue with the parsing.
 */
                    if (valList == (void *)NULL) break;

/*
 *                  If we've already parsed the expected number of values
 *                  there's an error.
 */
                    if (valCnt >= valsExpected) return(0);

/*
 *                  Save the value for the caller.
 */
                    switch (valType) {
                        case V_DBL:
                            ((double *)valList)[valCnt] = atof(token);
                            break;
                        case V_INT:
                            ((int *)valList)[valCnt] = atoi(token);
                            break;
                        case V_STRING:
                            if (valsExpected == 1) {
                                strcpy((char *)valList, token);
                            } else {
                                strcpy(((char **)(valList))[valCnt], token);
                            }
                            break;
                        default:
                            return(0);
                    }
                    break;

            }  /* end switch(tokenType) */

/*
 *          If we expected and found only a single data item, we're done.
 */
            if ((!isValList) && (valsExpected <= 1)) return(1);

        }  /* for (valCnt = 0; ... ) */
}


/*---------------------------------------------------------------------------
 *
 *      Function:     LookupParam
 *      Description:  Search the specified parameter list for a
 *                    binding matching the given token.
 *
 *      Arguments:
 *          token  Name of the parameter to look up in the paramater list
 *
 *      Returns:  -1 if no parameter matching <token> was found.  If
 *                a matching parameter is found and the matching parameter
 *                is an alias for another parameter, the index of the
 *                aliased parameter is returned.  Otherwise, the index
 *                of the matching parameter is returned.
 *
 *-------------------------------------------------------------------------*/
int LookupParam(ParamList_t *list, char *token)
{
        int  offset, index1, index2;
        char *name;
        void *addr1;

/*
 *      Check if the specified name matches any known parameter.
 */
        for (index1 = 0; index1 < list->paramCnt; index1++) {

            name = list->varList[index1].varName;

            for (offset = 0; ; offset++ ) {

                if ((name[offset] == 0) && (token[offset] == 0)) {
/*
 *                  If the specified identifer(token) is actually an
 *                  alias for another identifier, find the primary
 *                  identifier and return that index instead.
 */
                    if ((list->varList[index1].flags & VFLAG_ALIAS) != 0) {
                        addr1 = list->varList[index1].valList;
                        for (index2 = 0; index2 < list->paramCnt; index2++) {
                            if (index2 == index1) continue;
                            if (addr1 == list->varList[index2].valList) {
                                return(index2);
                            }
                        }
                    }
/*
 *                  If the identifier is not an alias for another 
 *                  identifier OR we cannot determine the identifier
 *                  being aliased (which should never happen), return
 *                  the index of the current identifier.
 */
                    return(index1);
                }
/*
 *              Parameter names are not case sensitive, so ignore
 *              any difference between upper/lower case
 */
                if (tolower(name[offset]) != tolower(token[offset])) {
                    break;
                }
            }
        }

        return(-1);
}
            

/*------------------------------------------------------------------------
 *
 *      Function:     GetNextToken
 *      Description:  Parse an input stream and return to the caller the
 *                    next token on the stream.
 *
 *                    NOTE: This parser is targeted toward the control
 *                    file input where data is specified in one of the
 *                    forms:
 *
 *                       identifier = numeric_value
 *                       identifier = "string_value"
 *                       identifier = 'string_value'
 *                    or identifier = [value_list]
 *
 *                    Value lists may span multiple lines so long as
 *                    individual values are delimited by white-space
 *                    or line-feeds.
 *
 *                    Any '#' not contained within quotes causes the
 *                    parser to skip the remainder of the current line.
 *
 *      Arguments:
 *          fp            file stream to parse
 *          token         location in which to return a token to the
 *                        caller.  Tokens returned are either identifiers
 *                        or values.  If the function finds delimiting
 *                        characters such as "=", "[" or "]", the
 *                        returned token will be a zero-length string but
 *                        the return value of the function will indicate
 *                        the type of delimiter found.
 *          maxTokenSize  number of bytes (including a NULL terminating
 *                        character) that can be stored into <token>.
 *
 *      Returns:  A value indicating what type of token or delimiting
 *                character was found.  Possible return values include:
 *
 *                    TOKEN_ERR             A parsing error was encountered
 *                    TOKEN_NULL            No token was found (i.e. EOF)
 *                    TOKEN_GENERIC         An identifier or value was
 *                                          returned to the caller in <token>
 *                    TOKEN_EQUAL           Found an "=" sign.  Nothing
 *                                          returned to caller in <token>
 *                    TOKEN_BEGIN_VAL_LIST  Found delimiter starting a
 *                                          value list.  Nothing returned
 *                                          to caller in <token>
 *                    TOKEN_END_VAL_LIST    Found delimiter terminating a
 *                                          value list.  Nothing returned
 *                                          to caller in <token>
 *
 *----------------------------------------------------------------------*/
int GetNextToken(FILE *fp, char *token, int maxTokenSize)
{
        int   inSingleQuote, inDblQuote;
        int   foundToken, tokenLen, tokenType, appendToToken, nextChar;

        inSingleQuote = 0;
        inDblQuote    = 0;
        foundToken    = 0;
        tokenLen      = 0;
        tokenType     = TOKEN_GENERIC;

        memset(token, 0, maxTokenSize);

        while (1) {

/*
 *          Read the next character from the iunput stream.  
 *          Start with the assumption that the next character is
 *          part of the token to be returned to the caller.
 *
 *          NOTE: In some cases, after reading a character, it
 *          is determined the character functions as an end-of-token
 *          delimiter and we do not want to consume that character
 *          yet.  In these cases, that single character is pushed
 *          back onto the input stream to be read later.
 */
            appendToToken = 1;
            nextChar = fgetc(fp);

            switch (nextChar) {
/*
 *              White space characters imbedded within quotes are
 *              considered part of the token.  Outside of quotes, they
 *              are ignored if they preceed the next token, and following
 *              a token they are treated as token delimiters.
 */
                case ' ':
                case '\t':
                    if (inSingleQuote || inDblQuote) break;
                    if (foundToken) return(tokenType);
                    appendToToken = 0;
                    break;
/*
 *              Comment characters imbedded within quotes are considered
 *              part of the token.  Outside of quotes they are treated
 *              as a token delimiter and anything following this character
 *              up to the next linefeed is ignored.
 */
                case '#':
                    if (inSingleQuote || inDblQuote) break;
                    if (foundToken) {
                        ungetc(nextChar, fp);
                        return(tokenType);
                    }
                    while ((nextChar != '\0') && (nextChar != '\n')) {
                        nextChar = fgetc(fp);
                    }
                    if (nextChar == '\n') ungetc(nextChar, fp);
                    appendToToken = 0;
                    break;
/*
 *              Single quotes imbedded within double quotes are considered
 *              part of the token.  Otherwise, quotes are treated as
 *              token delimiters and any character contained within a 
 *              a matched pair of quotes is treated as part of a token.
 *
 *              NOTE: Quoted strings may not span multiple lines.
 */
                case '\'':
                    if (inSingleQuote) return(tokenType);
                    if (inDblQuote) break;
                    if (foundToken) {
                        ungetc(nextChar, fp);
                        return(tokenType);
                    }
                    inSingleQuote = 1;
                    foundToken = 1;
                    appendToToken = 0;
                    break;
/*
 *              Double quotes imbedded within single quotes are considered
 *              part of the token.  Otherwise, quotes are treated as
 *              token delimiters and any character contained within a 
 *              a matched pair of quotes is treated as part of a token.
 *
 *              NOTE: Quoted strings may not span multiple lines.
 */
                case '\"':
                    if (inDblQuote) return(tokenType);
                    if (inSingleQuote) break;
                    if (foundToken) {
                        ungetc(nextChar, fp);
                        return(tokenType);
                    }
                    inDblQuote = 1;
                    foundToken = 1;
                    appendToToken = 0;
                    break;
/*
 *              '=' characters imbedded within quotes are considered
 *              part of the token.  Outside of quotes, they are
 *              treated as token delimiters.  The function will
 *              return a zero length string in <token> and a return
 *              code of TOKEN_EQUAL if this character is the first
 *              non-white-space character found.
 */
                case '=':
                    if (inSingleQuote || inDblQuote) break;
                    if (foundToken) {
                        ungetc(nextChar, fp);
                        return(tokenType);
                    }
                    return(TOKEN_EQUAL);
/*
 *              '[' characters imbedded within quotes are considered
 *              part of the token.  Outside of quotes, they are
 *              treated as token delimiters.  The function will
 *              return a zero length string in <token> and a return
 *              code of TOKEN_BEGIN_VAL_LIST if this character is
 *              the first non-white-space character found.
 */
                case '[':
                    if (inSingleQuote || inDblQuote) break;
                    if (foundToken) {
                        ungetc(nextChar, fp);
                        return(tokenType);
                    }
                    return(TOKEN_BEGIN_VAL_LIST);
/*
 *              ']' characters imbedded within quotes are considered
 *              part of the token.  Outside of quotes, they are
 *              treated as token delimiters.  The function will
 *              return a zero length string in <token> and a return
 *              code of TOKEN_END_VAL_LIST if this character is
 *              the first non-white-space character found.
 */
                case ']':
                    if (inSingleQuote || inDblQuote) break;
                    if (foundToken) {
                        ungetc(nextChar, fp);
                        return(tokenType);
                    }
                    return(TOKEN_END_VAL_LIST);
/*
 *              Line feeds found before the end of a quoted string
 *              are considered errors.  Outside of a quoted string,
 *              they will be treated as token delimiters.
 */
                case '\n':
                    if (inSingleQuote || inDblQuote) return(TOKEN_ERR);
                    if (foundToken) return(tokenType);
                    appendToToken = 0;
                    break;
/*
 *              NULL characters or an EOF found before the end of 
 *              a quoted string are considered errors.  Outside
 *              of quoted strings they are treated as token
 *              delimiters.  If either of these is found before
 *              a token, end-of-file is assumed and a return code
 *              of TOKEN_NULL is returned to the caller along with
 *              a zero length string in <token>.
 */
                case 0:
                case EOF:
                    if (inSingleQuote || inDblQuote) return(TOKEN_ERR);
                    if (foundToken) return(tokenType);
                    return(TOKEN_NULL);
/*
 *              Anything else is treated as part of the token to be
 *              returned to the caller.
 */
                default:
                    break;
            }

/*
 *          If this character is to be treated as part of the token
 *          to be returned to the caller, make sure there is enough
 *          space available and append the character to the token.
 */
            if (appendToToken) {
                if (++tokenLen > (maxTokenSize-1)) return(TOKEN_ERR);
                foundToken = 1;
                tokenType = TOKEN_GENERIC;
                *token = nextChar;
                token++;
            }
        }

/*
 *      Just for show; we should never hit this.
 */
        return(TOKEN_ERR);
}
