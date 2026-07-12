//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  jxf_parser_dh.c
 *  Date: 2013/01/22
 */

#include "jxf_euclid.h"

typedef struct _jxf_optionsNode jxf_OptionsNode;

struct _jxf_parser_dh
{
    jxf_OptionsNode *head;
    jxf_OptionsNode *tail;
};

struct _jxf_optionsNode
{
    char *name;
    char *value;
    jxf_OptionsNode *next;
};

static jxf_bool jxf_find( jxf_Parser_dh p, char *option, jxf_OptionsNode** ptr );
static void jxf_init_from_default_settings_private( jxf_Parser_dh p );

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhCreate"
void jxf_Parser_dhCreate( jxf_Parser_dh *p )
{
    JXF_START_FUNC_DH
    jxf_OptionsNode *ptr;
    struct _jxf_parser_dh *tmp = (struct _jxf_parser_dh *)JXF_MALLOC_DH(sizeof(struct _jxf_parser_dh)); JXF_CHECK_V_ERROR;
   *p = tmp;
    tmp->head = tmp->tail = (jxf_OptionsNode*)JXF_MALLOC_DH(sizeof(jxf_OptionsNode)); JXF_CHECK_V_ERROR;
    ptr = tmp->head;
    ptr->next = NULL;
    ptr->name  = (char *)JXF_MALLOC_DH(6*sizeof(char)); JXF_CHECK_V_ERROR;
    ptr->value = (char *)JXF_MALLOC_DH(6*sizeof(char)); JXF_CHECK_V_ERROR;
    strcpy(ptr->name, "JUNK");
    strcpy(ptr->value, "JUNK");
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhDestroy"
void jxf_Parser_dhDestroy( jxf_Parser_dh p )
{
    JXF_START_FUNC_DH
    jxf_OptionsNode *ptr2 = p->head, *ptr1 = ptr2;
    if (ptr1 != NULL)
    {
        do
        {
            ptr2 = ptr2->next;
            JXF_FREE_DH(ptr1->name);
            JXF_FREE_DH(ptr1->value);
            JXF_FREE_DH(ptr1);
            ptr1 = ptr2;
        } while (ptr1 != NULL);
    }
    JXF_FREE_DH(p);
    JXF_END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhUpdateFromFile"
void jxf_Parser_dhUpdateFromFile( jxf_Parser_dh p, char *filename )
{
    JXF_START_FUNC_DH_2
    char line[80], name[80], value[80];
    FILE *fp;
    
    if ((fp = fopen(filename, "r")) == NULL)
    {
        jxf_sprintf(jxf_msgBuf_dh, "can't open >>%s<< for reading", filename);
        JXF_SET_INFO(jxf_msgBuf_dh);
    }
    else
    {
        jxf_sprintf(jxf_msgBuf_dh, "updating parser from file: >>%s<<", filename);
        JXF_SET_INFO(jxf_msgBuf_dh);
        while (!feof(fp))
        {
            if (fgets(line, 80, fp) == NULL) break;
            if (line[0] != '#')
            {
                if (jxf_sscanf(line, "%s %s", name, value) != 2) break;
                jxf_Parser_dhInsert(p, name, value);
            }
        }
        fclose(fp);
    }
    JXF_END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhInit"
void jxf_Parser_dhInit( jxf_Parser_dh p, JXF_Int argc, char *argv[] )
{
    JXF_START_FUNC_DH_2
    JXF_Int j;
    
    jxf_init_from_default_settings_private(p); JXF_CHECK_V_ERROR;
    jxf_Parser_dhUpdateFromFile(p, "./database"); JXF_CHECK_V_ERROR;
    for (j = 1; j < argc; ++ j)
    {
        if (strcmp(argv[j],"-db_filename") == 0)
        {
            ++ j;
            if (j < argc)
            {
                jxf_Parser_dhUpdateFromFile(p, argv[j]); JXF_CHECK_V_ERROR;
            }
        }
    }
    JXF_Int i = 0;
    while (i < argc)
    {
        if (argv[i][0] == '-')
        {
            char value[] = {"1"};
            jxf_bool flag = jxf_false;
            
            if (i+1 < argc && argv[i+1][0] == '-' && argv[i+1][1] == '-')
            {
                flag = jxf_true;
            }
            if ((i+1 == argc || argv[i+1][0] == '-') && !flag)
            {
                jxf_Parser_dhInsert(p, argv[i], value);
            }
            else if (flag)
            {
                jxf_Parser_dhInsert(p, argv[i], argv[i+1]+1);
            }
            else
            {
                jxf_Parser_dhInsert(p, argv[i], argv[i+1]);
            }
        }
        ++ i;
    }
    JXF_END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhHasSwitch"
jxf_bool jxf_Parser_dhHasSwitch( jxf_Parser_dh p, char *s )
{
    JXF_START_FUNC_DH_2
    jxf_bool has_switch = jxf_false;
    jxf_OptionsNode *node;
    
    if (p != NULL && jxf_find(p,s,&node))
    {
        if (!strcmp(node->value, "0"))
        {
            has_switch = jxf_false;
        }
        else if (!strcmp(node->value, "jxf_false"))
        {
            has_switch = jxf_false;
        }
        else if (!strcmp(node->value, "False"))
        {
            has_switch = jxf_false;
        }
        else if (!strcmp(node->value, "FALSE"))
        {
            has_switch = jxf_false;
        }
        else
        {
            has_switch = jxf_true;
        }
    }
    JXF_END_FUNC_VAL_2(has_switch)
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhReadInt"
jxf_bool jxf_Parser_dhReadInt( jxf_Parser_dh p, char *in, JXF_Int *out )
{
    JXF_START_FUNC_DH_2
    jxf_bool has_switch = jxf_false;
    jxf_OptionsNode *node;
    
    if (p != NULL && jxf_find(p,in,&node))
    {
       *out = atoi(node->value);
        if (!strcmp(node->value, "0"))
        {
            has_switch = jxf_false;
        }
        else
        {
            has_switch = jxf_true;
        }
    }
    JXF_END_FUNC_VAL_2(has_switch)
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhReadDouble"
jxf_bool jxf_Parser_dhReadDouble( jxf_Parser_dh p, char *in, JXF_Real *out )
{
    JXF_START_FUNC_DH_2
    jxf_bool optionExists = jxf_false;
    jxf_OptionsNode *node;
    
    if (p != NULL && jxf_find(p,in,&node))
    {
       *out = atof(node->value);
        optionExists = jxf_true;
    }
    JXF_END_FUNC_VAL_2(optionExists)
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhReadString"
jxf_bool jxf_Parser_dhReadString( jxf_Parser_dh p, char *in, char **out )
{
    JXF_START_FUNC_DH_2
    jxf_bool optionExists = jxf_false;
    jxf_OptionsNode *node;
    
    if (p != NULL && jxf_find(p,in,&node))
    {
       *out = node->value;
        optionExists = jxf_true;
    }
    JXF_END_FUNC_VAL_2(optionExists)
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhPrint"
void jxf_Parser_dhPrint( jxf_Parser_dh p, FILE *fp, jxf_bool allPrint )
{
    JXF_START_FUNC_DH_2
    jxf_OptionsNode *ptr = p->head;
    
    if (fp == NULL) JXF_SET_V_ERROR("fp == NULL");
    if (jxf_myid_dh == 0 || allPrint)
    {
        jxf_fprintf(fp, "------------------------ registered options:\n");
        if (ptr == NULL)
        {
            jxf_fprintf(fp, "Parser object is invalid; nothing to print!\n");
        }
        else
        {
            ptr = ptr->next;
            while (ptr != NULL)
            {
                jxf_fprintf(fp, "   %s  %s\n", ptr->name, ptr->value);
                fflush(fp);
                ptr = ptr->next;
            }
        }
        jxf_fprintf(fp, "\n");
        fflush(fp);
    }
    JXF_END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "jxf_Parser_dhInsert"
void jxf_Parser_dhInsert( jxf_Parser_dh p, char *option, char *value )
{
    JXF_START_FUNC_DH_2
    jxf_OptionsNode *node;
    JXF_Int length;
    
    if (p == NULL) goto PARSER_NOT_INITED;
    if (jxf_find(p,option,&node))
    {
        JXF_Int length2 = strlen(node->value) + 1;
        
        length = strlen(value) + 1;
        if (length2 < length)
        {
            JXF_FREE_DH(node->value);
            node->value  = (char *)JXF_MALLOC_DH(length*sizeof(char)); JXF_CHECK_V_ERROR;
        }
        strcpy(node->value, value);
    }
    else
    {
        node = p->tail;
        p->tail = node->next = (jxf_OptionsNode *)JXF_MALLOC_DH(sizeof(jxf_OptionsNode)); JXF_CHECK_V_ERROR;
        node = node->next;
        length = strlen(option) + 1;
        node->name = (char *)JXF_MALLOC_DH(length*sizeof(char)); JXF_CHECK_V_ERROR;
        strcpy(node->name, option);
        length = strlen(value) + 1;
        node->value = (char *)JXF_MALLOC_DH(length*sizeof(char)); JXF_CHECK_V_ERROR;
        strcpy(node->value, value);
        node->next = NULL;
    }
    
PARSER_NOT_INITED: ;
    
    JXF_END_FUNC_DH_2
}

#undef __FUNC__
#define __FUNC__ "jxf_find"
jxf_bool jxf_find( jxf_Parser_dh p, char *option, jxf_OptionsNode **ptr )
{
    JXF_START_FUNC_DH_2
    jxf_OptionsNode *tmpPtr = p->head;
    jxf_bool foundit = jxf_false;
    
    while (tmpPtr != NULL)
    {
        if (strcmp(tmpPtr->name,option) == 0)
        {
            foundit = jxf_true;
           *ptr = tmpPtr;
            break;
        }
        tmpPtr = tmpPtr->next;
    }
    JXF_END_FUNC_VAL_2(foundit)
}

#undef __FUNC__
#define __FUNC__ "jxf_init_from_default_settings_private"
void jxf_init_from_default_settings_private( jxf_Parser_dh p )
{
    JXF_START_FUNC_DH_2
    jxf_Parser_dhInsert(p, "-sig_dh", "1"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-px", "1"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-py", "1"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-pz", "0"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-m", "4"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-xx_coeff", "-1.0"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-yy_coeff", "-1.0"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-zz_coeff", "-1.0"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-level", "1"); JXF_CHECK_V_ERROR;
    jxf_Parser_dhInsert(p, "-printStats", "0"); JXF_CHECK_V_ERROR;
    JXF_END_FUNC_DH_2
}
