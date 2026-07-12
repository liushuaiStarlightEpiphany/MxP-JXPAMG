#ifndef JXF_MV_HEADER
#include "jxf_mv.h"
#endif

#ifndef JXF_ENV_VARIABLE_H
#include "jxf_env_variable.h"
#endif

/*!
 * Gets the value of an environment variables.
 */
JXF_Int jxf_getenv(JXF_ENV_TYPE env_type)
{
    char *temp;
    switch(env_type)
    {
        case JXF_Bind_Type:
        temp = getenv("JXF_Bind_Type");
        if(temp)
        {
            if(strcmp(temp,"CORE") == 0)
            {
                return 2;
            }
            else if(strcmp(temp,"CPU") == 0)
            {
                return 1;
            }
            else if(strcmp(temp,"NODE") == 0)
            {
                return 0;
            }
            else
            {
                printf("\n >>> no such bind type \n");
                exit(0);
            }
        }
        else
        {
            return 2;
        }
        break;

        default:return 0;
    }
}