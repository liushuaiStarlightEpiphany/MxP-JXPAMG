#ifndef JX_MV_HEADER
#include "jx_mv.h"
#endif

#ifndef JX_ENV_VARIABLE_H
#include "jx_env_variable.h"
#endif

/*!
 * Gets the value of an environment variables.
 */
JX_Int jx_getenv(JX_ENV_TYPE env_type)
{
    char *temp;
    switch(env_type)
    {
        case Bind_Type:
        temp = getenv("Bind_Type");
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