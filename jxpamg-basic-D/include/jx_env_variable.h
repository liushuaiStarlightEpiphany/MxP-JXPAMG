#ifndef JX_ENV_VARIABLE_H
#define JX_ENV_VARIABLE_H


typedef int JX_Int;
/*!
 * \enum JX_ENV_TYPE
 */ 
typedef enum JX_ENV_TYPE{
    JX_ENV_HPCTOOLKIT_PROFILE_AMGLEVEL,
    Bind_Type,
}JX_ENV_TYPE;
/*!
 * Gets the value of an environment variables.
 */
JX_Int jx_getenv(JX_ENV_TYPE env_type);

#endif

