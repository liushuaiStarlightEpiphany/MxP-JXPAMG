#ifndef JXF_ENV_VARIABLE_H
#define JXF_ENV_VARIABLE_H


typedef int JXF_Int;
/*!
 * \enum JXF_ENV_TYPE
 */ 
typedef enum JXF_ENV_TYPE{
    JXF_ENV_HPCTOOLKIT_PROFILE_AMGLEVEL,
    JXF_Bind_Type,
}JXF_ENV_TYPE;
/*!
 * Gets the value of an environment variables.
 */
JXF_Int jxf_getenv(JXF_ENV_TYPE env_type);

#endif

