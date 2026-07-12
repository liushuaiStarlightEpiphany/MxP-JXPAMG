//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jx_util.h"
#include <stdarg.h>
#include <stdio.h>

/* these prototypes are missing by default for some compilers */
int vscanf( const char *format , va_list arg );
int vfscanf( FILE *stream , const char *format, va_list arg );
int vsscanf( const char *s , const char *format, va_list arg );

JX_Int
jx_new_format( const char *format, char **newformat_ptr )
{
   const char *fp;
   char       *newformat, *nfp;
   JX_Int   newformatlen;

   newformatlen = 2 * strlen(format) + 1;
   newformat = jx_TAlloc(char, newformatlen);

   nfp = newformat;
   for (fp = format; *fp != '\0'; )
   {
      if (*fp == '%')
      {
         *nfp = '%'; nfp ++; fp ++;

         if (*fp == '%')
         {
            *nfp = '%'; nfp ++; fp ++;
            continue;
         }

         while (*fp == 'l') fp ++;

         while (*fp != '\0' &&
                *fp != 'd' && *fp != 'i' &&
                *fp != 'f' && *fp != 'e' && *fp != 'E' &&
                *fp != 'g' && *fp != 'G' &&
                *fp != 'c' && *fp != 's' && *fp != 'p' &&
                *fp != 'n' && *fp != 'o' && *fp != 'u' &&
                *fp != 'x' && *fp != 'X' &&
                *fp != '%')
         {
            *nfp = *fp; nfp ++; fp ++;
         }

         if (*fp == '\0') break;

         switch (*fp)
         {
            case 'd':
            case 'i':
#if JX_USING_BIG_INT
              *nfp = 'l'; nfp ++;
              *nfp = 'l'; nfp ++;
#endif
              *nfp = *fp; nfp ++; fp ++;
               break;
            case 'f':
            case 'e':
            case 'E':
            case 'g':
            case 'G':
#if JX_USING_BIG_DOUBLE
              *nfp = 'L'; nfp ++;
#endif
              *nfp = *fp; nfp ++; fp ++;
               break;
            default:
              *nfp = *fp; nfp ++; fp ++;
               break;
         }
      }
      else
      {
         *nfp = *fp; nfp ++; fp ++;
      }
   }
   *nfp = '\0';

  *newformat_ptr = newformat;

   return 0;
}

JX_Int
jx_free_format( char *newformat )
{
   jx_TFree(newformat);
   return 0;
}

/* printf functions */

JX_Int
jx_printf( const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JX_Int ierr = 0;

   va_start(ap, format);
   jx_new_format(format, &newformat);
   ierr = vprintf(newformat, ap);
   jx_free_format(newformat);
   va_end(ap);

   return ierr;
}

JX_Int
jx_fprintf( FILE *stream, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JX_Int ierr = 0;

   va_start(ap, format);
   jx_new_format(format, &newformat);
   ierr = vfprintf(stream, newformat, ap);
   jx_free_format(newformat);
   va_end(ap);

   return ierr;
}

JX_Int
jx_sprintf( char *s, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JX_Int ierr = 0;

   va_start(ap, format);
   jx_new_format(format, &newformat);
   ierr = vsprintf(s, newformat, ap);
   jx_free_format(newformat);
   va_end(ap);

   return ierr;
}

/* scanf functions */

JX_Int
jx_scanf( const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JX_Int ierr = 0;

   va_start(ap, format);
   jx_new_format(format, &newformat);
   ierr = vscanf(newformat, ap);
   jx_free_format(newformat);
   va_end(ap);

   return ierr;
}

JX_Int
jx_fscanf( FILE *stream, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JX_Int ierr = 0;

   va_start(ap, format);
   jx_new_format(format, &newformat);
   ierr = vfscanf(stream, newformat, ap);
   jx_free_format(newformat);
   va_end(ap);

   return ierr;
}

JX_Int
jx_sscanf( char *s, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JX_Int ierr = 0;

   va_start(ap, format);
   jx_new_format(format, &newformat);
   ierr = vsscanf(s, newformat, ap);
   jx_free_format(newformat);
   va_end(ap);

   return ierr;
}
