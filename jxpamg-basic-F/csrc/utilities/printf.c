//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

#include "jxf_util.h"
#include <stdarg.h>
#include <stdio.h>

/* these prototypes are missing by default for some compilers */
int vscanf( const char *format , va_list arg );
int vfscanf( FILE *stream , const char *format, va_list arg );
int vsscanf( const char *s , const char *format, va_list arg );

JXF_Int
jxf_new_format( const char *format, char **newformat_ptr )
{
   const char *fp;
   char       *newformat, *nfp;
   JXF_Int   newformatlen;
   JXF_Int   foundpercent = 0;

   newformatlen = 2*strlen(format)+1; /* worst case is all %d's to %lld's */
   newformat = jxf_TAlloc(char, newformatlen);

   nfp = newformat;
   for (fp = format; *fp != '\0'; fp ++)
   {
      if (*fp == '%')
      {
         foundpercent = 1;
      }
      else if (foundpercent)
      {
         if (*fp == 'l')
         {
            fp ++; /* remove 'l' and maybe add it back in switch statement */
            if (*fp == 'l')
            {
               fp ++; /* remove second 'l' if present */
            }
         }
         switch (*fp)
         {
            case 'd':
            case 'i':
#if JXF_USING_BIG_INT
              *nfp = 'l'; nfp ++;
              *nfp = 'l'; nfp ++;
#endif
               foundpercent = 0; break;
            case 'f':
            case 'e':
            case 'E':
            case 'g':
            case 'G':
#if JXF_USING_BIG_DOUBLE
              *nfp = 'L'; nfp++; /* modify with 'L' */
#else
              *nfp = 'l'; nfp++; /* modify with 'l' (default is double) */
#endif
               foundpercent = 0; break;
            case 'c':
            case 'n':
            case 'o':
            case 'p':
            case 's':
            case 'u':
            case 'x':
            case 'X':
            case '%':
               foundpercent = 0; break;
         }
      }
     *nfp = *fp; nfp ++;
   }
  *nfp = *fp; nfp ++;

  *newformat_ptr = newformat;

   return 0;
}

JXF_Int
jxf_free_format( char *newformat )
{
   jxf_TFree(newformat);
   return 0;
}

/* printf functions */

JXF_Int
jxf_printf( const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JXF_Int ierr = 0;

   va_start(ap, format);
   jxf_new_format(format, &newformat);
   ierr = vprintf(newformat, ap);
   jxf_free_format(newformat);
   va_end(ap);

   return ierr;
}

JXF_Int
jxf_fprintf( FILE *stream, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JXF_Int ierr = 0;

   va_start(ap, format);
   jxf_new_format(format, &newformat);
   ierr = vfprintf(stream, newformat, ap);
   jxf_free_format(newformat);
   va_end(ap);

   return ierr;
}

JXF_Int
jxf_sprintf( char *s, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JXF_Int ierr = 0;

   va_start(ap, format);
   jxf_new_format(format, &newformat);
   ierr = vsprintf(s, newformat, ap);
   jxf_free_format(newformat);
   va_end(ap);

   return ierr;
}

/* scanf functions */

JXF_Int
jxf_scanf( const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JXF_Int ierr = 0;

   va_start(ap, format);
   jxf_new_format(format, &newformat);
   ierr = vscanf(newformat, ap);
   jxf_free_format(newformat);
   va_end(ap);

   return ierr;
}

JXF_Int
jxf_fscanf( FILE *stream, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JXF_Int ierr = 0;

   va_start(ap, format);
   jxf_new_format(format, &newformat);
   ierr = vfscanf(stream, newformat, ap);
   jxf_free_format(newformat);
   va_end(ap);

   return ierr;
}

JXF_Int
jxf_sscanf( char *s, const char *format, ... )
{
   va_list   ap;
   char     *newformat;
   JXF_Int ierr = 0;

   va_start(ap, format);
   jxf_new_format(format, &newformat);
   ierr = vsscanf(s, newformat, ap);
   jxf_free_format(newformat);
   va_end(ap);

   return ierr;
}
