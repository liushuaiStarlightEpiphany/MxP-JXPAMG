#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef int JX_Int;

#define jx_TAlloc(type, count) ((type*)malloc((count) * sizeof(type)))
#define jx_TFree(ptr) free(ptr)

JX_Int
jx_new_format( const char *format, char **newformat_ptr )
{
   const char *fp;
   char       *newformat, *nfp;
   JX_Int   newformatlen;
   JX_Int   foundpercent = 0;

   newformatlen = 2*strlen(format)+1;
   newformat = (char*)malloc(newformatlen);

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
            fp ++;
            if (*fp == 'l')
            {
               fp ++;
            }
         }
         switch (*fp)
         {
            case 'd':
            case 'i':
               foundpercent = 0; break;
            case 'f':
            case 'e':
            case 'E':
            case 'g':
            case 'G':
              *nfp = 'l'; nfp++;
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

int main() {
    const char *fmt = "[AHI-STAT] Level %d: N = %d, N-F = %d, N-EI = %d, r_EI_node = %.4f %%, r_EI_F = %.4f %%";
    char *newfmt = NULL;
    jx_new_format(fmt, &newfmt);
    printf("Original: [%s]\n", fmt);
    printf("New:      [%s]\n", newfmt);
    printf("Match:    %s\n", strcmp(fmt, newfmt) == 0 ? "YES" : "NO");
    if (strcmp(fmt, newfmt) != 0) {
        printf("Differences:\n");
        for (int i = 0; fmt[i] || newfmt[i]; i++) {
            if (fmt[i] != newfmt[i]) {
                printf("  pos %d: '%c'(%02x) vs '%c'(%02x)\n",
                       i, fmt[i] ? fmt[i] : '?', (unsigned char)fmt[i],
                       newfmt[i] ? newfmt[i] : '?', (unsigned char)newfmt[i]);
            }
        }
    }
    free(newfmt);
    return 0;
}
