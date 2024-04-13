#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int yydebug;
extern FILE *yyin, *yyout;

extern int MVL_LINNUM;

extern void yyerror(const char *s);
extern int yyparse ();


int main( int argc, const char* argv[] )
{
  char* infilename;

  /* missing parameter check */

  /* yydebug=1; */


  /* open files */
  /* input file */
  if (argc>1){
    infilename=(char*)malloc(sizeof(char)*strlen(argv[1])+1);
    strcpy(infilename, argv[1]);
    if ((yyin=fopen(argv[1],"r"))==NULL){
      yyerror("unable to open inputfile");
    }
    argc--;
    argv++;
  }else{
    infilename=" ";
  }
  /* output file */
  if (argc>1){
    if ((yyout=fopen(argv[1],"w"))==NULL){
      yyerror("could not open outputfile\n");
    }
  }

  yyparse();


  exit(0);

}
