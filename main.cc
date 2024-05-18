#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "1805082_SymbolTable.h"

extern int yydebug;
extern FILE *yyin, *yyout;
FILE *inputFile;

extern int MVL_LINNUM;

extern void yyerror(const char *s);
extern int yyparse ();

std::ofstream logFile, errFile, sym_tables;

extern int error_count, line_count;

extern SymbolTable table;

int main( int argc, const char* argv[] )
{
  /* missing parameter check */
    
  /* yydebug=1; */


  /* open files */
  if ((argc!=4) && (argc!=5)) {
      std::cout << "command: ./fparse input.c log.txt error.txt [output.c]" << std::endl;
      return 0;
  }

  if((inputFile=fopen(argv[1],"r"))==NULL) {
      printf("Cannot Open Input File.\n");
      exit(1);
  }

  logFile.open(argv[2]);
  errFile.open(argv[3]);

  if (argc==5){
    if ((yyout=fopen(argv[4],"w"))==NULL){
      yyerror("could not open outputfile\n");
    }
  }

  sym_tables.open("sym_tables.txt");

  yyin=inputFile;

  //yydebug = 1;  // Enable Bison's debug mode

  yyparse();

  table.exitScope();

  logFile << endl;
  logFile << "Total lines: " << line_count << endl;
  logFile << "Total errors: " << error_count << endl << endl;

  table.printCurrScopeTable(); // Print the current scope table

  logFile.close();
  sym_tables.close();
  errFile.close();
  fclose(yyin);

  return 0;

}
