#ifndef WRITER_H
#define WRITER_H


#include <iostream>
#include <string.h>
#include <cstdio>
#include <fstream>

using namespace std;

extern ofstream output;

int statementZone = 0;
int zonaPragma = 0;
int finalizeOK = 0;
int task = 0;
int deletePrePragma = 0;

char * _yytext=NULL;
char *linea = NULL;

extern void finSecuencial();
extern void MPITaskEnd();

void setYytext(char * yytext) {
    if (_yytext != NULL) {
        free(_yytext);  // Free the previously allocated memory
    }
    _yytext = yytext ? strdup(yytext) : NULL;  // Duplicate yytext_pre1
}

void updateText() {
    /*if(zonaPragma == 1 && deletePrePragma == 0) {
        delete[] linea;
        linea = NULL;
		printf("LIMPIEZA\n");
        deletePrePragma = 1;
    }*/
    if(strcmp(_yytext, "\n") != 0) {

        int size = (linea ? strlen(linea) : 0) + strlen(_yytext) + 1;
        char *addToLinea = new char[size];

        if (linea) {
            strcpy(addToLinea, linea);
            delete[] linea;
        } else {
            addToLinea[0] = '\0';
        }

        strcat(addToLinea, _yytext);
        linea = addToLinea;
    } else {    
        if(zonaPragma == 1 && finalizeOK == 0 && task == 0 ) {
            output << "\tif (__taskid == 0) {\n" << endl;
            zonaPragma = 0;
        }
        if(task != 0) {
            
            output << "\t" << (linea ? linea : "") << endl;
            if(strchr(linea, '}') != nullptr){
                MPITaskEnd();
                task = 0;
            }
            
        }
        else if(statementZone == 1 && zonaPragma == 0 && finalizeOK == 0) {
            output << (linea ? "\t" + string(linea) : "") << endl;
        }
        else{
            output << (linea ? linea : "") << endl;
        }
        delete[] linea;
        linea = NULL;
		deletePrePragma = 0;
    }
    
}

void lastLine(){
    output << (linea ? linea : "") << endl;
}

#endif
