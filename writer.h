#ifndef WRITER_H
#define WRITER_H


#include <iostream>
#include <string.h>
#include <cstdio>
#include <fstream>

using namespace std;

extern ofstream output;
int statementZone = 0;

char * _yytext=NULL;
char *linea = NULL;

void setYytext(char * yytext) {
    if (_yytext != NULL) {
        free(_yytext);  // Free the previously allocated memory
    }
    _yytext = yytext ? strdup(yytext) : NULL;  // Duplicate yytext_pre1
}

void updateText() {
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
        //if(statementZone == 1){
            //output << (linea ? "\t" + string(linea) : "") << endl;
        //}
        //else{
            output << (linea ? linea : "") << endl;
        //}
        delete[] linea;
        linea = NULL;
    }
}

#endif
