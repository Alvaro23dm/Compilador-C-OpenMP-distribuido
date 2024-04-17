%{
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include "y.tab.hh"
#include <FlexLexer.h>
#include "SymbolTable.h"


extern void yyerror(const char *);

static void count(void);
void log(const std::string& token, const std::string& value);
static void comment(void);
//static int check_type(void);
int line =1;
std::ofstream fich("./logs/tokens.txt", std::ios_base::out);

%}

D			[0-9]
L			[a-zA-Z_]
H			[a-fA-F0-9]
E			([Ee][+-]?{D}+)
P                       ([Pp][+-]?{D}+)
FS			(f|F|l|L)
IS                      ((u|U)|(u|U)?(l|L|ll|LL)|(l|L|ll|LL)(u|U))

%option nounput

%%
"/*"			{ comment(); }
"//"[^\n]*              { /* consume //-comment */ }


"auto"          { count(); log("AUTO", "-"); return(AUTO); }
"_Bool"         { count(); log("BOOL", "-"); return(BOOL); }
"break"         { count(); log("BREAK", "-"); return(BREAK); }
"case"          { count(); log("CASE", "-"); return(CASE); }
"char"          { count(); log("CHAR", "-"); return(CHAR); }
"_Complex"      { count(); log("COMPLEX", "-"); return(COMPLEX); }
"const"         { count(); log("CONST", "-"); return(CONST); }
"continue"      { count(); log("CONTINUE", "-"); return(CONTINUE); }
"default"       { count(); log("DEFAULT", "-"); return(DEFAULT); }
"do"            { count(); log("DO", "-"); return(DO); }
"double"        { count(); log("DOUBLE", "-"); return(DOUBLE); }
"else"          { count(); log("ELSE", "-"); return(ELSE); }
"enum"          { count(); log("ENUM", "-"); return(ENUM); }
"extern"        { count(); log("EXTERN", "-"); return(EXTERN); }
"float"         { count(); log("FLOAT", "-"); return(FLOAT); }
"for"           { count(); log("FOR", "-"); return(FOR); }
"goto"          { count(); log("GOTO", "-"); return(GOTO); }
"if"            { count(); log("IF", "-"); return(IF); }
"_Imaginary"    { count(); log("IMAGINARY", "-"); return(IMAGINARY); }
"inline"        { count(); log("INLINE", "-"); return(INLINE); }
"int"           { count(); log("INT", "-"); return(INT); }
"long"          { count(); log("LONG", "-"); return(LONG); }
"register"      { count(); log("REGISTER", "-"); return(REGISTER); }
"restrict"      { count(); log("RESTRICT", "-"); return(RESTRICT); }
"return"        { count(); log("RETURN", "-"); return(RETURN); }
"short"         { count(); log("SHORT", "-"); return(SHORT); }
"signed"        { count(); log("SIGNED", "-"); return(SIGNED); }
"sizeof"        { count(); log("SIZEOF", "-"); return(SIZEOF); }
"static"        { count(); log("STATIC", "-"); return(STATIC); }
"struct"        { count(); log("STRUCT", "-"); return(STRUCT); }
"switch"        { count(); log("SWITCH", "-"); return(SWITCH); }
"typedef"       { count(); log("TYPEDEF", "-"); return(TYPEDEF); }
"union"         { count(); log("UNION", "-"); return(UNION); }
"unsigned"      { count(); log("UNSIGNED", "-"); return(UNSIGNED); }
"void"          { count(); log("VOID", "-"); return(VOID); }
"volatile"      { count(); log("VOLATILE", "-"); return(VOLATILE); }
"while"         { count(); log("WHILE", "-"); return(WHILE); }

{L}({L}|{D})*		{ 
                        count(); 
                        log("IDENTIFIER", yytext); 
                        yylval.symPtr = new SymbolInfo(yytext, "IDENTIFIER");
                        yylval.symPtr -> setToPrint(yytext);
                        return IDENTIFIER;
                    }

0[xX]{H}+{IS}?		{ count(); log("CONSTANT", yytext); return(CONSTANT); }
0[0-7]*{IS}?		{ count(); log("CONSTANT", yytext); return(CONSTANT); }
[1-9]{D}*{IS}?		{ count(); log("CONSTANT", yytext); return(CONSTANT); }
L?'(\\.|[^\\'\n])+'	{ count(); log("CONSTANT", yytext); return(CONSTANT); }

{D}+{E}{FS}?		        { count(); log("CONSTANT", yytext); return(CONSTANT); }
{D}*"."{D}+{E}?{FS}?	    { count(); log("CONSTANT", yytext); return(CONSTANT); }
{D}+"."{D}*{E}?{FS}?	    { count(); log("CONSTANT", yytext); return(CONSTANT); }
0[xX]{H}+{P}{FS}?	        { count(); log("CONSTANT", yytext); return(CONSTANT); }
0[xX]{H}*"."{H}+{P}{FS}?    { count(); log("CONSTANT", yytext); return(CONSTANT); }
0[xX]{H}+"."{H}*{P}{FS}?    { count(); log("CONSTANT", yytext); return(CONSTANT); }


L?\"(\\.|[^\\"\n])*\"	{ count(); log("STRING_LITERAL", yytext); return(STRING_LITERAL); }
<<EOF>>         fich.close(); return 0;

"..."			{ count(); log("ELLIPSIS", "-"); return(ELLIPSIS); }
">>="			{ count(); log("RIGHT_ASSIGN", "-"); return(RIGHT_ASSIGN); }
"<<="			{ count(); log("LEFT_ASSIGN", "-"); return(LEFT_ASSIGN); }
"+="			{ count(); log("ADD_ASSIGN", "-"); return(ADD_ASSIGN); }
"-="			{ count(); log("SUB_ASSIGN", "-"); return(SUB_ASSIGN); }
"*="			{ count(); log("MUL_ASSIGN", "-"); return(MUL_ASSIGN); }
"/="			{ count(); log("DIV_ASSIGN", "-"); return(DIV_ASSIGN); }
"%="			{ count(); log("MOD_ASSIGN", "-"); return(MOD_ASSIGN); }
"&="			{ count(); log("AND_ASSIGN", "-"); return(AND_ASSIGN); }
"^="			{ count(); log("XOR_ASSIGN", "-"); return(XOR_ASSIGN); }
"|="			{ count(); log("OR_ASSIGN", "-"); return(OR_ASSIGN); }
">>"			{ count(); log("RIGHT_OP", "-"); return(RIGHT_OP); }
"<<"			{ count(); log("LEFT_OP", "-"); return(LEFT_OP); }
"++"			{ count(); log("INC_OP", "-"); return(INC_OP); }
"--"			{ count(); log("DEC_OP", "-"); return(DEC_OP); }
"->"			{ count(); log("PTR_OP", "-"); return(PTR_OP); }
"&&"			{ count(); log("AND_OP", "-"); return(AND_OP); }
"||"			{ count(); log("OR_OP", "-"); return(OR_OP); }
"<="			{ count(); log("LE_OP", "-"); return(LE_OP); }
">="			{ count(); log("GE_OP", "-"); return(GE_OP); }
"=="			{ count(); log("EQ_OP", "-"); return(EQ_OP); }
"!="			{ count(); log("NE_OP", "-"); return(NE_OP); }
";"			    { count(); log(";", "-"); return(';'); }
("{"|"<%")		{ count(); log("{", "-"); return('{'); }
("}"|"%>")		{ count(); log("}", "-"); return('}'); }
","			    { count(); log(",", "-"); return(','); }
":"			    { count(); log(":", "-"); return(':'); }
"="			    { count(); log("=", "-"); return('='); }
"("			    { count(); log("(", "-"); return('('); }
")"			    { count(); log(")", "-"); return(')'); }
("["|"<:")		{ count(); log("[", "-"); return('['); }
("]"|":>")		{ count(); log("]", "-"); return(']'); }
"."			    { count(); log(".", "-"); return('.'); }
"&"			    { count(); log("&", "-"); return('&'); }
"!"			    { count(); log("!", "-"); return('!'); }
"~"			    { count(); log("~", "-"); return('~'); }
"-"			    { count(); log("-", "-"); return('-'); }
"+"			    { count(); log("+", "-"); return('+'); }
"*"			    { count(); log("*", "-"); return('*'); }
"/"			    { count(); log("/", "-"); return('/'); }
"%"			    { count(); log("%", "-"); return('%'); }
"<"			    { count(); log("<", "-"); return('<'); }
">"			    { count(); log(">", "-"); return('>'); }
"^"			    { count(); log("^", "-"); return('^'); }
"|"			    { count(); log("|", "-"); return('|'); }
"?"			    { count(); log("?", "-"); return('?'); }

[ \t\v\n\f]		{ count(); }
.			{ /* Add code to complain about unmatched characters */ }

%%

int yywrap()
{
	return 1;
}


void comment()
{
	char c, prev = 0;

	std::cout << yytext << std::endl;

	while ((c = yyinput()) != 0) { // (EOF maps to 0)
        if (yyout != nullptr)
            fputc(c, yyout);
        if (c == '/' && prev == '*')
            return;
        prev = c;
    }
    yyerror("unterminated comment");
}

int column = 0;

void count() {
    if (yytext[0] == ' ') {
        //fprintf(fich, "Token encontrado en línea %d: ESPACIO\n",line);
    }
    else if (yytext[0] == '\t') {
        column += 8 - (column % 8);
        //fprintf(fich, "Token encontrado en línea %d: TABULADOR\n",line);
    }
    else if (yytext[0] == '\v') {
        //fprintf(fich, "Token encontrado en línea %d: TABULADOR VERTICAL\n",line);
    }
    else if (yytext[0] == '\n') {
        //fprintf(fich, "Token encontrado en línea %d: SALTO DE LÍNEA\n",line);
        column = 0;
        line ++;
    }
    else if (yytext[0] == '\f') {
        //fprintf(fich, "Token encontrado en línea %d: SALTO DE PÁGINA\n",line);
    }
    else{
        for (int i = 0; yytext[i] != '\0'; i++) {
            column++;
        }
        //fprintf(fich, "Token encontrado en línea %d: %s\n",line, yytext);
        // std::cout << yytext[i];
    }

}
/*
int check_type() {
    
    * pseudo code --- this is what it should check
    *
    *   if (yytext == type_name)
    *       return TYPE_NAME;
    *
    *   return IDENTIFIER;
    

    
    *   it actually will only return IDENTIFIER
    

    return IDENTIFIER;
}
*/
void log(const std::string& token, const std::string& value) {
     fich << "Token encontrado en línea " << line << ": <" << token << ", " << value << ">\n";
}