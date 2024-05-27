#ifndef TEMPLATE_H
#define TEMPLATE_H


#include <iostream>
#include <string>
#include <string.h>
#include <cstdio>
#include <fstream>
#include "1805082_SymbolTable.h"
#include <vector>

using namespace std;

vector<const char *> args;

extern ofstream output, errFile;
extern int statementZone;
extern int prePragma;
extern SymbolTable table;

//-------------------------------AUXILIARES--------------------
void setPrePragma(){
    prePragma = 1;
    output << "\t}" << endl;
}

void addArg(const char *arg){
    
    args.push_back(arg);
    /*else{
        printf("Argument \"%s\" at pragma not in Symbol Table\n", arg);
        exit(1);
    }*/
}

vector<string> extraerValores(string cadena) {
    vector<string> valores;
    
    // Encontrar la posición del carácter '[' y ']'
    size_t posIni = cadena.find('[');
    size_t posFin = cadena.find(']');
    
    // Si se encontraron '[' y ']'
    if (posIni != string::npos && posFin != string::npos) {
        // Extraer la subcadena entre '[' y ']'
        string subcadena = cadena.substr(posIni + 1, posFin - posIni - 1);
        
        // Encontrar la posición del carácter '*'
        size_t posAsterisco = subcadena.find('*');
        
        // Si se encontró '*'
        if (posAsterisco != string::npos) {
            // Extraer las subcadenas de 'count', 'm' y 'n'
            string countStr = cadena.substr(0, posIni);
            string mStr = subcadena.substr(0, posAsterisco);
            string nStr = subcadena.substr(posAsterisco + 1);

            // Agregar las subcadenas al vector
            valores.push_back(countStr);
            valores.push_back(mStr);
            valores.push_back(nStr);
        }
    }
    
    return valores;
}

int contarCorchetes(string str) {
    int contador = 0;
    for (char c : str) {
        if (c == '[') {
            contador++;
        }
    }
    return contador;
}

//-------------------------------PLANTILLAS--------------------
void MPIInit(){
    string init=    "\tint __taskid = -1, __numprocs = -1;\n\n"
                    "\tMPI_Init(&argc,&argv);\n"
                    "\tMPI_Comm_size(MPI_COMM_WORLD,&__numprocs);\n" 
                    "\tMPI_Comm_rank(MPI_COMM_WORLD,&__taskid);\n"
                    "\tMPI_Get_processor_name(processor_name,&namelen);\n\n"
                    "\tif (__taskid == 0) {\n";
    output << init << endl;
    statementZone = 1;

}

void MPIFinalize(){
    string finalize=    "\tMPI_Finalize();\n";
    output << finalize << endl;
}


void MPIAlloc(){
    if(prePragma==0){
        setPrePragma();
    }

    if(args.size() != 1){
        errFile << "Error: Alloc pragma must have 1 argument and it has " << args.size() << endl;
        exit(1);
    }
    
    string alloc = "";
    string arg1 = args.at(0);
    string variableName = arg1.substr(0, arg1.find('['));
    SymbolInfo *sim1 = table.getSymbolInfo(variableName);

    int nCorchetes = contarCorchetes(arg1);
    if(nCorchetes==2){
        alloc= "\t" + sim1->getVariableType() + " __" + arg1 + ";\n";
    }else if (nCorchetes==1){
        vector<string> items = extraerValores(arg1);
        alloc =  "\tif (__taskid != 0)\n"
                        "\t\t" + items.at(0) + "= ( " + sim1->getVariableType() + " * ) malloc ( " + items.at(1) + " * " + items.at(2) + " * sizeof (TYPE) );\n";
    }else{
        errFile << "Error: MPIAlloc argument has wrong format" << endl;
        exit(1);
    }
    output << alloc << endl;
    args.clear();
}

void MPIBroad(){
    if(prePragma==0){
        setPrePragma();
    }
    
    if(args.size() != 2){
        errFile << "Error: Broad pragma must have 2 argument and it has " << args.size() << endl;
        exit(1);
    }

    
    string arg1 = args.at(0);
    string arg2 = args.at(1);

    SymbolInfo *sim1 = table.getSymbolInfo(arg1);
    SymbolInfo *sim2 = table.getSymbolInfo(arg2);

    string broad =  "\tMPI_Bcast(&" + arg1 + ", 1, MPI_" + sim1->getVariableType() + ", 0, MPI_COMM_WORLD);\n"
                    "\tMPI_Bcast(" + arg2 + ", N, MPI_" + sim2->getVariableType()  + ", 0, MPI_COMM_WORLD);\n";

    output << broad << endl;
    args.clear();
}
#endif 
