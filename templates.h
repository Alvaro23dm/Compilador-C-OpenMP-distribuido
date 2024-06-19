#ifndef TEMPLATE_H
#define TEMPLATE_H


#include <iostream>
#include <string>
#include <string.h>
#include <cstdio>
#include <fstream>
#include "SymbolTable.h"
#include <vector>

using namespace std;

vector<const char *> args;

extern ofstream output, errFile;

extern int statementZone;
extern int zonaPragma;
extern int finalizeOK;
extern int contadorTask;

extern int chunk;
extern int task;

extern SymbolTable table;

//-------------------------------AUXILIARES--------------------
void finSecuencial(){
    zonaPragma = 1;
    output << "\t}" << endl;
}

void addArg(const char *arg){
    
    args.push_back(arg);
    /*else{
        printf("Argument \"%s\" at pragma not in Symbol Table\n", arg);
        exit(1);
    }*/
}

vector<string> extractValues(string str) {
    vector<string> values;
    
    // Find the position of '[' and ']'
    size_t startPos = str.find('[');
    size_t endPos = str.find(']');
    
    // If '[' and ']' are found
    if (startPos != string::npos && endPos != string::npos) {
        // Extract the substring between '[' and ']'
        string subStr = str.substr(startPos + 1, endPos - startPos - 1);
        
        // Find the position of '*'
        size_t asteriskPos = subStr.find('*');
        
        // If '*' is found
        if (asteriskPos != string::npos) {
            // Extract the substrings for 'count', 'm', and 'n'
            string countStr = str.substr(0, startPos);
            string mStr = subStr.substr(0, asteriskPos);
            string nStr = subStr.substr(asteriskPos + 1);

            // Add the substrings to the vector
            values.push_back(countStr);
            values.push_back(mStr);
            values.push_back(nStr);
        }
    }
    
    return values;
}

int countBrackets(string str) {
    int counter = 0;
    for (char c : str) {
        if (c == '[') {
            counter++;
        }
    }
    return counter;
}

vector<string> splitMatAIndices(const string& input) {
    size_t start_bracket1 = input.find('[');
    size_t end_bracket1 = input.find(']');
    size_t start_bracket2 = input.find('[', end_bracket1);
    size_t end_bracket2 = input.find(']', start_bracket2);

    string matA = input.substr(0, start_bracket1);
    string fA = input.substr(start_bracket1 + 1, end_bracket1 - start_bracket1 - 1);
    string cA = input.substr(start_bracket2 + 1, end_bracket2 - start_bracket2 - 1);

    return {matA, fA, cA};
}

void preConfPragma(){
    //statementZone = 1; //puede dar problemas
    if(zonaPragma==0){
        finSecuencial();
    }
}

vector<string> inTask(){
    int in = 0;
    vector<string> result;
    for(auto& arg : args){
        if(strcmp(arg, "IN") == 0){
            in = 1;
        }else if(strcmp(arg, "OUT") == 0){
            in = 0;
        }
        else if(in == 1){
            result.push_back(arg);
        }
    }
    return result;
}

vector<string> outTask(){
    int out = 0;
    vector<string> result;
    for(const auto& arg : args){
        if(strcmp(arg, "OUT") == 0){
            out = 1;
        }else if(strcmp(arg, "IN") == 0){
            out = 0;
        }
        else if(out == 1){
            result.push_back(arg);
        }
    }
    return result;
}

vector<string> argFormat(string arg) {
    vector<string> result;

    size_t startBracketPos = arg.find('[');
    size_t colonPos = arg.find(':');
    size_t endBracketPos = arg.find(']');

    if (startBracketPos != string::npos && colonPos != string::npos && endBracketPos != string::npos) {
        string part1 = arg.substr(0, startBracketPos);
        result.push_back(part1);

        string part2 = arg.substr(startBracketPos + 1, colonPos - startBracketPos - 1);
        result.push_back(part2);

        string part3 = arg.substr(colonPos + 1, endBracketPos - colonPos - 1);
        result.push_back(part3);
    } else {
        cerr << "Invalid format: " << arg << endl;
    }

    return result;
}

vector<string> scatterArgs(){
	vector<string> result = splitMatAIndices(args.at(0));
	result.push_back(args.at(1));
	for (size_t i = 0; i < result.size(); ++i) {
        cout << result[i] << endl;
    }
	return result;
}

string aMayuscula(string cadena) {
  for (string::size_type i = 0; i < cadena.length(); i++) cadena[i] = toupper(cadena[i]);
  return cadena;
}
string aMinuscula(string cadena) {
  for (string::size_type i = 0; i < cadena.length(); i++) cadena[i] = tolower(cadena[i]);
  return cadena;
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
    if(zonaPragma==0){
        finSecuencial();
    }
    finalizeOK = 1;
    string finalize=    "\tMPI_Finalize();\n";
    output << finalize << endl;
}


void MPIAlloc(){
    preConfPragma();

    if(args.size() != 1){
        errFile << "Error: Alloc pragma must have 1 argument and it has " << args.size() << endl;
        exit(1);
    }
    
    string alloc = "";
    string arg1 = args.at(0);
    string variableName = arg1.substr(0, arg1.find('['));
    SymbolInfo *sim1 = table.getSymbolInfo(variableName);

    int nCorchetes = countBrackets(arg1);
    if(nCorchetes==2){
        alloc= "\t" + sim1->getVariableType() + " __" + arg1 + ";\n";
    }else if (nCorchetes==1){
        vector<string> items = extractValues(arg1);
        alloc =  "\tif (__taskid != 0)\n"
                        "\t\t" + items.at(0) + "= ( " + aMinuscula(sim1->getVariableType()) + " * ) malloc ( " + items.at(1) + " * " + items.at(2) + " * sizeof (" + aMinuscula(sim1->getVariableType()) + ") );\n";
    }else{
        errFile << "Error: MPIAlloc argument has wrong format" << endl;
        exit(1);
    }
    output << alloc << endl;
    args.clear();
}

void MPIBroad(){
    preConfPragma();

	string broad = "";

    for(const auto& arg : args){
		SymbolInfo *sim = table.getSymbolInfo(arg);
		if(sim->isArray()){
			broad += "\tMPI_Bcast(" + string(arg) + ", " + sim->getSizeList() + ", MPI_" + sim->getVariableType()  + ", 0, MPI_COMM_WORLD);\n";
		}
		else{
			broad += "\tMPI_Bcast(&" + string(arg) + ", 1, MPI_" + sim->getVariableType() + ", 0, MPI_COMM_WORLD);\n";
		}
	}
    output << broad << endl;
    args.clear();
}

void MPIScatterChunk(){
    preConfPragma();
	
    vector<string> argsScatter = scatterArgs();
    
	string uno = argsScatter.at(0);//matA
	string dos = argsScatter.at(1);//fA
	string tres = argsScatter.at(2);//cA
	string cuatro = argsScatter.at(3);//cA

	SymbolInfo *sim1 = table.getSymbolInfo(uno);
    
    string scatter =    
        "\t{\n"
        "\t\tint __" + dos + "_chunk;\n"
        "\t\t" + aMinuscula(sim1->getVariableType()) + " **__" + uno + ";\n"
        "\t\tint *displsA = (int *)malloc(__numprocs * sizeof(int));\n"
        "\t\tint *countsA = (int *)malloc(__numprocs * sizeof(int));\n"
        "\t\tint offset;\n"
        "\n"
        "\t\tif (" + dos + " % __numprocs == 0)\n"
        "\t\t\t__" + dos + "_chunk = " + dos + " / __numprocs;\n"
        "\t\telse\n"
        "\t\t\t__" + dos + "_chunk = " + dos + " / __numprocs + 1;\n"
        "\t\tif (__taskid == 0) {\n"
        "\t\t\toffset = 0;\n"
        "\t\t\tfor (i = 0; i < __numprocs; i++) {\n"
        "\t\t\t\tcountsA[i] = __" + dos + "_chunk * " + cuatro + ";\n"
        "\t\t\t\tdisplsA[i] = offset;\n"
        "\t\t\t\toffset += __" + dos + "_chunk * " + cuatro + ";\n"
        "\t\t\t}\n"
        "\t\t\tcountsA[__numprocs - 1] = (" + dos + " - (__" + dos + "_chunk * (__numprocs - 1))) * " + cuatro + ";\n"
        "\t\t}\n"

        "\t\telse if (__taskid == __numprocs - 1)\n"
        "\t\t\tcountsA[__taskid] = (" + dos + " - (__" + dos + "_chunk * (__numprocs - 1))) * " + cuatro + ";\n"
        "\t\telse\n"
        "\t\t\tcountsA[__taskid] = __" + dos + "_chunk * " + cuatro + ";\n"
        "\t\tif (__taskid != 0) {\n"
        "\t\t\t__" + uno + " = MATRIX2D(sizeof(" + aMinuscula(sim1->getVariableType()) + "), " + dos + ", " + tres + ");\n"
        "\t\t}\n"
        "\t\t" + uno + " = MATRIX2D(sizeof(" + aMinuscula(sim1->getVariableType()) + "), __" + dos + "_chunk * __numprocs, " + tres + ");\n"
        "\n"

        "\t\t" + aMinuscula(sim1->getVariableType()) + " *" + uno + "aux = &" + uno + "[0][0];\n"
        "\t\t" + aMinuscula(sim1->getVariableType()) + " *__" + uno + "aux = &__" + uno + "[0][0];\n"
        "\n"
		
        "\t\tMPI_Scatterv(" + uno + "aux, countsA, displsA, MPI_" + sim1->getVariableType() + ", __" + uno + "aux + (__" + dos + "_chunk * " + tres + " * __taskid), countsA[__taskid], MPI_" + sim1->getVariableType() + ", 0, MPI_COMM_WORLD);\n"
        "\t}\n";
	output << scatter << endl;
	args.clear();
}

void MPITask(){
    preConfPragma();
    task = 1;
    string task =  
        "\tif (__taskid == " + to_string(contadorTask) + ") {\n";
    contadorTask++;

    vector<string> inArg = inTask();

    for(const auto& arg : inArg){
        vector<string> argItem = argFormat(arg);
        SymbolInfo *sim = table.getSymbolInfo(argItem.at(0));

        task += "\t\tMPI_Recv (&" + argItem.at(0) + "[" + argItem.at(1) + "]" +", " + argItem.at(2) + ", MPI_" +  sim->getVariableType() + ", prev, 0, MPI_COMM_WORLD, &__status);\n";
    }

    output << task;

}

void MPITaskEnd(){
    preConfPragma();
    string task =  "";

    vector<string> outArg = outTask();

    for(const auto& arg : outArg){
        vector<string> argItem = argFormat(arg);
        SymbolInfo *sim = table.getSymbolInfo(argItem.at(0));
        
        task += "\t\tMPI_Send (&" + argItem.at(0) + "[" + argItem.at(1) + "]" +", " + argItem.at(2) + ", MPI_" +  sim->getVariableType() + ", sig, 0, MPI_COMM_WORLD);\n";
    }
    task += "\t}";

    output << task << endl;
    args.clear();
}
#endif 
