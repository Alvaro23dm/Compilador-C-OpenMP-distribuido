#ifndef SYMBOLTABLE_H
#define SYMBOLTABLE_H

#include<iostream>
#include<string>
#include<vector>
#include <fstream>
#include <sstream>
using namespace std;

bool verbose = true;
stringstream outputStream;

string getOutputStream(){
    string temp = outputStream.str();
    outputStream.str(std::string());
    return temp;
}
uint32_t sdbm(string str)
    {
        uint32_t hash = 0;
       
        for (char c : str)
            hash = c + (hash << 6) + (hash << 16) - hash;

        return hash;
    }


class ForPrinting{
    string to_print;
    bool baseRecursion = true;
    bool startRecursion = true;
   

    public:

    int exp2_start_marker;
    int exp3_start_marker;
    int stmt1_start_marker;
    int stmt1_end_marker;
    int stmt2_start_marker;
    int stmt2_end_marker;
   
    ForPrinting(){

    }

    bool getBaseRecursion(){
        return baseRecursion;
    }

    void setBaseRecursion(bool b){
        baseRecursion = b;
    }

    bool getStartRecursion(){
        return startRecursion;
    }

    void setStartRecursion(bool b){
        startRecursion = b;
    }

    void setToPrint(string s){
        this->to_print = s;
    }

    string getToPrint(){
        return this->to_print;
    }

    ForPrinting(const ForPrinting &old_obj){
        to_print = old_obj.to_print;
    };
};

class SymbolInfo{
    string* Name;
    string* Type;
    SymbolInfo* next;
    vector<string> func_parameters;
    string func_return_type;
    string var_type;
    string to_print;
    string Asm;
    // bool isIndexedVar = false;
    int  stackLoc = -42069;
    int arrSize = 0;
    bool funcDefined = false;
    bool isGlobalArray = false;

    public:

    int exp2_start_marker;
    int exp3_start_marker;
    int stmt1_start_marker;
    int stmt1_end_marker;
    int stmt2_start_marker;
    int stmt2_end_marker;

    SymbolInfo(const SymbolInfo &old_obj){
        Name = new string(*old_obj.Name);
        Type = new string(*old_obj.Type);
        func_parameters = old_obj.func_parameters;
        func_return_type = old_obj.func_return_type;
        var_type = old_obj.var_type;
        to_print = old_obj.to_print;
        Asm =  old_obj.Asm;
        stackLoc = old_obj.stackLoc;
        arrSize = old_obj.arrSize;
        funcDefined = old_obj.funcDefined;
        isGlobalArray = old_obj.isGlobalArray;
        next = nullptr;
    };


    void setIsGlobalArray(bool b){
        isGlobalArray = b;
    }

    bool getIsGlobalArray(){
        return isGlobalArray;
    }
    SymbolInfo(string Name, string Type){
        this->Name = new string(Name);
        this->Type = new string(Type);
        this->next = nullptr;
        this->var_type = "undeclared type";
    } 
    

    //NEW CODE
    bool getFuncDefined(){
        return this->funcDefined;
    }

    //NEW CODE
    void setFuncDefined(bool b){
        this->funcDefined = b;
    }

    // bool getIsIndexedVar(){
    //     return isIndexedVar;
    // }
    // bool setIsIndexedVar(bool b){
    //     this->isIndexedVar = b;
    // }

    ~SymbolInfo(){
        delete Name;
        delete Type;
    } 

    void setStackLoc(int s){
        stackLoc = s;
    }

    int getStackLoc(){
        return stackLoc;
    }
    void setToPrint(string s){
        this->to_print = s;
    }

    string getToPrint(){
        return this->to_print;
    }

    void setAsm(string s){
        this->Asm = s;
    }

    string getAsm(){
        return this->Asm;
    }


    void setArrSize(int i){
        arrSize = i;
    }

    void setName(string s){
        this->Name = new string(s);
    }

    int getArrSize(){
        return arrSize;
    }

    string getType(){
        return *Type;
    }
    void setType(string Type){
        this->Type = new string(Type);
    }
    void setFuncParameters(vector<string> parameters){
        func_parameters = parameters;
    }

    vector<string> getFuncParameters(){
        return func_parameters;
    }

    void setVarType(string type){
        var_type = type;
    }

    string getVarType(){
        if (arrSize!=0)
            return var_type+" array";

        if(var_type.empty())
            return "Undefined type";
        return var_type;
    }

    string getElementType(){
        if (arrSize==0)
            return "Error: Not an array";
        return var_type;
    }
    string getVarType(bool debug){

        if (arrSize!=0)
            return var_type+" array";
            
        return var_type;
    }

    void setFuncReturnType(string type){
        func_return_type = type;
    }

    string getFuncReturnType(){
        return func_return_type;
    }






    void setNext(SymbolInfo* next){
        this->next = next;
    }

    SymbolInfo* getNext(){

        return this->next;
    }

    string getName(){
        return *Name;
    }

    void printSymbol(){
      
        outputStream<<"< "<<*Name<<" : "<<*Type<<", "<<getVarType()<<" > ";

    }

    

    void printEverything(){
        outputStream<< *Name <<endl;
        outputStream<< *Type <<endl;
        outputStream<< func_return_type <<endl;
        outputStream<< var_type <<endl;
        outputStream<< to_print <<endl;
        outputStream<< Asm <<endl;
        outputStream<< stackLoc <<endl;
        outputStream<< arrSize <<endl;
        outputStream<< funcDefined <<endl;
    }



};


class ScopeTable{

    ScopeTable* parentScope;
  
    int n;
    int numChildren;
    string id;

    public:

    vector<SymbolInfo* > table;

    string getID(){
        return this->id;
    }

    int getNumChildren(){
        return this->numChildren;
    }

    ScopeTable(int n, string parentID, int order){
        this->numChildren = 0;
        this->n = n;
        this->id = parentID+to_string(order);
        for(int i=0;i<n;i++){
            table.push_back(nullptr);
        }
    }

    void incrementNumChildren(){
        this->numChildren ++;
    }

    ~ScopeTable(){
        
        for (std::vector<SymbolInfo*>::iterator i = table.begin(); i != table.end(); ++i) {

        SymbolInfo* p = *i;

        while(p!=nullptr){
                SymbolInfo* next = p->getNext();
                delete p;
                p = next;
            }

        }

        table.clear();
        // delete[] table;
     }




    void setParentScope(ScopeTable* parentScope){
        this->parentScope = parentScope;
    }

    ScopeTable* getParentScope(){
        return this->parentScope;
    }

    bool insertSymbol(string Name, string Type){
        outputStream.str(std::string());
        int k = sdbm(Name)%this->n;
       

        int offset = 0;
        
        
        if (this->table[k] == nullptr)
            this->table[k] = new SymbolInfo(Name, Type);
        
        else{
        
            SymbolInfo* current = this->table[k];
            offset = 1;

            while (current->getNext()!= nullptr){
                if(current->getName().compare(Name)==0){
                    if(verbose){
                    current->printSymbol();
                    outputStream<<" already exists in current ScopeTable"<<endl;
                    }
                    return false;
                }

                offset++;
                current = current->getNext();
            }

            if(current->getName().compare(Name)==0){
                if(verbose){
                        current->printSymbol();
                        outputStream<<" already exists in current ScopeTable"<<endl;
                }
                    return false;
                }


            current->setNext(new SymbolInfo(Name, Type));
        }
        
        string to_print = "Inserted  <" + Name +" , " + Type + "> in ScopeTable# "+this->getID()+" at position "+to_string(k)+", "+ to_string(offset);
        if (verbose)
            outputStream<< to_print <<endl;
        return true;
            
    }
    
    bool deleteSymbol(string Name){

        int k = sdbm(Name)%this->n;
        int offset = 0;
        //list empty
        if (this->table[k] == nullptr){
            if (verbose)
                outputStream<<Name<<" not found"<<endl;
            return false;
        }
    

        SymbolInfo* current = this->table[k];

        //first Item is to be removed
        if (current->getName().compare(Name)==0){
            this->table[k] = current->getNext();
            delete current; 
            if (verbose)
             outputStream<<"Deleted Entry " << k<<", "<< offset <<" from current ScopeTable"<<endl;
            return true;
        }


        bool found_flag = false;
        SymbolInfo* before = nullptr;
        

        offset = 0;
        while (current != nullptr ){
            
            if (current->getName().compare(Name)==0){
                found_flag = true;
                break;
            }
            
            before = current;
            current = current->getNext();
            offset++;
        }

        if (found_flag){
            before->setNext(current->getNext());
            delete current; 

            if (verbose)
             outputStream<<"Deleted Entry " << k<<", "<< offset <<" from current ScopeTable"<<endl;
        }
        else{
            if (verbose)
             outputStream<<Name<<" not found"<<endl;
        }

        return found_flag;

    }
    
    SymbolInfo* lookupSymbol(string Name){
        
        int k = sdbm(Name)%this->n;
       
        int offset = 0;
        
        SymbolInfo* current = this->table[k];

        while (current != nullptr ){

            if (current->getName().compare(Name)==0)
                break;
            
            current = current->getNext();
            offset++;
        }
        
        if (current!=nullptr)
            if (verbose)
             outputStream<<"Found in ScopeTable# "+this->getID()+" at position "+to_string(k)+", "+ to_string(offset)<<endl;
        // else
        //     outputStream<<"Not found in Current Scope Table"<<endl;
        return current;


    }

    

    void printOneRow(int k){
        if (this->table[k] == nullptr)
            return;

        SymbolInfo* current = this->table[k];

        while (current != nullptr ){
         
            current -> printSymbol(); 
            outputStream<<" ";
            current = current->getNext();
        }

        

        
    }
    void printScopeTable(){

        outputStream<<"ScopeTable# "<<this->getID()<<endl;

        for (int i=0;i<n;i++){
            if (this->table[i] == nullptr)
                continue;;
            outputStream<<i<<"--> ";
            printOneRow(i);
            outputStream<<endl;

        }

        outputStream<<endl;
    }
};


class SymbolTable{
    ScopeTable* current;
    int n; //SCOPETABLE SIZE
    int globalScopeCount;

    public:

    SymbolTable(int n){
        this->n = n;
        this->globalScopeCount = 1;
        current = new ScopeTable(n, "", globalScopeCount);
        current->setParentScope(nullptr);
    }


    void enterScope(){

        string parentID;
        int order;
        if (current != nullptr){
            parentID = current->getID()+".";
            current->incrementNumChildren();
            order = current->getNumChildren();
        }  
        else{
            parentID = "";
            this->globalScopeCount++;
            order = this->globalScopeCount;  
        }

        ScopeTable* temp = current;
       

        current = new ScopeTable(n, parentID, order);
    
        current->setParentScope(temp);
      
        if (verbose)
             outputStream<<"New ScopeTable with id "+current->getID() + " created"<<endl;
    }

    void printCurrentScopeID(){
        if (current == nullptr)
            return;

        outputStream<<current->getID()<<endl;
    }
    string getCurrentScopeID(){
        if (current == nullptr)
            return "";

        return current->getID();
    }


    void exitScope(){
        if (current == nullptr){
            outputStream<<"No current Scope. Cannot exit."<<endl;
            return;
        }

        ScopeTable* temp = current->getParentScope();
        if (verbose)
             outputStream<<"ScopeTable with id "+current->getID()+" removed"<<endl;
        delete current; 
        current = temp;
    }

    bool insert(string Name, string Type){
        
        if (current == nullptr){
            outputStream<<"No current Scope. Cannot Insert."<<endl;
            return false;
        }
        return current->insertSymbol(Name,Type);
    }

    bool remove(string Name){
        if (current == nullptr){
            outputStream<<"No current Scope. Cannot Remove."<<endl;
            return false;
        }

        bool flag = current->deleteSymbol(Name);
     
        return flag;
    }

    SymbolInfo* lookup(string Name){
        if (current == nullptr){
            outputStream<<"No current Scope. Cannot Lookup."<<endl;
            return nullptr;
        }

        ScopeTable* toSearch = current;

        while (toSearch != nullptr){
            SymbolInfo* returnVal = toSearch->lookupSymbol(Name);
            if (returnVal!=nullptr)
                return returnVal;
            toSearch = toSearch->getParentScope();
           
        }

        outputStream<<"Not found"<<endl;
        return nullptr;
    }

    void printCurrentScopeTable(){
        outputStream.str(std::string());
        if (current == nullptr){
            outputStream<<"No current Scope. Cannot Print."<<endl;
            return;
        }
        current->printScopeTable();
    }
    

    void recurseTablesForPrinting(ScopeTable* x){
        if (x==nullptr)
            return;
        
        x->printScopeTable();
        recurseTablesForPrinting(x->getParentScope());

    };


    void printAllScopeTable(){
        outputStream.str(std::string());
 
        if (current == nullptr){
            outputStream<<"No current Scope. Cannot Print."<<endl;
            return;
        }
        ScopeTable* x = current;
        recurseTablesForPrinting(x);
    }





};


vector<string> stringSplit(string str){
    stringstream ss(str);
    string temp;
    vector<string> words;


    while(ss){
        ss>>temp;
        words.push_back(temp);
    }
    return words;
}
/*
int main(){

    string str;
    int n;
 
    
    ifstream inputFile("1805053input.txt"); 

    inputFile >> n;
    SymbolTable table(n);
    

   while (getline (inputFile, str))  {
        
        vector<string> words = stringSplit(str);

        if (words[0].compare("I")==0){
            string name = words[1];
            string value = words[2];
            table.insert(name,value);
        }
        else if (words[0].compare("L")==0){
            string temp = words[1];
            table.lookup(temp);
        }
        else if (words[0].compare("D")==0){
            string temp = words[1];
            table.remove(temp);
        }
        else if (words[0].compare("P")==0){
            if (words[1].compare("A")==0)
                table.printAllScopeTable();
            else if (words[1].compare("C")==0)
                table.printCurrentScopeTable();
        }
        else if (words[0].compare("S")==0){
            table.enterScope();
        }
        else if (words[0].compare("E")==0){
            table.exitScope();
        }


        
       
       
    }
    inputFile.close();

}
*/
#endif