//Input for ReSoC: read all input data from input file(s)
//Created by Jun Li, Dimue company, Wuhan.
//Created: 2017-6-8
//Last modified: 2017-6-26

#ifndef RESOC_ULT_INPUTMANAGER_H
#define RESOC_ULT_INPUTMANAGER_H
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

class InputManager
{
public:
    InputManager();
    ~InputManager();
    bool readInputFile(string a_fileName);
    map<string, vector<string> > * getInputDataBlock();

private:
    map <string, vector<string> > m_inputDataBlock_;
//    string m_Buffer_;

};



#endif