#include <vector>
#include <fstream>
#include "Framework/Util/UtilLog.h"
#include "Framework/Util/UtilString.h"
#include "InputManager.h"

InputManager::InputManager()
{

}

InputManager::~InputManager()
{

}

map<string, vector<string> > * InputManager::getInputDataBlock()
{
    int inputDataSize = (int)m_inputDataBlock_.size();
    if (inputDataSize > 0)
    {
        return &m_inputDataBlock_;
    }
    else
    {
		LOG(ERROR) << "No Data! Please check!" << endl;
        exit(-1);
    }
    
}

bool InputManager::readInputFile(string a_fileName)
{
	ifstream inFile;
    inFile.open(a_fileName.c_str(), fstream::in);
    string templeWord;
    string wholeLine;
    vector<string> lineString;
    vector<string> linesOfDataBlock;
    int indexDataBolck = 0;
	int l_numberZones = 0;
    if (inFile.is_open())
    {
        while (!inFile.eof())
        {
            getline(inFile, wholeLine);
            stringstream stringin(wholeLine);
            lineString.clear();
            while (stringin>>templeWord)
            {
                lineString.push_back(templeWord);
            }
            if (lineString.size() > 0 && lineString[0].substr(0, 2) == "--")
            {
                continue;
            }

            //NO contents!
            if (lineString.size() == 0)
            {
                continue;
            }

            if (lineString.size() >= 2 && lineString[1] == "BEGIN")
            {
				string l_keyString = lineString[0];
				l_keyString.erase(l_keyString.end()-1);
                if (lineString[0]=="PVT")
                {
    //                dataBlock.dataBlockName = "PVT";
                    linesOfDataBlock.clear();
                    indexDataBolck = 5;
                    continue;
                }
			    else if (lineString[0] == "RELPERM" || l_keyString == "RELPERM")
			    {
				    //dataBlock.dataBlockName = "PROJ";
					l_numberZones += 1;
				    linesOfDataBlock.clear();
				    indexDataBolck = 9;
				    continue;
			    }
			    else if (lineString[0] == "CAPPRESS" || l_keyString == "CAPPRESS")
			    {
				    //dataBlock.dataBlockName = "CAPPRESS";
				    linesOfDataBlock.clear();
				    indexDataBolck = 10;
				    continue;
			    }
                else if (lineString[0] == "FRAC_CAPPRESS")
                {
                    //dataBlock.dataBlockName = "FRAC_CAPPRESS";
                    linesOfDataBlock.clear();
                    indexDataBolck = 12;
                    continue;
                }
                else if (lineString[0] == "PROJ")
                {
    //                dataBlock.dataBlockName = "PROJ";
                    linesOfDataBlock.clear();
                    indexDataBolck = 1;
                    continue;
                }

                else if (lineString[0] == "GEO")
                {
					//dataBlock.dataBlockName = "GEO";
                    linesOfDataBlock.clear();
                    indexDataBolck = 2;
                    continue;
                }
                else if (lineString[0] == "ROCK")
                {
                    //dataBlock.dataBlockName = "GEO";
                    linesOfDataBlock.clear();
                    indexDataBolck = 11;
                    continue;
                }
                else if (lineString[0] == "PRO")
                {
     //               dataBlock.dataBlockName = "PRO";
                    linesOfDataBlock.clear();
                    indexDataBolck = 3;
                    continue;
                }

                else if (lineString[0] == "SAT")
                {
     //               dataBlock.dataBlockName = "SAT";
                    linesOfDataBlock.clear();
                    indexDataBolck = 4;
                    continue;
                }

                else if (lineString[0] == "INIDATA")
                {
      //              dataBlock.dataBlockName = "INI";
                    linesOfDataBlock.clear();
                    indexDataBolck = 6;
                    continue;
                }

                else if (lineString[0] == "WELL")
                {
    //                dataBlock.dataBlockName = "WELL";
                    linesOfDataBlock.clear();
                    indexDataBolck = 7;
                    continue;
                }

                else if (lineString[0] == "SCHEDULE")
                {
    //                dataBlock.dataBlockName = "SCHEDULE";
                    linesOfDataBlock.clear();
                    indexDataBolck = 8;
                    continue;
                }

                else if (lineString[0] == "GCM")
                {
                    linesOfDataBlock.clear();
                    indexDataBolck = 11;
                    continue;
                }

                else if (lineString[0] == "THERMAL")
                {
                    linesOfDataBlock.clear();
                    indexDataBolck = 13;
                    continue;
                }

                else
                {
                    LOG(ERROR) << "Wrong keywords" << endl;
                    return false;
                }
            }

            if (lineString.size() >= 2 && lineString[1] == "END")
            {
				string l_keyString = lineString[0];
				l_keyString.erase(l_keyString.end() - 1);
                if (lineString[0] == "PVT")
                {
                    if (indexDataBolck != 5)
                    {
                        LOG(ERROR) << "Wrong inputs for PVT datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("PVT", linesOfDataBlock));
                    }
                    continue;
                }

			    else if (lineString[0] == "RELPERM" || l_keyString == "RELPERM")
			    {
				    if (indexDataBolck != 9)
				    {
                        LOG(ERROR) << "Wrong inpts for RELPERM datablock!" << endl;
					    return false;
				    }
				    else
				    {
					    indexDataBolck = 0;
					    m_inputDataBlock_.insert(map<string, vector<string> >::value_type(lineString[0], linesOfDataBlock));
				    }
				    continue;
			    }
			    else if (lineString[0] == "CAPPRESS" || l_keyString == "CAPPRESS")
			    {
				    if (indexDataBolck != 10)
				    {
                        LOG(ERROR) << "Wrong inpts for CAPPRESS datablock!" << endl;
					    return false;
				    }
				    else
				    {
					    indexDataBolck = 0;
					    m_inputDataBlock_.insert(map<string, vector<string> >::value_type(lineString[0], linesOfDataBlock));
				    }
				    continue;
			    }
                else if (lineString[0] == "FRAC_CAPPRESS")
                {
                    if (indexDataBolck != 12)
                    {
                        LOG(ERROR) << "Wrong inpts for FRAC_CAPPRESS datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("FRAC_CAPPRESS", linesOfDataBlock));
                    }
                    continue;
                }
                else if (lineString[0] == "PROJ")
                {
                    if (indexDataBolck != 1)
                    {
                        LOG(ERROR) << "Wrong inpts for PROJ datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("PROJ", linesOfDataBlock));
                    }
                    continue;
                }

                else if (lineString[0] == "GEO")
                {
                    if (indexDataBolck != 2)
                    {
                        LOG(ERROR) << "Wrong inputs for GEO datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("GEO", linesOfDataBlock));
                    }
                    continue;
                }
                else if (lineString[0] == "ROCK")
                {
                    if (indexDataBolck != 11)
                    {
                        LOG(ERROR) << "Wrong inputs for ROCK datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("ROCK", linesOfDataBlock));
                    }
                    continue;
                }
                else if (lineString[0] == "PRO")
                {
                    if (indexDataBolck != 3)
                    {
                        LOG(ERROR) << "Wrong inputs for PRO datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("PRO", linesOfDataBlock));
                    }
                    continue;
                }

                else if (lineString[0] == "SAT")
                {
                    if (indexDataBolck != 4)
                    {
                        LOG(ERROR) << "Wrong inputs for GEO datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("SAT", linesOfDataBlock));
                    }
                    continue;
                }

                else if (lineString[0] == "INIDATA")
                {
                    if (indexDataBolck != 6)
                    {
                        LOG(ERROR) << "Wrong inputs for INI datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("INIDATA", linesOfDataBlock));
                    }
                    continue;
                }

                else if (lineString[0] == "WELL")
                {
                    if (indexDataBolck != 7)
                    {
                        LOG(ERROR) << "Wrong inputs for WELL datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("WELL", linesOfDataBlock));
                    }
                    continue;
                }

                else if (lineString[0] == "SCHEDULE")
                {
                    if (indexDataBolck != 8)
                    {
                        LOG(ERROR) << "Wrong inputs for SCHEDULE datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("SCHEDULE", linesOfDataBlock));
                    }
                    continue;
                }

                else if (lineString[0] == "GCM")
                {
                    if (indexDataBolck != 11)
                    {
                        LOG(ERROR) << "Wrong inputs for GCM datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("GCM", linesOfDataBlock));
                    }
                }
                else if (lineString[0] == "GCM")
                {
                    if (indexDataBolck != 13)
                    {
                        LOG(ERROR) << "Wrong inputs for THERMAL datablock!" << endl;
                        return false;
                    }
                    else
                    {
                        indexDataBolck = 0;
                        m_inputDataBlock_.insert(map<string, vector<string> >::value_type("THERMAL", linesOfDataBlock));
                    }
                }
                else
                {
                    LOG(ERROR) << "Wrong keywords" << endl;
                    return false;
                }
            }
            wholeLine = UtilString::trim(wholeLine);
            linesOfDataBlock.push_back(wholeLine);
       
         }
		 linesOfDataBlock.clear();
		 linesOfDataBlock.push_back(static_cast<ostringstream*>(&(ostringstream() << l_numberZones))->str());
		 m_inputDataBlock_.insert(map<string, vector<string> >::value_type("NUMZONES", linesOfDataBlock));
         inFile.close();
         return true;
    }
    else
    {
        LOG(ERROR) << "file can not be opened" << endl;
        return false;
    }
    
}
