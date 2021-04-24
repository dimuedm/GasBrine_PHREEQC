
#include "GeoChemInputs.h"


GeoChemInputs::GeoChemInputs()
{
    m_isKinetics = false;

    /*
    0 - first colume is time;
    1 - first colume is temperature;
    2 - first colume is pressure;
    */
    m_printType = -1; 

    m_timeListIn_ = false;
    m_redox_ = false;
    m_printInitialEquil_ = false;

    m_gasMoleNumber_ = 0.0;
}

GeoChemInputs::~GeoChemInputs()
{
}

bool GeoChemInputs::readInputs(map<string, vector<string> >* a_dataBlock)
{
    int inputDataBlockSize = (int)a_dataBlock->size();
    map <string, vector<string> >::iterator it;
    vector<string>* v_GCMData;

    if (inputDataBlockSize <= 0)
    {
        cerr << "No data! Please check!" << endl;
        return false;
    }
    else
    {
        it = a_dataBlock->find("GCM");

        if (it == a_dataBlock->end())
        {
            cerr << "There is no GCM data block!" << endl;
            return false;
        }
        else
        {
            v_GCMData = &(it->second);
            if (m_readInputs_(v_GCMData))
            {
                return true;
            }
            else
            {
                cerr << "Read GCM data error! Please check!"<<endl;
                return false;
            }
        }
    }
}

bool GeoChemInputs::m_readInputs_(vector<string>* a_v_GCMData)
{
    int iGCMDataSize = (int)a_v_GCMData->size();
    vector< vector<double> > doubleData;
    vector< vector<int> > intData;
    bool databaseExist = false;
 
    bool readInSuccess;
    int rowNumber;
    int columNumber;
  
    bool phrqcInputFileExist = 0;
    if (iGCMDataSize <= 0)
    {
        cerr << "ERROR: GCM: NO GCM data!" << endl;
        return false;
    }

    for (int i = 0; i < iGCMDataSize; i++)
    {

        if ((*a_v_GCMData)[i] == "DATABASE")
        {
            vector<vector<string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber != 1)
            {
                cerr << "ERROR: GCM: Error database file clarification!" << endl;
                return false;
            }
            m_databaseFile_ =  UtilString::trim(stringData[0][0]);

            databaseExist = 1;
        }

        if ((*a_v_GCMData)[i] == "TIME")
        {
            if (m_timeListIn_)
            {
                cerr << "ERROR: time list existing already! Check the inputs!" << endl;
                return false;
            }
            m_timeListIn_ = true;
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error simulation time data inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = doubleData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error simulation time data inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_simTime_.push_back(doubleData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRQCINPUT")
        {
            vector<vector<string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber != 1)
            {
                cerr << "ERROR: GCM: Error geochemistry input file clarification!" << endl;
                return false;
            }

            m_phrqc_inputFile_ = UtilString::trim(stringData[0][0]);
            phrqcInputFileExist = 1;
        }
       
        if ((*a_v_GCMData)[i] == "GASMOLEFRACTION")
        {
            vector<vector<string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);
            rowNumber = (int)stringData.size();
            for (int j = 0; j < rowNumber; j++)
            {
                if (stringData[j].size() < 2)
                {
                    cerr << "ERROR: GCM: imcomplete mineral molar volume information! Please check input file!" << endl;
                }
                string l_gasName = stringData[j][0];
                float l_gasMoleFrac;
                istringstream ss(stringData[j][1]);
                ss >> l_gasMoleFrac;
                m_gasMoleFraction_[l_gasName] = l_gasMoleFrac;
            }
        }

        if ((*a_v_GCMData)[i] == "GASMOLENUMBERS")
        {
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            if (doubleData.size() <= 0)
            {
                cerr << "ERROR: GCM: Error gas mole number inputs! Please check!" << endl;
                exit(-1);
                return false;
            }
            m_gasMoleNumber_ = doubleData[0][0];

            if (m_gasMoleNumber_ < 0)
            {
                cout << "WARNING: GCM: negative gas mole number input. Set to 0.0. Please check!" << endl;
            }
        }

        if ((*a_v_GCMData)[i] == "TEMPERATURE_PRESSURE_LIST")
        {
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error temperature-pressure list data inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = doubleData[k].size();
                if (columNumber != 2)
                {
                    cerr << "ERROR: GCM: error input for temperature-pressure list! Check the inputs!" << endl;
                    return false;
                }

                pair<double, double> l_temp_press;
                l_temp_press.first = doubleData[k][0];
                l_temp_press.second = doubleData[k][1];
                m_temperature_pressureList_.push_back(l_temp_press);
            }
        }

        if ((*a_v_GCMData)[i] == "GASMOLEFRACTION_LIST")
        {
            vector<string> componentName;
            i++;
            stringstream split;
            string tempWord;
            split.str((*a_v_GCMData)[i]);
            while (split >> tempWord)
            {
                componentName.push_back(tempWord);
            }
            m_gasComponentName_ = componentName;

            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error gas mole fraction list data inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0; k < rowNumber; k++)
            {
                columNumber = doubleData[k].size();
                if (columNumber != componentName.size())
                {
                    cerr<< "ERROR: GCM: Error gas mole fraction list data inputs! Please check!" << endl;
                    return false;
                }
                m_gasMoleFractionList_.push_back(doubleData[k]);
            }
        }

        if ((*a_v_GCMData)[i] == "PRESSURE_LIST") //unit: bar
        {
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error pressure list data inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = doubleData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error simulation time data inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_pressureList_.push_back(doubleData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "TEMPERATURE_LIST") //Unit: K
        {
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error temperature list data inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = doubleData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error temperature list data inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_temperatureList_.push_back(doubleData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "TIME_STEPS")
        {
            if (m_timeListIn_)
            {
                cerr << "ERROR: time list existing already! Check the inputs!" << endl;
                return false;
            }
            m_timeListIn_ = true;
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error temperature list data inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = doubleData[k].size();
                if (columNumber != 2)
                {
                    cerr << "ERROR: GCM: There should be 3 data for time step inputs! Please check!" << endl;
                    return false;
                }

                double l_endTime = doubleData[k][0];
                int l_steps = (int) doubleData[k][1];
                if (l_steps <= 0)
                {
                    cerr << "ERROR: negative time steps! Please check the inputs!" << endl;
                    return false;
                }

                double l_iniTime;
                if (k == 0)
                {
                    l_iniTime = 0.0;
                }
                else
                {
                    l_iniTime = doubleData[k - 1][0];
                }
                
                if (l_iniTime > l_endTime)
                {
                    cerr << "ERROR: time list and steps input error! Please check the inputs!" << endl;
                    return false;
                }

                double l_deltaTime = (l_endTime - l_iniTime) / l_steps;

                if (k == 0)
                {
                    m_simTime_.push_back(l_iniTime);
                }
                for (int j = 0;j < l_steps;j++)
                {
                    m_simTime_.push_back(l_iniTime + (j + 1) * l_deltaTime);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_MOLALITY_WATER")
        {
            vector< vector< string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error component molality for print inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = stringData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error component molality  for print inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_printWaterMolality_.push_back(stringData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_GASMOLEFRACTION")
        {
            vector< vector< string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error gas mole fraction for print inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = stringData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error gas mole fraction for print inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_printGasMolefraction_.push_back(stringData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_GASMOLENUMBER")
        {
            vector< vector< string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error gas mole number for print inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = stringData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error gas mole number for print inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_printGasMoleNumbers_.push_back(stringData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_MASTER")
        {
            vector< vector< string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error master species for print inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = stringData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error master species for print inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_printMaster_.push_back(stringData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_SECONDARY")
        {
            vector< vector< string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error secondary species for print inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = stringData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error secondary species for print inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_printSecondary_.push_back(stringData[k][j]);
                }
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_MINERAL")
        {
            vector< vector< string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = (int)stringData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error master species for print inputs! Please check!" << endl;
                return false;
            }

            for (int k = 0;k < rowNumber;k++)
            {
                columNumber = stringData[k].size();
                if (columNumber <= 0)
                {
                    cerr << "ERROR: GCM: Error master species for print inputs! Please check!" << endl;
                    return false;
                }
                for (int j = 0;j < columNumber;j++)
                {
                    m_printMineral_.push_back(stringData[k][j]);
                }
            }
        }
        if ((*a_v_GCMData)[i] == "REDOX")
        {
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);
            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error redox inputs! Please check!" << endl;
                return false;
            }
            if (doubleData[0][0] > 0)
            {
                m_redox_ = true;
            }
            else
            {
                m_redox_ = false;
            }
        }

        if ((*a_v_GCMData)[i] == "PRINT_INITIAL")
        {
            vector< vector< double>> doubleData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);
            rowNumber = (int)doubleData.size();
            if (rowNumber <= 0)
            {
                cerr << "ERROR: GCM: Error PRINT_INITIAL inputs! Please check!" << endl;
                return false;
            }
            if (doubleData[0][0] > 0)
            {
                m_printInitialEquil_ = true;
            }
            else
            {
                m_printInitialEquil_ = false;
            }
        }
    }

    if (!databaseExist)
    {
        cout << "WARNING: GCM: PHREEQC database file is not clarified, and the default file name is used!"<<endl;
        m_databaseFile_ = "phreeqc.dat";
    }

    if (!phrqcInputFileExist)
    {
        cout << "WARNING: GCM: PHREEQC input file is not clarified";
        m_phrqc_inputFile_ =  "gcm.inp";
    }

    bool readPhrqcInputFileSuccess = m_readPhrqcInputFile_(m_phrqc_inputFile_);
    if (!readPhrqcInputFileSuccess)
    {
        cerr<<"ERROR: GCM: Fail to read PHREEQC input file!"<<endl;
        return false;
    }

    if (m_simTime_.size() > 1)
    {
        m_isKinetics = true;
        if ((m_pressureList_.size() != 1) || (m_temperatureList_.size() != 1))
        {
            cerr << "ERROR: the number of pressure or temperature inputs is not proper! Please check the inputs!" << endl;
            exit(-1);
            return false;
        }
        m_printType = 0;
    }
    else
    {
        if (((m_pressureList_.size() == 0) || (m_temperatureList_.size() == 0)) && (m_temperature_pressureList_.size()==0) )
        {
            cerr << "ERROR: no proper inputs of pressure or temperature! Please check the inputs!" << endl;
            exit(-1);
            return false;
        }

        if ((m_pressureList_.size() > 1) && (m_temperatureList_.size() > 1))
        {
            cerr << "ERROR: too many inputs for temperature or pressure! Please check the inputs" << endl;
            exit(-1);
            return false;
        }

        if (m_gasMoleFractionList_.size() > 1)
        {
            if (m_temperatureList_.size() > 1 || m_pressureList_.size() > 1 || m_temperature_pressureList_.size() > 1)
            {
                cerr << "ERROR: too many list inputs! Please check the input file!" << endl;
                exit(-1);
                return false;
            }
            else if (m_gasMoleFraction_.size() > 0)
            {
                cerr << "ERROR: too many list inputs! Please check the input file!" << endl;
                exit(-1);
                return false;
            }
        }

        if (m_pressureList_.size() > 1)
        {
            m_printType = 2;
        }
        else if (m_temperatureList_.size() > 1)
        {
            m_printType = 1;
        }
        else
        {
            m_printType = -1;
        }
    }

    if (m_temperature_pressureList_.size() > 0)
    {
        if (m_temperatureList_.size() > 0 || m_pressureList_.size() > 0)
        {
            cerr << "ERROR: too many temperature/pressure inputs!" << endl;
            exit(-1);
            return false;
        }
    }
    return true;
}

bool GeoChemInputs::m_readPhrqcInputFile_(string& a_phrqc_input)
{
    string fold_inputFile = a_phrqc_input;
    //   fold_inputFile = REAL_FILE_PATH(fold_inputFile);

    ifstream phrqcInputs(fold_inputFile.c_str());
    string wholeLine;
    string tempWord;
    string dataBlockString;
    vector<string> lineString;
    //    int bigNumber = 9999;

    int solutionIndex = -1;
    int equilPhaseIndex = -1;
    int kineticIndex = -1;
    int inDataBlock = 0;
    int dataIndex;

    if (phrqcInputs.is_open())
    {
        while (!phrqcInputs.eof())
        {
            getline(phrqcInputs, wholeLine);
            wholeLine = UtilString::trim(wholeLine);
            stringstream stringin(wholeLine);
            lineString.clear();
            while (stringin >> tempWord)
            {
                lineString.push_back(UtilString::trim(tempWord));
            }

            if (lineString.size() == 0)
            {
                continue;
            }

            if (lineString[0] == "SOLUTION")
            {
                if (inDataBlock == 0)
                {
                    dataBlockString = wholeLine + "\n";
                }

                if (lineString.size() < 2)
                {
                    cerr << "ERROR << GCM: a solution index (an integer) should be provided!" << endl;
                    return false;
                }

                stringstream solutionIndexRead(lineString[1]);
                solutionIndexRead >> solutionIndex;
                if (inDataBlock > 0)
                {
                    /*    if ((solutionIndex != equilPhaseIndex)&&(solutionIndex!=kineticIndex))
                        {
                            logError logerror;
                            logerror.LOGERROR( "GCM: solution index and equilibrium phase index are not equal!");
                            return false;
                        }
                        */

                    dataBlockString += wholeLine + "\n";
                }
                else if (inDataBlock == 0)
                {
                    dataIndex = solutionIndex;
                }

                if (inDataBlock == 0)
                {
                    inDataBlock = 1;
                }
            }
            else if (lineString[0] == "EQUILIBRIUM_PHASES")
            {
                if (inDataBlock == 0)
                {
                    dataBlockString = wholeLine;
                }
                if (lineString.size() < 2)
                {
                    cerr << "ERROR: GCM: a solution index (an integer) should be provided!" << endl;
                    return false;
                }

                stringstream equiliPhaseIndexRead(lineString[1]);
                equiliPhaseIndexRead >> equilPhaseIndex;
                if (inDataBlock > 0)
                {
                    /*
                    if ((solutionIndex != equilPhaseIndex)&& (equilPhaseIndex!=kineticIndex))
                     {
                         logError logerror;
                         logerror.LOGERROR("GCM: the equlibrium phase index should be equal to that of the two others!");
                         return false;
                     }
                     */

                    dataBlockString += wholeLine + "\n";
                }
                else if (inDataBlock == 0)
                {
                    dataIndex = equilPhaseIndex;
                }

                if (inDataBlock == 0)
                {
                    inDataBlock = 1;
                }
            }
            else if (lineString[0] == "KINETICS")
            {
                if (inDataBlock == 0)
                {
                    dataBlockString = wholeLine;
                }
                if (lineString.size() < 2)
                {
                    cerr << "ERROR: GCM: a solution index (an integer) should be provided!" << endl;
                    return false;
                }

                stringstream kineticIndexRead(lineString[1]);
                kineticIndexRead >> kineticIndex;
                if (inDataBlock > 0)
                {
                    /*
                    if ((kineticIndex != equilPhaseIndex) && (kineticIndex!=solutionIndex))
                    {
                        logError logerror;
                        logerror.LOGERROR("GCM: the kinetic index should be equal to that of the two others!");
                        return false;
                    }
                    */
                    dataBlockString += wholeLine + "\n";
                }
                else if (inDataBlock == 0)
                {
                    dataIndex = kineticIndex;
                }

                if (inDataBlock == 0)
                {
                    inDataBlock = 1;
                }
            }
            else if (lineString[0] == "END")
            {
                inDataBlock = 0;
                dataBlockString += wholeLine + "\n";
                m_phrqcInputBlock_[dataIndex] = dataBlockString;
            }
            else
            {
                if (inDataBlock == 0)
                {
                    continue;
                }
                dataBlockString += wholeLine + "\n";
                //               stringstream dataBlockStringStream(dataBlockString);
            }
        }
        phrqcInputs.close();
        return true;
    }
    else
    {
        cerr << "ERROR: GCM: CANNOT open file " + fold_inputFile << endl;
        return false;
    }
}

string* GeoChemInputs::m_getDataBaseFile()
{
    return &m_databaseFile_;
}

string* GeoChemInputs::m_getPhrqcInputFile()
{
    return &m_phrqc_inputFile_;
}

list<double>* GeoChemInputs::m_getSimTime()
{
    return &m_simTime_;
}

list<double>* GeoChemInputs::m_getPressure()
{
    return &m_pressureList_;
}

list<double>* GeoChemInputs::m_getTemperature()
{
    return &m_temperatureList_;
}

vector<string>* GeoChemInputs::m_getMasterPrint()
{
    return &m_printMaster_;
}

vector<string>* GeoChemInputs::m_getSecondaryPrint()
{
    return &m_printSecondary_;
}

vector<string>* GeoChemInputs::m_getMineralPrint()
{
    return &m_printMineral_;
}

map<string, double>* GeoChemInputs::m_getGasMoleFraction()
{
    return &m_gasMoleFraction_;
}

double  GeoChemInputs::m_getGasMoleNumber()
{
    return m_gasMoleNumber_;
}

vector<string>* GeoChemInputs::m_waterMolalityPrint()
{
    return &m_printWaterMolality_;
}

vector<string>* GeoChemInputs::m_gasMoleFractionPrint()
{
    return &m_printGasMolefraction_;
}

vector<string>* GeoChemInputs::m_gasMoleNumbers()
{
    return &m_printGasMoleNumbers_;
}

list<pair<double, double>>* GeoChemInputs::m_getTempPressList()
{
    return &m_temperature_pressureList_;
}

bool GeoChemInputs::m_getRedox()
{
    return m_redox_;
}

list < vector<double> >* GeoChemInputs::m_getMoleFractionList()
{
    return &m_gasMoleFractionList_;
}

vector<string>* GeoChemInputs::m_getGasComponentName()
{
    return &m_gasComponentName_;
}

bool GeoChemInputs::m_getPrintInitial()
{
    return m_printInitialEquil_;
}