#include "GeoChemEnvironment.h"
#include "logError.h"

#include "Phreeqc_ReSoC.h"
#include "logError.h"
using namespace std;

GeoChemEnvironment::GeoChemEnvironment()
{
}

GeoChemEnvironment::~GeoChemEnvironment()
{
    /*for (map<int, istream*>::iterator it = m_phrqcIstram_.begin(); it != m_phrqcIstram_.end(); it++)
    {
        delete it->second;
    }*/
    if (m_databaseStream_->goodbit)
    {
        delete m_databaseStream_;
    }
}

bool GeoChemEnvironment::
m_readInputs_(map <string, vector<string> > * a_pinputDataBlock_)
{

    int inputDataBlockSize = a_pinputDataBlock_->size();
    map <string, vector<string> >::iterator it;
    vector<string> *v_GCMData;

    if (inputDataBlockSize <= 0)
    {
        logError logerror;
        logerror.LOGERROR("Wrong size for input data!" );
        return false;
    }
    else
    {
        it = a_pinputDataBlock_->find("GCM");

        if (it == a_pinputDataBlock_->end())
        {
            LOG(INFO) << "There is no GCM data block!" ;
            return false;
        }
        else
        {
            v_GCMData = &(it->second);
            if (m_readGCMDataBlock_(v_GCMData))
            {
                return true;
            }
            else
            {
                logError logerror;
                logerror.LOGERROR("Read GCM data error! Please check!" );
                return false;
            }
        }
    }
}

bool GeoChemEnvironment::
m_readGCMDataBlock_(vector<string> *a_v_GCMData)
{
    int iGCMDataSize = a_v_GCMData->size();
    vector< vector<double> > doubleData;
    vector< vector<int> > intData;
    vector<zone> v_tempZone;
    bool readInSuccess;
    int rowNumber;
    int columNumber;
    bool zoneindexExist=0, zoneExist=0, databaseExist=0;
    bool phrqcInputFileExist = 0;
    if (iGCMDataSize <= 0)
    {
        logError logerror;
        logerror.LOGERROR("GCM: NO GCM data!");
        return false;
    }

    for (int i = 0; i < iGCMDataSize; i++)
    {
        if ((*a_v_GCMData)[i] == "ZONEINDEX")
        {
            intData.clear();
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, intData);
            rowNumber = intData.size();

            if (!readInSuccess || rowNumber <= 0)
            {
                logError logerror;
                logerror.LOGERROR("GCM: Fail to read ZONEINDEX!");
                return false;
            }            

            for (int j = 0; j < rowNumber; j++)
            {
                columNumber = intData[j].size();
                if (columNumber <= 0)
                {
                    logError logerror;
                    logerror.LOGERROR("GCM: Less data than needed!");
                    return false;
                }

                for (int k = 0; k < columNumber; k++)
                {
                    m_zoneIndex_.push_back(intData[j][k]);
                }
            }

            m_defaultZoneIndex_ = m_zoneIndex_[m_zoneIndex_.size() - 1];
            zoneindexExist = 1;
        }

        if ((*a_v_GCMData)[i] == "ZONE")
        {
            doubleData.clear();
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, doubleData);

            rowNumber = doubleData.size();
            if (rowNumber <= 0 || !readInSuccess)
            {
                logError logerror;
                logerror.LOGERROR("GCM: Fail to read ZONE!");
                return false;
            }
            m_zoneNumber_ = rowNumber;
            for (int j = 0; j < rowNumber; j++)
            {
                columNumber = doubleData[j].size();
                if (columNumber != 6)
                {
                    logError logerror;
                    logerror.LOGERROR("GCM: NO enough zone coordinate data!");
                    return false;
                }

                zone tempZone;
                tempZone.xmin = doubleData[j][0];
                tempZone.xmax = doubleData[j][1];
                tempZone.ymin = doubleData[j][2];
                tempZone.ymax = doubleData[j][3];
                tempZone.zmin = doubleData[j][4];
                tempZone.zmax = doubleData[j][5];

                v_tempZone.push_back(tempZone);
            }
            zoneExist = 1;
        }

        if ((*a_v_GCMData)[i] == "DATABASE")
        {
            vector<vector<string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = stringData.size();
            if (rowNumber != 1)
            {
                logError logerror;
                logerror.LOGERROR("GCM: Error database file clarification!");
                return false;
            }
            m_databaseFile_ = "Input/" + UtilString::trim(stringData[0][0]);

            databaseExist = 1;
        }

        if ((*a_v_GCMData)[i] == "INJSOLUTION")
        {
            intData.clear();
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, intData);
            rowNumber = intData.size();
            if (rowNumber <= 0)
            {
                logError logerror;
                logerror.LOGERROR("GCM: Error injection solution reading!");
                return false;
            }
            int l_numberOfinjectionSolution = 0;
            for (int j = 0; j < rowNumber; j++)
            {
                int columNumber = intData[j].size();
                for (int k = 0; k < columNumber; k++)
                {
                    m_injectedSolution_.push_back(intData[j][k]);
                }              
            }
        }

        if ((*a_v_GCMData)[i] == "PRQCINPUT")
        {
            vector<vector<string>> stringData;
            readInSuccess = UtilIO::hardDataReadIn(i, a_v_GCMData, stringData);

            rowNumber = stringData.size();
            if (rowNumber != 1)
            {
                logError logerror;
                logerror.LOGERROR("GCM: Error geochemistry input file clarification!");
                return false;
            }

            m_phrqc_inputFile_ = UtilString::trim(stringData[0][0]);
            phrqcInputFileExist = 1;
        }
    }

    if (!databaseExist)
    {
        LOG(INFO) << "GCM: PHREEQC database file is not clarified, and the default file name is used!";
        m_databaseFile_ = "Input/phreeqc.dat";
    }

    if (!phrqcInputFileExist)
    {
        LOG(INFO) << "GCM: PHREEQC input file is not clarified";
        m_phrqc_inputFile_ = "Input/gcm.inp";
    }

    bool readPhrqcInputFileSuccess = m_readPhrqcInputFile_(m_phrqc_inputFile_);
    if (!readPhrqcInputFileSuccess)
    {
        logError logerror;
        logerror.LOGERROR("GCM: Fail to read PHREEQC input file!");
        return false;
    }

    if (!zoneindexExist)
    {
        logError logerror;
        logerror.LOGERROR("GCM: zone index is not provided!");
        return false;
    }

    if (zoneExist)
    {
        for (int i = 0; i < m_zoneNumber_; i++)
        {
            if ( i >= (m_zoneIndex_.size()-1) )
            {
                break;
            }
            m_mp_Zones_[m_zoneIndex_[i]] = v_tempZone[i];
        }
    }
    
    map<int, zone>::iterator zoneIter = m_mp_Zones_.begin();
    for (; zoneIter != m_mp_Zones_.end(); zoneIter++)
    {
        map<int, string>::iterator it = m_phrqcInputBlock_.find(zoneIter->first);
        if (it == m_phrqcInputBlock_.end())
        {
            logError logerror;
            logerror.LOGERROR("GCM: NO related PHREEQC inputs for the zone " + to_string(zoneIter->first) + "!");
            return false;
        }        
    }
}

vector<int> * GeoChemEnvironment::
m_getInjectedSolution()
{
    return &m_injectedSolution_;
}

bool GeoChemEnvironment::
m_readPhrqcInputFile_(string &a_phrqc_input)
{
    string fold_inputFile = "Input/" + a_phrqc_input;
    fold_inputFile = REAL_FILE_PATH(fold_inputFile);

    ifstream phrqcInputs( fold_inputFile.c_str() );
    string wholeLine;
    string tempWord;
    string dataBlockString;
    vector<string> lineString;
//    int bigNumber = 9999;
    
    int solutionIndex;
    int equilPhaseIndex;
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
                    dataBlockString = wholeLine+"\n";
                }
                
                if (lineString.size() < 2)
                {
                    logError logerror;
                    logerror.LOGERROR( "GCM: a solution index (an integer) should be provided!");
                    return false;
                }

                stringstream solutionIndexRead(lineString[1]);
                solutionIndexRead >> solutionIndex;
                if (inDataBlock > 0)
                {
                    if (solutionIndex != equilPhaseIndex)
                    {
                        logError logerror;
                        logerror.LOGERROR( "GCM: solution index and equilibrium phase index are not equal!");
                        return false;
                    }

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
                    logError logerror;
                    logerror.LOGERROR("GCM: a solution index (an integer) should be provided!");
                    return false;
                }

                stringstream equiliPhaseIndexRead(lineString[1]);
                equiliPhaseIndexRead >> equilPhaseIndex;
                if (inDataBlock > 0)
                {
                    if (solutionIndex != equilPhaseIndex)
                    {
                        logError logerror;
                        logerror.LOGERROR("GCM: the solution index and equilibrium phase index are not equal!");
                        return false;
                    }

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
        
    }
    else
    {
        logError logerror;
        logerror.LOGERROR("GCM: CANNOT open file " + fold_inputFile);
        return false;
    }

//     for (map<int, string>::iterator it = m_phrqcInputBlock_.begin(); it != m_phrqcInputBlock_.end(); it++)
//    {
//        string tempFileName = "phrqcInput_Zone"+ to_string(it->first);
//        ofstream tempOF(tempFileName);
//
//        tempOF << it->second;
//
////        m_phrqcIstram_[it->first] = new ifstream(tempFileName, ios_base::in);
//    }
     m_databaseFile_ = UtilFile::getFilePath(m_databaseFile_);
     m_databaseStream_ = new ifstream(m_databaseFile_, ios_base::in);


    return true;
}

bool GeoChemEnvironment::
m_GeochemEvironmentSetup(map <string, vector<string> > * a_pinputDataBlock_)
{
    return m_readInputs_(a_pinputDataBlock_);
}

//map<int, istream *>* GeoChemEnvironment::
//m_getPhrqcIstream()
//{
//    return &m_phrqcIstram_;
//}
//
istream ** GeoChemEnvironment::
m_getDatabaseStream()
{
    return &m_databaseStream_;
}

map<int, string> *GeoChemEnvironment::
m_getPrqcInputData()
{
    return &m_phrqcInputBlock_;
}

string *GeoChemEnvironment::
m_getDatabaseFileName()
{
    return &m_databaseFile_;
}

map<int, zone> *GeoChemEnvironment::
m_getZones()
{
    return &m_mp_Zones_;
}

int GeoChemEnvironment::
m_getDefaultZoneIndex()
{
    return m_defaultZoneIndex_;
}