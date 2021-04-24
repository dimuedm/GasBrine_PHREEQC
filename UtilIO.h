/*

Tool for input data process
Created by Jun Li, Dimue company, Wuhan, China
Created: 2018-1-3
Modified: 2018-1-3

*/

#ifndef INPUT_UTILIO_H_
#define INPUT_UTILIO_H_

#include <vector>
#include <typeinfo>
#include "Framework/Util/UtilString.h"

using namespace std;
class UtilIO
{
public:

    template<typename T>
    static bool hardDataReadIn(int &iVectorPosition, vector<string> *v_PVTData, vector< vector<T> > &data)
    {
        bool keywordBlockEnd = false;
        stringstream split;
        T word;
        int dataSize = (int)((*v_PVTData).size());
        vector<T> lineString;
        while (iVectorPosition < dataSize)
        {
            ++iVectorPosition;
            if ((*v_PVTData)[iVectorPosition].back() == '/')
            {
                (*v_PVTData)[iVectorPosition].erase((*v_PVTData)[iVectorPosition].end() - 1);
                keywordBlockEnd = true;
            }
            split.clear();
            split.str((*v_PVTData)[iVectorPosition]);

            while (split >> word)
            {                
                lineString.push_back(word);
            }

            if (!lineString.empty())
            {
                data.push_back(lineString);
            }
            lineString.clear();
            if (keywordBlockEnd)
            {
                break;
            }
        }
        if (keywordBlockEnd)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

	//Function added by Raheel Ahmed (Dimue Tech): 11-2019
	template<typename T>
	static bool getLineDataOfAFile(ifstream &a_inFile, vector<T> &a_lineData)
	{
		bool l_readSuccess = false;
		a_lineData.clear();
		T tempWord;
		string wholeLine;
		getline(a_inFile, wholeLine);
		stringstream stringIn(wholeLine);
		while (stringIn >> tempWord)
		{
			a_lineData.push_back(tempWord);
		}

		if (a_lineData.size() > 0)
		{
			l_readSuccess = true;
		}

		return l_readSuccess;
	}


	//Function added by Raheel Ahmed (Dimue Tech): 11-2019
	static vector<vector<string> > getModelsData(vector<string> * a_pdataBlock, string a_key1, string a_key2)
	{
		vector<vector<string> > l_modelsData;
		vector<string > l_dataString;
		int l_keyFoundCounter = 0;
		for (int i = 0; i < (int)a_pdataBlock->size(); i++)
		{
			if (a_pdataBlock->at(i) == a_key1)
			{
				if (l_keyFoundCounter > 0)
				{
					l_modelsData.push_back(l_dataString);
					l_dataString.clear();
				}
				//i++;
				l_keyFoundCounter++;
				//continue;
			}
			else if (a_pdataBlock->at(i) == a_key2)
			{
				l_modelsData.push_back(l_dataString);
				l_dataString.clear();
			}
			l_dataString.push_back(a_pdataBlock->at(i));
		}
		l_modelsData.push_back(l_dataString);
		return l_modelsData;
	}

	//Function added by Raheel Ahmed (Dimue Tech): 11-2019
	static vector<vector<string> > getModelsData(vector<string>* a_pdataBlock, string a_key1, string a_key2, string a_key3)
	{
		vector<vector<string> > l_modelsData;
		vector<string > l_dataString;
		int l_keyFoundCounter = 0;
		for (int i = 0; i < (int)a_pdataBlock->size(); i++)
		{
			if (a_pdataBlock->at(i) == a_key1)
			{
				if (l_keyFoundCounter > 0)
				{
					l_modelsData.push_back(l_dataString);
					l_dataString.clear();
				}
				//i++;
				l_keyFoundCounter++;
				//continue;
			}
			else if (a_pdataBlock->at(i) == a_key2 || a_pdataBlock->at(i) == a_key3)
			{
				l_modelsData.push_back(l_dataString);
				l_dataString.clear();
			}
			l_dataString.push_back(a_pdataBlock->at(i));
		}
		l_modelsData.push_back(l_dataString);
		return l_modelsData;
	}
private:

};

//UtilIO::UtilIO()
//{
//}
//
//UtilIO::~UtilIO()
//{
//}


//template<typename T>
//bool UtilIO::hardDataReadIn(int &iVectorPosition, vector<string> *v_PVTData, vector< vector<T> > &data)
//{
//    bool keywordBlockEnd = false;
//    stringstream split;
//    T word;
//    int dataSize = (*v_PVTData).size();
//    vector<T> lineString;
//    while (iVectorPosition < dataSize)
//    {
//        ++iVectorPosition;
//        if ((*v_PVTData)[iVectorPosition].back() == '/')
//        {
//            (*v_PVTData)[iVectorPosition].erase((*v_PVTData)[iVectorPosition].end() - 1);
//            keywordBlockEnd = true;
//        }
//        split.clear();
//        split.str((*v_PVTData)[iVectorPosition]);
//
//        while (split >> word)
//        {
//            lineString.push_back(word);
//        }
//
//        if (!lineString.empty())
//        {
//            data.push_back(lineString);
//        }
//        lineString.clear();
//        if (keywordBlockEnd)
//        {
//            break;
//        }
//    }
//    if (keywordBlockEnd)
//    {
//        return true;
//    }
//    else
//    {
//        return false;
//    }
//}


#endif
