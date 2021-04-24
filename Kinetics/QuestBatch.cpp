
#include "QuestBatch.h"

QuestBatch::QuestBatch()
{
}

QuestBatch::~QuestBatch()
{
}

void QuestBatch::initialize(conditionChange &conditionIn, ofstream &a_outputs)
{
    condition = conditionIn;
    condition.gasPhaseIn = false;
    outPuts = &a_outputs;
}

bool QuestBatch::questBatchModel(string inputFile, map<string, double> &gasMoleNumber)
{
    //initialize the solution with the mineral composition
   gasWater.initializePhreeqc(inputFile);
   condition.masterChange = gasWater.finalWater.masterSpecies_molality;
   condition.gasPhaseIn = false;
   condition.mass_of_H2O = 1.0;
   map<string, minerals>::iterator it = gasWater.finalMineral.begin();
   for (; it != gasWater.finalMineral.end(); it++)
   {
       condition.mineralsUpdate[it->first] = it->second.final_mole;
   }

  // updateSolubility(gasMoleNumber);
   double totalTime = 90.0* 3153600.0;
   double deltaTime = 3153600.0;
   double currentTime = 0;
   
   reaction(0, gasWater, gasMoleNumber);
   resultPrint(*outPuts, currentTime);
   
   int i_step = 0;
   while (currentTime < totalTime)
   {
       i_step++;
       reaction(deltaTime, gasWater, gasMoleNumber);
       currentTime += deltaTime;
       resultPrint(*outPuts, currentTime);
   }
   currentTime = totalTime;
   totalTime += 10.0 * 31536000.0;
   deltaTime = 31536000.0;
 //  currentTime = totalTime - deltaTime;
   while (currentTime < totalTime)
   {
       reaction(deltaTime, gasWater, gasMoleNumber);
       currentTime += deltaTime;
       resultPrint(*outPuts, currentTime);
   }

   currentTime = totalTime;
   totalTime += 100.0 * 315360000.0;
   deltaTime = 315360000.0;
   while (currentTime < totalTime)
   {
       reaction(deltaTime, gasWater, gasMoleNumber);
       currentTime += deltaTime;
       resultPrint(*outPuts, currentTime);
   }

   outPuts->close();

   return OK;
}

void QuestBatch::reaction(double deltaTime, GasIniWaterEquilibria &gasWater, map<string, double> &gasMoleNumber)
{
    // Calculate CO2 solubility in water (current water, and use phreeqc to calcualte CO2 activity in water)
   // gasWater.calPhaseEquilibria(gasMoleNumber, )
    updateSolubility(gasMoleNumber);

    condition.deltaTime = deltaTime;
   // condition.masterChange = gasWater.finalWater.masterSpecies_molality;
    gasWater.updateReaction(condition);

    condition.masterChange = gasWater.finalWater.masterSpecies_molality;
    //condition.mineralsUpdate = gasWater.finalMineral
    map<string, minerals>::iterator it = gasWater.finalMineral.begin();
    for (; it != gasWater.finalMineral.end(); it++)
    {
        condition.mineralsUpdate[it->first] = it->second.final_mole;
    }
    
}


void QuestBatch::updateSolubility()
{
    LDBLE atm_1bar = 0.9869233;
    double pressureBar = 212.0; //condition.pressure/atm_1bar;
    double temperaureK = condition.temperature+273.15;

    condition.masterChange = gasWater.finalWater.masterSpecies_molality;
}

void QuestBatch::updateSolubility(map<string, double> &gasMoleNumber)
{
    LDBLE atm_1bar = 0.9869233;
    double pressureBar = 212.0; //condition.pressure/atm_1bar;
    double temperaureK = condition.temperature + 273.15;

    condition.masterChange = gasWater.finalWater.masterSpecies_molality;
    gasWater.calPhaseEquilibria(gasMoleNumber, condition.masterChange, temperaureK, pressureBar);
    double CO2solubility = gasWater.xAqueous[0] / gasWater.xAqueous[1] * 55.508;

    map<string, double> ::iterator it = condition.masterChange.find("C(4)");
    if (it != condition.masterChange.end())
    {
        it->second = CO2solubility;
    }
    else
    {
        it = condition.masterChange.find("C(+4)");
        if (it != condition.masterChange.end())
        {
            it->second = CO2solubility;
        }
        else
        {
            it = condition.masterChange.find("C");
            if (it != condition.masterChange.end())
            {
                it->second = CO2solubility;
            }
            else
            {
                condition.masterChange["C"] = CO2solubility;
            }
        }

    }
}

void QuestBatch::resultPrint(ofstream &outPuts, double currentTime)
{
    if (currentTime < 1.e-5)
    {
        outPuts << "TIME " << "PH ";
        map<string, double>::iterator it = gasWater.finalWater.masterSpecies_molality.begin();
        for (; it != gasWater.finalWater.masterSpecies_molality.end(); it++)
        {
            outPuts << it->first << " ";
        }

        map<string, minerals>::iterator it_mineral = gasWater.finalMineral.begin();
        for (; it_mineral != gasWater.finalMineral.end(); it_mineral++)
        {
            outPuts << it_mineral->first << " ";
        }
        map<string, kineticPhases>::iterator it_kmineral = gasWater.finalKineticMineral.begin();
        for (; it_kmineral != gasWater.finalKineticMineral.end(); it_kmineral++)
        {
            outPuts << it_kmineral->first << " ";
        }
        outPuts << endl;
   //     return;
    }
   // gasWater.finalWater.
    double pH = gasWater.finalWater.ph;
    outPuts << currentTime << " " << pH << " ";

 //   int numberOfMasterSpecies = gasWater.finalWater.masterSpecies_molality.size();
    map<string, double>::iterator it = gasWater.finalWater.masterSpecies_molality.begin();
    for (; it!= gasWater.finalWater.masterSpecies_molality.end(); it++)
    {
        outPuts << it->second << " ";
    }

    map<string, minerals>::iterator it_mineral = gasWater.finalMineral.begin();
    for (; it_mineral != gasWater.finalMineral.end(); it_mineral++)
    {
        outPuts << it_mineral->second.final_mole << " ";
    }

    map<string, kineticPhases> ::iterator it_kmineral = gasWater.finalKineticMineral.begin();
    for (; it_kmineral != gasWater.finalKineticMineral.end(); it_kmineral++)
    {
        outPuts << it_kmineral->second.final_mole << " ";
    }
    outPuts << endl;
}