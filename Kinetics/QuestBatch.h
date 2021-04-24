#pragma once
#include "../MultiPhaseEquilibria/GasIniWaterEquilibria.h"

class QuestBatch
{
public:
    QuestBatch();
    ~QuestBatch();

public:
    GasIniWaterEquilibria gasWater;

    void initialize(conditionChange &condition_in, ofstream &a_outputs);
    conditionChange condition;

    bool questBatchModel(string inputFile, map<string, double> &gasMoleNumber);
    void resultPrint(ofstream &a_out, double currentTime);

private:

    void reaction(double time, GasIniWaterEquilibria &gasWater, map<string, double> &gasMoleNumber);
    void updateSolubility(map<string, double> &gasMoleNumber);
    void updateSolubility();
    ofstream *outPuts;

};

