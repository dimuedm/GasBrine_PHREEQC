/*

Try to reproduce calcite solubility in CO2-H2O fluids
Created by Jun Li, Dimue company
Created: 201-5-20

*/

#ifndef ELLISREPRODUCE_H_
#define ELLISREPRODUCE_H_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
class EllisReproduce
{
public:
    EllisReproduce();
    ~EllisReproduce();

    void initialize();
    void calcCalciteSolubility(double temperature, double pressure, double mNaCl, double &calciteMolality);
private:
    GasWaterEquilibria gasWaterEquili;
    Phreeqc phreeqc;
};


#endif