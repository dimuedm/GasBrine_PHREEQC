/*

Try to reproduce magnesite solubility in CO2-H2O fluids
Reference: Pascale B¨¦n¨¦zeth, Giuseppe D. Saldi, Jean-Louis Dandurand, Jacques Schott, Chemical Geology 286 (2011) 21¨C31
Created by Jun Li, Dimue company
Created: 201-5-21

*/
#ifndef MAGNESITESOLUBILITY_H_
#define MAGNESITESOLUBILITY_H_

#include "../MultiPhaseEquilibria/GasWaterEquilibria.h"
class MagnesiteSolubility
{
public:
    MagnesiteSolubility();
    ~MagnesiteSolubility();

    void initialize();
    void calcMagnesiteSolubility(double temperature, double pressure, double mNaCl, double CO2moleNumber, double &magMolality, double &Mg2Molality, double &pH);
    void calcMagnesiteSolubilityNoCO2(double temperature, double pressure, double mNaCl, double &magMolality);
private:
    GasWaterEquilibria gasWaterEquili;
    Phreeqc phreeqc;
};


#endif