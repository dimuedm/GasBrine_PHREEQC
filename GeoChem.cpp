#include "Phreeqc.h"
#include "global_structures.h"
#include "GeoChem.h";


using namespace std;

GeoChem::GeoChem()
{
    geoChemExist = 1;
}
GeoChem::~GeoChem()
{
    if (m_phreeqc_)
    {
        delete m_phreeqc_;
    }
}

GeoChem::GeoChem(const GeoChem &obj)
{
    m_phreeqc_ = new Phreeqc();
    *m_phreeqc_ = *(obj.m_phreeqc_);
}

GeoChem&
GeoChem::
 operator = (GeoChem &obj1)
{
    /*if (this->m_phreeqc_)
    {
       delete  this->m_phreeqc_;
    }*/

    this->m_phreeqc_ = new Phreeqc();
    *(this->m_phreeqc_) = *(obj1.m_phreeqc_);
    this->m_phreeqc_->ReSoCMode = 1;
    this->m_phreeqc_->ReSoCMode_DEBUG = 0;
    this->m_phreeqc_->ReSoC_timeStepCalc = 0;
  
    return *this;
}

void GeoChem::m_assign(Phreeqc *a_phreeqc)
{
    //if (!m_phreeqc_)
    //{
        m_phreeqc_ = new Phreeqc();
    //}

    *m_phreeqc_ = *a_phreeqc;

    m_phreeqc_->ReSoCMode = 1;
    m_phreeqc_->ReSoCMode_DEBUG = 0;
    m_phreeqc_->ReSoC_timeStepCalc = 0;
}

Phreeqc* GeoChem::m_getPhreeqc()
{
    return m_phreeqc_;
}

void GeoChem::save_solution()
{
    m_phreeqc_->xsolution_save(-1);
}

//int GeoChem::m_getMasterNumber()
//{
//    int* l_numberOfMasterSpecies = new int();
//    *l_numberOfMasterSpecies = m_phreeqc_->count_
//}

std::map<int, std::string> *
GeoChem::getMasterSpecies()
{
    return m_phreeqc_->get_masterSpecies_aq();
}

std::map<std::string, double>*
GeoChem::getMasterSpeciesMolality()
{
    return m_phreeqc_->get_masterMolality_aq();
}

void GeoChem::print_solution(int n_user)
{
    m_phreeqc_->print_solution(n_user);
}

void GeoChem::initialize_database(istream *a_databaseStream)
{
    //if (!m_phreeqc_)
    //{
        m_phreeqc_ = new Phreeqc();
    //}

    m_phreeqc_->ReSoCMode = 1;
    m_phreeqc_->ReSoCMode_DEBUG = 0;
    m_phreeqc_->ReSoC_timeStepCalc = 0;

    int errors = m_phreeqc_->do_initialize();

    m_phreeqc_->Get_phrq_io()->push_istream(a_databaseStream);
    errors = m_phreeqc_->read_database();
    /*if (a_databaseStream)
    {
        delete a_databaseStream;
    }*/
    m_phreeqc_->Get_phrq_io()->clear_istream();
}

void GeoChem::initialize(string *a_inputString)
{
    istream *input_cookie = NULL;

//    m_phreeqc_->Get_phrq_io()->clear_istream();

    input_cookie = new stringstream(*a_inputString);

    m_phreeqc_->Get_phrq_io()->push_istream(input_cookie);
    m_phreeqc_->run_simulations();

    /* if (db_cookie)
    {
    delete db_cookie;
    }*/

    /*if (input_cookie)
    {
        delete input_cookie;
    }*/
    m_phreeqc_->Get_phrq_io()->clear_istream();
}

//void GeoChem::initialize(istream *a_inputStream, istream *a_databaseStream)
//void GeoChem::initialize(string *a_inputString, string *a_databaseFileName)
void GeoChem::initialize(string *a_inputString, istream *a_databaseStream)
{
//    std::istream *db_cookie = NULL;
    std::istream *input_cookie = NULL;
    int errors;

//    m_phreeqcInputFile_ = "Ini-solution-for-solubility-test.dat";
//    std::string outFileName;
//    std::string databaseFileName = "pitzer.dat";
    m_phreeqc_ = new Phreeqc();
    m_phreeqc_->ReSoCMode = 1;
    m_phreeqc_->ReSoCMode_DEBUG = 0;
    m_phreeqc_->ReSoC_timeStepCalc = 0;

    /*string a_inputFileName = "test";
    errors = m_phreeqc_->process_file_names(&db_cookie, &input_cookie, &a_inputFileName, &outFileName, a_databaseFileName, TRUE);*/

    /*ifstream fin(*a_databaseFileName);
    if (!fin)
    {
        cout << "Fail to open database file!" << endl;
    }*/


//    db_cookie = new ifstream(*a_databaseFileName, std::ios_base::in);
    input_cookie = new stringstream(*a_inputString);


    {
        std::string outFileName;
        std::ofstream * local_output_stream;
        int outFileNameSize = outFileName.size();
        if (outFileNameSize <= 0)

        {
            outFileName =  "ReSoC_Phreeqc.out";
        }

        local_output_stream = m_phreeqc_->open_output_stream(&outFileName, std::ios_base::out, false);
        m_phreeqc_->Get_phrq_io()->Set_output_ostream(local_output_stream);
    }

    errors = m_phreeqc_->do_initialize();

    m_phreeqc_->Get_phrq_io()->push_istream(a_databaseStream);
    errors = m_phreeqc_->read_database();
    m_phreeqc_->Get_phrq_io()->clear_istream();

    m_phreeqc_->Get_phrq_io()->push_istream(input_cookie);
    m_phreeqc_->run_simulations();

   /* if (db_cookie)
    {
        delete db_cookie;
    }*/

    if (input_cookie)
    {
        delete input_cookie;
    }
}

void GeoChem::update_timeStep(conditionChange &a_newCondition)
{
    m_phreeqc_->run_simulation_timeStep(a_newCondition);
    m_phreeqc_->results_for_ReSoC(m_gas_, m_water_, m_mineral_);
}

waterSolution_results * GeoChem::getWater()
{
    return &m_water_;
}

std::map<std::string, minerals> *GeoChem::getMineral()
{
    return &m_mineral_;
}


