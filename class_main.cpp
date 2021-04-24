#include "Phreeqc.h"

#include "NameDouble.h"
#include "Solution.h"
#include "Reaction.h"
#include "PPassemblage.h"
#include "Exchange.h"
#include "Surface.h"
#include "GasPhase.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
//#include "logError.h"
//#include <sys/signal.h>
//#include <fenv.h>
/* ----------------------------------------------------------------------
 *   MAIN
 * ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- */
std::ifstream * Phreeqc::
open_input_stream(std:: string *default_name, std::ios_base::openmode mode, bool batch)
/* ---------------------------------------------------------------------- */
{
//	char name[MAX_LENGTH];
    std::string name;
	std::ifstream *new_stream;
//	int l;
#ifdef ERROR_OSTREAM
	std::ostream * error_ostream_save = phrq_io->Get_error_ostream();
#else
	FILE * error_file_save = phrq_io->Get_error_file();
#endif

	for (;;)
	{
/*
 *   Get file name
 */
//		strcpy(name, default_name);
        name = *default_name;
		if (!batch )
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
/*			screen_msg(sformatf("%s\n", query));
			if ((*default_name)[0] != '\0')
			{
				screen_msg(sformatf("Default: %s\n", default_name));
			}*/                               // Li Jun comments the paragraph;
//			char *s_ptr = fgets(name, MAX_LENGTH, stdin);
            int nameSize = name.size();
            if (nameSize <= 0)
            {
                //logError log_error;
                //log_error.LOGERROR("PHREEQC: File name error!");

                return NULL;
            }
		}
/*
 *   Open existing file to read
 */
		new_stream = new std::ifstream(name, mode);
		if (new_stream == NULL || !new_stream->is_open())
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
//			error_string = sformatf( "\nERROR: Cannot open file, %s.\n", name); //Li Jun comments this line. 2017-8-6.
            /*logError log_error;
            log_error.LOGERROR("PHREEQC: Cannot open file!");*/
 //           std::cerr << "ERROR: Cannot open file " << name << std::endl;
 //           exit(-1);

//			screen_msg(error_string);  //Li Jun comments this line. 2017-8-6.
#ifdef NPP
			error_msg(sformatf( "\nERROR: Cannot open file, %s.\n       Please check, and give the correct, full path + name.\n", name), STOP);
			break;
#endif
			error_flush();
			batch = FALSE;
			continue;		
		}
		break;
	}
//	strncpy(default_name, name, MAX_LENGTH);
	if (!batch )
	{
		//phrq_io->Set_error_ostream(error_file_save);
#ifdef ERROR_OSTREAM
		phrq_io->Set_error_ostream(error_ostream_save);
#else
		phrq_io->Set_error_file(error_file_save);
#endif
	}
	return (new_stream);
}
/* ---------------------------------------------------------------------- */
std::ofstream * Phreeqc::
open_output_stream(std::string *default_name, std::ios_base::openmode mode, bool batch)
/* ---------------------------------------------------------------------- */
{
//	char name[MAX_LENGTH];   // Li Jun comments the line, 2017-8-3
    std::string name;        // Define name as string instead of char*, Li Jun: 2017-8-3  
	std::ofstream *new_stream;
//	int l;
#ifdef ERROR_OSTREAM
	std::ostream * error_ostream_save = phrq_io->Get_error_ostream();
#else
	FILE * error_file_save = phrq_io->Get_error_file();
#endif
	
	for (;;)
	{
/*
 *   Get file name
 */
//		strcpy(name, default_name);
        name = *default_name;
		if (!batch )
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			
/*			screen_msg(sformatf("%s\n", query));
			if (default_name[0] != '\0')
			{
				screen_msg(sformatf("Default: %s\n", default_name));
			}
			char *s_ptr = fgets(name, MAX_LENGTH, stdin);
			if (s_ptr == NULL)
			{
			    std::cerr << "Failed defining name." << std::endl;
			}
			l = (int) strlen(name);
			name[l - 1] = '\0';
			if (name[0] == '\0')
			{
				strcpy(name, default_name);
			}*/                                     //Li Jun comments the paragragh, 2017-8-3
		}
/*
 *   Open existing file to read
 */ 
        
        int nameSize = name.size(); // Li Jun added, 2017-8-3
        if (nameSize <= 0)
        {
            std::cerr << "File name size error!" << std::endl;
            exit(-1);
        }       //Li Jun added this paragragh, 2017-8-3

		new_stream = new std::ofstream(name, mode);
     

/*		if (new_stream == NULL || !new_stream->is_open())
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			error_string = sformatf( "\nERROR: Cannot open file, %s.\n", name);
			screen_msg(error_string);
			error_flush();
			batch = FALSE;
			continue;		
		}*/
		break;
	}
//	strncpy(default_name, name, MAX_LENGTH);    //Li Jun comments the line. 2017-8-3
	if (!batch )
	{
#ifdef ERROR_OSTREAM
		phrq_io->Set_error_ostream(error_ostream_save);
#else
		phrq_io->Set_error_file(error_file_save);
#endif
	}
	return (new_stream);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
std::ofstream * Phreeqc::
open_output_file(char *query, char *default_name, std::ios_base::openmode mode, bool batch)
/* ---------------------------------------------------------------------- */
{
	char name[MAX_LENGTH];
	std::ofstream *new_stream;
	int l;
#ifdef ERROR_OSTREAM
		std::ostream * error_ostream_save = phrq_io->Get_error_ostream();
#else
		FILE * error_file_save = phrq_io->Get_error_file();
#endif
	

	for (;;)
	{
/*
 *   Get file name
 */
		strcpy(name, default_name);
		if (!batch )
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			screen_msg(sformatf("%s\n", query));
			if (default_name[0] != '\0')
			{
				screen_msg(sformatf("Default: %s\n", default_name));
			}
			char *s_ptr = fgets(name, MAX_LENGTH, stdin);
			if (s_ptr == NULL)
			{
			    std::cerr << "Failed defining name." << std::endl;
			}
			l = (int) strlen(name);
			name[l - 1] = '\0';
			if (name[0] == '\0')
			{
				strcpy(name, default_name);
			}
		}
/*
 *   Open existing file to read
 */
		new_stream = new std::ofstream(name, mode);
		if (new_stream == NULL || !new_stream->is_open())
		{
#ifdef ERROR_OSTREAM
			phrq_io->Set_error_ostream(&std::cerr);
#else
			phrq_io->Set_error_file(stderr);
#endif
			error_string = sformatf( "\nERROR: Cannot open file, %s.\n", name);
			screen_msg(error_string);
			error_flush();
			batch = FALSE;
			continue;		
		}
		break;
	}
	strncpy(default_name, name, MAX_LENGTH);
	if (!batch )
	{
#ifdef ERROR_OSTREAM
		phrq_io->Set_error_ostream(error_ostream_save);
#else
		phrq_io->Set_error_file(error_file_save);
#endif
	}
	return (new_stream);
}
#endif
