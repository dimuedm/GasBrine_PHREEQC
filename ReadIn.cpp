
#include <iostream>
#include "Phreeqc.h"
#include "global_structures.h"

// Li Jun changed: first two arguments were changed.
int Phreeqc::
process_file_names(std::istream **db_cookie,
std::istream **input_cookie, std::string *inputFileName, std::string *outFileName, std::string *databaseFileName, int log)
/* ---------------------------------------------------------------------- */
{
    int l;
    char token[2 * MAX_LENGTH];// , default_name[2 * MAX_LENGTH];
//    char query[2 * MAX_LENGTH];
//    char in_file[2 * MAX_LENGTH], out_file[2 * MAX_LENGTH], db_file[2 * MAX_LENGTH];    //Li Jun coomments the lines. 2017-8-3
    char *env_ptr;
    char *ptr;
    /*
    *   Prepare error handling
    */
    try {
        
        /*
        *   Prep for get_line
        */
        max_line = MAX_LINE;
        space((void **)((void *)&line), INIT, &max_line, sizeof(char));
        space((void **)((void *)&line_save), INIT, &max_line, sizeof(char));
        /*
        *   Open error ostream
        */
       
            phrq_io->error_open(NULL);

        /*
        *   Open user-input file
        */
//        strcpy_s(query, "Name of input file?"); //Li Jun comments this line - 2017-8-3
        std::ifstream * local_input_stream = NULL;
//        if (argc <= 1)      //Li Jun comments this line
//        {                   //Li Jun comments this line
 //           default_name[0] = '\0';
            local_input_stream = open_input_stream(inputFileName, std::ios_base::in, false);
//        }   //Li Jun comments this line
//        else   //Li Jun comments this line
//        {      //Li Jun comments this line
//            strcpy(default_name, argv[1]);     //Li Jun comments this line
//            local_input_stream = open_input_stream(query, default_name, std::ios_base::in, true);  //Li Jun comments this line
//        }   //Li Jun comments this line
//        screen_msg(sformatf("Input file: %s\n\n", default_name));  //Li Jun comments this line 2017-8-3
//        strcpy(in_file, default_name);                             //Li Jun comments this line 2017-8-3
        /*
        *   Open file for output
        */
  //      strcpy(query, "Name of output file?");  //Li Jun comments this line 2017-8-3
 //       ptr = default_name;                     //Li Jun comments this line 2017-8-3
 //       copy_token(token, &ptr, &l);            //Li Jun comments this line 2017-8-3
 //       strcat(token, ".out");                  //Li Jun comments this line 2017-8-3
            if (ReSoCMode_DEBUG)
            {
                std::ofstream * local_output_stream;
                //        if (argc <= 1)             //Li Jun comments this line
                //        {                          //Li Jun comments this line
                int outFileNameSize = (int) outFileName->size();
                if (outFileNameSize <= 0)

                {
                    *outFileName = *inputFileName + ".out";
                }

                local_output_stream = open_output_stream(outFileName, std::ios_base::out, false);
                //        }                            //Li Jun comments this line
                //        else if (argc == 2)         //Li Jun comments this line
                //        {                           //Li Jun comments this line
                //            local_output_stream = open_output_stream(query, token, std::ios_base::out, true); //Li Jun comments this line
                //        }  //Li Jun comments this line   
                //        else if (argc >= 3)  //Li Jun comments this line
                //        {  //Li Jun comments this line
                //            strcpy(token, argv[2]);  //Li Jun comments this line
                //            local_output_stream = open_output_stream(query, token, std::ios_base::out, true);  //Li Jun comments this line
                //        }  //Li Jun comments this line
                //        screen_msg(sformatf("Output file: %s\n\n", token)); //Li Jun comments this line. 2017-8-3.
                //        strcpy(out_file, token); //Li Jun comments this line. 2017-8-3.
                phrq_io->Set_output_ostream(local_output_stream);
            }   /// Li Jun modified. 2017-12-10.
        /*
        *   Open log file
        */
            if (ReSoCMode_DEBUG)
            {
                if (log == TRUE)
                {
                    if (!phrq_io->log_open("phreeqc.log"))
                    {
                        error_msg("Cannot open log file, phreeqc.log.", STOP);
                    }
                }
            } // Li Jun modified at 2017-12-10,
        /*
        *  Read input file for DATABASE keyword
        */
        if (local_input_stream->is_open())
        {
            phrq_io->push_istream(local_input_stream);
            if (get_line() == KEYWORD)
            {
                ptr = line;
                copy_token(token, &ptr, &l);
                if (strcmp_nocase(token, "database") == 0)
                {
#ifdef PHREEQ98
                    user_database = string_duplicate(prefix_database_dir(ptr));
#else
                    user_database = string_duplicate(ptr);
#endif
                    if (string_trim(user_database) == EMPTY)
                    {
                        warning_msg("DATABASE file name is missing; default database will be used.");
                        user_database = (char *)free_check_null(user_database);
                    }
                }
            }
            phrq_io->pop_istream();
        }
        else
        {
            delete local_input_stream;
            //error_string = sformatf("Error opening file, %s.", in_file); //Li Jun comments the line. 2017-8-3.
            std::cerr << "Error opening input file " << *inputFileName << std::endl; // Li Jun added the line. 2017-8-3.
//            error_msg(error_string, STOP); //Li Jun comments the line. 2017-8-3.
            exit(-1);
        }

        /*
        *   Open data base
        */
 //       strcpy(query, "Name of database file?");  //Li Jun comments the line. 2017-8-6.
        env_ptr = getenv("PHREEQC_DATABASE");
        if (user_database != NULL)
        {
            strcpy(token, user_database);
        }
        else if (env_ptr != NULL)
        {
            strcpy(token, env_ptr);
        }
        else
        {
            strcpy(token, default_data_base);
        }

        std::ifstream * local_database_file;
 //       if (argc <= 1) //Li Jun comments this line
//        {
        std::string dataBaseFileName = token;
            local_database_file = open_input_stream( &dataBaseFileName, std::ios_base::in, false);            
//        }
//        else if (argc < 4) //Li Jun comments this line
//        {
//            local_database_file = open_input_stream(query, token, std::ios_base::in, true); //Li Jun comments this line
//        }
//        else if (argc >= 4) //Li Jun comments this line
//        {
//            if (user_database == NULL) //Li Jun comments this line
//            {
//                strcpy(token, argv[3]); //Li Jun comments this line
//            }
//            else
//            {
//#ifndef PHREEQCI_GUI //Li Jun comments this line
//                warning_msg("Database file from DATABASE keyword is used; command line argument ignored."); //Li Jun comments this line
//#endif
//            } //Li Jun comments this line
//            local_database_file = open_input_stream(query, token, std::ios_base::in, true); //Li Jun comments this line
//        } //Li Jun comments this line
        local_database_file->close();
        delete local_database_file;
//        screen_msg(sformatf("Database file: %s\n\n", token)); //Li Jun comments the line. 2017-8-6.
//        strcpy(db_file, token);                               //Li Jun comments the line. 2017-8-6.

 /*       output_msg(sformatf("   Input file: %s\n", in_file));
        output_msg(sformatf("  Output file: %s\n", out_file));
        output_msg(sformatf("Database file: %s\n\n", token)); */   //Li Jun comments the line. 2017-8-6.
        /*
        *   local cleanup
        */
        user_database = (char *)free_check_null(user_database);
        line = (char *)free_check_null(line);
        line_save = (char *)free_check_null(line_save);

        *db_cookie = new std::ifstream(dataBaseFileName, std::ios_base::in);
        *input_cookie = new std::ifstream(*inputFileName, std::ios_base::in);
    }
    catch (const PhreeqcStop&)
    {
        return get_input_errors();
    }
    return 0;
}