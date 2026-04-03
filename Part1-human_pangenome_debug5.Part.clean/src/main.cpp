//
//  main.cpp
//  human_pangenome
//
//  Created by kongshuang on 03/07/2022.
//


#include <iostream>
#include <getopt.h>
#include <vector>
#include <string.h>
#include <signal.h>
#include <exception>



#include "config.h"
#include "logger.h"
#include "parallel.h"
#include "endmapping.h"

#include "rmcont.h"

bool parseArgs(int argc, char**argv, std::string& inputdir, std::string& inputPath ,std::string& outputdir, std::string&configfile)
{
    auto printUsage = []()
    {
        std::cerr << "Usage: pan-genome pipeline \n"
                  << "--indir: input directory which contain sam file \n"
                  << "--outdir output directory\n"
		  << "--inPath input path file(can not set with indir simultaneously) \n"
                  << "--config(optional) configuration file\n";
    };
    int optionIndex = 0;
    static option longOptions[] =
    {
        {"indir",required_argument,0,0},
        {"outdir",required_argument,0,0},
        {"config",required_argument,0,0},
	{"inPath",required_argument,0,0},
        {0,0,0,0}
    };
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex))!= -1)
    {
        switch(opt)
        {
            case 0:
                if(!strcmp(longOptions[optionIndex].name, "indir"))
                    inputdir = optarg;
                else if (!strcmp(longOptions[optionIndex].name, "outdir"))
                    outputdir = optarg;
		else if (!strcmp(longOptions[optionIndex].name, "inPath"))
		    inputPath = optarg;		
                else if (!strcmp(longOptions[optionIndex].name, "config"))
                    configfile = optarg;
                break;

            case 'h':
                printUsage();
                exit(0);
        }
    }
    if(inputdir.empty() || outputdir.empty())
    {
        printUsage();
        return false;
    }
    return true;
}




int main(int argc, char** argv){

    signal(SIGSEGV, segfaultHandler);
    std::set_terminate(exceptionHandler);


    std::string inputdir = "None";
    std::string inputPath = "None";
    std::string outputdir;
    std::string configfile = "/src/ConfigFile" ;
    //std::string suffix = "sam";
    bool debugging = false;
    if(!parseArgs(argc, argv, inputdir, inputPath, outputdir, configfile)) return 1;
    Config::load(configfile); 
    
    if(outputdir.back() == '/') outputdir.pop_back();
    std::string logFile = outputdir + "/pan_genome.log";
    if(!demo_exists(outputdir)){
    	fs::create_directories(outputdir);
    }
    Logger::get().setDebugging(debugging);
    if (!logFile.empty()) Logger::get().setOutputFile(logFile);
    if(inputdir.back() == '/') inputdir.pop_back();
    Logger::get().info() << "Get input direcotry: \n" << inputdir;
    Logger::get().info() << "Get output direcotry: \n" << outputdir;
    Logger::get().info() << "Get config file: \n" << configfile;
    Logger::get().info() << "Get input path file: \n" << inputPath;
    std::string suffix = Config::get("suffix");
    FileContainer* files;

  //  FileContainer files(inputdir, suffix);
  //  files.getfilenames();
    if(inputdir == "None"){
    	files = new FileContainer(inputPath,suffix,true);
        files->getfilenamePath();

    }else{
        files = new FileContainer(inputdir,suffix,false);
        files->getfilenames();
    }

    auto infostruct = getinfo(*files);
    Logger::get().info() << "Loading files-- " << infostruct.filenumbers << " files have been found";
    Logger::get().info() << "Loading files-- " << infostruct.totalsize << "  totol size of files";
    Read2CleanCon r2c(infostruct.filenames, outputdir);
    r2c.build();
    if(Config::get("Vrmt") == "slow"){
        r2c.rmcont();
    }
    else {
        r2c.rmcontFast();
    };
    size_t Steps = (size_t)std::atoi(Config::get("Steps").c_str());
    if(Steps == 0 || Steps == 1 || Steps == 2){
        Endmapping EdMp(r2c);
        EdMp.buildMapping();
        EdMp.Endmap();
    }

    
    delete files;

 
    return 0;
}





