//
//  rmcont.hpp
//  human_pangenome
//
//  Created by kongshuang on 03/07/2022.
//

#ifndef rmcont_hpp
#define rmcont_hpp

//control the pipe flow by the building structure.

/* utils.h
#include <vector>
#include <algorithm>
#include <sstream>
#include <execinfo.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <unistd.h> config.h
#include <unordered_map>
#include <fstream>
#include "logger.h"
#include "utils.h"

#include <iostream> logger.h
#include <fstream>
#include <ctime>
#include <stdexcept>
*/
#include "parallel.h"
#include "cuckoohash_map.hh"
#include "stdc++.h"
#include "config.h"
#include "utils.h"
#include <utility>
#include <stdlib.h>
#include <zlib.h>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <regex>
#include <sstream> // string to size_t in curl
#include <functional>  //std function 

struct buildfiles
{
    buildfiles(std::string& fullpath,std::string& dirpath, std::string& filename):
    fullpath(fullpath), dirpath(dirpath), filename(filename)
    {}
    void getaligfiles(std::string& outpath){
        unalignedreads.first = outpath + "/unaligned_first_" + filename + ".fastq";
        unalignedreads.second = outpath + "/unaligned_second_" + filename + ".fastq";
        matealignedreads.first = outpath + "/matealigned_first_" + filename + ".fastq";
        matealignedreads.second = outpath + "/matealigned_second_" + filename + ".fastq";
        linkagesam = outpath + "/linkage_" + filename + ".sam";
        contigfileDir = outpath + "/contigs_" + filename ;
        contigfile = contigfileDir + "/final.contigs.fa";
        outputpath = outpath;
	linkageSamPopins = contigfileDir + "/non_ref.bam";
        
    }
    std::string fullpath;
    std::string dirpath;
    std::string filename;
    std::pair<std::string, std::string> unalignedreads;
    std::pair<std::string, std::string> matealignedreads;
    std::string linkagesam;
    std::string contigfile;
    std::string contigfileDir;
    std::string outputpath;
    std::string linkageSamPopins; 
};

//typedef cuckoohash_map<std::string, std::vector<ContigInfo>>  MapAcVecContigInfo; // accession number -- vector of contiginfo structure
//typedef std::unordered_map<std::string, ContigInfo> MapHdContigInfo; //unordered map header(contig header) ---- ContigInfo

struct ContigInfo
{
    ContigInfo():
    isnovel(false),haveClassInfo(false), haveAcinfo(false), haveBlastAc(false) {}
    uint32_t checkAcClinfo(){
        return 0;
    };
    std::string contigname;
    std::string accession;
    std::string contigSequence;
    bool haveAcinfo;
    bool haveBlastAc; 
    bool haveClassInfo; 
    bool isnovel; 
    size_t taxid;
    std::string classinfo;
    double blastIdentity; 
    double balstLength;
};


//typedef cuckoohash_map<std::string, std::vector<ContigInfo>>  MapAcVecContigInfo; // accession number -- vector of contiginfo structure
typedef std::unordered_map<std::string, ContigInfo> MapHdContigInfo; //unordered map header(contig header) ---- ContigInfo 

struct buildMapinfo{
    buildMapinfo(std::string& cleancontig,buildfiles& buildfile,std::string InOutput, std::string& Infilename):
    cleanContig(cleancontig),buildFile(buildfile), output(InOutput),filename(Infilename){}
    
    std::string cleanContig;
    std::string contigIndex;
    std::string ReadToContigSam;
    std::string ReadToContigTxt;
    std::string mateRegion;
    buildfiles buildFile;
    std::string output;
    std::string filename;
    std::pair<std::string, std::string> matealignedreads;
    std::string linkageSam;
    std::string passMate;
    std::string sortMateRegion;
    void buildingmap(){ // output/sample1BowIndex
        contigIndex = output +  filename + "BowIndex";
        matealignedreads.first = buildFile.contigfileDir + "/single.fastq";
	matealignedreads.second = buildFile.contigfileDir + "/paired.2.fastq";
        ReadToContigSam = output + filename + "ReadToContig.sam" ;
        ReadToContigTxt = output + filename + "ReadToContig.txt" ;
        linkageSam = buildFile.linkageSamPopins;
        mateRegion = output + filename + "MateRegion";
        passMate = output + filename + "PassMate";
        mateRegion = output + filename + "MateRegion";
        sortMateRegion = mateRegion + "Sort";
    }
};
typedef cuckoohash_map<std::string, buildMapinfo> MapHdBuildMap;



//raw read to clean contigs
class Read2CleanCon{
public:
    class ParseException : public std::runtime_error 
    {
    public:
        ParseException(const std::string & what):
            std::runtime_error(what)
        {}
    };
    Read2CleanCon(std::vector<std::string>& samFilenames, std::string& outpath):
    _samFilenames(samFilenames), _outpath(outpath)
    {}
    // build  unmapped pair reads files and one mate mapped read files and one mate mapped sam files. And assemble contigs.
    // please all the aditional software is available in your Linux system path
    void build();
    void rmcont();
    void loadAccession();
    void loadLineage();
    void loadBlastfile(MapHdContigInfo& ContigInfos,const std::string& blastpath);
    size_t searchNcbi(std::string& accessionId);
    void writecleanContig(MapHdContigInfo& ContigInfos, std::string& outputname, std::string& outUnmap,const bool isFast);
    void loadContig(MapHdContigInfo& ContigInfos, const std::string& contigpath);
    void rmcontFast();
    /*
    struct buildfiles
    {
        buildfiles(std::string& fullpath,std::string& dirpath, std::string& filename):
        fullpath(fullpath), dirpath(dirpath), filename(filename)
        {}
        void getaligfiles(std::string& outpath){
            unalignedreads.first = outpath + "/unaligned_first_" + filename + ".fastq";
            unalignedreads.second = outpath + "/unaligned_second_" + filename + ".fastq";
            matealignedreads.first = outpath + "/matealigned_first_" + filename + ".fastq";
            matealignedreads.second = outpath + "/matealigned_second_" + filename + ".fastq";
            linkagesam = outpath + "/linkage_" + filename + ".sam";
            contigfile = outpath + "/contigs_" + filename + ".fasta";
        }
        std::string fullpath;
        std::string dirpath;
        std::string filename;
        std::pair<std::string, std::string> unalignedreads;
        std::pair<std::string, std::string> matealignedreads;
        std::string linkagesam;
        std::string contigfile;
    }; */
    struct buildrmtfiles  
    {
        buildrmtfiles(std::string& Incontigpath, std::string& Infilename, buildfiles buildfiles):
        contigPath(Incontigpath),fileName(Infilename),firstStageBuildfile(buildfiles)
        {}
        void buildingpath(std::string& outpath)
        {
            blastRes = outpath + "/blast_" + fileName + ".bl";
            cleanContig = outpath + "/novelContig_" + fileName + ".fa";
            cleanHuman = outpath + "/humanContig_" + fileName + ".fa";
            unmapRead = outpath + "/unaligned_" + fileName + ".fa";
        }
        std::string fileName;
        std::string contigPath;
        std::string blastRes;
        std::string cleanContig;
        std::string cleanHuman;
        std::string unmapRead;
        buildfiles firstStageBuildfile;
    };
    void buildpath();
    void buildContigclass(MapHdContigInfo& ContigInforst/* MapAcVecContigInfo& AcContigInfo */);
    
    MapHdBuildMap buildEndmap() ;

private:
    size_t _novelLength;
    size_t _novlCount;
    //cuckoohash_map<std::string, std::vector<ContigClass>> _buildclass; 
    cuckoohash_map<std::string, buildfiles> _buildfiles;
    cuckoohash_map<std::string, uint32_t > _buildfails;
    cuckoohash_map<std::string, buildrmtfiles> _buildrmtfiles;
    cuckoohash_map<std::string, size_t> _accession2id;
    cuckoohash_map<size_t, std::string> _id2lineage;
    cuckoohash_map<std::string , size_t> _contigclassfail;
    //cuckoohash_map<std::string, std::vector<std::string>> _contigfiles;
    std::vector<std::string> _samFilenames;
    std::string _outpath;
    
};


#endif /* rmcont_hpp */

