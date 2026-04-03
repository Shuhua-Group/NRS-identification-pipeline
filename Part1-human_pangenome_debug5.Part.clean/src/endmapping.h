//
//  endmapping.hpp
//  human_pangenome
//
//  Created by kongshuang on 03/07/2022.
//

#ifndef endmapping_hpp
#define endmapping_hpp

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
/* rmcont.h
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
*/
// typedef std::unordered_map<std::string, ContigInfo> MapHdContigInfo;

#include "utils.h"
#include "rmcont.h"

const static std::unordered_map<std::string, int> CHRLENGTHMAP ={
    {"chr1" , 248956422},
    {"chr2" , 242193529},
    {"chr3" , 198295559},
    {"chr4" , 190214555},
    {"chr5" , 181538259},
    {"chr6" , 170805979},
    {"chr7" , 159345973},
    {"chr8" , 145138636},
    {"chr9" , 138394717},
    {"chr10" , 133797422},
    {"chr11" , 135086622},
    {"chr12" , 133275309},
    {"chr13" , 114364328},
    {"chr14" , 107043718},
    {"chr15" , 101991189},
    {"chr16" , 90338345},
    {"chr17" , 83257441},
    {"chr18" , 80373285},
    {"chr19" , 58617616},
    {"chr20" , 64444167},
    {"chr21" , 46709983},
    {"chr22" , 50818468},
    {"chrX" , 156040895},
    {"chrY" , 57227415},
};

typedef std::vector<std::pair<int,int>> deltaMatch; 

typedef std::pair<std::string,std::string> endRegions;

struct ContigLine{
    ContigLine(std::string IncontigName, std::string Inchr, uint32_t InChrStart, uint32_t InchrEnd, std::string InLinkRead, uint32_t IncontigStart, uint32_t IncontigEnd, std::string InstrandContig, std::string InstrandRef, bool InisRight):
    LinecontigName(IncontigName),chr(Inchr), chrStart(InChrStart), chrEnd(InchrEnd), LinkRead(InLinkRead),
    contigStart(IncontigStart), contigEnd(IncontigEnd), strandContig(InstrandContig),strandRef(InstrandRef),isRight(InisRight){}
    ContigLine(){
        LinecontigName = "NONE";
    }
    std::string LinecontigName;
    std::string chr;
    uint32_t chrStart;
    uint32_t chrEnd;
    std::string LinkRead;
    uint32_t contigStart;
    uint32_t contigEnd;
    std::string strandContig; 
    std::string strandRef;
    bool isRight;
};

struct matchInfo{
    matchInfo(const int& InsertPoint, const bool& Reverse):
    InsertPoint(InsertPoint),Reverse(Reverse){}
    matchInfo(){}
    
    void setRef(const int& InRefLen, const int& InRefStart, const int& InRefEnd){
        RefLen = InRefLen;
        RefStart = InRefStart;
        RefEnd = InRefEnd;

    }
    void setRef2(const std::string& InrefName, const int& InRefSTART, const int&InRefEND){
        refName = InrefName;
        RefSTART = InRefSTART;
        RefEND = InRefEND;
    }
    void setContig(const int& IncontigIns, const int& IncontigStart,const int& IncontigEnd){
        contigMatchIns = IncontigIns;
        contigStart = IncontigStart;
        contigEnd = IncontigEnd;
    }
    int contigMatchIns;
    int InsertPoint;
    bool Reverse;
    int RefLen;
    int RefStart; 
    int RefEnd; 
    //int contigLen;
    int contigStart;
    int contigEnd;
    std::string refName;
    
    int RefSTART; 
    int RefEND;
};

struct ContigMapInfo{
    ContigMapInfo(std::string& IncontigName, uint32_t Inlength, std::vector<ContigLine> InContigLinkMap, std::string InoutPath):
    ContigName(IncontigName),exRightRegion(false),exLeftRegion(false),exBothRegion(false), contigLength(Inlength), contigLinkMap(InContigLinkMap), outPath(InoutPath), leftMatch(false), rightMatch(false),BothMatch(false){}
    std::string ContigName;
    endRegions leftRegion;
    endRegions rightRegion;
    bool exRightRegion;
    bool exLeftRegion;
    bool exBothRegion;
    
    uint32_t contigLength;
    std::string contigEndInfo;
    
    std::string outPath;
    
    std::string rightDelta; 
    std::string filtedRightDelta;
    std::string leftDelta;
    std::string filtedLeftDelta;
    
    bool leftMatch;
    matchInfo leftMatchInfo;
    bool rightMatch;
    matchInfo rightMatchInfo;
    
    std::pair<std::vector<int>, std::vector<int>> leftInsertSitua;
    std::pair<std::vector<int>, std::vector<int>> rightInsertSitua;
    
    bool BothMatch;
    
    void ConstructRegion(uint32_t chrRegion, uint32_t contigRegion, std::string& cleanContig); 
    std::pair<std::vector<int>,  std::vector<int>> getRegionMatch(const bool& Right);

    std::vector<ContigLine> contigLinkMap ; 
};

typedef std::unordered_map<std::string, ContigMapInfo> MpContigMapinfo; 



class Endmapping
{
public:
    class ParseException : public std::runtime_error 
    {
    public:
        ParseException(const std::string & what):
            std::runtime_error(what)
        {}
    };
    Endmapping(Read2CleanCon rd2con): 
    _rd2con(rd2con){}
    void buildMapping();
    void Endmap();
    MpContigMapinfo contigMapping(std::string& mateRegionTxt,int32_t contigTailRange, int32_t chrRange,int32_t mapContigTail, std::string& output, std::string& cleanContig);
    void writeoutput(std::vector<std::string>& v1,std::string& outPath);
    //void writeBothEnd(std::string& outPath);
private:
    MapHdBuildMap _buildmap;
    //cuckoohash_map<std::string, buildrmtfiles> _buildmappingfiles;
    Read2CleanCon _rd2con;
};

/*
struct CountValue{
    CountValue(uint32_t InchrStart, int Incount, bool InisRight):
    count(Incount){
        isRight.push_back(InisRight);
        chrStart.push_back(InchrStart);
    }
    std::vector<bool> isRight;
    std::vector<uint32_t> chrStart;
    int count;
};  
*/


#endif /* endmapping_hpp */
