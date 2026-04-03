//
//  utils.h
//  human_pangenome
//
//  Created by kongshuang on 03/07/2022.
//

#ifndef utils_h
#define utils_h

#include <vector>
#include <algorithm>
#include <sstream>
#include <execinfo.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>  // for opendir
#include <iterator>
#include <cmath> // for ceil in file size recognise
#include <filesystem> // for file size , exist check
#include <fstream>  // for file size check
#include <cstdint> // check file exist
#include <string> // to_string
#include <curl/curl.h> // curl related function
#include <stdexcept> // throw runtime error in endmapping , read delta file 
#include <limits.h>
#include <zlib.h>
#include <regex> // regex_search
#include <tuple> // for get largest chr; 



#include "logger.h"


namespace fs = std::filesystem;
using std::filesystem::directory_iterator;

template<class T>
void vecRemove(std::vector<T>& v, T val)
{
    v.erase(std::remove(v.begin(), v.end(), val), v.end());
}

template<typename T>
T quantile(const std::vector<T>& vec, int percent)
{
    if (vec.empty()) return 0;
    //NOTE: there's a bug in libstdc++ nth_element,
    //that sometimes leads to a segfault. This is why
    //we have this inefficient impleemntation here
    //std::nth_element(vec.begin(), vec.begin() + vec.size() / 2,
    //                 vec.end());
    auto sortedVec = vec;
    std::sort(sortedVec.begin(), sortedVec.end());
    size_t targetId = std::min(vec.size() * (size_t)percent / 100,
                               vec.size() - 1);
    return sortedVec[targetId];
}

template<typename T>
T median(const std::vector<T>& vec)
{
    return quantile(vec, 50);
}

inline std::vector<std::string>
splitString(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) elems.push_back(item);
    return elems;
}

inline bool fileExists(const std::string& path)
{
    std::ifstream fin(path);
    return fin.good();
}

inline void segfaultHandler(int signal __attribute__((unused)))
{
    void *stackArray[20];
    size_t size = backtrace(stackArray, 10);
    Logger::get().error() << "Segmentation fault! Backtrace:";
    char** backtrace = backtrace_symbols(stackArray, size);
    for (size_t i = 0; i < size; ++i)
    {
        Logger::get().error() << "\t" << backtrace[i];
    }
    abort();
}

inline void exceptionHandler()
{
    static bool triedThrow = false;
    try
    {
        if (!triedThrow)
        {
            triedThrow = true;
            throw;
        }
    }
    catch (const std::exception &e)
    {
        Logger::get().error() << "Caught unhandled exception: " << e.what();
    }
    catch (...) {}

    void *stackArray[20];
    size_t size = backtrace(stackArray, 10);
    char** backtrace = backtrace_symbols(stackArray, size);
    for (size_t i = 0; i < size; ++i)
    {
        Logger::get().error() << "\t" << backtrace[i];
    }
    abort();
}

/*
struct pathhandle{
    pathhandle(std::string& fullpath,std::string& dirpath, std::string& filename):
    fullpath(fullpath), dirpath(dirpath), filename(filename)
    {}
    std::string fullpath;
    std::string dirpath;
    std::string filename;
};

inline pathhandle processpath(const std::string& strings){
    size_t dotPos = strings.rfind(".");
    size_t slash = strings.find_last_of("/");
    if(dotPos == std::string::npos || slash == std::string::npos) return strings;
    std::string dirpath = strings.substr(0, slash);
    std::string filename = std::string(&strings[slash+1], &strings[dotPos]);
    pathhandle paths(strings, dirpath, filename);
    return paths;
}
*/
struct FilesInfo{
    FilesInfo(std::vector<std::string>& filenames, uint16_t& filenumbers,std::string& totalsize):
    filenames(filenames),filenumbers(filenumbers),totalsize(totalsize)
    {}
    std::vector<std::string> filenames;
    uint16_t filenumbers;
    std::string totalsize;
};



class FileContainer
{
public:
    FileContainer(std::string& Input, std::string& suffix, const bool isPath):
    _suffix(suffix), _filecount(0), _isPath(isPath)
    {
        if(_isPath)
        {
            _pathFilePath = Input;
        }else{
            _dirname = Input;
        }
    }
    
    void getfilenamePath(){
        unsigned long long int totFilesize{};
        int i{};
        size_t BUF_SIZE = 32 * 1024 * 1024;
        char* rawBuffer = new char[BUF_SIZE];
        auto* fd = gzopen(_pathFilePath.c_str(),"rb");
        if(!fd)
        {
            delete[] rawBuffer;
            Logger::get().error()<< " can not open pathfile" << "\n";
	    throw std::runtime_error("Can't open " + _pathFilePath);
        }
        int lineNo = 1;
        std::string nextLine;
        try
        {
            while(!gzeof(fd))
            {
                for(;;)
                {  
                    char* read = gzgets(fd, rawBuffer, BUF_SIZE);
                    if(!read) break;
                    nextLine += read;
                    if(nextLine.empty()) break;
                    if(nextLine.back() == '\n')
                    {
                        nextLine.pop_back();
                        break;
                    }
                }
                if(nextLine.empty()) continue;
                if(nextLine.back() == '\r') nextLine.pop_back();
                if(nextLine[0] == '#')
                {
                    nextLine.clear();
                    lineNo++;
                    continue;
                }
                else
                {
                    _filenames.push_back(nextLine);
                    totFilesize += fs::file_size(nextLine.c_str());
                    _filecount++;
                }
                ++lineNo;
                nextLine.clear();
            }
        }
        catch(...)
        {
            gzclose(fd);
            Logger::get().error()<< " can not open pathfile" << "\n";
	    throw std::runtime_error("SomeThing Wrong during handle input files ");	    
        }
        delete[] rawBuffer;
        gzclose(fd);
     	for (; totFilesize >= 1024.; totFilesize /= 1024., ++i) { };
        totFilesize = std::ceil(totFilesize * 10.) / 10.;
     	 _totalsize = std::to_string(totFilesize) + "BKMGTPE"[i] ;
        
    }
    
    void getfilenames(){      
        std::string dirpath;
        size_t slash = _dirname.find_last_of("/");
        if (slash != std::string::npos) dirpath = _dirname.substr(0, slash + 1); //we had removed the last character if it is  '/'
        unsigned long long int totFilesize{};
        int i{};
        for(const auto& file : directory_iterator(_dirname)){
                std::string filename = file.path();
                size_t dotPos = filename.rfind(".");
                if(dotPos == std::string::npos) continue;
                std::string suffix = filename.substr(dotPos + 1);
                if(suffix == _suffix) {
                    totFilesize += fs::file_size(filename);
                    _filecount++;
                    _filenames.push_back(filename);
                }
        }
        
        for (; totFilesize >= 1024.; totFilesize /= 1024., ++i) { };
        totFilesize = std::ceil(totFilesize * 10.) / 10.;
        _totalsize = std::to_string(totFilesize) + "BKMGTPE"[i] ;
    }
private:
    bool _isPath;
    std::string _pathFilePath;
    std::string _dirname;
    std::string _suffix;
    std::string _totalsize;
    uint16_t _filecount;
    std::vector<std::string> _filenames;
    friend  FilesInfo getinfo(FileContainer);
};






/*
class FileContainer
{
public:
    FileContainer(std::string& dirname, std::string& suffix):
    _dirname(dirname), _suffix(suffix), _filecount(0)
    {}
    void getfilenames(){      
        std::string dirpath;
        //DIR *dpdf = opendir(_dirname.c_str());
        //struct dirent *epdf;
        size_t slash = _dirname.find_last_of("/");
        if (slash != std::string::npos) dirpath = _dirname.substr(0, slash + 1); //we had removed the last character if it is  '/'
        unsigned long long int totFilesize{};
        int i{};
        //if(dpdf != NULL){
          //  while ( epdf = readdir(dpdf) ){
        for(const auto& file : directory_iterator(_dirname)){
                std::string filename = file.path();
                size_t dotPos = filename.rfind(".");
                if(dotPos == std::string::npos) continue;
                std::string suffix = filename.substr(dotPos + 1);
                if(suffix == _suffix) {
                    //struct stat stat_buff;
                    //std::string filepath;
                    //long filesize;
                    //filepath = dirpath + "/" + filename;
                    //int rc = stat(filepath.c_str(), &stat_buff);
                    //rc == 0 ? filesize=stat_buff.st_size : -1;
                    //_totalsize += filesize;
                    totFilesize += fs::file_size(filename);
                    _filecount++;
                    _filenames.push_back(filename);
                }
        }
        
        for (; totFilesize >= 1024.; totFilesize /= 1024., ++i) { };
        totFilesize = std::ceil(totFilesize * 10.) / 10.;
        _totalsize = std::to_string(totFilesize) + "BKMGTPE"[i] ;
    }
private:
    std::string _dirname;
    std::string _suffix;
    std::string _totalsize;
    uint16_t _filecount;
    std::vector<std::string> _filenames;
    friend  FilesInfo getinfo(FileContainer);
};

*/

inline FilesInfo getinfo(FileContainer t){
    FilesInfo filsinformation(t._filenames, t._filecount, t._totalsize);
        return filsinformation;
}

template <typename Out>
inline void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


inline bool demo_exists(const fs::path& p, fs::file_status s = fs::file_status{}) 
{
    //std::cout << p;
    if(fs::status_known(s) ? fs::exists(s) : fs::exists(p))
        return true;
    else
        return false;
}

inline size_t curl_to_string(char* ptr, size_t size, size_t nmemb, void* data) 
{
        std::string* str = (std::string*)data;
        int x;
        for (x = 0; x < size * nmemb; ++x)
        {
                (*str) += ptr[x];
        }

        return size*nmemb;
}

inline std::string curlGetHtmlSource(std::string& link) 
{
    CURL* curl;
    CURLcode res;
    std::string html_txt;
    curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, &link[0]);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, true);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_to_string);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &html_txt);

        /* Perform the request, res will get the return code */
        res = curl_easy_perform(curl);

        /* Check for errors */
        if (res != CURLE_OK)
        {
            fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
                curl_easy_cleanup(curl);
                throw std::runtime_error("Can't get html source");
        }

        /* always cleanup */
        curl_easy_cleanup(curl);
    }
    return html_txt;
}

template <class T>
std::tuple<std::string, bool, float, float> getLargestChr(T& Inmap){
    float sum = 0;
    float max = 0 ;
    std::string maxChr;
    for(const auto& i : Inmap.lock_table()){
        sum = sum + i.second;
	    Logger::get().info()<< "i. second is " << i.second << "\n";
        if (i.second > max){
            max = i.second;
            maxChr = i.first;
        }
    }
    if(sum == 0){
        return std::make_tuple(maxChr, false, sum, max);
    }
    else{
        return std::make_tuple(maxChr, true, sum, max);
    }
}


inline float getClust(std::vector<int>& chrVec){
    std::sort(chrVec.begin(), chrVec.end());
    float median = chrVec[chrVec.size() / 2];
    float clusterNum = 0;
    std::vector<int> deleteVec;
    Logger::get().info() << " the size of input is "<<chrVec.size() << "\n";
    for(auto i : chrVec)
    {
        if((i <= median+2000) && (i >= median-2000))
        {
            clusterNum++ ;
	    Logger::get().info() << "get cluster num " << i <<"\n";
	    Logger::get().info() << "cluster num is " << clusterNum << "\n";
        }
        else{
            deleteVec.push_back(i);
	    Logger::get().info() << "remove un cluster value "<< i << "\n";
        }
    }
    for(auto i:deleteVec){
        if(std::find(chrVec.begin(), chrVec.end(), i) != chrVec.end()){
            chrVec.erase(std::remove(chrVec.begin(), chrVec.end(), i), chrVec.end());
        }
    }
    
    for(auto i:chrVec){
        std::cout<<i << "\t";
    }
    return clusterNum;
}


template<class T>
int getExtreamDeltaMatch(T& v1, const bool& largest, const bool& first){
    int index = 0;
    int reIndex = 0;
    int SmallestValue = INT_MIN ;
    int LargestValue = INT_MAX ;
    for(const auto& i : v1){
        if(largest){
            if(first){
                if(i.first >= SmallestValue){
                    SmallestValue = i.first;
                    reIndex = index;
                }
            }
            else{
                if(i.second >= SmallestValue){
                    SmallestValue = i.second;
                    reIndex = index;
                }
            }
        }
        else{
            if(first){
                if(i.first <= LargestValue){
                    LargestValue = i.first;
                    reIndex = index;
                }
            }
            else{
                if(i.second <= LargestValue){
                    LargestValue = i.second;
                    reIndex = index;
                }
            }
            
        }
        index++;
        
    }
    return reIndex;
}





template<class T1>
bool check_key(T1 m, std::string key)
{
	if(m.count(key) == 0)
		return false;
	return true;
}

inline void FastaLenSelect(std::string InFaFile, std::string OutFaFile){
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(InFaFile.c_str(),"rb");
    FILE* fout = fopen(OutFaFile.c_str(),"w");
    if (!fout) throw std::runtime_error("Can't open " + OutFaFile);
    
    if(!fd)
    {
        delete[] rawBuffer;
    }
    int lineNo = 1;
    std::string nextLine;
    std::string header;
    std::string sequence;
    bool keep = true;
    try
    {
        while(!gzeof(fd))
        {
            for(;;)
            {
                char* read = gzgets(fd, rawBuffer, BUF_SIZE);
                if(!read) break;
                nextLine += read;
                if(nextLine.empty()) break;
                if(nextLine.back() == '\n')
                {
                    nextLine.pop_back();
                    break;
                }
            }
            if(nextLine.empty()) continue;
            if(nextLine.back() == '\r') nextLine.pop_back();

            if(nextLine[0] == '>' )
            {
                if(!header.empty())
                {
                    if (sequence.empty()) break;
                    header.clear();
                    sequence.clear();
                }
                header = nextLine.erase(0,1);
                std::regex const e{".*LN:i:(\\d+)\\b.*"};
                std::smatch matches;
                if(std::regex_search(header, matches, e)){
                    int length = std::stoi(matches[1].str());
                    if(length >= 500)
                    {
                        keep = true;
                    }else{
                        keep = false;
                    }
                }
            }
            else
            {
                std::copy(nextLine.begin(), nextLine.end(),
                          std::back_inserter(sequence));
                if(keep == true){
                    header = ">k" + header + "\n";
                    sequence = sequence + "\n";
                    fwrite(header.data(), sizeof(header.data()[0]),
                           header.size(),fout);
                    fwrite(sequence.data(), sizeof(sequence.data()[0]),
                           sequence.size(),fout);
                }
                header.clear();
                sequence.clear();
            }
            ++lineNo;
            nextLine.clear();
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw std::runtime_error(" something wrong withing the read write process  " );

    }
    delete[] rawBuffer;
    gzclose(fd);
    fclose(fout);
}

inline std::vector<std::string> getFaHead(std::string InFaFile){
    std::vector<std::string> Headers;
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(InFaFile.c_str(),"rb");
    if(!fd)
    {
        delete[] rawBuffer;
        throw std::runtime_error("Can't open " + InFaFile);
    }
    int lineNo = 1;
    std::string nextLine;
    std::string header;
    std::string sequence;
    bool keep = true;
    try
    {
        while(!gzeof(fd))
        {
            for(;;)
            {
                char* read = gzgets(fd, rawBuffer, BUF_SIZE);
                if(!read) break;
                nextLine += read;
                if(nextLine.empty()) break;
                if(nextLine.back() == '\n')
                {
                    nextLine.pop_back();
                    break;
                }
            }
            if(nextLine.empty()) continue;
            if(nextLine.back() == '\r') nextLine.pop_back();

            if(nextLine[0] == '>' )
            {
                if(!header.empty())
                {
                    if (sequence.empty()) break;
                    header.clear();
                    sequence.clear();
                }
                header = nextLine.erase(0,1);
                Headers.push_back(header);
            }
            else
            {
                header.clear();
                sequence.clear();
            }
            ++lineNo;
            nextLine.clear();
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw std::runtime_error(" something wrong withing the read write process  " );

    }
    delete[] rawBuffer;
    gzclose(fd);
    return Headers;
}


#endif /* utils_h */



