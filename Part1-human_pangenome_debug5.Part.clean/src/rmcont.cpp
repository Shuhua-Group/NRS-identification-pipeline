//
//  rmcont.cpp
//  human_pangenome
//
//  Created by kongshuang on 03/07/2022.
//

#include "rmcont.h"
#include "logger.h"

#include <functional>
#include <regex>


// the file will stored in the related container and generate corresonding linux command output
void Read2CleanCon::build()
{
    Logger::get().info() << "Start building aligned/unaligned files for " << _samFilenames.size()<<" sam files";
    
    if(!demo_exists(_outpath)){
        fs::create_directories(_outpath);
    }
    this->buildpath();
    
    std::function<void(const std::string& filepath)> buildfiles =
    [this](const std::string& filepath)
    {
        std::string IniEnv = Config::get("IniEnv");
        std::string AssembleTool = Config::get("Assembletool");
        auto pathstruct = _buildfiles.find(filepath); 
        //samfile = filename.substr(slash+1, filename.length())



        uint32_t statvalue5 = 0;
        uint32_t stateSum = 0;

        if(AssembleTool == "megahit")
        {                
            std::vector<std::string> megahitfiles = {pathstruct.unalignedreads.first,pathstruct.unalignedreads.second,pathstruct.matealignedreads.first, pathstruct.matealignedreads.second};

            std::string megaInputfiles = "";
            
            for(auto& i : megahitfiles){
                if(fs::file_size(i) != 0){
                    megaInputfiles = megaInputfiles + "," + i ;
                }
            }
            uint32_t statvalue1 = system((IniEnv + " samtools fastq -f 12 " + filepath + " -1 " + \
                        pathstruct.unalignedreads.first + " -2 " + pathstruct.unalignedreads.second ).c_str());
            uint32_t statvalue2 = system((IniEnv + " samtools fastq -f 68 -F 8 " + filepath + " > " +  pathstruct.matealignedreads.first).c_str());
            uint32_t statvalue3 = system((IniEnv + " samtools fastq -f 132 -F 8 " + filepath + " > " +  pathstruct.matealignedreads.second).c_str());
            uint32_t statvalue4 = system((IniEnv + " samtools view -h -f 8 -F 4 " + filepath+ " > " +  pathstruct.linkagesam).c_str());
            statvalue5 = system((IniEnv  +" megahit -r" + megaInputfiles + " -o " + pathstruct.contigfileDir).c_str());
            stateSum = statvalue1 + statvalue2 + statvalue3 + statvalue4 + statvalue5;
        }
        else
        {
            std::string Popins = Config::get("PopinsPath");
            statvalue5 = system((Popins + " assemble -s " + pathstruct.contigfileDir + " " + filepath).c_str());
            std::string InputFile = pathstruct.contigfileDir + "/assembly_final.contigs.fa";
            FastaLenSelect(InputFile, pathstruct.contigfile);
            stateSum = statvalue5;		
	    }

        if(stateSum != 0){
            _buildfiles.erase(filepath);
           _buildfails.insert(pathstruct.filename, stateSum);
        }// remove files from the next step analysis.
    };
    // for bam files
    /*
    std::function<void(const std::string& filepath)> buildfilesBam =
    [this](const std::string& filepath)
    {
        auto pathstruct = _buildfiles.find(filepath);
        //samfile = filename.substr(slash+1, filename.length())
        uint32_t statvalue1 = system(("samtools fastq -f 12 " + filepath + " -1 " +  \
                    pathstruct.unalignedreads.first + " -2 " + pathstruct.unalignedreads.first ).c_str());
        uint32_t statvalue2 = system(("samtools fastq -f 68 -F 8 " + filepath + " > " +  pathstruct.matealignedreads.first).c_str());
        uint32_t statvalue3 = system(("samtools fastq -f 132 -F 8 " + filepath + " > " +  pathstruct.matealignedreads.second).c_str());
        uint32_t statvalue4 = system(("samtools view -h -f 8 -F 4 " + filepath+ " > " +  pathstruct.linkagesam).c_str());
        uint32_t statvalue5 = system(("megahit -r " + pathstruct.unalignedreads.first + " " + \
                     " " + pathstruct.unalignedreads.first + " " + pathstruct.matealignedreads.first + \
                     " " + pathstruct.matealignedreads.second + " " + "-o " + pathstruct.contigfile ).c_str());
        statvalue1  = statvalue1 + statvalue2 + statvalue3 + statvalue4 + statvalue5;
        _buildfails.insert(pathstruct.filename, statvalue1);
        if(statvalue1!=0) _buildfiles.erase(filepath); // remove files from the next step analysis.
    };
     */
    
    size_t threadnum = (size_t)std::atoi(Config::get("buildThread").c_str());
    size_t Steps = (size_t)std::atoi(Config::get("Steps").c_str());
    if(Steps == 0){
        processInParallel(_samFilenames, buildfiles, threadnum, true);
        uint16_t failnumber{};
        for (const auto& filename : _buildfails.lock_table())
        {
            if(filename.second != 0)
            {
                Logger::get().warning() << filename.first << " failed for build, please check file states";
                failnumber ++ ;
            }
        }
        Logger::get().warning() << "There are totally "<<failnumber << " files failed for build";
        }
}




void Read2CleanCon::buildpath(){

    Logger::get().info()<<"build paths for "<< _samFilenames.size() << " files" <<"\n";
    //std::string outBuildings = _outpath + "/hitContig";
    //fs::create_directories(outBuildings);
    for(auto fullpath : _samFilenames)
    {
        
        std::string dirpath;
        //std::string filename; 
        size_t slash = fullpath.find_last_of("/");
        size_t dotPos = fullpath.rfind(".");
        //dirpath = filename.substr(0,slash);
        dirpath = fullpath.substr(0,slash);
        std::string filename(fullpath.c_str() + slash + 1, fullpath.c_str() + dotPos);
        
        std::string outBuildings = _outpath + "/pan_" + filename;
        fs::create_directories(outBuildings);
        
        buildfiles buildedfile(fullpath, dirpath, filename);
        buildedfile.getaligfiles(outBuildings);
        _buildfiles.insert(fullpath, buildedfile);
    }

}



void Read2CleanCon::rmcontFast()
{
    std::vector<std::string> Contigfiles;
    for (const auto& buildingfiles : _buildfiles.lock_table())  
    {
        std::string contigfile = buildingfiles.second.contigfile ;
        std::string filename = buildingfiles.second.filename;
        buildfiles bfs = buildingfiles.second;
        //buildrmtfiles buildrmt(buildingfiles.second.contigfile, buildingfiles.second.filename);
        buildrmtfiles buildrmt(contigfile, filename, bfs);
        std::string outpath = _outpath + "/pan_" + filename;
        buildrmt.buildingpath(outpath);
        _buildrmtfiles.insert(contigfile, buildrmt);
        Contigfiles.push_back(contigfile);
    }
    std::function<void(const std::string& contigfilepath)> blastrmt =
    [this](const std::string& contigfilepath)
    {
        buildrmtfiles curbuild = _buildrmtfiles.find(contigfilepath);
        std::string Inienv = Config::get("IniEnv");
    // blast parallel
                
        
        uint32_t statvalue1 = system((Inienv + " blastn -query " + contigfilepath + " -db " + Config::get("fastBlastindex") \
                                       + " -out " + curbuild.blastRes + " -evalue " +  Config::get("evalue") \
                                       + " -outfmt " + Config::get("outfmt") + " -max_target_seqs " + Config::get("maxTargetSeqs") \
                                       + " -num_threads " + Config::get("blastNThread")).c_str());
        MapHdContigInfo contigs; // string -> Contiginfo // out the scope and then it will die! the memory will release
        this->loadContig(contigs,curbuild.contigPath); //
        loadBlastfile(contigs, /*AcContigInfo*/ curbuild.blastRes); 

        writecleanContig(contigs, curbuild.cleanContig, curbuild.unmapRead, true);
 	uint32_t statvalue3 = system((Inienv + " samtools faidx " + curbuild.cleanContig ).c_str());       
    };
    size_t threadnum = (size_t)std::atoi(Config::get("blastThread").c_str());
    size_t Steps = (size_t)std::atoi(Config::get("Steps").c_str());
    if(Steps == 0 || Steps == 1){
        processInParallel(Contigfiles, blastrmt, threadnum, true);
    }
    
    
}



void Read2CleanCon::rmcont() 
{
    std::vector<std::string> Contigfiles;
    for (const auto& buildingfiles : _buildfiles.lock_table())  
    {
        std::string contigfile = buildingfiles.second.contigfile ;
        std::string filename = buildingfiles.second.filename;
        buildfiles bfs = buildingfiles.second;
        //buildrmtfiles buildrmt(buildingfiles.second.contigfile, buildingfiles.second.filename);
        buildrmtfiles buildrmt(contigfile, filename, bfs);
        std::string outpath = _outpath + "/pan_" + filename;
        buildrmt.buildingpath(outpath);
        _buildrmtfiles.insert(contigfile, buildrmt);
        Contigfiles.push_back(contigfile);
    }
    // load accessiontoid
    size_t Steps = (size_t)std::atoi(Config::get("Steps").c_str());
    if(Steps == 0 || Steps == 1){
        this->loadAccession();
        // load id to class file

        this->loadLineage();

        //parallel plast
    }

    // parallel convert accession to id , then transfer the id to class information.
    std::function<void(const std::string& contigfilepath)> blastrmt =
    [this](const std::string& contigfilepath)
    {
        buildrmtfiles curbuild = _buildrmtfiles.find(contigfilepath);
    // blast parallel
        
        uint32_t statvalue2 = system(("module load blast; blastn -query " + contigfilepath + " -db " + Config::get("blastindex") \
                                        + " -out " + curbuild.blastRes + " -evalue " +  Config::get("evalue") \
                                        + " -outfmt " + Config::get("outfmt") + " -max_target_seqs " + Config::get("maxTargetSeqs") \
                                        + " -num_threads " + Config::get("blastNThread")).c_str());
        //build vector of contigclass for all contigs in this contig file.
        // Read in contig file --> unordered map(with in this scope)
        // read in blast file --> unordered map(with in this scope)
        // match between contig and blast ; match between blast and accession; match between accession and lineage
        // every step update the fail function and remove contig which failed in the searching process.
        // or we could put these things into a class
        MapHdContigInfo contigs; // string -> Contiginfo // out the scope and then it will die! the memory will release
        
        this->loadContig(contigs,curbuild.contigPath); 

        //MapAcVecContigInfo AcContigInfo;
        loadBlastfile(contigs, /*AcContigInfo*/ curbuild.blastRes); 
        buildContigclass(contigs); 
        writecleanContig(contigs, curbuild.cleanContig, curbuild.unmapRead, false);
    };
    size_t threadnum = (size_t)std::atoi(Config::get("endmapThread").c_str());
    
    if(Steps == 0 || Steps == 1){
        processInParallel(Contigfiles, blastrmt, threadnum, true);
    }
    


}

void Read2CleanCon::buildContigclass(MapHdContigInfo& ContigInforst/* MapAcVecContigInfo& AcContigInfo */) // for different contig file: build the info vector
{
    /*for(auto& vecContiginfo : AcContigInfo.lock_table())
    {
        if(_accession2id.contains(vecContiginfo.first))
        {
            for(auto i : vecContiginfo.second)
            {
                i.haveAcinfo = true;
            }
        }
    }*/
    // new features !!!!!!! add NcbiKeyWord which is used to select specific species
    std::string keyword = Config::get("NcbiKeyWord");
    std::vector<std::string> keyWords = split(keyword, '-');
    

    for(auto& StContiginfo : ContigInforst)  //st : structure
    {
        if(_accession2id.contains(StContiginfo.second.accession))
        {
            StContiginfo.second.haveAcinfo = true;
            StContiginfo.second.taxid = _accession2id.find(StContiginfo.second.accession);
        }
        else
        { // searching from internet
//            size_t searchID = searchNcbi(StContiginfo.second.accession);
            size_t searchID = 0 ; 
            if(searchID != 0){
                StContiginfo.second.taxid = searchID;
                StContiginfo.second.haveAcinfo = true;
            }
        }
        if(StContiginfo.second.haveAcinfo == true){
            if(_id2lineage.contains(StContiginfo.second.taxid))
            {
                StContiginfo.second.classinfo = _id2lineage.find(StContiginfo.second.taxid);
                StContiginfo.second.haveClassInfo = true;
                /*
                if (std::regex_match(StContiginfo.second.classinfo.c_str(), std::regex("genus:Homo")) |
                    std::regex_match(StContiginfo.second.classinfo.c_str(), std::regex("order:Primates")))
                {
                    StContiginfo.second.isnovel = true;
                }
                */
               for(auto keyword:keyWords)
               {
                    if(std::regex_match(StContiginfo.second.classinfo.c_str(), std::regex(".*:" + keyword + ";.*")))
                    {
                        StContiginfo.second.isnovel = true;
                    }
               }

            }
        }
    }
}


void Read2CleanCon::loadAccession()
{

    std::string accessionfile = Config::get("accession");
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(accessionfile.c_str(),"rb");
    if(!fd)
    {
        delete[] rawBuffer;
        throw ParseException("Can't open accession2taxid file "); 
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
            if(nextLine[0] == '#' || lineNo == 1 )
            {
                nextLine.clear();
                lineNo++;
                continue;
            }
            else
            {
                std::vector<std::string> x = split(nextLine, '\t');
                _accession2id.insert(x[0], std::atoi(x[2].c_str()));
            }
            ++lineNo;
            nextLine.clear();
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw ParseException("read in accession file error, please check again.");
    }
    delete[] rawBuffer;
    gzclose(fd);
}



void Read2CleanCon::loadLineage()
{
    std::string accessionfile = Config::get("lineage");
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(accessionfile.c_str(),"rb");
    if(!fd)
    {
        delete[] rawBuffer;
        throw ParseException("Can't open accession2taxid file "); 
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
            if(lineNo == 1 )
            {
                nextLine.clear();
                lineNo++;
                continue;
            }
            else
            {
                std::regex e1("\\t\\|\\t"); 
                std::string re1 = std::regex_replace(nextLine.c_str(), e1, "-");
                std::regex e2("\\t\\|");
                std::string re2 = std::regex_replace(re1.c_str(), e2, "-");
                //std::cout<< re2 <<"\n";
                std::vector<std::string> x = split(re2, '-');
                std::string res;
                res = "name:" + x[1] + "; species:" + x[2] + "; " + "genus:" + x[3] + "; " + "family:" + x[4] +
                        "; " + "order:" + x[5] + "; " + "class:" + x[6] + "; " + "phylum:" +
                        x[7] + "; " + "kingdom:" + x[8] + "; " + "superkingdom:" + x[9] + "; ";
                _id2lineage.insert(std::atoi(x[0].c_str()), res);

            }
            ++lineNo;
            nextLine.clear();
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw ParseException("read in accession file error, please check again.");
    }
    delete[] rawBuffer;
    gzclose(fd);
}



void Read2CleanCon::loadContig(MapHdContigInfo& ContigInfos, const std::string& contigpath)
{

    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(contigpath.c_str(),"rb");
    //std::vector<ContigInfo> v1;
    if(!fd)
    {
        delete[] rawBuffer;
        throw ParseException("Can't open contig files"); 
    }
    int lineNo = 1;
    std::string nextLine;
    std::string header;
    std::string sequence;
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
                    if (sequence.empty()) throw ParseException("empty sequence error"); // when something throw
                    header.clear();
                    sequence.clear();
                }
                header = nextLine.erase(0,1);
            }
            else
            {
                std::copy(nextLine.begin(), nextLine.end(),
                          std::back_inserter(sequence));
                
                ContigInfo c1;
                std::vector<std::string> x1 = split(header, ' ');
                
                c1.contigname = x1[0]; // remove larger than sigh
                
                c1.contigSequence = sequence;
                //v1.push_back(c1);
                ContigInfos.insert({c1.contigname, c1});
                sequence.clear();
                header.clear();
            }
            ++lineNo;
            nextLine.clear();
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw ParseException("read in accession file error, please check again.");
    }
    delete[] rawBuffer;
    gzclose(fd);

}
//void loadBlastfile(contigs, AcContigInfo, curbuild.contigPath);

void Read2CleanCon::loadBlastfile(MapHdContigInfo& ContigInfos, /* MapAcVecContigInfo& AcContigInfo */ const std::string& blastpath)
{
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(blastpath.c_str(),"rb");

    if(!fd)
    {
        delete[] rawBuffer;
        throw ParseException("Can't open contig files "); 
    }
    int lineNo = 1;
    std::string nextLine;
    std::string header;
    std::string sequence;
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

            if(nextLine[0] == '#' )
            {
                nextLine.clear();
                lineNo++;
                continue;
            }
            else
            {
                std::vector<std::string> x1 = split(nextLine, '\t'); 
                std::vector<std::string> x2 = split(x1[1], '.'); 
                //ContigInfos.at(x1[0]).haveBlastAc = true;
                ContigInfos.at(x1[0]).accession = x2[0] ;
//                std::vector<ContigInfo>    v1 = {ContigInfos.at(x1[0])}
//                AcContigInfo.upsert(x2[0],
//                    [](std::vector<ContigInfo>& v){v.push_back(ContigInfos.at(x1[0]))}, v1);
                double identity = std::atof(x1[2].c_str());
                double length =std::atof(x1[3].c_str());
                ContigInfos.at(x1[0]).blastIdentity = identity;
                ContigInfos.at(x1[0]).balstLength = length;
                if(identity >= std::atof(Config::get("blastIdentity").c_str()) & length >= std::atof(Config::get("blastLength").c_str()) ){
                    ContigInfos.at(x1[0]).haveBlastAc = true;
                }
            }
            ++lineNo;
            nextLine.clear();
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw ParseException("read in accession file error, please check again."); 
    }
    delete[] rawBuffer;
    gzclose(fd);

}

void Read2CleanCon::writecleanContig(MapHdContigInfo& ContigInfos, std::string& outputname, std::string& outUnmap, const bool isFast)
{
    FILE* fout = fopen(outputname.c_str(), "w");
    if (!fout) throw std::runtime_error("Can't open " + outputname);
    
    FILE* foutUnmap = fopen(outUnmap.c_str(), "w");
    if (!foutUnmap) throw std::runtime_error("Can't open " + outUnmap);
    
    
    Logger::get().debug()<< "Writing clean Contigs";
    for(const auto& Stcontig : ContigInfos)
    {
        bool writeInfo = false;
        if(!isFast){
            writeInfo = Stcontig.second.haveAcinfo & Stcontig.second.haveBlastAc & Stcontig.second.haveClassInfo & Stcontig.second.isnovel ;
        }
        else{
            writeInfo = Stcontig.second.haveBlastAc;
        }
        if(writeInfo){
            std::string header = '>' + Stcontig.second.contigname + "\n";
            std::string sequence = Stcontig.second.contigSequence + "\n";
            fwrite(header.data(), sizeof(header.data()[0]),
                   header.size(),fout);
            fwrite(sequence.data(), sizeof(sequence.data()[0]),
                   sequence.size(),fout);
        }
        else{
            std::string header = '>' + Stcontig.second.contigname + "\n";
            std::string sequence = Stcontig.second.contigSequence + "\n";
            fwrite(header.data(), sizeof(header.data()[0]),
                   header.size(),foutUnmap);
            fwrite(sequence.data(), sizeof(sequence.data()[0]),
                   sequence.size(),foutUnmap);
        }

    }
    fclose(fout);
    fclose(foutUnmap);
}


size_t Read2CleanCon::searchNcbi(std::string& accessionId)
{
//    下面的代码因为网络访问问题而暂时搁置(a song called huahai, for those broken hearts)
//    std::string web = "https://www.ncbi.nlm.nih.gov/nuccore/" + accessionId;
    std::string web = "https://www.mulanci.org/lyric/sl114574/";
    std::string httpText;
    httpText = curlGetHtmlSource(web);
    std::regex const e{".*ORGANISM=(\\d+)&amp.*"};
    std::smatch sm;
    if(std::regex_search(httpText.cbegin(),httpText.cend(),sm,e)){
       // std::cout << sm[1] << "\n";
        std::istringstream iss(sm[1]); // convert string to  size_t and return
        size_t accessionID;
        iss >> accessionID;
        return accessionID;
    }
    else{
        return 0;
    };

}

MapHdBuildMap Read2CleanCon::buildEndmap()
{
    MapHdBuildMap buildedmap;
    for (const auto& buildingrmt : _buildrmtfiles.lock_table())
    {
        std::string filename = buildingrmt.second.firstStageBuildfile.filename;
        std::string outpath = _outpath + "/pan_" + filename + "/";
        std::string cleanContig = buildingrmt.second.cleanContig;
        buildfiles bdf = buildingrmt.second.firstStageBuildfile;
        buildMapinfo buildmap(cleanContig, bdf, outpath, filename); 
        buildmap.buildingmap();
        buildedmap.insert(cleanContig, buildmap);
    }
    return buildedmap;
}
