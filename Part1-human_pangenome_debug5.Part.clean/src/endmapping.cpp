//
//  endmapping.cpp
//  human_pangenome
//
//  Created by kongshuang on 03/07/2022.
//

//typedef std::unordered_map<std::string, ContigInfo> MapHdContigInfo;

#include "endmapping.h"

void Endmapping::buildMapping()
{
    _buildmap = _rd2con.buildEndmap();
}
// clean contig 
void Endmapping::Endmap()
{
    std::vector<std::string> EndmapContigs;
    for(auto i : _buildmap.lock_table()){
        EndmapContigs.push_back(i.first);
    }
    int32_t chrRange = std::atoi(Config::get("ChrMappingRange").c_str());
    int32_t contigTailRange = std::atoi(Config::get("ContigTail").c_str());
    int32_t mapContigTail = std::atoi(Config::get("MapContigTail").c_str());
    
    std::function<void(const std::string& ContigPath)> Endmapping =
    [this, chrRange, contigTailRange, mapContigTail](const std::string& ContigPath){
        buildMapinfo curBuild = this->_buildmap.find(ContigPath);
        
        auto Headers = getFaHead(curBuild.cleanContig);
        Logger::get().info() <<"curBuild.matealignedreads.first is " << curBuild.matealignedreads.first << "\n";
        Logger::get().info() <<"curBuild.matealignedreads.second is " << curBuild.matealignedreads.second << "\n";
        Logger::get().info() <<"Linkage sam is " << curBuild.linkageSam << "\n";
 
        
        
        std::string IniEnv = Config::get("IniEnv");
        uint32_t state1 = system((IniEnv + "bowtie2-build " + curBuild.cleanContig + " " + curBuild.contigIndex ).c_str());
        
        
        uint32_t state2 = system((IniEnv + "bowtie2 -x " + curBuild.contigIndex + " -U " + \
                                  curBuild.matealignedreads.first + " -S " + \
                                  curBuild.ReadToContigSam ).c_str());
       
    	Logger::get().info() <<"curBuild.matealignedreads.first is " << curBuild.matealignedreads.first << "\n";
    	Logger::get().info() <<"curBuild.matealignedreads.second is " << curBuild.matealignedreads.second << "\n";


        std::string bashCmd1 = (IniEnv + "samtools view -h -F 2304 " + curBuild.ReadToContigSam + \
                            " | samtools sort -n -O bam | bedtools bamtobed -i stdin "
                            " | awk '{OFS=\"\\t\"} {print $4,$1,$6,$2,$3}'"
                            " | sed -e 's/\\/[1-2]//g' "
                            " | awk 'NR==FNR{a[$1]=$2;next}{if(a[$2]!=\"\"){print $0\"\\t\"a[$2]}}' <(samtools view -H " + curBuild.ReadToContigSam + \
                            " | " + " awk -F '[\\t:]' '{if($1==\"@SQ\")print $3,$5}') - | sort > " + curBuild.ReadToContigTxt);
        //std::cout << test << "\n";
            
        
        std::string bashHead = ("echo #!/bin/bash\n");

        
        {
        FILE* fout = fopen((curBuild.ReadToContigTxt + ".sh" ).c_str(), "w");
        fwrite(bashHead.data(), sizeof(bashHead.data()[0]),bashHead.size(),fout);
        fwrite(bashCmd1.data(), sizeof(bashCmd1.data()[0]),bashCmd1.size(),fout);
        fclose(fout);
        }
        /*
        uint32_t state3 = system((  IniEnv + "samtools view -h -F 2304 " + curBuild.ReadToContigSam + \
                                  " | samtools sort -n -O bam | bedtools bamtobed -i stdin " \
                                  " | awk '{OFS=\"\\t\"} {print $4,$1,$6,$2,$3}'" \
                                  " | sed -e 's/\\/[1-2]//g' " \
                                  " | awk 'NR==FNR{a[$1]=$2;next}{if(a[$2]!=\"\"){print $0\"\\t\"a[$2]}}' <(samtools view -H " + curBuild.ReadToContigSam + \
                                  " | " + " awk -F '[\\t:]' '{if($1==\"@SQ\")print $3,$5}') - | sort > " + curBuild.ReadToContigTxt).c_str());
         */
        uint32_t state3 = system(("bash " + curBuild.ReadToContigTxt + ".sh").c_str());
        Logger::get().info() <<"Linkage sam is " << curBuild.linkageSam << "\n";
        std::string bashCmd2 = (  IniEnv + " samtools view -H  " + curBuild.linkageSam + \
                             " | cat - <(awk 'NR==FNR{ a[$1]; next }$1 in a{ print $0 ; delete a[$1]; next }' " + curBuild.ReadToContigTxt + \
                             " <( samtools view " + curBuild.linkageSam + " ))" + \
                             " | samtools sort -n -O bam "
                             " | bedtools bamtobed -i stdin "
                             " | awk '{OFS=\"\\t\"}{print $4,$1,$6,$2,$3}' "
                             " | sed -e \"s/\\/[1-2]//g\" | sort > " + curBuild.passMate );
        {
            FILE* fout = fopen((curBuild.passMate + ".sh" ).c_str(), "w");
            fwrite(bashHead.data(), sizeof(bashHead.data()[0]),bashHead.size(),fout);
            fwrite(bashCmd2.data(), sizeof(bashCmd2.data()[0]),bashCmd2.size(),fout);
            fclose(fout);
        }
        uint32_t state4 = system(("bash " + curBuild.passMate + ".sh").c_str());
        /*
        uint32_t state4 = system((  IniEnv + " samtools view -H  " + curBuild.linkageSam + \
                                  " | cat - <(awk 'NR==FNR{ a[$1]; next }$1 in a{ print $0 ; delete a[$1]; next }' " + curBuild.ReadToContigTxt + \
                                  " <( samtools view " + curBuild.linkageSam + " ))" + \
                                  " | samtools sort -n -O bam "
                                  " | bedtools bamtobed -i stdin "
                                  " | awk '{OFS=\"\\t\"}{print $4,$1,$6,$2,$3}' "
                                  " | sed -e \"s/\\/[1-2]//g\" | sort > " + curBuild.passMate ).c_str());
         */
        uint32_t state5 = system(("join -j 1 " + curBuild.ReadToContigTxt + " " + curBuild.passMate + " > " + curBuild.mateRegion).c_str());
        //std::string sortMateRegion = curBuild.mateRegion + "Sort";
        system(("sort -V -k2 " + curBuild.mateRegion + " > " + curBuild.sortMateRegion).c_str());
        
//        uint32_t judge = state1 + state2 + state3 + state4 + state5 ;
        //panGenomeDebug2
        //uint32_t judge = 0;
        std::string outPutRegion = curBuild.output + "contigRegion/";
        fs::create_directory(outPutRegion);
	    int judge =0 ; 
        if(judge == 0){
            MpContigMapinfo curMpConMapinfo; // unordered map ; contig name
            curMpConMapinfo = this->contigMapping(curBuild.sortMateRegion,contigTailRange,chrRange, mapContigTail, outPutRegion, curBuild.cleanContig);

            std::vector<std::string> outStreamLeft;
            std::vector<std::string> outStreamRight;
            std::vector<std::string> outStreamBoth;
            std::vector<std::string> outStreamUnmap;
            
            for(auto& contigs : Headers)
            {
                if(curMpConMapinfo.count(contigs) != 0)
                {
                    for(auto& i : curMpConMapinfo)
                    {
                        if(i.first == contigs)
                        {
                            if(i.second.leftMatch && i.second.rightMatch)
                            {
                				Logger::get().info()<< "line 133 contig " << i.first << "\n";
                				Logger::get().info()<< "line 133 contig " << i.first << "\n";
                                if(i.second.leftMatchInfo.Reverse == i.second.rightMatchInfo.Reverse && (i.second.leftMatchInfo.refName == i.second.rightMatchInfo.refName))
                                {
                                    if(i.second.leftMatchInfo.Reverse)
                                    {
                                        int leftMapSize = i.second.leftInsertSitua.first.size();
                                        int rightMapSize = i.second.rightInsertSitua.first.size();
                    					Logger::get().info()<< "contig " << i.first << " and leftMapSize" << leftMapSize << "\n";
                    					Logger::get().info()<< "contig " << i.first << " and rightMapSize" << rightMapSize << "\n";
                                        int distance = INT_MAX;
                                        int suitaLeftInsertion=0;
                                        int suitaRightInsertion=0;
                                        int leftContigIns = 0;
                                        int rightContigIns = 0;
                                        for(int leftMap=0; leftMap<leftMapSize; leftMap++){
                    					    for(int rightMap =0; rightMap<rightMapSize; rightMap++){
                        						Logger::get().info()<< " we have left insertion "<< i.second.leftInsertSitua.first[leftMap] << "\n";
                        						Logger::get().info()<< " we have right insertion "<< i.second.rightInsertSitua.first[rightMap] << "\n";
                                                int curDist = i.second.leftInsertSitua.first[leftMap]-i.second.rightInsertSitua.first[rightMap];
                                                if(curDist>=0 && curDist<=distance && curDist<=1000){
                                                    distance = curDist;
                                                    suitaLeftInsertion = i.second.leftInsertSitua.first[leftMap];
                                                    suitaRightInsertion = i.second.rightInsertSitua.first[rightMap];
                                                    leftContigIns = i.second.leftInsertSitua.second[leftMap];
                                                    rightContigIns = i.second.rightInsertSitua.second[rightMap];
                                                }
                                            }
                                        }
					 
                                        if(suitaRightInsertion != 0)
                                        {
                                            i.second.BothMatch = true;
                                            std::string outLine = i.first + "\t" + i.second.leftMatchInfo.refName + "\t" +
                                            std::to_string(suitaLeftInsertion) + "\t" +
                                            std::to_string(leftContigIns) + "\t" +
                                            std::to_string(suitaRightInsertion) + "\t" +
                                            std::to_string(rightContigIns) + "\t" +
                                            std::to_string(i.second.leftMatchInfo.Reverse) + "\t" +
                                            std::to_string(i.second.contigLength) + "\n";
                                            outStreamBoth.push_back(outLine);
                                        }
                                    }
                                    else
                                    {
                                        int leftMapSize = i.second.leftInsertSitua.first.size();
                                        int rightMapSize = i.second.rightInsertSitua.first.size();
                                        int distance = INT_MAX;
                                        int suitaLeftInsertion=0;
                                        int suitaRightInsertion=0;
					int leftContigIns = 0;
					int rightContigIns = 0;
                                        for(int leftMap=0; leftMap<leftMapSize; leftMap++){
                                            for(int rightMap =0; rightMap<rightMapSize; rightMap++){
                                                int curDist = i.second.rightInsertSitua.first[rightMap] - i.second.leftInsertSitua.first[leftMap];
                                                if(curDist>=0 && curDist<=distance && curDist<=1000 ){
                                                    distance = curDist;
                                                    suitaLeftInsertion = i.second.leftInsertSitua.first[leftMap];
                                                    suitaRightInsertion = i.second.rightInsertSitua.first[rightMap];
                                                    leftContigIns = i.second.leftInsertSitua.second[leftMap];
                                                    rightContigIns = i.second.rightInsertSitua.second[rightMap];
                                                }
                                            }
                                        }
                                        if(suitaRightInsertion != 0)
                                        {
                                            i.second.BothMatch = true;
                                            std::string outLine = i.first + "\t" + i.second.leftMatchInfo.refName + "\t" +
                                                std::to_string(suitaLeftInsertion) + "\t" +
                                                std::to_string(leftContigIns) + "\t" + 
                                                std::to_string(suitaRightInsertion) + "\t" +
                                                std::to_string(rightContigIns) + "\t" + 
                                                std::to_string(i.second.leftMatchInfo.Reverse) + "\t" +
                                                std::to_string(i.second.contigLength) + "\n";
                                            outStreamBoth.push_back(outLine);
                                        }
                                    }
                                }
                            }
                            else if(i.second.leftMatch)
                            {

                                std::string outLine = i.first + "\t" + i.second.leftMatchInfo.refName + "\t" +
                                    std::to_string(i.second.leftMatchInfo.InsertPoint) + "\t" +
                                    std::to_string(i.second.leftMatchInfo.Reverse) + "\t" +
                                    std::to_string(i.second.contigLength) + "\t" +
                                    std::to_string(CHRLENGTHMAP.at(i.second.leftMatchInfo.refName)) + "\t" +
                                    std::to_string(i.second.leftMatchInfo.RefSTART + i.second.leftMatchInfo.RefStart) + "\t" +
                                    std::to_string(i.second.leftMatchInfo.RefSTART + i.second.leftMatchInfo.RefEnd) + "\t" +
                                    std::to_string(i.second.leftMatchInfo.contigMatchIns) + "\n";
                                outStreamLeft.push_back(outLine);

                            }
                            else if(i.second.rightMatch)
                            {
                                std::string outLine = i.first + "\t" + i.second.rightMatchInfo.refName + "\t" +
                                    std::to_string(i.second.rightMatchInfo.InsertPoint) + "\t" +
                                    std::to_string(i.second.rightMatchInfo.Reverse) + "\t" +
                                    std::to_string(i.second.contigLength) + "\t" +
                                    std::to_string(CHRLENGTHMAP.at(i.second.rightMatchInfo.refName)) + "\t" +
                                    std::to_string(i.second.rightMatchInfo.RefSTART + i.second.rightMatchInfo.RefStart) + "\t" +
                                    std::to_string(i.second.rightMatchInfo.RefSTART + i.second.rightMatchInfo.RefEnd) + "\t" +
                                    std::to_string(i.second.rightMatchInfo.contigMatchIns) + "\n";
                                outStreamRight.push_back(outLine);
                            }
                            else
                            {
                                std::string outLine = i.first + "\n";
                                outStreamUnmap.push_back(outLine);
                            }
                        }
                    }
                }else
                {
                    std::string outLine = contigs + "\n";
                    outStreamUnmap.push_back(outLine);
                }
            }
                
        
            std::string outLeft = outPutRegion + "leftEndMapping" ;
            std::string outRight = outPutRegion + "rightEndMapping" ;
            std::string outBoth = outPutRegion + "BothEndMapping" ;
            std::string outUnmap = outPutRegion + "unmappedContig";
            
            writeoutput(outStreamLeft, outLeft);
            writeoutput(outStreamRight, outRight);
            writeoutput(outStreamBoth, outBoth);
            writeoutput(outStreamUnmap, outUnmap);

        }
        else{
            Logger::get().warning() << ContigPath << " build doesn't success, please check your pipeline related files" << "\n";
        };
        
          
    };
    size_t threadnum = (size_t)std::atoi(Config::get("blastThread").c_str());
    processInParallel(EndmapContigs ,Endmapping , threadnum, true);
    
}


MpContigMapinfo Endmapping::contigMapping(std::string &mateRegionTxt, int32_t contigTailRange,int32_t chrRange, int32_t mapContigTail, std::string& output, std::string& cleanContig){
    MpContigMapinfo contigsinfo; // contig head
    //std::string accessionfile = Config::get("accession");
    
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(mateRegionTxt.c_str(),"rb");
    if(!fd)
    {
        delete[] rawBuffer;
        throw ParseException("Can't open mate RegionTxt file ");
    }
    int lineNo = 1;
    std::string nextLine;
    std::string contigName;
    int32_t length;
    std::vector<ContigLine> lineRecord;
    try
    {
        while(!gzeof(fd))
        {
            for(;;)
            {  // gzget 
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
                std::vector<std::string> x = split(nextLine, ' ');
                if(lineNo == 1){
                    length = std::atoi(x[5].c_str());
                    if(std::atoi(x[4].c_str()) > (length - contigTailRange ) && check_key(CHRLENGTHMAP, x[6])){
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7], true);
                        lineRecord.push_back(CtL);
                        contigName = x[1];
                    }
                    else if(std::atoi(x[3].c_str()) < contigTailRange && check_key(CHRLENGTHMAP, x[6]) ){
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7], false);
                        lineRecord.push_back(CtL);
                        contigName = x[1];
                    }
                    else{
                        contigName = x[1];
                    }
                }
                else if (x[1] != contigName && lineRecord.size() != 0){
                    ContigMapInfo curContigInfo(contigName, length, lineRecord, output);
                    curContigInfo.ConstructRegion(chrRange, mapContigTail, cleanContig);
                    contigsinfo.insert({contigName,curContigInfo});
                    lineRecord.clear();
                    if(std::atoi(x[4].c_str()) > (length - contigTailRange) && check_key(CHRLENGTHMAP, x[6]))
                    {
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7],true);
                        lineRecord.push_back(CtL);
                    }
                    else if (std::atoi(x[3].c_str()) < contigTailRange && check_key(CHRLENGTHMAP, x[6])){
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7],false);
                        lineRecord.push_back(CtL);
                    }
                    length = std::atoi(x[5].c_str());
                    contigName = x[1];
                }
                else if(x[1] == contigName){
                    if(std::atoi(x[4].c_str()) > (length - contigTailRange) && check_key(CHRLENGTHMAP, x[6]))
                    {
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7],true);
                        lineRecord.push_back(CtL);
                    } // contigName and length not going to change
                    else if (std::atoi(x[3].c_str()) < contigTailRange && check_key(CHRLENGTHMAP, x[6]) ){
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7],false);
                        lineRecord.push_back(CtL);
                    }
                }
                else{ //  contig 
                    ContigMapInfo curContigInfo(contigName, length, lineRecord, output);
                    contigsinfo.insert({contigName,curContigInfo});
                    lineRecord.clear();
                    length = std::atoi(x[5].c_str());
                    contigName = x[1];
                    if(std::atoi(x[4].c_str()) > (length - contigTailRange) && check_key(CHRLENGTHMAP, x[6]))
                    {
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7],true);
                        lineRecord.push_back(CtL);
                    } // contigName and length not going to change
                    else if (std::atoi(x[3].c_str()) < contigTailRange && check_key(CHRLENGTHMAP, x[6]) ){
                        ContigLine CtL(x[1], x[6], std::atoi(x[8].c_str()), std::atoi(x[9].c_str()),x[0],
                                       std::atoi(x[3].c_str()), std::atoi(x[4].c_str()), x[2],x[7],false);
                        lineRecord.push_back(CtL);
                    }
                    
                }
                
            }
            ++lineNo;
            nextLine.clear();
        }
        if(lineRecord.size() != 0){
            ContigMapInfo curContigInfo(contigName,length, lineRecord, output);
            curContigInfo.ConstructRegion(chrRange, mapContigTail, cleanContig);
            contigsinfo.insert({contigName,curContigInfo});
        }
    }
    catch(...)
    {
        gzclose(fd);
        throw ParseException(("error. contigmapping " + mateRegionTxt ).c_str());
    }
    delete[] rawBuffer;
    gzclose(fd);
    
    return contigsinfo;
}



void ContigMapInfo::ConstructRegion(uint32_t chrRegion, uint32_t contigRegion, std::string& cleanContig){
    //cuckoohash_map<std::string,int> chrCount;

    
    std::string maxChrLeft;
    std::string maxChrRight;
    cuckoohash_map<std::string, int> leftMap;
    cuckoohash_map<std::string, int> rightMap;
    std::vector<int> leftChrStart;
    std::vector<int> rightChrStart;
    
    for(const auto& i : contigLinkMap){
        if(i.isRight == true){
            rightMap.upsert(i.chr,[](int& num){num++;},1);
        }
        else{
            leftMap.upsert(i.chr, [](int& num){num++;},1);
        }
    };
    std::string nucmerL = Config::get("nucmerL");
    std::string nucmerB = Config::get("nucmerB");
    std::string nucmerC = Config::get("nucmerC");
    
    std::string inienv = Config::get("IniEnv");
    
    auto maxLeft = getLargestChr(leftMap); // chr match。
    if(std::get<1>(maxLeft)){
        maxChrLeft = std::get<0>(maxLeft);
        float sum = std::get<2>(maxLeft);
        float max = std::get<3>(maxLeft);

        for(const auto& i : contigLinkMap){
            if(i.chr == maxChrLeft && i.isRight==false){
                leftChrStart.push_back(i.chrStart);
            }
        }
        uint32_t contigStart = 1 ;
        uint32_t contigEnd = contigRegion > contigLength ? contigLength : contigRegion;

        float clustnum = getClust(leftChrStart);
        bool judgeNum = true ? (clustnum/sum >= 0.95) : false;
       	
	Logger::get().info()<< "left contig " << ContigName << " get; sum is  "<<  sum << " judgeNum is " << judgeNum << "clusteNum is "<< clustnum<<"\n"; 
        if((*max_element(leftChrStart.begin(),leftChrStart.end()) - *min_element(leftChrStart.begin(),leftChrStart.end()))<= chrRegion && judgeNum && leftChrStart.size()>=2){ // 2000
            exLeftRegion = true;
            int LeftmedianNum = median(leftChrStart);
            int LeftLeftMost = (LeftmedianNum - 1500 < 0) ? 1 : LeftmedianNum - 1499;
            int LeftRightMost = (LeftmedianNum + 1500 > CHRLENGTHMAP.at(maxChrLeft)) ? CHRLENGTHMAP.at(maxChrLeft) : LeftmedianNum + 1500;
            std::string G38path = Config::get("Grch38Noalt");
            std::string outLeftChr = outPath + "LeftChrRegion" + ContigName + ".fa";
            
            system((inienv + " samtools faidx " + G38path + " " + maxChrLeft + ":" + std::to_string(LeftLeftMost) + "-" + std::to_string(LeftRightMost) + " > " + outLeftChr).c_str()); // output 
            
            std::string outLeftContig = outPath + "LeftContigRegion" + ContigName + ".fa";
            system((inienv + " samtools faidx " + cleanContig + " " + ContigName + ":" + std::to_string(contigStart) + "-" + std::to_string(contigEnd) + " > " + outLeftContig).c_str()); // output 
            leftRegion = {outLeftChr, outLeftContig};
            leftDelta = outPath + "LeftDelta" + ContigName + ".Delta";
            uint32_t state1 = system((inienv + " nucmer --maxmatch -l " + nucmerL + " -b " + nucmerB + " -c " + nucmerC + " -p " + leftDelta + " " +
                    outLeftChr + " " + outLeftContig ).c_str());
            filtedLeftDelta = outPath + "filtedLeftDelta" + ContigName + ".delta";
            uint32_t state2 = system((inienv + " delta-filter -q -r -o 0 -g " + leftDelta + ".delta" + " > " + filtedLeftDelta).c_str());
            if(true){
                leftInsertSitua = getRegionMatch(false);
            }
            else{
                Logger::get().error() <<  "ERROR: delta or nucmer wrong for " << ContigName << "\n";
            }//
            
        }
    };
    
    auto maxRight = getLargestChr(rightMap); // chr match
    if(std::get<1>(maxRight)){
        maxChrRight = std::get<0>(maxRight);
        float sum = std::get<2>(maxRight);
        float max = std::get<3>(maxRight);

        for(const auto& i : contigLinkMap){
            if(i.chr == maxChrRight && i.isRight==true){
                rightChrStart.push_back(i.chrStart);
            }
        }
        
        uint32_t contigEnd = contigLength;
        uint32_t contigStart = contigLength - contigRegion > 0 ?  (contigLength-contigRegion + 1) : 0;

        float clustnum = getClust(rightChrStart);
        bool judgeNum = true ? (clustnum/sum >= 0.95) : false;

         Logger::get().info()<< "right contig " <<ContigName <<" get; sum is  "<<  sum << " judgeNum is " << judgeNum << "clusteNum is "<< clustnum<<"\n";
        if((*max_element(rightChrStart.begin(),rightChrStart.end()) - *min_element(rightChrStart.begin(),rightChrStart.end()))<= chrRegion && judgeNum && rightChrStart.size() >= 2){ 
            exRightRegion = true;
            int RightmedianNum = median(rightChrStart);
            int RightLeftMost = (RightmedianNum - 1500 < 0) ? 0 : RightmedianNum - 1499;
            int RightRightMost = (RightmedianNum + 1500 > CHRLENGTHMAP.at(maxChrRight)) ? CHRLENGTHMAP.at(maxChrRight) : RightmedianNum + 1500;
            std::string G38path = Config::get("Grch38Noalt");
            std::string outRightChr = outPath + "RightChrRegion" + ContigName + ".fa";
            
            system((inienv + " samtools faidx " + G38path + " " + maxChrRight + ":" + std::to_string(RightLeftMost) + "-" + std::to_string(RightRightMost) + " > " + outRightChr).c_str()); 
            
            std::string outRightContig = outPath + "RightContigRegion" + ContigName + ".fa";
            system((inienv + " samtools faidx " + cleanContig + " " + ContigName + ":" + std::to_string(contigStart) + "-" + std::to_string(contigEnd) + " > " + outRightContig).c_str()); 
            leftRegion = {outRightChr, outRightContig};
            
            rightDelta = outPath + "RightDelta" + ContigName + ".Delta";
            
            std::string test  = (inienv + " nucmer --maxmatch -l " + nucmerL + " -b " + nucmerB + " -c " + nucmerC + " -p " + rightDelta + " " + \
                                 outRightChr + " " + outRightContig );
 //           std::cout << test << "\n";
            
            uint32_t state1 = system((inienv + " nucmer --maxmatch -l " + nucmerL + " -b " + nucmerB + " -c " + nucmerC + " -p " + rightDelta + " " + \
                    outRightChr + " " + outRightContig ).c_str());
            filtedRightDelta = outPath + "filtedRightDelta" + ContigName + ".delta";
            
            std::string test_filter = (inienv + " delta-filter -q -r -o 0 -g " + rightDelta + ".delta" + " > " + filtedRightDelta);
            
            uint32_t state2 = system((inienv + " delta-filter -q -r -o 0 -g " + rightDelta + ".delta" + " > " + filtedRightDelta).c_str());
            if(true){
                rightInsertSitua  = getRegionMatch(true);
            }
            else{
                Logger::get().error() <<  "ERROR: delta or nucmer wrong for " << ContigName << "\n";
            }//
        }


    };


    

}



std::pair<std::vector<int>, std::vector<int>> ContigMapInfo::getRegionMatch(const bool& Right){
    std::string deltaFile;
    if(Right){
         deltaFile = filtedRightDelta;
    }
    else{
         deltaFile = filtedLeftDelta;
    }
    size_t BUF_SIZE = 32 * 1024 * 1024;
    char* rawBuffer = new char[BUF_SIZE];
    auto* fd = gzopen(deltaFile.c_str(),"rb");
    if(!fd)
    {
        delete[] rawBuffer;
        throw std::runtime_error("ERROR: can't open delta file "); 
    }
    int lineNo = 1;
    std::string nextLine;
    int refStart = 0;
    int refEnd = 0;
    int contigMatchStart = 0;
    int contigMatchEnd = 0;
    std::string ref;
    deltaMatch deltaMatchRef;
    deltaMatch deltaMatchContig;
    std::vector<bool> reverse;
    std::vector<int> insertions;
    std::vector<int> seqInsertions;
    try
    {
        while(!gzeof(fd))
        {
            for(;;)
            {  //
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
            if(nextLine[0] == '#' || lineNo == 1 || lineNo ==2 )
            {
                nextLine.clear();
                lineNo++;
                continue;
            }
            else
            {
                std::vector<std::string> x = split(nextLine, ' ');
                if(x.size() == 4){
                    std::vector<std::string> x1 = split(x[0], ':');
                    ref = x1[0].erase(0,1);
                    std::vector<std::string> x2 = split(x1[1], '-');
                    refStart = std::atoi(x2[0].c_str());
                    refEnd = std::atoi(x2[1].c_str());
                    std::vector<std::string> x3 = split(x[1], ':');
                    std::vector<std::string> xContig = split(x3[1], '-');
                    contigMatchStart = std::atoi(xContig[0].c_str());
                    contigMatchEnd = std::atoi(xContig[1].c_str());
                }
                else if(x.size() == 7)
                {
                    //rightMatch.push
                    int endDist;
                    int EndDist = std::atoi(Config::get("EndDist").c_str());
                    if(Right)
                    {
                        endDist = 150 - std::max(std::atoi(x[2].c_str()),std::atoi(x[3].c_str()));
                    }
                    else
                    {
                        endDist = std::min(std::atoi(x[2].c_str()),std::atoi(x[3].c_str()));
                    }
                    if(endDist <= EndDist)
                    {
                        std::pair<int,int> pRef = {std::atoi(x[0].c_str()),std::atoi(x[1].c_str())};
                        std::pair<int,int> pContig = {std::atoi(x[2].c_str()),std::atoi(x[3].c_str())};
                        deltaMatchRef.push_back(pRef);
                        deltaMatchContig.push_back(pContig);
                        if(std::atoi(x[2].c_str())<std::atoi(x[3].c_str()))
                        {
                            reverse.push_back(false);
                        }
                        else{
                            reverse.push_back(true);
                        }
                        
                    }
                    
                }
            }
            ++lineNo;
            nextLine.clear();
        }
        std::set<int> s;
        unsigned size = reverse.size();
        for( unsigned i = 0; i < size; ++i ) s.insert( reverse[i] );

        int vecSize = reverse.size();
        Logger::get().info() << "reverse vector size is "<< vecSize << "\n";
        if(reverse.size() != 0 && s.size()==1){
            if(!Right){
                leftMatch = true; 

                if(!reverse[0]){
                    Logger::get().info()<< "checkpoin1"<<"\n";
                    int minIndex = getExtreamDeltaMatch(deltaMatchRef,true,true);
                    Logger::get().info()<< "checkpoin2"<<"\n";
                    matchInfo leftMatch((refStart + deltaMatchRef[minIndex].second - 1 ),false);
                    leftMatch.setRef((deltaMatchRef[minIndex].second - deltaMatchRef[minIndex].first), deltaMatchRef[minIndex].first, deltaMatchRef[minIndex].second);
                    leftMatch.setRef2(ref, refStart, refEnd);
                    leftMatch.setContig((deltaMatchContig[minIndex].second + 1), deltaMatchContig[minIndex].first, deltaMatchContig[minIndex].second);
                    leftMatchInfo =leftMatch;
                    for(int i=0; i<vecSize; i++){
                        insertions.push_back((refStart + deltaMatchRef[i].second - 1));
                        seqInsertions.push_back((deltaMatchContig[minIndex].second + 1));
                        Logger::get().info()<<" left no_reverse insertion point is "<< (refStart + deltaMatchRef[i].second - 1 )<< "\n";
                    }
                }
                else{

                    int maxIndex = getExtreamDeltaMatch(deltaMatchRef,false,false);
                    matchInfo leftMatch((refStart + deltaMatchRef[maxIndex].first -1),true);
                    leftMatch.setRef((deltaMatchRef[maxIndex].second - deltaMatchRef[maxIndex].first), deltaMatchRef[maxIndex].first, deltaMatchRef[maxIndex].second);
                    leftMatch.setRef2(ref, refStart, refEnd);
                    leftMatch.setContig((deltaMatchContig[maxIndex].first +1 ), deltaMatchContig[maxIndex].first, deltaMatchContig[maxIndex].second);
                    leftMatchInfo = leftMatch;
                    for(int i=0; i<vecSize; i++){
                        Logger::get().info()<< "left reverse insertion point get " << (refStart + deltaMatchRef[i].first -1) << "\n";
                        insertions.push_back((refStart + deltaMatchRef[i].first -1));
                        seqInsertions.push_back((deltaMatchContig[maxIndex].first +1 ));
                    }
                }
            }
            else {
                rightMatch = true;
                if(!reverse[0]){
                    
                    int maxIndex = getExtreamDeltaMatch(deltaMatchRef,false,false);
                    matchInfo rightMatch((refStart + deltaMatchRef[maxIndex].first -1), false);
                    rightMatch.setContig((contigMatchStart + deltaMatchContig[maxIndex].first -2), deltaMatchContig[maxIndex].first, deltaMatchContig[maxIndex].second);
                    rightMatch.setRef((deltaMatchRef[maxIndex].second-deltaMatchRef[maxIndex].first), deltaMatchRef[maxIndex].first, deltaMatchRef[maxIndex].second);
                    rightMatch.setRef2(ref, refStart, refEnd);
                    rightMatchInfo = rightMatch ;
                    for(int i=0; i<vecSize; i++){
                        insertions.push_back((refStart + deltaMatchRef[i].first -1));
                        seqInsertions.push_back((contigMatchStart + deltaMatchContig[maxIndex].first -2));
			Logger::get().info()<<" right no_reverse insertion point is "<< (refStart + deltaMatchRef[i].first - 1 )<< "\n";		
                    }
                }
                else{
                    int minIndex = getExtreamDeltaMatch(deltaMatchRef, true, true);
                    matchInfo rightMatch((refStart + deltaMatchRef[minIndex].second -1 ),true);
                    rightMatch.setContig((contigMatchStart + deltaMatchContig[minIndex].second -2 ), deltaMatchContig[minIndex].first, deltaMatchContig[minIndex].second);
                    rightMatch.setRef((deltaMatchRef[minIndex].second-deltaMatchRef[minIndex].first), deltaMatchRef[minIndex].first, deltaMatchRef[minIndex].second);
                    rightMatch.setRef2(ref, refStart, refEnd);
                    rightMatchInfo = rightMatch ;
                    for(int i=0; i<vecSize; i++){
                        insertions.push_back((refStart + deltaMatchRef[i].second - 1 ));
                        seqInsertions.push_back((contigMatchStart + deltaMatchContig[minIndex].second -2 ));
			Logger::get().info()<<" right reverse insertion point is "<< (refStart + deltaMatchRef[i].second - 1 )<< "\n";

                    }
                };
            }
        }

    }
    catch(...)
    {
        gzclose(fd);
        throw std::runtime_error(("error get region match." + deltaFile ).c_str());
    }
    std::pair<std::vector<int>, std::vector<int>> insertSituation = {insertions, seqInsertions};
    delete[] rawBuffer;
    gzclose(fd);
    return insertSituation;
}

void Endmapping::writeoutput(std::vector<std::string>& v1,std::string& outPath){
    FILE* fout = fopen(outPath.c_str(), "w");
    if (!fout) throw std::runtime_error("Can't open " + outPath);

    Logger::get().debug()<< "Writing Endmapping result";
    for(auto& i : v1){
        fwrite(i.data(), sizeof(i.data()[0]),
               i.size(),fout);
    }
    fclose(fout);
}



