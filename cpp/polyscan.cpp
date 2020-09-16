#include <iostream>
#include <sstream>
#include <bitset>
#include  <map>
#include <unordered_map>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <numeric>
#include <cmath>
#include "utilities.h"
#include "polyscan.h"
#include "bamreader.h"
#include "param.h"
#include "sample.h"
#include "distribution.h"
//#include "train.h"

extern Param paramd;
extern std::map<std::string, int> SitesSupport;
extern bit8_t alphabet[];
extern bit8_t rev_alphabet[];
extern char uhomo_code[];
extern char homo_code[];
extern Sample sample;
//extern Train train;

PolyScan::PolyScan() {
    homosBuffer.reserve(paramd.bufSize);
    totalSites.reserve(paramd.bufSize);
}

PolyScan::~PolyScan() {
    totalSites.clear();
}

// eliminates a character from the input string
void PolyScan::eliminate(const char ch, std::string &str) {
    size_t eliminateCharPos = str.find(ch);
    while (eliminateCharPos != std::string::npos) {
        str.erase(eliminateCharPos, 1);
        eliminateCharPos = str.find(ch);
    }
}

// Parse one region 
bool PolyScan::ParseOneRegion(const std::string &regionString) {
    size_t separatorPos = regionString.find(":");
    bool correctParse = false;
    bool m_endDefined = false;
    bool m_startDefined = false;
    int m_start = -1;
    int m_end = -1;
    std::string m_targetChromosomeName;
    // found a separator
    if (separatorPos != std::string::npos) {
        m_targetChromosomeName = regionString.substr(0, separatorPos);
        std::string coordinates = regionString.substr(separatorPos + 1);
        // removes the ',' in 1,000 or 1,000,000 that users may add
        // for readability but wreak havoc with atoi
        eliminate(',', coordinates);
        size_t startEndSeparatorPos = coordinates.find("-");
        // there are two coordinates
        if (startEndSeparatorPos != std::string::npos) {
            std::string secondPositionStr = coordinates.substr(
                    startEndSeparatorPos + 1);
            m_end = atoi(secondPositionStr.c_str());
            m_endDefined = true;
        }

        m_start = atoi(coordinates.c_str());
        m_startDefined = true;
        if (m_start < 0 || (m_endDefined && (m_end < m_start))) {
            correctParse = false;
        } else {
            correctParse = true;
        }
    }
        // no separator found
    else {
        m_targetChromosomeName = regionString;
        correctParse = true;
    }
    // assign values
    region_one.chr = m_targetChromosomeName;
    region_one.start = m_start;
    region_one.end = m_end;

    return correctParse;
}

// Loading bed regions
void PolyScan::LoadBeds(std::ifstream &fin) {
    std::string chr;
    std::string line;
    std::string tempChr = "";
    int start;
    int stop;
    int i = -1;
    while (getline(fin, line)) {
        std::stringstream linestream(line);
        linestream >> chr;
        linestream >> start;
        linestream >> stop;
        /*
         std::cout<< chr <<"\t"
         << start <<"\t"SitesSupport
         << stop << "\n";
         */
        if (chr == tempChr) {
            BedRegion tempBedRegion;
            tempBedRegion.start = start;
            tempBedRegion.end = stop;
            beds[i].regions_list.push_back(tempBedRegion);
        } else {
            ++i;
            BedChr tempBedChr;
            tempBedChr.chr = chr;
            beds.push_back(tempBedChr);
            // load mapping
            chrMaptoIndex.insert(std::pair<std::string, bit16_t>(chr, i));
            BedRegion tempBedRegion;
            tempBedRegion.start = start;
            tempBedRegion.end = stop;
            beds[i].regions_list.push_back(tempBedRegion);
            tempChr = chr;
        }
        linestream.clear();
        linestream.str("");
    }
}

// loading bam list
// only load one bam file
void PolyScan::LoadBams(const std::string &bam1, const std::string &bam2) {
    BamPairs t_bampair;
    t_bampair.sName = "sample_name";

    if (bam1.find(".bam") != std::string::npos) {
        t_bampair.normal_bam = bam1;
    } else {
        std::cerr << "please provide valid format normal bam file ! \n";
        exit(0);
    }
    if (bam2.find(".bam") != std::string::npos) {
        t_bampair.tumor_bam = bam2;
    } else {
        std::cerr << "please provide valid format tumor bam file ! \n";
        exit(0);
    }

    // loading
    totalBamPairs.push_back(t_bampair);
    totalBamPairsNum++;
}

// loading bam list
// load tumor bam file only
void PolyScan::LoadBam(const std::string &bam) {
    BamTumors t_bamtumor;
    t_bamtumor.sName = "sample_name";
//    std::cout<<bam<<"\n";
    if (bam.find(".bam") != std::string::npos) {
        t_bamtumor.tumor_bam = bam;
    } else {
        std::cerr << "please provide valid format tumor bam file ! \n";
        exit(0);
    }        //add by yelab

    // loading
    totalBamTumors.push_back(t_bamtumor);
    totalBamTumorsNum++;
}

//add by yelab
void PolyScan::LoadBamn(const std::string &bam, const std::string &Name) {
    BamNormals t_bamnormal;

//   std::cout<<bam<<"\n";
//   std::cout<<bam.find(".bam")<<"\n";
//   std::cout<<std::string::npos<<"\n";
//	if (bam.find(".bam") != std::string::npos) {
    if (bam.find(".bam") != std::string::npos) {
        t_bamnormal.sName = Name;
//		t_bamnormal.normal_bam = abs_path(bam);
        t_bamnormal.normal_bam = bam;

    } else {
        std::cerr << "please provide valid format normal bam file ! \n";
        exit(0);
    }

    // loading
    totalBamNormals.push_back(t_bamnormal);
    totalBamNormalsNum++;

}

// read and load sites
bool PolyScan::LoadHomosAndMicrosates(std::ifstream &fin) {
    std::string chr;
    std::string bases;
    std::string fbases;
    std::string ebases;
    std::string line;
    std::string tChr = "";
    // count total loading sites
    //
    totalHomosites = 0;
    int loc;
    bit8_t siteLength;
    bit16_t tsiteLength;
    bit16_t siteBinary;
    bit16_t siteRepeats;
    bit16_t frontF;
    bit16_t tailF;
    double thres;

    int j = 0;
    BedChr tbedChr;
    BedRegion tbedRegion;
    bit16_t tIndex;

    // skip title
    getline(fin, line);
    while (getline(fin, line)) {
        std::stringstream linestream(line);
        linestream >> chr;
        linestream >> loc;
        linestream >> tsiteLength;
        linestream >> siteBinary;
        linestream >> siteRepeats;
        linestream >> frontF;
        linestream >> tailF;
        //xxx
        linestream >> bases;
        linestream >> fbases;
        linestream >> ebases;
        if (!paramd.hard) {
            linestream >> thres;
        } else {
            thres = paramd.hunterThreshold;
        }

        // filtering
        if (tsiteLength > 1 && paramd.HomoOnly == 1)
            continue;
        if (tsiteLength == 1 && paramd.MicrosateOnly == 1)
            continue;
        if (tsiteLength == 1
            && ((siteRepeats < paramd.MininalHomoForDis)
                || (siteRepeats > paramd.MaxHomoSize)))
            continue;
        if (tsiteLength > 1
            && ((siteRepeats < paramd.MinMicrosateForDis)
                || (siteRepeats > paramd.MaxMicrosateForDis)))
            continue;

        siteLength = tsiteLength & 255;
        // defined one region
        if (ifUserDefinedRegion) {
            if (chr != region_one.chr) {
                continue;
            } else {
                if ((loc < region_one.start)
                    || ((loc + siteLength * siteRepeats) > region_one.end)) {
                    continue;
                }
            }
        }
        // bed filtering
        if (ifUserDefinedBed) {
            // new chr
            if (tChr != chr) {
                j = 0;
                if (chrMaptoIndex.count(chr) > 0) {
                    tbedChr = beds[chrMaptoIndex[chr]];
                    tChr = tbedChr.chr;
                    tbedRegion = tbedChr.regions_list[j];
                } else {
                    continue;
                }
            }
            // filtering
            /*
             if (loc < tbedRegion.start) continue;
             if (loc > tbedRegion.end) {
             for (j; j<tbedChr.regions_list.size(); j++) {
             tbedRegion = tbedChr.regions_list[j];
             if (loc < tbedRegion.end) { break; }
             }
             if (j >= tbedChr.regions_list.size()) continue;
             }
             if ((loc + siteLength * siteRepeats) > tbedRegion.end) continue;
             */
            // filtering
            if ((loc + siteLength * siteRepeats) < tbedRegion.start)
                continue;
            if (loc > tbedRegion.end) {
                for (j; j < tbedChr.regions_list.size(); j++) {
                    tbedRegion = tbedChr.regions_list[j];
                    if (loc <= tbedRegion.end) {
                        break;
                    }
                }
                if ((loc + siteLength * siteRepeats) < tbedRegion.start)
                    continue;
                if (j >= tbedChr.regions_list.size())
                    continue;
            }
        }

        // load sites
        //HomoSite *toneSite = new HomoSite;
        HomoSite toneSite;
        toneSite.chr = chr;
        toneSite.location = loc;
        toneSite.typeLen = siteLength;
        toneSite.homoType = siteBinary;
        toneSite.length = siteRepeats;
        toneSite.frontKmer = frontF;
        toneSite.endKmer = tailF;
        toneSite.bases = bases;
        toneSite.fbases = fbases;
        toneSite.ebases = ebases;
        if (paramd.pro) {
            toneSite.thres = thres;
//        	std::cout<<"dfjjf"<<"\n";
        }
//        std::cout<<toneSite.thres<<"\n";
//        std::cout<<toneSite.thres<<"\n";

        toneSite.lowcut =
                ((loc - MAX_READ_LENGTH) > 0) ? (loc - MAX_READ_LENGTH) : 0;
        toneSite.highcut = loc + MAX_READ_LENGTH;

        totalSites.push_back(toneSite);
        totalHomosites++;

        linestream.clear();
        linestream.str("");

    } // end while

    if (totalHomosites != 0)
        return true;
    return false;
}

// bed regions ?
void PolyScan::BedFilterorNot() {
    if (beds.size() > 0)
        ifUserDefinedBed = true;
}

// test sites loading
void PolyScan::TestHomos() {
    for (unsigned long i = 0; i < totalHomosites; i++) {
        HomoSite *toneSite = &totalSites[i];
        std::cout << toneSite->chr << "\t" << toneSite->location << "\t"
                  << int(toneSite->typeLen) << "\t" << toneSite->homoType << "\t"
                  << toneSite->length << "\t" << toneSite->frontKmer << "\t"
                  << toneSite->endKmer << "\t" << sizeof(*toneSite) << "\n";
    }
}

// split windows
void PolyScan::SplitWindows() {
    Window oneW;
    HomoSite *second;
    HomoSite *first = &totalSites[0];

    oneW._start = first->location;
    oneW._end = oneW._start;
    oneW._chr = first->chr;
    oneW._startSite = oneW._endSite = &totalSites[0];
    for (int i = 1; i < totalHomosites; i++) {
        first = &totalSites[i];
        if ((first->chr == oneW._chr)
            && (first->location - oneW._start) < paramd.windowSize) {
            continue;
        }
        oneW._end = totalSites[i - 1].location + MAX_SPAN_SIZE;
        oneW._endSite = &totalSites[i - 1];
        oneW._siteCount = oneW._endSite - oneW._startSite + 1;
        // record one window
        oneW.ChangeStart();
        totalWindows.push_back(oneW);
        totalWindowsNum++;
        oneW._start = first->location;
        oneW._end = oneW._start;
        oneW._chr = first->chr;
        oneW._startSite = oneW._endSite = &totalSites[i];
    }
    oneW._end = first->location + MAX_SPAN_SIZE;
    oneW._endSite = first;
    oneW._siteCount = oneW._endSite - oneW._startSite + 1;
    // record this window
    oneW.ChangeStart();
    totalWindows.push_back(oneW);
    totalWindowsNum++;
}

// test windows
void PolyScan::TestWindows() {
    Window *oneW;
    for (int i = 0; i < totalWindowsNum; i++) {
        oneW = &totalWindows[i];
        std::cout << oneW->_chr << "\t" << oneW->_start << "\t"
                  << oneW->_siteCount << "\t" << oneW->_startSite->chr << "\t"
                  << oneW->_startSite->location << "\t" << oneW->_endSite->chr
                  << "\t" << oneW->_endSite->location << "\n";
    }
}

// initial distribution
void PolyScan::InithializeDistributions() {
    for (int i = 0; i < totalWindowsNum; i++) {
        totalWindows[i].InitialDisW();
    }
}

// release distribution
void PolyScan::releaseDistributions() {
    for (int i = 0; i < totalWindowsNum; i++) {
        totalWindows[i].ClearDis();
    }
}

// output distribution
void PolyScan::outputDistributions() {
    for (int i = 0; i < totalWindowsNum; i++) {
        totalWindows[i].OutputDisW();
    }
}

// get distribution 
void PolyScan::GetHomoDistribution(Sample &oneSample,
                                   const std::string &prefix) {
    oneSample.iniOutput(prefix);
    std::vector <SPLIT_READ> readsInWindow;
    for (int i = 0; i < totalWindowsNum; i++) {
        totalWindows[i].InitialDisW();
        totalWindows[i].GetDistribution(readsInWindow);
        totalWindows[i].DisGenotypingW(oneSample);
        totalWindows[i].PouroutDisW(oneSample);
        totalWindows[i].ClearDis();
        readsInWindow.clear();
        std::cout << "window: " << i << " done...:" << totalWindows[i]._chr
                  << ":" << totalWindows[i]._start << "-" << totalWindows[i]._end
                  << std::endl;
    }
    // FDR
    oneSample.calculateFDR();
    oneSample.pourOutSomaticFDR();
    // MSI score
    oneSample.pourOutMsiScore();
    oneSample.closeOutStream();
    oneSample.VerboseInfo();

}

// for tumor only input
void PolyScan::GetHomoTumorDistribution(Sample &oneSample,
                                        const std::string &prefix) {
    oneSample.iniTumorDisOutput(prefix);
    std::vector <SPLIT_READ> readsInWindow;
    for (int i = 0; i < totalWindowsNum; i++) {
        totalWindows[i].InitialTumorDisW();
        totalWindows[i].GetTumorDistribution(readsInWindow);
        totalWindows[i].PourTumoroutDisW(oneSample);
        totalWindows[i].PouroutTumorSomatic(oneSample);
        totalWindows[i].ClearTumorDis();
        readsInWindow.clear();
        std::cout << "window: " << i << " done...:" << totalWindows[i]._chr
                  << ":" << totalWindows[i]._start << "-" << totalWindows[i]._end
                  << std::endl;
    }
    oneSample.pourOutMsiScore();
    oneSample.closeOutStream();
    oneSample.VerboseInfo();

}

// add by YeLab
void PolyScan::GetHunterTumorDistribution(Sample &oneSample,
                                          const std::string &prefix) {
    oneSample.hunterIniTumorDisOutput(prefix);
    std::vector <SPLIT_READ> readsInWindow;
    for (int i = 0; i < totalWindowsNum; i++) {
        totalWindows[i].InitialTumorDisW();
        totalWindows[i].GetTumorDistribution(readsInWindow);
        if (0) {

            totalWindows[i].ClearTumorDis();
            break;
        }
        totalWindows[i].PourTumoroutDisW(oneSample);
        totalWindows[i].PouroutTumorSomaticH(oneSample);
        totalWindows[i].ClearTumorDis();
        readsInWindow.clear();
        std::cout << "window: " << i << " done...:" << totalWindows[i]._chr
                  << ":" << totalWindows[i]._start << "-" << totalWindows[i]._end
                  << std::endl;
    }
    oneSample.pourOutMsiScore();
    oneSample.closeOutStream();
    oneSample.VerboseInfo();

}

std::string getRel(std::string str, std::string pattern) {
    std::string::size_type pos;
    std::vector <std::string> result;

    str += pattern; //扩展字符串以方便操作
    int size = str.size();

    for (int i = 0; i < size; i++) {
        pos = str.find(pattern, i);
        if (pos < size) {
            std::string s = str.substr(i, pos - i);
            result.push_back(s);
            i = pos + pattern.size() - 1;
        }
    }
    return result[result.size() - 1];
}

void PolyScan::GetNormalDistrubution(Sample &oneSample,
                                     const std::string &prefix) { //add by yelab

    totalBamTumorsNum = 1;
//    BamTumors t_bamtumor;
//    t_bamtumor.sName = "sample_name";
//    t_bamtumor.tumor_bam = "";
//    totalBamTumors.push_back(t_bamtumor);
    std::ofstream fout;
    fout.open(prefix + getRel(paramd.homoFile, "/") + "_baseline");
    if (!fout) {
        std::cerr << "EOORR: please check your output path!\n<<\n";
    }
    fout << "chromosome" << "\t" << "location" << "\t" << "repeat_unit_length"
         << "\t" << "repeat_unit_binary" << "\t" << "repeat_times" << "\t"
         << "left_flank_binary " << "\t" << "right_flank_binary" << "\t"
         << "repeat_unit_bases" << "\t" << "left_flank_bases" << "\t"
         << "right_flank_bases" << "\t"
         //		<< "covReads"          <<"\t"
         << "threshold" << "\t" << "supportSamples" << "\n";

    omp_set_num_threads(paramd.numberThreads);
#pragma omp parallel
    for (unsigned short i = 0; i < totalBamNormalsNum; i++) {
        unsigned short j = i;
//		const char* per=(prefix+"/"+totalBamNormals[j].sName).c_str();
//		int isCreate = mkdir(per,00755);
        std::cout << "Process the " << j + 1 << " case : "
                  << totalBamNormals[j].sName << " "
                  << totalBamNormals[j].normal_bam << "\n";
//		std::cout<<j<<"\t"<<totalBamNormals[j].sName<<"\n";
//		train.trainIniNormalDisOutput(prefix+"/"+totalBamNormals[j].sName+"/"+totalBamNormals[j].sName);
        oneSample.trainIniNormalDisOutput(
                prefix + "/detail/" + totalBamNormals[j].sName);
        std::vector <SPLIT_READ> readsInWindow;
//        totalBamTumors[0].tumor_bam = totalBamNormals[j].normal_bam;

        for (int i = 0; i < totalWindowsNum; i++) {
            totalWindows[i].InitialTumorDisW();

//				totalWindows[i].GetTumorDistribution(readsInWindow);
            if (!totalBamNormals[j].normal_bam.empty()) {
                // extract reads


                totalWindows[i].LoadReads(readsInWindow, totalBamNormals[j].normal_bam.c_str());


                totalWindows[i].ScanReads(readsInWindow, 0, true);
//					std::cout<<j<<"\t"<<totalBamNormals[j].sName<<"\n";
                readsInWindow.clear();
            }

//			if (0) {
//
//				totalWindows[i].ClearTumorDis();
//				break;
//			}
//				totalWindows[i].PourTumoroutDisW(oneSample);
            totalWindows[i].PouroutTumorSomaticH(oneSample);
            totalWindows[i].ClearTumorDis();
            readsInWindow.clear();
            std::cout << "    Process " << j + 1 << " "
                      << totalBamNormals[j].sName << " " << "window: " << i
                      << " done...:" << totalWindows[i]._chr << ":"
                      << totalWindows[i]._start << "-" << totalWindows[i]._end
                      << std::endl;
        }
//					oneSample.pourOutMsiScore();
        oneSample.closeOutStreamTrain();
//					oneSample.VerboseInfo();
    }
//#pragma omp atomic;
    std::cout << "Build baseline for microsatellites ..." << std::endl;

    std::map<std::string, int>::iterator it;
    std::map<std::string, double> threshold;
    std::unordered_map <std::string, std::vector<double>> FinalSites;
    std::unordered_map < std::string, std::vector < double >> ::iterator
    iter;
    int supportNumCutOff = (int) totalBamNormalsNum * paramd.sampleRatio; // 2 is parameter
    std::vector<double> buf;
    std::vector<double> tmp;
    for (it = SitesSupport.begin(); it != SitesSupport.end(); ++it) {

//	    std::cout << (it->first) << " => " << (it->first)<< (it->second )<< '\n';
        if ((it->second) >= supportNumCutOff) {
            FinalSites[it->first] = buf;
        }
    }

    std::string chr;
    std::string location;
    std::string left_flank_bases;
    std::string right_flank_bases;
    std::string repeat_times;
    std::string rub;
    std::string rfb;
    std::string cov;
    std::string site;
    double pro_u;
    double pro_v;
    for (unsigned short j = 0; j < totalBamNormalsNum; j++) {
//		std::cout<<prefix+"/detail/"+totalBamNormals[j].sName<<"\n";

        std::ifstream fin;
        std::string oneLine = "";
        std::istringstream linestream;
        fin.open(
                (prefix + "/detail/" + totalBamNormals[j].sName + "_all").c_str());
//		if(!fin)
//		{
//			std::cout<<"open fail!\n";
//		}
        while (getline(fin, oneLine)) {

            linestream.str(oneLine);
            linestream >> chr >> location >> left_flank_bases
                       >> right_flank_bases >> repeat_times >> rub >> rfb >> cov
                       >> pro_u;
            linestream.clear();
            site = chr + "_" + location;
            if (FinalSites.find(site) != FinalSites.end())
                FinalSites[site].push_back(pro_u);

//			buf=FinalSites[];
//			buf

//			std::cout<<chr<<"  11111  "<<pro_u<<"\n";
        }
        fin.close();
    }

    for (iter = FinalSites.begin(); iter != FinalSites.end(); iter++) {
        tmp = iter->second;
//    	std::cout<<tmp.size()<<"\n";
        double sum = std::accumulate(std::begin(tmp), std::end(tmp), 0.0);
        double mean = sum / tmp.size();

        double accum = 0.0;
        std::for_each(std::begin(tmp), std::end(tmp), [&](const double d) {
            accum += (d - mean) * (d - mean);
        });
        double stdev = std::sqrt(accum / (tmp.size() - 1));
        threshold[iter->first] = mean + 3 * stdev;
    }
//    std::cout<<threshold.size()<<"fjjfjfsjs"<<"\t"<<"\n";

//    std::cout<<totalHomosites<<"\n";
    for (int k = 0; k < totalHomosites; k++) {
        site = totalSites[k].chr + "_" + std::to_string(totalSites[k].location);
//    	std::cout<<totalSites[k].chr<<"\t"
//				 <<totalSites[k].location<<"\t"<<"\n";
        if (threshold.find(site) != threshold.end()) {
//    		std::cout<<totalHomosites<<"\n";
            fout << totalSites[k].chr << "\t" << totalSites[k].location << "\t"
                 << std::to_string(totalSites[k].typeLen) << "\t"
                 << totalSites[k].homoType << "\t" << totalSites[k].length
                 << "\t" << totalSites[k].frontKmer << "\t"
                 << totalSites[k].endKmer << "\t" << totalSites[k].bases
                 << "\t" << totalSites[k].fbases << "\t"
                 << totalSites[k].ebases << "\t" << threshold[site] << "\t"
                 << SitesSupport[site] << "\n";
            "\n";
//    	    	std::cout<<"hhh"<<"\n";
        }
    }
    fout.close();
    std::cout << "The micorsatellites file  with baseline is in: "
              << prefix + getRel(paramd.homoFile, "/") + "_baseline\n"
              << std::endl;
}
