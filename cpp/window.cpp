
#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>

#include "window.h"
#include "param.h"
#include "polyscan.h"


#include <fstream>
extern Param paramd;
extern std::map<std::string, int> SitesSupport;
extern PolyScan polyscan;
extern char homo_code[];
extern char uhomo_code[];
extern bit8_t alphabet[];

Window::Window() :
		_start(0), _end(0), _chr(""), _siteCount(0), _startSite(NULL), _endSite(
				NULL) {
	//xxxx
	//
}
;

Window::~Window() {
	//xxx
}
;

void Window::InitialDisW() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->InitialDis();
	}
}
;

void Window::InitialTumorDisW() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->InitialTumorDis();
	}
}
;

void Window::InitialNormalDisW() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->InitialNormalDis();
	}
}

void Window::ClearDis() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->ReleaseMemory();
	}
}
;

void Window::ClearTumorDis() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->ReleaseTumorMemory();
	}
}
;

void Window::ClearNormalDis() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->ReleaseNormalMemory();
	}
}

void Window::OutputDisW() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->OutputDis();
	}
}
;

void Window::OutputTumorDisW() {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;

		p->OutputTumorDis();
	}
}
;

void Window::PouroutDisW(Sample &oneSample) {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		if (paramd.outputzeroDis) {
			if ((p->outDislabel(oneSample) )) {
				//        std::cout << p->outDislabel(oneSample) << "\n";
				p->PouroutDis(oneSample);
			}

		} else {
			p->PouroutDis(oneSample);

		}

	}
}
;

void Window::PourTumoroutDisW(Sample &oneSample) {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;


//        std::cout<<p->thres<<"\n";
		if (paramd.outputzeroDis) {
			if (p->outDislabelTumorOnly(oneSample)) {
				p->PourTumoroutDis(oneSample);
			}
		} else {
			p->PourTumoroutDis(oneSample);

		}
	}
}
;

void Window::DisGenotypingW(Sample &oneSample) {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->DisGenotyping(oneSample);
	}
}
;

void Window::PouroutTumorSomatic(Sample &oneSample) {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;

		p->DisTumorSomatic(oneSample);
	}
}
;

// add by YeLab ,for hunter
void Window::PouroutTumorSomaticH(Sample &oneSample) {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
		p->HunterDisTumorSomatic(oneSample);
//		std::cout<<p->thres<<"hhshdsh"<<"\n";
	}
}
;

void Window::PouroutNormalH(Sample &oneSample) {
	HomoSite *p = NULL;
	for (unsigned short i = 0; i < _siteCount; i++) {
		p = _startSite + i;
//		std::cout<<_startSite<<"\n";
		p->HunterTrain(oneSample);
	}
}
;

// change start
void Window::ChangeStart() {
	if ((_start - MAX_SPAN_SIZE) < 0) {
		_start = 0;
	} else {
		_start -= MAX_SPAN_SIZE;
	}
}

void Window::GetDistribution(std::vector<SPLIT_READ> &readsInWindow) {
	for (unsigned short j = 0; j < polyscan.totalBamPairsNum; j++) {
		// normal
		if (!polyscan.totalBamPairs[j].normal_bam.empty()) {
			// extract reads
			LoadReads(readsInWindow,
					polyscan.totalBamPairs[j].normal_bam.c_str(),polyscan.refPath);
			ScanReads(readsInWindow, j, false);
			readsInWindow.clear();
		}
		// tumor
		if (!polyscan.totalBamPairs[j].tumor_bam.empty()) {
			// extract reads
			LoadReads(readsInWindow,
					polyscan.totalBamPairs[j].tumor_bam.c_str(),polyscan.refPath);
			ScanReads(readsInWindow, j, true);
			readsInWindow.clear();
		}
	}
}

void Window::GetTumorDistribution(std::vector<SPLIT_READ> &readsInWindow) {

	for (unsigned short j = 0; j < polyscan.totalBamTumorsNum; j++) {
		// tumor
		if (!polyscan.totalBamTumors[j].tumor_bam.empty()) {
			// extract reads

			//添加额外参数refPath，支持读取cram文件
			LoadReads(readsInWindow,
					polyscan.totalBamTumors[j].tumor_bam.c_str(),polyscan.refPath);

			ScanReads(readsInWindow, j, true);
			readsInWindow.clear();
		}
	}
//	out.flush();
}

void Window::GetNormalDistribution(std::vector<SPLIT_READ> &readsInWindow) { //add by yelab
	for (unsigned short j = 0; j < polyscan.totalBamTumorsNum; j++) {
		// normal
		if (!polyscan.totalBamNormals[j].normal_bam.empty()) {
			// extract reads
			LoadReads(readsInWindow,
					polyscan.totalBamNormals[j].normal_bam.c_str(),polyscan.refPath);
			ScanReads(readsInWindow, j, true);
			readsInWindow.clear();
		}
	}
}

void Window::LoadReads(std::vector<SPLIT_READ> &readsInWindow,
		const std::string bam,std::string ref) {
	std::string tag = "";
	if (!bam.empty()) {
		// extract reads
//        std::cout<<bam.c_str()<<"\t"<<"bam_path"<<std::endl;

        if(ref.empty())
        {
            ReadInBamReads(bam.c_str(), _chr, _start, _end, readsInWindow, tag);
//            std::cout<<readsInWindow.size()<<std::endl;
        }
        else
        {
            ReadInCramReads(bam.c_str(), ref.c_str(),_chr, _start, _end, readsInWindow, tag);
//            std::cout<<readsInWindow.size()<<std::endl;
        }
	}
}

void Window::ScanReads(const std::vector<SPLIT_READ> &readsInWindow,
		unsigned short bamIndex, bool isTumor) {

	// openmp parallel
	omp_set_num_threads(paramd.numberThreads);
#pragma omp parallel for
	for (unsigned short i = 0; i < _siteCount; i++) {
		HomoSite *p = _startSite + i;
		unsigned long tsize = readsInWindow.size();
		for (unsigned long j = 0; j < tsize; j++) {
//		std::cout <<j << " "
//			  << readsInWindow[j].Mapped << " "
//			 << p->lowcut << " "
//			<< readsInWindow[j].MatchedRelPos << " " 
//			<<p->highcut << " "
//			 << "\n";
			if (readsInWindow[j].Mapped) {
				if ((readsInWindow[j].MatchedRelPos < p->lowcut)
						|| (readsInWindow[j].MatchedRelPos > p->highcut))
					continue;
			}
			unsigned short tCount = DoOneRead(readsInWindow[j].ReadSeq, p);
			if (isTumor){

		//	std::cout << readsInWindow[j].ReadSeq <<" "   << tCount <<" " << p->fbases << "\n";
			}
			if ((tCount > 0) && (tCount < paramd.s_dispots)) {
				if (isTumor) {
					p->tumorDis[bamIndex][tCount - 1]++;
				} else {
					p->normalDis[bamIndex][tCount - 1]++;
				}
			} else {
				// don't scan reverse if mapped
				if (readsInWindow[j].Mapped)
					continue;
				// reverse
				std::string tStr = readsInWindow[j].ReadSeq;
				ReverseComplement(tStr);
				tCount = DoOneRead(tStr, p);
				if ((tCount > 0) && (tCount < paramd.s_dispots)) {
					if (isTumor) {
						p->tumorDis[bamIndex][tCount - 1]++;
					} else {
						p->normalDis[bamIndex][tCount - 1]++;
					}
				}
			}
		}
	}
}

unsigned short Window::DoOneRead(const std::string &oneRead,
		const HomoSite *p) {
	std::string::size_type startPos = 0;
	unsigned short count = 0;
	while (std::string::npos != (startPos = oneRead.find(p->fbases, startPos))) {
		//std::cout<<startPos<<std::endl;
		count = 0;
		std::string::size_type tstart0 = startPos + p->fbases.length();
		std::string::size_type tstart = tstart0;
		while (tstart0 == (tstart = oneRead.find(p->bases, tstart))) {
			count++;
			tstart += p->bases.length();
			tstart0 = tstart;
		}
		// if get one 
		// std::cout << paramd.MininalHomoForDis << " "<< count <<"\n";
//		if ((p->typeLen == 1) && (count >= paramd.MininalHomoForDis)
//				|| (p->typeLen > 1) && (count >= paramd.MinMicrosate)) 
		if (count>0)
		{
			tstart = tstart0;
			if (tstart == (tstart0 = oneRead.find(p->ebases, tstart0))) {
				return count;
			}
		}
		startPos++;
	}
	return 0;
}

// reverse complement string
void Window::ReverseComplement(std::string &theWord) {
	char tempChar;
	unsigned int t_uint0;
	unsigned int t_uint1;
	for (int i = 0; i < theWord.length() / 2; i++) {
		tempChar = theWord[i];
		t_uint0 = alphabet[tempChar];
		t_uint1 = alphabet[theWord[theWord.length() - i - 1]];
		theWord[i] = uhomo_code[t_uint1];
		theWord[theWord.length() - i - 1] = uhomo_code[t_uint0];
	}
	if (theWord.length() % 2) {
		tempChar = theWord[theWord.length() / 2];
		t_uint0 = alphabet[tempChar];
		theWord[theWord.length() / 2] = uhomo_code[t_uint0];
	}
}

