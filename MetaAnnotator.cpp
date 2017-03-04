/*
 * (1) use scores of contigs to predict taxon of contigs
 * and (2) use taxon of contigs to predict taxon of clusters
 */
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <omp.h>

#include "Methods.h"
#include "LearnLamda.h"
using namespace std;
/////////////////////////////////////////////////////////////////////
//parameters
//static const unsigned UNASSIGNED  = ~(0U);
//const double ScoreThresh = -1.0;
int ReadLen = 0;
int MetaKmerLen = 5;
int ClusterSize = 0;
int MaxSpecies = 200;
int MinSpecies = 2;
int StrLen = 50;
int StrLenLowCover = 25;
double MC3_Thresh = 0.97;
int OrphanLen = 400;

int Num_Thread = 0;
int CtgLenThresh = 500;
int AlignThresh = 0;
/////////////////////////////////////////////////////////////////////
//global variables
KmerNode** KmerMap;
ContigsClass Ctgs;
ReadsClass Reads;
KmerNodeAloc NodePool;
USet uset;
NCBI_nodes_dmp NodesDmp;
vector<ClustTaxoInfoClass> TaxoOfClust;
vector<double> ClustCoverage;

AccTester acc_tester;
////////////////////////////////////////////////////////////////////////////
/////////for test
struct SpInfoOfCtgs
{
	int length;
	int NReads;

	int matchedLength;
	int spTaxId;
	int refLevel;
	int spFId;
	SpInfoOfCtgs():length(0),matchedLength(0),spTaxId(0),refLevel(0),spFId(0),NReads(0){}
};
vector<SpInfoOfCtgs> spInfoOfCtgs;
////////////////////////////////////////////////////////////////////////////
void usage()
{
	cerr << "usage:\t MetaAnnotator contig.fa read.fa bwt.idx nodes.dmp bwt.idx.mtx nrbwt.idx [options]" << endl;
	cerr << endl;
	exit(-1);
}

void printtime(string str = "")
{
	static time_t rawtime;
	static tm* timeinfo;
	cerr << str << endl;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	cerr<<asctime(timeinfo);
}
//////////////////////////////
void init()
{
	{
		unsigned long long u64Exp = ((((1ULL << 50)%HASHSIZE)<<14))%HASHSIZE;
		unsigned long long tExp = 1;
		for(int i=0;i<10;++i)
		{
			U64HASH[i] = tExp;
			tExp *= u64Exp;
			tExp %= HASHSIZE;
		}
	}
	KmerMap = new KmerNode*[HASHSIZE];
	for(unsigned i=0;i<HASHSIZE;++i)
		KmerMap[i] = NULL;
}

void input(int argc, char* argv[])
{
	//Important! readlength should belong to (50,128]. Default value = 75.
//	AddParameter("ReadLen", &ReadLen, INTEGER);
	//We use q-mer distribution in Phase 3. Here you can set qmer = 4 or 5. Default value = 4.
	AddParameter("qmer", &MetaKmerLen, INTEGER);
	//You can fix the #. of species. Default value = 0.
	//In this case, the program will predict the #. of species.
	AddParameter("Species", &ClusterSize, INTEGER);
	//We use binary search to predict #. of species, this is max value for the first iteration.
	//Set this number properly will improve the performance and speed up the program. Default = 600.
	AddParameter("MaxSpecies", &MaxSpecies, INTEGER);
	//We use binary search to predict #. of species, this is min value for the first iteration.
	//Set this number properly will improve the performance and speed up the program. Default = 2.
	AddParameter("MinSpecies", &MinSpecies, INTEGER);
	//wmer in phase 1, default = 50;
	AddParameter("wmer", &StrLen, INTEGER);
	//wmer for low coverage read , default = 25; if a smaller value is needed, mapinto() should be modifed to mappint more 16-mers.
	AddParameter("lowwmer", &StrLenLowCover, INTEGER);
	//Threshold for further-merge in MetaCluster 3.0. Default valule = 0.9. If you set this value larger, program will conceptually predict more species.
	AddParameter("MetaC3Threshold", &MC3_Thresh, FLOAT);
	//groups of size smaller than Orphan16merThresh will be considered as orphans
	AddParameter("Orphan", &OrphanLen, INTEGER);
	//
	//////////////////////////////////////////////////////
	//get read length
	{
		string path = argv[2];
		const int MAXLINE = Ctgs.MAXLINE;
		char* Buf = new char[MAXLINE];
		ReadLen = 0;
		ifstream ifs(path.c_str());
		if(ifs.fail()){
			cerr << "File open failed:\t"<<path<<endl;
			exit(-1);
		}
		ifs.getline(Buf,MAXLINE);
		ifs.getline(Buf,MAXLINE);
		while(!ifs.eof() && Buf[0]!='>'){
			ReadLen += strlen(Buf);
			ifs.getline(Buf,MAXLINE);
		}
		ifs.close();
		delete[]Buf;
	}
	//////////////////////////////////////////////////////////////////////////
	//set minimum length of contigs.
	AddParameter("CtgLenThresh", &CtgLenThresh, INTEGER);
	AlignThresh = ReadLen*0.95;
	AddParameter("AlignThresh", &AlignThresh, INTEGER);

	ProcessParameters(argc,argv);
	/////////////////////////////////////////////////////////////////
	cerr << dec << "ReadLen:\t " << ReadLen << endl;
	cerr << "CtgLenThresh:\t " << CtgLenThresh << endl;
	cerr << "AlignThresh:\t " << AlignThresh << endl;
	cerr << "MC3_Thresh:\t" << MC3_Thresh << endl;
	cerr << "TN_WEIGHT:\t" << TN_WEIGHT << endl;
	/////////////////////////////////////////////////////////////////
	if(INTEST){
		printtime("before loading contigs. ");
		system("ps ux");
	}
	Ctgs.init(argv[1], CtgLenThresh,KmerMap, NodePool);
	if(INTEST){
		cerr << "# of Contigs: " << Ctgs.CtgNum << "\tNodeNum:" << NodePool.getNodeNum()<< endl;
		printtime("Contigs loaded. ");
	}
	NodePool.fixCtgNum();
	//Reads.init(argv[2], ReadLen, AlignThresh, KmerMap, NodePool, Ctgs, acc_tester);
	if(INTEST){
		cerr << "NodeNum:" << NodePool.getNodeNum()<< endl;
		printtime("initializing uset: ");
		system("ps ux");
	}
	{
		int CtgNum = Ctgs.CtgNum;
		unsigned* ctglen = new unsigned[CtgNum];
		for(int i=0;i<CtgNum;++i)
			ctglen[i] = Ctgs.contigs[i]->length;
		//uset.init(Ctgs.CtgNum + Reads.ReadNum, CtgNum, ctglen);
		////////////////////////////////////////////
		long long totalCtgLen = 0;
		for(int i=0;i<CtgNum;++i)
			totalCtgLen += ctglen[i];
		if(ClusterSize<=0)
			ClusterSize = totalCtgLen/300000;
		cerr << "estimated #. clusters: " << ClusterSize << endl;
		////////////////////////////////////////////
		delete[]ctglen;
	}
	delete[] KmerMap;
	KmerMap = NULL;
	///////////////////////////////////////////////////////////
	if(INTEST){
		cerr << dec << "Total Reads:\t"<< Reads.TotalNum << endl;
		cerr << "Unmaped Reads:\t" << Reads.ReadNum << endl;
		cerr << "AlignThresh:\t" << Reads.AlignThresh << endl;
		printtime("initializing acc_tester:\t");
		system("ps ux");
		acc_tester.init(Reads.TotalNum, Reads.MatchId, Reads.ReadNum,Reads.NewIdToOldId);
	}
//	NodePool.shrinkSize(NULL,2);
}

bool tCtgInfoSort(const vector<int>& t1,const vector<int>&t2)
{
	if(t1[5]==t2[5])
	{
		if(t1[6]==t2[6])
			return t1[1]>t2[1];
		else return t1[6]<t2[6];
	}
	return t1[5]<t2[5];
}

 const int L7[] = {4,8,12,17,21,24,28};
const int NL = 7;
void printReMtx(int Evaluate[10][10][10]){
		cerr << "print cluster: " << endl;
		for(int l=0;l<10;++l){
			long long sum = 0;
			for(int i=0;i<=NL;++i)
				for(int j=0;j<=NL;++j)
					sum += Evaluate[l][i][j];
			if(sum==0)continue;
			cout << "print reference level:" << l << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << Evaluate[l][i][j] << '\t';
				cout << endl;
			}
			cout << endl<< endl;

			int lower = 0,high=0,tinc=0,corr=Evaluate[l][l][l];
			for(int i=0;i<l;++i)
				for(int j=0;j<=l;++j)
					lower += Evaluate[l][i][j];
			for(int i=l+1;i<NL;++i)
				high += Evaluate[l][i][i];
			for(int i=0;i<NL;++i)
				for(int j=0;j<NL;++j)
					tinc += Evaluate[l][i][j];
			tinc -= (lower+high+corr);
			cerr << tinc << '\t' << lower << '\t' << high << '\t' << corr << '\t' <<sum << endl;
		}
}

void calCtgScore(vector<map<pair<int,int>,vector<double> > >&taxid_score,char* argv[]){
        BWTs bwts(argv[3]);
        //vector<map<int,int> >totalhit(Ctgs.CtgNum);
   #pragma omp parallel for schedule(dynamic)
        for(int i=0;i<Ctgs.CtgNum;i++){
                 calCtg(bwts,Ctgs.contigs[i]->str,taxid_score[i],Ctgs.impbool[i]->strbool);
		}

	bwts.clear();
  
        ofstream ctgout(argv[7]);
	for (int i=0;i<Ctgs.CtgNum;++i){
            map<int,double> genomescore;
            ctgout<<i<<" "<<Ctgs.contigs[i]->str.length()<<" "<<bvToString(Ctgs.impbool[i]->strbool)<<endl;
	    	
            for(map<pair<int,int>,vector<double> >::const_iterator itr = taxid_score[i].begin();itr!=taxid_score[i].end();++itr)
               if (itr->second[0] >= 5){
                  ctgout<<itr->first.first<<" "<<itr->second[0]<<" ";
                  //genomescore[itr->first.first>>16] += itr->second[0];
               }
            /*
            vector<pair<int,double> > sortscore;
            for(map<int,double>::iterator itr1 = genomescore.begin();itr1 != genomescore.end();itr1++)
                  sortscore.push_back(make_pair(itr1->first,itr1->second));
            sort(sortscore.begin(),sortscore.end(),cmpsecond2);
            for(int i=0;i<sortscore.size() && i<20;i++){
                  if(sortscore[i].second < sortscore[0].second/20)
                       break;
                  ctgout<<sortscore[i].first<<" "<<sortscore[i].second<<" ";
            } 
            */
            ctgout<<endl;
	}
	
	ctgout.close();

}
int ClustNumEst(vector<map<pair<int,int>,vector<double> > >&taxid_score){

        ofstream fout("ctgscore100.txt");
	set<int>cand;
        for(int i=0;i<Ctgs.CtgNum;i++){

		        vector<pair<int,int> >singlehit;
                        for(map<pair<int,int>,vector<double> >::iterator itr = taxid_score[i].begin();itr!=taxid_score[i].end();itr++){
                          singlehit.push_back(make_pair(itr->first.first,(int)itr->second[0]));
                        }
			sort(singlehit.begin(),singlehit.end(),cmpsecond);
		        map<int,int>genome_hit;
			int sum=0;
			for(int j=0;j<singlehit.size() && j<10;j++){
				 genome_hit[singlehit[j].first >> 16] += singlehit[j].second;
                                 sum += singlehit[j].second;
			}



		        for(map<int,int>::iterator itr = genome_hit.begin();itr!= genome_hit.end();itr++){
                                if(itr->second  > sum * 0.5)
                                      cand.insert(itr->first);


			        fout<<itr->first<<" "<<itr->second<<" ";

                        }
			fout<<endl;


        }
        fout.close();
        cout<<"candsize: "<<cand.size()<<endl;
        return cand.size();
        //fout.close();
};
void anaCluster(MCPara& mcpara, MetaCluster& metacluster,vector<unsigned>& ReadNumInCtg,vector<map<pair<int,int>,vector<double> > >&taxid_score,char* argv[])
{
        //read info file
        map<int,int> gnum_to_tax;
        string bwtidx(argv[3]);
        bwtidx += ".info";
        ifstream idxinfo(bwtidx.c_str());
        assert(!idxinfo.fail());
        int gnum = 0;
        int ncbigid,gtaxid,glen;
        string path;
        while(!idxinfo.eof()){
                idxinfo>>ncbigid>>gtaxid>>glen>>path;
                if(idxinfo.eof())break;
                gnum_to_tax[gnum++] = gtaxid;
        }
        idxinfo.close();


	NodesDmp.init(argv[4]);
		//////////////////////////////////////////////////////////////////
/*	if(INTEST){
		///////////////////////////////////////////////
		//initial Forbid
		spInfoOfCtgs.resize(Ctgs.CtgNum);
		for(int i=0;i<Reads.TotalNum;++i)
			if(Reads.MatchId[i]>=0)
				++spInfoOfCtgs[Reads.MatchId[i]].NReads;
		ifstream ifs("/home/ywang/DB/list/No_Taxid.list");
		assert(!ifs.fail());
		const int MAXFile = 5000;
		int Fid2Taxid[MAXFile];
		for(int i=0;i<MAXFile;++i)
			Fid2Taxid[i] = 0;
		while(!ifs.eof()){
			int tn=-1,ttax;
			ifs >> tn >> ttax;
			if(tn>=0)Fid2Taxid[tn]=ttax;
		}
		ifs.close();
		//////////////////////////////////////////////////////////////////
		for(int i=0;i<Ctgs.CtgNum;++i){
			char tbuf[100];
			sscanf(Ctgs.info[i].c_str(),"%s\t%d/%d\t%d", tbuf,&spInfoOfCtgs[i].matchedLength,&spInfoOfCtgs[i].length,&spInfoOfCtgs[i].spFId);
			spInfoOfCtgs[i].spTaxId = Fid2Taxid[spInfoOfCtgs[i].spFId];
		}
		map<int,int>Fid2Lid;
		for(int i=0;i<Ctgs.CtgNum;++i)
			spInfoOfCtgs[i].refLevel = Fid2Lid[spInfoOfCtgs[i].spFId];
	}*/

	int VCtgNum = metacluster.Size;
	int ClustNum = 0;
	for(int i=0;i<VCtgNum;++i)
		ClustNum = max(metacluster.best[i],ClustNum);
	++ClustNum;

	ClustNum = Ctgs.CtgNum;
	vector<int>toClustId(Ctgs.CtgNum,0);
	for(int i=0;i<Ctgs.CtgNum;++i){
		//toClustId[i] = metacluster.best[mcpara.getNewId(uset.find(i))];
		toClustId[i] = i;
	}

	vector<vector<int> >CtgIdInCluster(ClustNum);
	for(int i=0;i<Ctgs.CtgNum;++i)
		CtgIdInCluster[toClustId[i]].push_back(i);
	TaxoOfClust.resize(ClustNum);
	vector<int>ClustLength(ClustNum);
	for(int i=0;i<CtgIdInCluster.size();++i)
		for(vector<int>::const_iterator itr2 = CtgIdInCluster[i].begin();itr2!=CtgIdInCluster[i].end();++itr2)
			ClustLength[i] += Ctgs.contigs[*itr2]->str.length();
	///////////////////////////////////////////////////////////////////////
	{
		ClustCoverage.resize(ClustNum);
		for(int i=0;i<CtgIdInCluster.size();++i){
			double sumN = 0,sumL = 0;
			for(vector<int>::const_iterator itr=CtgIdInCluster[i].begin();itr!=CtgIdInCluster[i].end();++itr){
				sumN += Ctgs.SumCoverage[*itr];
				sumL += Ctgs.SumBase[*itr];
			}
			ClustCoverage[i] = sumN/sumL;
		}
	}
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////



        /*
        map<string,int>freqmap;
        ifstream freqkmer("freqkmer-string-gene.txt");
        string kmer20;
        int kmerge;
        while(!freqkmer.eof()){
            freqkmer>>kmer20>>kmerge;
            freqmap[kmer20]=kmerge;

        }
        freqkmer.close();
        vector<pair<int,int> >freq_score(Ctgs.CtgNum);
        */



	/*
        {
		if(INTEST){printtime("Before initiating bwts. ");system("ps ux");}
		BWTs bwts(argv[3]);
		if(INTEST){printtime("Before Annotating contigs:");system("ps ux");}
	#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<Ctgs.CtgNum;++i){
			calTaxoForCtg(bwts, Ctgs.contigs[i]->str, taxid_score[i]);
			//calfreqkmer(Ctgs.contigs[i]->str,freqmap,freq_score[i]);
                }
		bwts.clear();
	}
        */



        /*
        ofstream hit("ctghit.txt");
        for(int i=0;i<Ctgs.CtgNum;i++){
             hit<<i<<" "<<toClustId[i]<<" ";
             vector<pair<int,int> > singlehit;
             for(map<int,int>::iterator itr=totalhit[i].begin();itr!=totalhit[i].end();itr++)
                  singlehit.push_back(make_pair(itr->first,itr->second));
             sort(singlehit.begin(),singlehit.end(),cmpsecond);
             for(int j=0;j<singlehit.size() && j<15;j++){
                  hit<<singlehit[j].first<<" "<<singlehit[j].second<<" ";

             }
             hit<<endl;
        }
        hit.close();
        */
        /*
        ofstream clust_kmer("clust_kmer.txt");

        for(int i=0;i<CtgIdInCluster.size();++i){
                int maxsp = -1;
                int maxnum = 0;
                map<int,int>clustsp;


                for(vector<int>::const_iterator itr1 = CtgIdInCluster[i].begin(); itr1!=CtgIdInCluster[i].end();++itr1){
                                int hitsp = freq_score[*itr1].first;
                                int hitnum = freq_score[*itr1].second;
                                if(hitnum > 0){

                                         if(clustsp.count(hitsp) == 0)
                                              clustsp[hitsp] = hitnum;
                                         else
                                              clustsp[hitsp] += hitnum;


                                }
                }

                for(map<int,int>::iterator itr2 = clustsp.begin();itr2!=clustsp.end();itr2++)
                                if(itr2->second > maxnum){
                                         maxsp = itr2->first;
                                         maxnum ++;
                                }
                clust_kmer<<i<<" "<<maxsp<<" "<<maxnum<<endl;
                cout<<i<<" "<<maxsp<<" "<<maxnum<<endl;

        }


        clust_kmer.close();
        return;
        */
	///////////////////////////////////////////////////
	////cal lamda
	if(INTEST){printtime("Before calculate lamda:");system("ps ux");}
	const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
	const int KLCA=(LCAK-MINLAMDAK)/STEPLAMDAK;
	map<int,vector<double> >SidLamda[kN];
	map<int,vector<double> >SidPrec[kN];
	map<int,vector<pair<double,double> > >SidNormal[kN];
	map<int,vector<int> >SidNormalPrior[kN];
	{
		const int MAX=1000000;
		char* Buf;
		Buf = new char[MAX];

		vector<int>Tid;
		char tBuf[1000];
		ifstream ifs((string(argv[3])+".info").c_str());
		assert(!ifs.fail());
		while(!ifs.eof()){
			ifs.getline(Buf,MAX);
			int gid,tid=-1,len=0;
			sscanf(Buf,"%d\t%d\t%d\t%s",&gid,&tid,&len,tBuf);
			if(tid==-1)break;
			Tid.push_back(tid);
		}
		ifs.close();
		delete[]Buf;
		int N = Tid.size();
		cerr << "# seq: " << N << endl;
		vector<vector<int> >Taxon(N);
		for(int i=0;i<N;++i)
			Taxon[i] = NodesDmp.getTaxo(Tid[i]);

		vector<vector<double> >Mtx[kN];
		for(int l=0;l<kN;++l){
			Mtx[l].resize(N);
			for(int i=0;i<N;++i)
				Mtx[l][i].resize(N,0);
		}
		ifs.open(argv[5]);
		assert(!ifs.fail());
		for(int k = 0;k<kN;++k){
			for(int i=0;i<N;++i)
				for(int j=0;j<N;++j)
					ifs >> Mtx[k][i][j];
		}
		ifs.close();
                /*
		for(int i=0;i<kN;++i){
			int kmer = MINLAMDAK+STEPLAMDAK*i;
			cerr << "load lamda for kmer:" << kmer << endl;
			getLamda(Taxon,Mtx[i],SidLamda[i],SidPrec[i],SidNormal[i],SidNormalPrior[i]);
		}
                */
	}

	/*
	for(int i=0;i<kN;++i){
		int kmer = MINLAMDAK+STEPLAMDAK*i;
		int j=strlen(argv[5]);
		for(;j>=0 && (argv[5][j]<='0' || argv[5][j]>'9');--j);
		argv[5][j] = i+'0'+MINLAMDAK/STEPLAMDAK;
		cerr << "load lamda for kmer:" << argv[5] << endl;
		getLamda(NodesDmp,(string(argv[3])+".info"),argv[5],SidLamda[i],SidPrec[i],SidNormal[i],SidNormalPrior[i]);
	}*/
	if(INTEST){printtime("Before annotating clusters:");system("ps ux");}
	///////////////////////////////////////////////////
	for(int i=0;i<kN;++i)
		cerr << "Lamda #: " << SidLamda[i].size() << endl;
	///////////////////////////////////////////////////////////////////////
	//getNRScore
	bool isNRIdxValid = true;
	{
		ifstream ifs(argv[6]);
		if(ifs.fail()){
			cerr << "NRIdx file open failed! Only complete genome database will be used!" << endl;
			isNRIdxValid = false;
		}
		ifs.close();
	}

	vector<map<int,double> >ClustTidScoreNr(CtgIdInCluster.size());
	vector<int> MaxClustTidByBoth(CtgIdInCluster.size());
	vector<int> MaxClustSidByBoth(CtgIdInCluster.size());

	if(isNRIdxValid){
	{
		vector<map<int,double> >CtgTidScore;
		getProtSearch(Ctgs, argv[6], CtgTidScore);
		////////////////////////////////////////////////////////////////////////
		for(int i=0;i<CtgTidScore.size();++i){
			int clusterId = toClustId[i];
			for(map<int,double>::const_iterator itr1=CtgTidScore[i].begin();itr1 !=CtgTidScore[i].end();++itr1){
				vector<int> allTaxon = NodesDmp.getTaxo(itr1->first);
				if(allTaxon[L7[6]]!=2759){
					/*
					int non0 = 0;
					while(non0 < 7 && allTaxon[L7[non0]] <= 0)++non0;
					if(non0<7)
					ClustTidScoreNr[clusterId][allTaxon[L7[non0]]] += itr1->second;
					*/
					ClustTidScoreNr[clusterId][allTaxon[L7[0]]] += itr1->second;
				}
			}
		}
	}
	{
		vector<map<int,double> >CtgTidScoreCG(taxid_score.size());
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<taxid_score.size();++i){
			for(map<pair<int,int>,vector<double> >::const_iterator itr1=taxid_score[i].begin();itr1 != taxid_score[i].end();++itr1)
				CtgTidScoreCG[i][itr1->first.second] = max(CtgTidScoreCG[i][itr1->first.second],itr1->second[KLCA]);
		}

		vector<map<int,double> >ClustTidScoreCG(CtgIdInCluster.size());
		for(int i=0;i<CtgTidScoreCG.size();++i){
			int clusterId = toClustId[i];
			for(map<int,double>::const_iterator itr1=CtgTidScoreCG[i].begin();itr1 !=CtgTidScoreCG[i].end();++itr1){
				vector<int> allTaxon = NodesDmp.getTaxo(itr1->first);
				if(allTaxon[L7[6]]!=2759)
					ClustTidScoreCG[clusterId][allTaxon[L7[0]]] += itr1->second;
			//	ClustTidScoreCG[clusterId][itr1->first] += itr1->second;
			}
		}

		//////////////////////////////////////////////////////////////////
		//test
		{
			cout << "clust max score: " << endl;
			for(int i=0;i<ClustTidScoreNr.size();++i){
				int maxNrTid = 0; double maxNrScore = 0;
				for(map<int,double>::const_iterator itr1=ClustTidScoreNr[i].begin();itr1!=ClustTidScoreNr[i].end();++itr1){
					if(itr1->second > maxNrScore){
						maxNrTid = itr1->first;
						maxNrScore = itr1->second;
					}
				}
				int maxCGTid = 0; double maxCGScore = 0;
				for(map<int,double>::const_iterator itr1=ClustTidScoreCG[i].begin();itr1!=ClustTidScoreCG[i].end();++itr1){
					if(itr1->second > maxCGScore){
						maxCGTid = itr1->first;
						maxCGScore = itr1->second;
					}
				}
				cout << maxNrTid << '\t' << maxCGTid << '\t' << maxNrScore << '\t' << maxCGScore << endl;
			}
		}
		//////////////////////////////////////////////////////////////////
		vector<map<int,double> >ClustTidScoreBoth(CtgIdInCluster.size());
		for(int i=0;i<ClustTidScoreBoth.size();++i){
			for(map<int,double>::const_iterator itr1=ClustTidScoreNr[i].begin();itr1!=ClustTidScoreNr[i].end();++itr1)
				if(itr1->first > 0)
				ClustTidScoreBoth[i][itr1->first] +=  (itr1->second);
			for(map<int,double>::const_iterator itr1=ClustTidScoreCG[i].begin();itr1!=ClustTidScoreCG[i].end();++itr1)
				if(itr1->first > 0)
		//		ClustTidScoreBoth[i][itr1->first] =  max(ClustTidScoreBoth[i][itr1->first],itr1->second);
				ClustTidScoreBoth[i][itr1->first] +=  (itr1->second);
		}
		//////////////////////////////////////////////////////////////////
#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<ClustTidScoreBoth.size();++i){
			int maxtid = 0; double maxscore = 0;
			for(map<int,double>::const_iterator itr2=ClustTidScoreBoth[i].begin();itr2!=ClustTidScoreBoth[i].end();++itr2){
				if(itr2->first > 0 && itr2->second > maxscore){
					maxscore = itr2->second;
					maxtid = itr2->first;
				}
			}
			MaxClustTidByBoth[i] = maxtid;
		}
	}
	{
		vector<map<pair<int,int>,double> >ClustSidTidScoreCG(CtgIdInCluster.size());
//#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<taxid_score.size();++i){
			int clusterId = toClustId[i];
			for(map<pair<int,int>,vector<double> >::const_iterator itr1=taxid_score[i].begin();itr1 != taxid_score[i].end();++itr1)
				ClustSidTidScoreCG[clusterId][itr1->first] += itr1->second[KLCA];
		}
		vector<double> ClustMaxScore(CtgIdInCluster.size());
		vector<int> ClustMaxLevel(CtgIdInCluster.size(),8);
		for(int i=0;i<ClustSidTidScoreCG.size();++i){
			vector<int> maxTaxon = NodesDmp.getTaxo(MaxClustTidByBoth[i]);
			for(map<pair<int,int>,double>::const_iterator itr1=ClustSidTidScoreCG[i].begin();itr1 != ClustSidTidScoreCG[i].end();++itr1){
				vector<int> curTaxon = NodesDmp.getTaxo(itr1->first.second);
				int matchlevel = 0;
				for(;matchlevel < 7 && (maxTaxon[L7[matchlevel]]<=0 || maxTaxon[L7[matchlevel]]!=curTaxon[L7[matchlevel]]);++matchlevel);
				if(matchlevel < 7 && matchlevel < ClustMaxLevel[i]){
					ClustMaxLevel[i] = matchlevel;
					MaxClustSidByBoth[i] = itr1->first.first;
					ClustMaxScore[i] = itr1->second;
				}else if(matchlevel == ClustMaxLevel[i]){
					if(ClustMaxScore[i] < itr1->second){
						ClustMaxScore[i] = itr1->second;
						MaxClustSidByBoth[i] = itr1->first.first;
					}
				}
			}
		}
		cerr << "Prot test 3." << endl;
	}
	}
	///////////////////////////////////////////////////
	//anotate clusters
	//////////////annotate clusters from score

	set<int>ignore;
        /*
        ifstream ignlist("ignore.txt");
        while(!ignlist.eof()){
                int gnum;
                ignlist>>gnum;
                ignore.insert(gnum);

        }
        ignlist.close();
	*/
	map<int,int>ClustMaxTid;
	map<int,int>ClustMaxScore;
	double P_CLUST = 0.5;
        ofstream fout("ana.txt");
        ofstream fout2("ana.txt.max");
        ofstream ctgscore("ctg.score");
	for(int i=0;i<CtgIdInCluster.size();++i){
		int maxSidB = MaxClustSidByBoth[i];
		vector<double> maxScore(kN,0);
		vector<int> maxSid(kN,-1);
		vector<int> maxTid(kN,-1);
		map<int,double> ClustTaxid_LCA_Score;
                map<int,double> ClustCtgScore;
		int curkLCA = -1;
		for(int j=0;j<kN;++j){
			map<int,double> ClustTaxidScore;

			for(vector<int>::const_iterator itr1 = CtgIdInCluster[i].begin(); itr1!=CtgIdInCluster[i].end();++itr1){
                                ctgscore<<*itr1<<" "<<Ctgs.contigs[*itr1]->str.length();
				for(map<pair<int,int>,vector<double> >::const_iterator itr = taxid_score[*itr1].begin();itr!=taxid_score[*itr1].end();++itr){
                                                 //if(ignore.count(itr->first.first >> 16) != 0)
                                                      //continue;
                                                 ClustTaxidScore[itr->first.first] += itr->second[j];
                                                 if(ClustCtgScore[itr->first.first] < itr->second[j])
                                                      ClustCtgScore[itr->first.first] = itr->second[j];
                                                 if(itr->second[j] >= 5)
                                                      ctgscore<<" "<<itr->first.first<<" "<<itr->second[j];
                                }
                                ctgscore<<endl;

                        }
			for(map<int,double>::const_iterator itr = ClustTaxidScore.begin();itr!=ClustTaxidScore.end();++itr){
			//	if(itr->second > maxScore[j]){
				if((isNRIdxValid && itr->first==maxSidB)||(!isNRIdxValid && itr->second>maxScore[j])){
					maxScore[j] = itr->second;
					maxSid[j] = itr->first;
					maxTid[j] = itr->first;
				}
			}
                        fout<<i<<" "<<ClustLength[i]<<" "<<maxSid[j]<<" "<<maxTid[j]<<" "<<maxScore[j]<<endl;
                        fout2<<i<<" "<<ClustLength[i]<<" "<<maxSid[j]<<" "<<maxTid[j]<<" "<<maxScore[j]<<endl;
			if(true){
				ClustTaxid_LCA_Score = ClustTaxidScore;
				curkLCA = j;
			}
		}
                curkLCA = 0;
		//ClustMaxScore[i] = maxScore[0];
		int maxkLCA = -1;
		for(int j=0;j<kN && maxSid[j]>=0;++j)
			maxkLCA = j;
		set<int>cand_taxo;
                vector<pair<int, double> > cand_score;
                vector<pair<int, double> > max_score;
		for(map<int,double>::const_iterator itr = ClustTaxid_LCA_Score.begin();itr!=ClustTaxid_LCA_Score.end();++itr){
			if(itr->second >= 0){

                                //int gnum = itr->first.first >> 16;




				//cand_taxo.insert(gnum_to_tax[gnum]);
                                cand_score.push_back(make_pair(itr->first,itr->second));



                        }
                        //cand_score.push_back(make_pair(itr->first.first,itr->second));
                }
                sort(cand_score.begin(),cand_score.end(),cmpsecond2);
                fout<<cand_score.size()<<" ";

                for(int k=0;k<10000 && k<cand_score.size();k++){
                        fout<<cand_score[k].first<<" "<<cand_score[k].second<<" ";
                }
                fout<<endl;

                for(map<int,double>::const_iterator itr = ClustCtgScore.begin();itr!= ClustCtgScore.end();itr++){
                        if(itr->second > 0){
                               max_score.push_back(make_pair(itr->first,itr->second));
                        }

                }
                sort(max_score.begin(),max_score.end(),cmpsecond2);

                fout2<<max_score.size()<<" ";
                for(int k=0;k<10000 && k<max_score.size();k++){
                        fout2<<max_score[k].first<<" "<<max_score[k].second<<" ";

                }
                fout2<<endl;

                continue;

		vector<vector<int> >V;
		for(set<int>::const_iterator itr=cand_taxo.begin();itr!=cand_taxo.end();++itr)
			V.push_back(NodesDmp.getTaxo(*itr));
		if(V.size()==0){
			cout << "cluster info:\t";
			cout << i << "\t" << 0 << '\t' << 9 << endl;
			continue;
		}
		int taxon = 0,level = 99;
		for(int k=0;k<NL;++k){
			bool allSame = true;
			for(int j=1;j<V.size() && allSame;++j){
				if(V[j][L7[k]] != V[0][L7[k]])
					allSame = false;
			}
			if(allSame){
				taxon = V[0][L7[k]];
				level = k;
				break;
			}
		}
		int firstLCA = taxon;
		//////////////////////////////////////
		//use lamda to decide taxon level.
		double finalscore = maxScore[curkLCA]/ClustLength[i];
                /*
		if(SidLamda[0].find(maxSid[curkLCA])!=SidLamda[0].end()){
			for(int k=level;k<NL-1;++k){
				int curKmerid = k<=maxkLCA?(maxkLCA-k):0;
				double candScore = maxScore[curKmerid]/ClustLength[i];
				if(candScore >= SidLamda[curKmerid][maxSid[curkLCA]][k])
					break;
				else{
					taxon = V[0][L7[k+1]];
					level = k+1;
					finalscore = candScore;
				}
			}
		}else cerr << "Sid not found." << i << '\t' << curkLCA << '\t' << maxSid[curkLCA] << endl;
                */
                /*
                if(SidLamda[0].find(maxSid[curkLCA])!=SidLamda[0].end()){
                        for(int k=level;k<NL-1;++k){
                                int curKmerid = 0;
                                double candScore = maxScore[curKmerid]/ClustLength[i];
                                if(candScore >= SidNormal[curKmerid][maxSid[curkLCA]][k].first){
                                        if(candScore >= SidNormal[curKmerid][maxSid[curkLCA]][k-1].first * 0.8){
                                               taxon = V[0][L7[k-1]];
                                               level = k-1;
                                               finalscore = candScore;
                                        }
                                        break;

                                }
                                else{
                                        taxon = V[0][L7[k+1]];
                                        level = k+1;
                                        finalscore = candScore;
                                }
                        }
                }else cerr << "Sid not found." << i << '\t' << curkLCA << '\t' << maxSid[curkLCA] << endl;
                */
		/////////////////////////////////////
		TaxoOfClust[i].taxon = taxon;
		TaxoOfClust[i].level = level;
//		TaxoOfClust[i].confidence = getTaxonConfidence(SidNormal[curkLCA][maxSid[curkLCA]],level,finalscore);
		vector<double> TaxonConfidence = getTaxonConfidence(SidNormal[curkLCA][maxSid[curkLCA]],SidNormalPrior[curkLCA][maxSid[curkLCA]],level,finalscore);
		for(int j=0;j<7;++j){
			TaxoOfClust[i].taxo[L7[j]].score = TaxonConfidence[j];
			TaxoOfClust[i].taxo[L7[j]].taxid = V[0][L7[j]];
		}
		TaxoOfClust[i].taxo[TaxoOfClust[i].taxo.size()-1].score = TaxonConfidence[TaxonConfidence.size()-1];
//				getTaxonConfidence(SidNormal[curkLCA][maxSid[curkLCA]],j,finalscore);
//		cout << i << "\tcluster taxon:\t" << taxon << '\t' << level << endl;
		if(INTEST){
			cout << "cluster: " << i <<'\t'<<ClustLength[i] <<'\t' << maxScore[curkLCA] <<'\t'<<maxSid[curkLCA]<<'\t'<<maxTid[curkLCA] << endl;
			vector<int> maxTaxon = NodesDmp.getTaxo(maxTid[curkLCA]);
			for(int j=0;j<NL;++j)
				cout << maxTaxon[L7[j]]<<'\t';
			cout << endl;
			vector<int> predTaxon = NodesDmp.getTaxo(taxon);
			for(int j=0;j<NL;++j)
				cout << predTaxon[L7[j]]<<'\t';
			cout << endl;
			ClustMaxTid[i] = taxon;
//			ClustMaxTid[i] = firstLCA;
//			ClustMaxTid[i] = MaxClustTidByBoth[i];
		}
	}
        fout.close();
        fout2.close();
        ctgscore.close();
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	if(INTEST)
	{
		cerr << "P_CLUST = " << P_CLUST << endl;
		cerr << "cal taxonomy for clusters" << endl;
		int Evaluate[10][10][10];
		int CorNearest[10][10][10];
		int IncNearest[10][10][10];
		int CorNNearest[10][10][10];
		int IncNNearest[10][10][10];
		for(int i=0;i<10;++i)
			for(int j=0;j<10;++j)
				for(int k=0;k<10;++k){
					Evaluate[i][j][k] = 0;
					CorNearest[i][j][k] = IncNearest[i][j][k] = 0;
					CorNNearest[i][j][k] = IncNNearest[i][j][k] = 0;
				}
		for(int i=0;i<CtgIdInCluster.size();++i){
			int maxTid = ClustMaxTid[i];
			vector<int> maxTaxon = NodesDmp.getTaxo(maxTid);

			int taxon =  TaxoOfClust[i].taxon;
			int level =  TaxoOfClust[i].level;
	//		if(taxon <= 0 || level >=6)continue;//140317
			vector<int> predTaxon = NodesDmp.getTaxo(taxon);
			int predLevel = 0;
			for(;predLevel<NL;++predLevel)
				if(predTaxon[L7[predLevel]]>0)
					break;

			for(int i1=0;i1<CtgIdInCluster[i].size();++i1){
				vector<int> correctTaxon = NodesDmp.getTaxo(spInfoOfCtgs[CtgIdInCluster[i][i1]].spTaxId);
				int correctLevel = 0;
				for(;correctLevel<NL;++correctLevel)
					if(predTaxon[L7[correctLevel]]==correctTaxon[L7[correctLevel]] && predTaxon[L7[correctLevel]])
						break;
				Evaluate[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
				/////////////////////////////////////////////////
				int correctmaxLevel = 0;
				for(;correctmaxLevel<NL;++correctmaxLevel)
					if(maxTaxon[L7[correctmaxLevel]]==correctTaxon[L7[correctmaxLevel]] && maxTaxon[L7[correctmaxLevel]])
						break;
				if(correctmaxLevel<=spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel){
					CorNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
					++CorNNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel];
				}
				else{
					IncNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
					++IncNNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel];
				}
				/////////////////////////////////////////////////////////////
				cerr<<"confidence of ctg " << CtgIdInCluster[i][i1]<<"\t"<<taxon<<"\t"<<level<<"\t" << spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel<<"\t"<<predLevel<<"\t"<<correctLevel<<"\t"<<spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads<<"\t";
				cerr << ClustLength[i] << "\t" << ClustMaxScore[i] << "\t";
				cerr << ClustMaxScore[i]/(double) ClustLength[i] << "\t";
				cerr << TaxoOfClust[i].taxo[L7[level]].score<<":"<<TaxoOfClust[i].taxo[L7[level]].taxid<<"\t";
				for(int j=0;j<7;++j)
					cerr << TaxoOfClust[i].taxo[L7[j]].score<<":"<<TaxoOfClust[i].taxo[L7[j]].taxid<<"\t";
				cerr << endl;
			}
		}
		cout << "cg evaluate: " << endl;
		cerr << "cg eva cluster: " << endl;
		for(int l=0;l<10;++l){
			long long sum = 0;
			for(int i=0;i<=NL;++i)
				for(int j=0;j<=NL;++j)
					sum += Evaluate[l][i][j];
			if(sum==0)continue;
			cout << "cg reference level:" << l << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << Evaluate[l][i][j] << '\t';
				cout << endl;
			}
			cout << endl<< endl;

			int lower = 0,high=0,tinc=0,corr=Evaluate[l][l][l];
			for(int i=0;i<l;++i)
				for(int j=0;j<=l;++j)
					lower += Evaluate[l][i][j];
			for(int i=l+1;i<NL;++i)
				high += Evaluate[l][i][i];
			for(int i=0;i<NL;++i)
				for(int j=0;j<NL;++j)
					tinc += Evaluate[l][i][j];
			tinc -= (lower+high+corr);
			cerr << tinc << '\t' << lower << '\t' << high << '\t' << corr << '\t' <<sum << endl;
		}
	}
	///////////////////////////////////////////////////////////////////////
	//predict by NR database
	/*
	map<int,vector<double> >NRSidLamda;
	map<int,vector<double> >NRSidPrec;
	map<int,vector<pair<double,double> > >NRSidNormal;
	map<int,vector<int> >NRSidNormalPrior;
	readNRLamda(argv[7],NRSidLamda,NRSidPrec,NRSidNormal,NRSidNormalPrior);
	cerr << "in lamda: " << NRSidLamda.size() << '\t'<< (NRSidLamda.begin()->first) << endl;
	*/
//	getNRLamda(NodesDmp,(string(argv[3])+".info.species"),argv[7],NRSidLamda,NRSidPrec,NRSidNormal,NRSidNormalPrior);
	///////////////////////////////////////////////////////////////////////
	/*
	const double PNR_CLUST = 0.01;
	//modified
#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<CtgIdInCluster.size();++i){
		{
			int tlevel = 0, taxon = TaxoOfClust[i].taxon;
			vector<int> allTaxon = NodesDmp.getTaxo(taxon);
			for(;tlevel<7 && allTaxon[L7[tlevel]]!=taxon;++tlevel);
			if(taxon>0 && tlevel<=3)continue;
			else{
				cout << "With NR, processing Cluster " << i << "\t" << CtgIdInCluster[i].size() << " : ";
				for(int j=0;j<CtgIdInCluster[i].size();++j)
					cout << CtgIdInCluster[i][j] << '\t';
				cout << endl;
			}
		}
		int maxtid = 0; double maxscore = 0;
		for(map<int,double>::const_iterator itr2=ClustTidScoreNr[i].begin();itr2!=ClustTidScoreNr[i].end();++itr2){
			if(itr2->first > 0 && itr2->second > maxscore){
				maxscore = itr2->second;
				maxtid = itr2->first;
			}
		}
		//TODO
		{
			for(map<int,double>::const_iterator itr2=ClustTidScoreNr[i].begin();itr2!=ClustTidScoreNr[i].end();++itr2){
				if(itr2->first == MaxClustTidByBoth[i]){
					maxscore = itr2->second;
					maxtid = itr2->first;
				}
			}
		}
	//	ClustMaxTid[i] = maxtid;
		set<int>cand_taxo;
		for(map<int,double>::const_iterator itr=ClustTidScoreNr[i].begin();itr!=ClustTidScoreNr[i].end();++itr){
			if(itr->second >= maxscore*(1.0-PNR_CLUST))
				cand_taxo.insert(itr->first);
		}
		vector<vector<int> >V;
		for(set<int>::const_iterator itr=cand_taxo.begin();itr!=cand_taxo.end();++itr)
			V.push_back(NodesDmp.getTaxo(*itr));
		if(V.size()==0)
			continue;
		int taxon = 0,level = 99;
		for(int k=0;k<NL;++k){
			bool allSame = true;
			for(int j=1;j<V.size() && allSame;++j){
				if(V[j][L7[k]] != V[0][L7[k]])
					allSame = false;
			}
			if(allSame){
				taxon = V[0][L7[k]];
				level = k;
				break;
			}
		}
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		//use lamda to decide taxon level.
		double finalscore = maxscore/ClustLength[i];
		if(NRSidLamda.find(maxtid)!=NRSidLamda.end()){
			for(int k=level;k<NL-1;++k){
				if(finalscore >= NRSidLamda[maxtid][k])
					break;
				else{
					taxon = V[0][L7[k+1]];
					level = k+1;
				}
			}
		}
		if(taxon>0 && level<7){
			TaxoOfClust[i].taxon = taxon;
			TaxoOfClust[i].level = level;
		}
	}
	if(INTEST)
	{
		cerr << "P_CLUST = " << P_CLUST << endl;
		cerr << "PNR_CLUST = " << PNR_CLUST << endl;
		cerr << "cal taxonomy for clusters" << endl;
		int Evaluate[10][10][10];
		int CorNearest[10][10][10];
		int IncNearest[10][10][10];
		int CorNNearest[10][10][10];
		int IncNNearest[10][10][10];
		for(int i=0;i<10;++i)
			for(int j=0;j<10;++j)
				for(int k=0;k<10;++k){
					Evaluate[i][j][k] = 0;
					CorNearest[i][j][k] = IncNearest[i][j][k] = 0;
					CorNNearest[i][j][k] = IncNNearest[i][j][k] = 0;
				}
		vector<int>unassigned(4);
		for(int i=0;i<CtgIdInCluster.size();++i){
			int maxTid = ClustMaxTid[i];
			vector<int> maxTaxon = NodesDmp.getTaxo(maxTid);

			int taxon =  TaxoOfClust[i].taxon;
			int level =  TaxoOfClust[i].level;
			if(taxon <= 0 || level >=6){
				for(int i1=0;i1<CtgIdInCluster[i].size();++i1)
					unassigned[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
				continue;//140317
			}
			vector<int> predTaxon = NodesDmp.getTaxo(taxon);
			int predLevel = 0;
			for(;predLevel<NL;++predLevel)
				if(predTaxon[L7[predLevel]]>0)
					break;

			for(int i1=0;i1<CtgIdInCluster[i].size();++i1){
				vector<int> correctTaxon = NodesDmp.getTaxo(spInfoOfCtgs[CtgIdInCluster[i][i1]].spTaxId);
				int correctLevel = 0;
				for(;correctLevel<NL;++correctLevel)
					if(predTaxon[L7[correctLevel]]==correctTaxon[L7[correctLevel]] && predTaxon[L7[correctLevel]])
						break;
				Evaluate[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
				/////////////////////////////////////////////////
				int correctmaxLevel = 0;
				for(;correctmaxLevel<NL;++correctmaxLevel)
					if(maxTaxon[L7[correctmaxLevel]]==correctTaxon[L7[correctmaxLevel]] && maxTaxon[L7[correctmaxLevel]])
						break;
				if(correctmaxLevel<=spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel){
					CorNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
					++CorNNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel];
				}
				else{
					IncNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel]+= spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads;
					++IncNNearest[spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel][predLevel][correctLevel];
				}
				/////////////////////////////////////////////////////////////
				cerr<<"confidence of ctg " << CtgIdInCluster[i][i1]<<"\t"<<taxon<<"\t"<<level<<"\t" << spInfoOfCtgs[CtgIdInCluster[i][i1]].refLevel<<"\t"<<predLevel<<"\t"<<correctLevel<<"\t"<<spInfoOfCtgs[CtgIdInCluster[i][i1]].NReads<<"\t";
				cerr << TaxoOfClust[i].taxo[L7[level]].score<<":"<<TaxoOfClust[i].taxo[L7[level]].taxid<<"\t";
				for(int j=0;j<7;++j)
					cerr << TaxoOfClust[i].taxo[L7[j]].score<<":"<<TaxoOfClust[i].taxo[L7[j]].taxid<<"\t";
				cerr << endl;
			}
		}
		for(int i=0;i<4;++i)
			cerr << "unassigned of level " << i << ":\t" << unassigned[i] << endl;
		cout << "evaluate: " << endl;
		cerr << "eva cluster: " << endl;
		for(int l=0;l<10;++l){
			long long sum = 0;
			for(int i=0;i<=NL;++i)
				for(int j=0;j<=NL;++j)
					sum += Evaluate[l][i][j];
			if(sum==0)continue;
			cout << "reference level:" << l << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << Evaluate[l][i][j] << '\t';
				cout << endl;
			}
			cout << endl<< endl;

			int lower = 0,high=0,tinc=0,corr=Evaluate[l][l][l];
			for(int i=0;i<l;++i)
				for(int j=0;j<=l;++j)
					lower += Evaluate[l][i][j];
			for(int i=l+1;i<NL;++i)
				high += Evaluate[l][i][i];
			for(int i=0;i<NL;++i)
				for(int j=0;j<NL;++j)
					tinc += Evaluate[l][i][j];
			tinc -= (lower+high+corr);
			cerr << tinc << '\t' << lower << '\t' << high << '\t' << corr << '\t' <<sum << endl;
		}
		cout << "eva max: " << endl;
		for(int l=0;l<10;++l){
			long long sum = 0;
			for(int i=0;i<=NL;++i)
				for(int j=0;j<=NL;++j)
					sum += (CorNearest[l][i][j]+IncNearest[l][i][j]);
			if(sum==0)continue;
			cout << "reference level:" << l << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << CorNearest[l][i][j]<<'/'<<IncNearest[l][i][j]<< '\t';
				cout << endl;
			}
			cout << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << CorNearest[l][i][j]/(double)(CorNearest[l][i][j]+IncNearest[l][i][j])<< '\t';
				cout << endl;
			}
			cout << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << CorNNearest[l][i][j]<<'/'<<IncNNearest[l][i][j]<< '\t';
				cout << endl;
			}
			cout << endl;
			for(int i=0;i<=NL;++i){
				for(int j=0;j<=NL;++j)
					cout << CorNNearest[l][i][j]/(double)(CorNNearest[l][i][j]+IncNNearest[l][i][j])<< '\t';
				cout << endl;
			}
			cout << endl;
			cout << endl<< endl;
		}
		printReMtx(CorNearest);
		printReMtx(IncNearest);
		printReMtx(CorNNearest);
		printReMtx(IncNNearest);
	}*/
}

void printVirtualContigs(map<int,vector<int> >& cluster_ctg,string filename)
{
	const int insertionN = 100;
	char buf[1000];
	for(int i=0;i<insertionN;++i)
		buf[i]='N';
	buf[insertionN] = 0;
	string insertion(buf);

	ofstream ofs(filename.c_str());
	assert(!ofs.fail());
	for(map<int,vector<int> >::const_iterator itr = cluster_ctg.begin();itr!=cluster_ctg.end();++itr)
	{
		ofs << Ctgs.info[*(itr->second.begin())] <<'\t';
		for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
			ofs << *itr2 << '_';
		ofs << endl;
		for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
		{
			if(itr2 != itr->second.begin())
			   	ofs << insertion;
			ofs << (Ctgs.contigs[*itr2]->str);
		}
		ofs << endl;
	}
	ofs.close();
}

void printClusterIdOfCtgs(string filename)
{
	if(INTEST){
		map<int,vector<int> > cluster_ctg;
		for(int i=0;i<Ctgs.CtgNum;++i)
			cluster_ctg[uset.find(i)].push_back(i);
		printVirtualContigs(cluster_ctg,filename);
	}
}

void dblearn(char* argv[]){
        NodesDmp.init(argv[4]);
        const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
	const int KLCA=(LCAK-MINLAMDAK)/STEPLAMDAK;
	map<int,vector<double> >SidLamda[kN];
	map<int,vector<double> >SidPrec[kN];
	map<int,vector<pair<double,double> > >SidNormal[kN];
	map<int,vector<int> >SidNormalPrior[kN];
	{
		const int MAX=1000000;
		char* Buf;
		Buf = new char[MAX];

		vector<int>Tid;
		char tBuf[1000];
		ifstream ifs((string(argv[3])+".info").c_str());
		assert(!ifs.fail());
		while(!ifs.eof()){
			ifs.getline(Buf,MAX);
			int gid,tid=-1,len=0;
			sscanf(Buf,"%d\t%d\t%d\t%s",&gid,&tid,&len,tBuf);
			if(tid==-1)break;
			Tid.push_back(tid);
		}
		ifs.close();
		delete[]Buf;
		int N = Tid.size();
		cerr << "# seq: " << N << endl;
		vector<vector<int> >Taxon(N);
		for(int i=0;i<N;++i)
			Taxon[i] = NodesDmp.getTaxo(Tid[i]);

		vector<vector<double> >Mtx[kN];
		for(int l=0;l<kN;++l){
			Mtx[l].resize(N);
			for(int i=0;i<N;++i)
				Mtx[l][i].resize(N,0);
		}
		ifs.open(argv[5]);
		assert(!ifs.fail());
		for(int k = 0;k<kN;++k){
			for(int i=0;i<N;++i)
				for(int j=0;j<N;++j)
					ifs >> Mtx[k][i][j];
		}
		ifs.close();
		for(int i=0;i<kN;++i){
			int kmer = MINLAMDAK+STEPLAMDAK*i;
			cerr << "load lamda for kmer:" << kmer << endl;
			getLamda(Taxon,Mtx[i],SidLamda[i],SidPrec[i],SidNormal[i],SidNormalPrior[i]);
		}
	}


}
int main(int argc, char* argv[])
{
	if(argc < 7)
		usage();
        
        init();
	input(argc,argv);
	if(INTEST)system("ps ux");
       
        vector<map<pair<int,int>,vector<double> > > taxid_score(Ctgs.CtgNum);
        calCtgScore(taxid_score,argv);
        NodePool.clear();
        return 0;        


        string mode(argv[6]);
        
	
        if(mode == "group")
           Groupctg(uset,argv[7]);
	NodePool.clear();
	if(INTEST)acc_tester.calAcc(uset);
        if(mode == "map"){
          string filename(argv[7]);
          outputReads(Reads,filename);
         return 0;
        }
        vector<int>newbin;
	printtime("Main: before MetaCluster. ");
	if(INTEST)system("ps ux");
	////////////////////////////////////////
//	if(INtEST)printClusterIdOfCtgs((string(argv[1])+".step1.contig.fa").c_str());
	////////////////////////////////////////
        //ClusterSize = 600;
	MCPara mcpara(MetaKmerLen, Ctgs, uset, acc_tester);
	MetaCluster metacluster(mcpara.KmerLen, mcpara.Size, mcpara.ReverKmerDistri, ClusterSize, MaxSpecies, MinSpecies, mcpara.GenoNum , mcpara.Component);
	if(INTEST)system("ps ux");

        //vector<map<pair<int,int>,vector<double> > > taxid_score(Ctgs.CtgNum);
        calCtgScore(taxid_score,argv);
        //ClusterSize = ClustNumEst(taxid_score);

        if(ClusterSize < 100) ClusterSize = 100;
        if(ClusterSize > 600) ClusterSize = 600;
        //ClusterSize = 600;

	ClusterSize ? (metacluster.muiltkmeans(10,ClusterSize,newbin)) : (metacluster.iterMeta(10,MC3_Thresh));
	//metacluster.iterMeta(10,MC3_Thresh);
        metacluster.clear();//best[]&Component[][] are not cleared.
	if(INTEST)system("ps ux");

        outputCluster(Reads.MatchId,Reads.TotalNum, uset,metacluster.best, mcpara,argv[2]);
	printtime("Main: before anaCluster. ");



	anaCluster(mcpara, metacluster, Reads.ReadNumInCtg,taxid_score,argv);
	printtime("Main: before output. ");
        return 0;
	outputResultReads(ClustCoverage,Reads.MatchId, Reads.TotalNum, uset,metacluster.best, mcpara,argv[2], TaxoOfClust);
	//outputResultContigs(ClustCoverage,Ctgs,metacluster.Size, uset, metacluster.best,  mcpara,argv[1], TaxoOfClust);
	if(INTEST){
		system("ps ux");
		printtime("Clustering finished. ");
		acc_tester.getPreSen4Other(metacluster.getComp(),metacluster.Size,metacluster.best);
		map<int,vector<int> > cluster_ctg;
		for(int i=0;i<Ctgs.CtgNum;++i)
			cluster_ctg[metacluster.best[mcpara.getNewId(uset.find(i))]].push_back(i);
		cout << "cluster id after step3: " << Ctgs.CtgNum << '\t' << cluster_ctg.size()<<endl;
		for(map<int,vector<int> >::const_iterator itr = cluster_ctg.begin();itr!=cluster_ctg.end();++itr)
		{
			cout << (itr->first) << '\t' << (itr->second.size()) << ":\t";
			for(vector<int>::const_iterator itr2 = itr->second.begin();itr2!=itr->second.end();++itr2)
				cout << *itr2 << '\t';
			cout << endl;
		}
		cout << "<end> cluster id after step3: " << Ctgs.CtgNum<<endl;
//		printVirtualContigs(cluster_ctg,(string(argv[1])+".step3.contig.fa"));
	}
	cerr << "Finished." << endl;
	return 0;
}
