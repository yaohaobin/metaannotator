#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <bitset>
#include "Methods.h"

using namespace std;

bool cmpReads(const KmerNode* p1, const KmerNode* p2)
{
	return (p1->VSize - p1->CtgNum) < (p2->VSize - p2->CtgNum);

}
bool cmpsecond(const pair<int,int>& p1,const pair<int,int>& p2){
        return p1.second > p2.second;

}

bool cmpsecond2(const pair<int,double>& p1,const pair<int,double>& p2){
        return p1.second > p2.second;


}
void calTaxoForCtg(const BWTs& bwts, const string& ctgstr, map<pair<int,int>,vector<double> >& ctg_sp_score,vector<bool>& bv  ){
  const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
  vector<int>substrLength;
  for(int lamdaK = MINLAMDAK;lamdaK <= MAXLAMDAK;lamdaK += STEPLAMDAK)
	substrLength.push_back(lamdaK);
  int Nsubstr = substrLength.size();
  int SHORTEST_PARA_KMER = substrLength[0];
  int LONGEST_PARA_KMER = substrLength[Nsubstr-1];
  int startidx = 0;
  int fragidx = 0;
 
  while(startidx < ctgstr.length()){
    int nextlen = 1000;
	bool imp = 0;
	if(ctgstr.length() - startidx < nextlen)
	    nextlen = ctgstr.length() - startidx;
	string str = ctgstr.substr(startidx,nextlen);	
	map<pair<int,int>,vector<double> > sp_score;
	set<pair<string,string> >kmers;
	{
		string rstr = str;
		for(int i=0;i<str.length();++i)
		{
			char ch = str[str.length()-1-i];
			switch(ch)
			{
				case 'A':rstr[i]='T';break;
				case 'C':rstr[i]='G';break;
				case 'G':rstr[i]='C';break;
				case 'T':rstr[i]='A';break;
				default:rstr[i]='N';
			}
		}
		int strend = (str.length()-(LONGEST_PARA_KMER-1));
		for(int i=0;i<strend;++i)
		{
			string k1 = str.substr(i,LONGEST_PARA_KMER);
			string k2 = rstr.substr(str.length()-LONGEST_PARA_KMER-i,LONGEST_PARA_KMER);
			if(k1>k2)swap(k1,k2);
			kmers.insert(make_pair(k1,k2));
		}
		int strend2 = (str.length()-(SHORTEST_PARA_KMER-1));
		for(int i=strend;i<strend2;++i)
		{
			string k1 = str.substr(i);
			string k2 = rstr.substr(0,str.length()-i);
			if(k1>k2)swap(k1,k2);
			kmers.insert(make_pair(k1,k2));
		}
		for(int i=0;i<LONGEST_PARA_KMER-SHORTEST_PARA_KMER;++i){
			string k1 = str.substr(0,LONGEST_PARA_KMER-1-i);
			string k2 = rstr.substr(str.length()-(LONGEST_PARA_KMER-1-i));
			if(k1>k2)swap(k1,k2);
			kmers.insert(make_pair(k1,k2));
		}
	}
	for(set<pair<string,string> >::const_iterator itr = kmers.begin();itr!=kmers.end();++itr){
               
		vector<set<pair<int,int> > >gross_cand_taxid;
		bwts.substrSearch(itr->first, substrLength, gross_cand_taxid);
		vector<set<pair<int,int> > >gross_cand_taxid2;
		bwts.substrSearch(itr->second, substrLength, gross_cand_taxid2);
		for(int i=0;i<Nsubstr;++i){
			for(set<pair<int,int> >::const_iterator itr2=gross_cand_taxid2[i].begin();itr2!=gross_cand_taxid2[i].end();++itr2)
				gross_cand_taxid[i].insert(*itr2);
		}
		for(int i=0;i<Nsubstr;++i){
			set<pair<int,int> >cand_taxid;
			for(set<pair<int,int> >::const_iterator itr2 = gross_cand_taxid[i].begin();itr2!=gross_cand_taxid[i].end();++itr2){
				cand_taxid.insert(*itr2);
				//genomehit[itr2->first]++;
			}
			if(cand_taxid.size()==0)continue;
			double unit = 1.0;
                      
			for(set<pair<int,int> >::const_iterator itr2 = cand_taxid.begin();itr2!=cand_taxid.end();++itr2){
                                //cout<<itr2->first<<" "<<itr2->second<<endl;
                                
				if(sp_score[*itr2].size()==0)
					sp_score[*itr2].resize(kN);
				sp_score[*itr2][i] += unit;
                                
			}
		}
	}
	
	for(map<pair<int,int>,vector<double> >::iterator itr = sp_score.begin();itr!=sp_score.end();itr++){
	    if(ctg_sp_score[itr->first].size() == 0)
		    ctg_sp_score[itr->first].resize(kN);
	    if(itr->second[0] > ctg_sp_score[itr->first][0])
		    ctg_sp_score[itr->first][0] = itr->second[0];
	    if(!imp && itr->second[0]>=5){
		    bv[fragidx] = 1;
		    imp = 1;
	    }		   
	}	    
	startidx += nextlen / 2;
	fragidx++;
  }
}



void calfreqkmer(const string& str,map<string,int>& freq,pair<int,int>& freq_score){
        freq_score.first = -1; 
        freq_score.second = 0;
        const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
        vector<int>substrLength;
        for(int lamdaK = MINLAMDAK;lamdaK <= MAXLAMDAK;lamdaK += STEPLAMDAK)
                substrLength.push_back(lamdaK);
        int Nsubstr = substrLength.size();
        int SHORTEST_PARA_KMER = substrLength[0];
        int LONGEST_PARA_KMER = substrLength[Nsubstr-1];
        map<int,int>freqge;
        int maxge = -1;
       
        {
                string rstr = str;
                for(int i=0;i<str.length();++i)
                {
                        char ch = str[str.length()-1-i];
                        switch(ch)
                        {
                                case 'A':rstr[i]='T';break;
                                case 'C':rstr[i]='G';break;
                                case 'G':rstr[i]='C';break;
                                case 'T':rstr[i]='A';break;
                                default:rstr[i]='N';
                        }
                }
                int strend = (str.length()-(LONGEST_PARA_KMER-1));
                for(int i=0;i<strend;++i)
                {
                        string k1 = str.substr(i,LONGEST_PARA_KMER);
                        string k2 = rstr.substr(str.length()-LONGEST_PARA_KMER-i,LONGEST_PARA_KMER);
                        if(k1>k2)swap(k1,k2);
                        if(freq.count(k1) !=0){
                             int ge = freq[k1];
                             if(freqge.count(ge) == 0)
                                     freqge[ge] = 1;
                             else
                                     freqge[ge] ++;

                        }
                        if(freq.count(k2) !=0){
                             int ge = freq[k2];
                             if(freqge.count(ge) == 0)
                                     freqge[ge] = 1;
                             else
                                     freqge[ge] ++;

                        }
                }
                int strend2 = (str.length()-(SHORTEST_PARA_KMER-1));
                for(int i=strend;i<strend2;++i)
                {
                        string k1 = str.substr(i);
                        string k2 = rstr.substr(0,str.length()-i);
                        if(k1>k2)swap(k1,k2);
                        if(freq.count(k1) !=0){
                             int ge = freq[k1];
                             if(freqge.count(ge) == 0)
                                     freqge[ge] = 1;
                             else
                                     freqge[ge] ++;

                        }
                        if(freq.count(k2) !=0){
                             int ge = freq[k2];
                             if(freqge.count(ge) == 0)
                                     freqge[ge] = 1;
                             else
                                     freqge[ge] ++;

                        }

                }
                for(int i=0;i<LONGEST_PARA_KMER-SHORTEST_PARA_KMER;++i){
                        string k1 = str.substr(0,LONGEST_PARA_KMER-1-i);
                        string k2 = rstr.substr(str.length()-(LONGEST_PARA_KMER-1-i));
                        if(k1>k2)swap(k1,k2);
                        if(freq.count(k1) !=0){
                             int ge = freq[k1];
                             if(freqge.count(ge) == 0)
                                     freqge[ge] = 1;
                             else
                                     freqge[ge] ++;

                        }
                        if(freq.count(k2) !=0){
                             int ge = freq[k1];
                             if(freqge.count(ge) == 0)
                                     freqge[ge] = 1;
                             else
                                     freqge[ge] ++;

                        }

                }
                int maxnum=-1;
                for(map<int,int>::iterator itr=freqge.begin();itr!=freqge.end();itr++){
                        if(itr->second > maxnum){
                             maxge = itr->first;
                             maxnum = itr->second;
                       }                             
                }
                freq_score.first = maxge;
                freq_score.second = maxnum;
        }


}
void calTaxoForCtg(const BWTs& bwts, const string& str, map<pair<int,int>,vector<double> >& sp_score){
	const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
	vector<int>substrLength;
	for(int lamdaK = MINLAMDAK;lamdaK <= MAXLAMDAK;lamdaK += STEPLAMDAK)
		substrLength.push_back(lamdaK);
	int Nsubstr = substrLength.size();
	int SHORTEST_PARA_KMER = substrLength[0];
	int LONGEST_PARA_KMER = substrLength[Nsubstr-1];

	set<pair<string,string> >kmers;
	{
		string rstr = str;
		for(int i=0;i<str.length();++i)
		{
			char ch = str[str.length()-1-i];
			switch(ch)
			{
				case 'A':rstr[i]='T';break;
				case 'C':rstr[i]='G';break;
				case 'G':rstr[i]='C';break;
				case 'T':rstr[i]='A';break;
				default:rstr[i]='N';
			}
		}
		int strend = (str.length()-(LONGEST_PARA_KMER-1));
		for(int i=0;i<strend;++i)
		{
			string k1 = str.substr(i,LONGEST_PARA_KMER);
			string k2 = rstr.substr(str.length()-LONGEST_PARA_KMER-i,LONGEST_PARA_KMER);
			if(k1>k2)swap(k1,k2);
			kmers.insert(make_pair(k1,k2));
		}
		int strend2 = (str.length()-(SHORTEST_PARA_KMER-1));
		for(int i=strend;i<strend2;++i)
		{
			string k1 = str.substr(i);
			string k2 = rstr.substr(0,str.length()-i);
			if(k1>k2)swap(k1,k2);
			kmers.insert(make_pair(k1,k2));
		}
		for(int i=0;i<LONGEST_PARA_KMER-SHORTEST_PARA_KMER;++i){
			string k1 = str.substr(0,LONGEST_PARA_KMER-1-i);
			string k2 = rstr.substr(str.length()-(LONGEST_PARA_KMER-1-i));
			if(k1>k2)swap(k1,k2);
			kmers.insert(make_pair(k1,k2));
		}
	}
	for(set<pair<string,string> >::const_iterator itr = kmers.begin();itr!=kmers.end();++itr){
		vector<set<pair<int,int> > >gross_cand_taxid;
		bwts.substrSearch(itr->first, substrLength, gross_cand_taxid);
		vector<set<pair<int,int> > >gross_cand_taxid2;
		bwts.substrSearch(itr->second, substrLength, gross_cand_taxid2);
		for(int i=0;i<Nsubstr;++i){
			for(set<pair<int,int> >::const_iterator itr2=gross_cand_taxid2[i].begin();itr2!=gross_cand_taxid2[i].end();++itr2)
				gross_cand_taxid[i].insert(*itr2);
		}
		
		for(int i=0;i<Nsubstr;++i){
			vector<pair<int,int> >cand_taxid;
			for(set<pair<int,int> >::const_iterator itr2 = gross_cand_taxid[i].begin();itr2!=gross_cand_taxid[i].end();++itr2)
				cand_taxid.push_back(*itr2);
		
			if(cand_taxid.size()==0)continue;
			double unit = 1.0/cand_taxid.size();
			for(vector<pair<int,int> >::const_iterator itr2 = cand_taxid.begin();itr2!=cand_taxid.end();++itr2){
				if(sp_score[*itr2].size()==0)
					sp_score[*itr2].resize(kN);
				sp_score[*itr2][i] += unit;
			}
		}
	}
}

void calTaxoForCtg(const BWTs& bwts, const string& str, map<pair<int,int>,vector<double> >& sp_score,int scope,map<int,int>& tomap){
        const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
        vector<int>substrLength;
        for(int lamdaK = MINLAMDAK;lamdaK <= MAXLAMDAK;lamdaK += STEPLAMDAK)
                substrLength.push_back(lamdaK);
        int Nsubstr = substrLength.size();
        int SHORTEST_PARA_KMER = substrLength[0];
        int LONGEST_PARA_KMER = substrLength[Nsubstr-1];

        set<pair<string,string> >kmers;
        {
                string rstr = str;
                for(int i=0;i<str.length();++i)
                {
                        char ch = str[str.length()-1-i];
                        switch(ch)
                        {
                                case 'A':rstr[i]='T';break;
                                case 'C':rstr[i]='G';break;
                                case 'G':rstr[i]='C';break;
                                case 'T':rstr[i]='A';break;
                                default:rstr[i]='N';
                        }
                }
                int strend = (str.length()-(LONGEST_PARA_KMER-1));
                for(int i=0;i<strend;++i)
                {
                        string k1 = str.substr(i,LONGEST_PARA_KMER);
                        string k2 = rstr.substr(str.length()-LONGEST_PARA_KMER-i,LONGEST_PARA_KMER);
                        if(k1>k2)swap(k1,k2);
                        kmers.insert(make_pair(k1,k2));
                }
                int strend2 = (str.length()-(SHORTEST_PARA_KMER-1));
                for(int i=strend;i<strend2;++i)
                {
                        string k1 = str.substr(i);
                        string k2 = rstr.substr(0,str.length()-i);
                        if(k1>k2)swap(k1,k2);
                        kmers.insert(make_pair(k1,k2));
                }
                for(int i=0;i<LONGEST_PARA_KMER-SHORTEST_PARA_KMER;++i){
                        string k1 = str.substr(0,LONGEST_PARA_KMER-1-i);
                        string k2 = rstr.substr(str.length()-(LONGEST_PARA_KMER-1-i));
                        if(k1>k2)swap(k1,k2);
                        kmers.insert(make_pair(k1,k2));
                }
        }
        for(set<pair<string,string> >::const_iterator itr = kmers.begin();itr!=kmers.end();++itr){
                vector<set<pair<int,int> > >gross_cand_taxid;
                bwts.substrSearch(itr->first, substrLength, gross_cand_taxid);
                vector<set<pair<int,int> > >gross_cand_taxid2;
                bwts.substrSearch(itr->second, substrLength, gross_cand_taxid2);
                for(int i=0;i<Nsubstr;++i){
                        for(set<pair<int,int> >::const_iterator itr2=gross_cand_taxid2[i].begin();itr2!=gross_cand_taxid2[i].end();++itr2)
                                gross_cand_taxid[i].insert(*itr2);
                }

                for(int i=0;i<Nsubstr;++i){
                        vector<pair<int,int> >cand_taxid;
                        for(set<pair<int,int> >::const_iterator itr2 = gross_cand_taxid[i].begin();itr2!=gross_cand_taxid[i].end();++itr2)
                                if(tomap[itr2->second] == scope)
                                   cand_taxid.push_back(*itr2);

                        if(cand_taxid.size()==0)continue;
                        double unit = 1.0/cand_taxid.size();
                        for(vector<pair<int,int> >::const_iterator itr2 = cand_taxid.begin();itr2!=cand_taxid.end();++itr2){
                                if(sp_score[*itr2].size()==0)
                                        sp_score[*itr2].resize(kN);
                                sp_score[*itr2][i] += unit;
                        }
                }
        }
}




/*
void calTaxoForCtg(const BWTs& bwts, const string& str, map<pair<int,int>,vector<double> >& sp_score){
	const int kN = (MAXLAMDAK-MINLAMDAK)/STEPLAMDAK + 1;
	for(int lamdaK = MINLAMDAK;lamdaK <= MAXLAMDAK;lamdaK += STEPLAMDAK){
		const int lamdaID=(lamdaK-MINLAMDAK)/STEPLAMDAK;
		set<pair<string,string> >kmers;
		{
			string rstr = str;
			for(int i=0;i<str.length();++i)
			{
				char ch = str[str.length()-1-i];
				switch(ch)
				{
					case 'A':rstr[i]='T';break;
					case 'C':rstr[i]='G';break;
					case 'G':rstr[i]='C';break;
					case 'T':rstr[i]='A';break;
					default:rstr[i]='N';
				}
			}
			int strend = (str.length()-(lamdaK-1));
			for(int i=0;i<strend;++i)
			{
				string k1 = str.substr(i,lamdaK);
				string k2 = rstr.substr(str.length()-lamdaK-i,lamdaK);
				if(k1>k2)swap(k1,k2);
				kmers.insert(make_pair(k1,k2));
			}
		}
		for(set<pair<string,string> >::const_iterator itr = kmers.begin();itr!=kmers.end();++itr)
		{
			set<pair<int,int> > gross_cand_taxid;
			bwts.search(itr->first, gross_cand_taxid);
			set<pair<int,int> > gross_cand_taxid2;
			bwts.search(itr->second, gross_cand_taxid2);
			for(set<pair<int,int> >::const_iterator itr2=gross_cand_taxid2.begin();itr2!=gross_cand_taxid2.end();++itr2)
				gross_cand_taxid.insert(*itr2);
			//////////////////////////////////////////////
			vector<pair<int,int> >cand_taxid;
			for(set<pair<int,int> >::const_iterator itr2 = gross_cand_taxid.begin();itr2!=gross_cand_taxid.end();++itr2)
				cand_taxid.push_back(*itr2);
			//////////////////////////////////////////////
			if(cand_taxid.size()==0)continue;
			double unit = 1.0/cand_taxid.size();
			for(vector<pair<int,int> >::const_iterator itr2 = cand_taxid.begin();itr2!=cand_taxid.end();++itr2){
				if(sp_score[*itr2].size()==0)
					sp_score[*itr2].resize(kN);
				sp_score[*itr2][lamdaID] += unit;
			}
		}
	}
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compa_read(int ReadLen,const int &position1,const int &position2,int &match,int &mismatch,const ULLN &read1,const ULLN &read2)
{
	if(position1 > position2)
	{
		int k=position1-position2;
		ULLN differ= (read2>>(2*k));
//		differ ^= (read1 & AndMask[ReadLen-k]);
		differ ^= read1; differ.keeptopk(ReadLen-k);
		mismatch = differ.non0base();
		match = ReadLen-k-mismatch;
		//Note: can change == to ^, so that can return a value;
	}
	else if(position1 < position2)
	{
		int k=position2-position1;
		ULLN differ= (read1>>(2*k));
//		differ ^= (read2 & AndMask[ReadLen-k]);
		differ ^= read2; differ.keeptopk(ReadLen-k);
		mismatch = differ.non0base();
		match = ReadLen-k-mismatch;
		//Note: can change == to ^, so that can return a value;
	}
	else
	{
		ULLN differ= read1;
		differ ^= read2;
		mismatch = differ.non0base();
		match = ReadLen-mismatch;
	}
}

void Groupctg(USet& uset,string filename){
    ifstream fin(filename.c_str());
	int from,to;
	while(!fin.eof()){
	     fin>>from>>to;
		 
		 uset.Union(from,to);
	}
	fin.close();
}

void MergeCtgs(const ReadsClass& Reads, const ContigsClass& Ctgs, USet& uset,string groupfile){
     ifstream group(groupfile.c_str());
	 int rid, type,newtype;
	 int cid1,cid2;
	 while(!group.eof()){
         group>>rid>>type>>newtype;
         cid1 = Reads.MatchId[rid];
		 cid2 = Reads.MatchId[type];
     	 if(cid1 >= 0 && cid2 >= 0){
		     if(uset.find(cid1) != uset.find(cid2)){
			     uset.Union(cid1,cid2);
                             cout<<rid<<" "<<cid1<<" "<<cid2<<endl;
                     }
		 
		 }
	
	 
	 
	 
	 }
     group.close();
}
void MergeReads(const KmerNodeAloc& NodePool, const ReadsClass& Reads, USet& uset, int CtgNum)
{
	////////////////////const variables
//	const int StrLen = 50;
	const int LARGELENGTH = 8000;
	const int MINSCORE = 50;
	const int FragLen = 1000;
	const double Cdf_Cl = 1.65;

	////////////////////const variables
	int confidence[PARA_READ+1];
	{
		double r = 0.01;
		double p = 1+r*r-5*r/3;
		double var = p*(1-p);
		for(int i=0;i<=PARA_READ;++i)
			confidence[i] = i*p-sqrt(var*i)*Cdf_Cl;
	}
	if(INTEST)
	{
		cerr << "confidence levels: " << endl;
		for(int i=0;i<=PARA_READ;++i)
			cerr << i << ":\t" << confidence[i] << endl;
	}

	////////////////////sort KmerNodes
	int NodeNum = NodePool.getNodeNum();
	KmerNode** SortedNodes = new KmerNode*[NodeNum];
	for(int i=0; i<NodeNum; ++i)
		SortedNodes[i] = NodePool.getRef(i);
	sort(SortedNodes, SortedNodes+NodeNum, cmpReads);
	////////////////////merge nodes only with reads
	omp_lock_t type_lock;
	omp_init_lock(&type_lock);
	int ReadLen = Reads.ReadLen;
#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<NodeNum; ++i)
	{
		NodeIDPosi* vect = SortedNodes[i]->myvector;
		int ctgnum = SortedNodes[i]->CtgNum;
		int sub_size = SortedNodes[i]->VSize - ctgnum;
		int* indexes = new int[sub_size];
		int* positions = new int[sub_size];
		bool* isRev = new bool[sub_size];
		for(int j=0;j<sub_size;++j)
		{
			indexes[j] = vect[j+ctgnum].id;
			positions[j] = vect[j+ctgnum].posi >> 2;
			isRev[j] = vect[j+ctgnum].isRev();
		}
	
		ULLN read1,read2;
		for(int i2=0;i2<sub_size;++i2)
		{
			int indexes1 = indexes[i2];
			int uidx1 = indexes1 + CtgNum;

			if(isRev[i2])	read1 = (*Reads.R2[indexes1]);
			else	read1 = (*Reads.R1[indexes1]);

	/*		int j=i+1;
			while(j<sub_size && positions[j]-positions[i2]<=ReadLen-StrLen)
				++j;*/
			for(int j=i2+1;j<sub_size;++j)
			{
				int indexes2 = indexes[j];
				int uidx2 = indexes2 + CtgNum;
				if(uset.find(uidx1)==uset.find(uidx2))
					continue;
			/*	if(positions[j]-positions[i2]>ReadLen-StrLenLowCover)
					break;*/
				/////////////////////////////////////////////////////////////////////////////////////
	//			if((uset.getReadNum(uidx1)<LARGESIZE && uset.getReadNum(uidx2)<LARGESIZE )||(uset.getReadNum(uidx1)<Frag) || uset.getReadNum(uidx2)<Frag)
				if((uset.getCtgLen(indexes1)<LARGELENGTH && uset.getCtgLen(indexes2)<LARGELENGTH )||(uset.getCtgLen(indexes1)<FragLen) || uset.getCtgLen(indexes2)<FragLen)
				{
					if(isRev[j])	read2 = (*Reads.R2[indexes2]);
					else read2 = (*Reads.R1[indexes2]);
					int match,mismatch;
					compa_read(ReadLen, positions[i2],positions[j],match,mismatch,read1,read2);
					assert(match >= 32);
				/////////////////bin small ones
			//		if(match-StrLenLowCover>=0 && (match-StrLenLowCover >= confidence[match+mismatch-StrLenLowCover]))
					if(match>=MINSCORE && (match-MINSCORE >= confidence[match+mismatch-MINSCORE]))
					{
						omp_set_lock(&type_lock);
						uset.Union(uidx1,uidx2);
						omp_unset_lock(&type_lock);
					}
				}
			}
		}
		delete[]indexes;delete[]positions;delete[]isRev;
	}
	omp_destroy_lock(&type_lock);
	delete[]SortedNodes;
	///////////////////////////////////////////////////////////////////////
}

void MergeCtgReads(const KmerNodeAloc& NodePool, const ReadsClass& Reads, const ContigsClass& Ctgs, USet& uset)
{
	////////////////////const variables
//	const int StrLenLowCover = 32;
//	const int StrLen = 50;
	const int LARGELENGTH = 8000;
	const int FragLen = 1000;
	const int MINSCORE = 50;
	const double Cdf_Cl = 1.65;

	////////////////////const variables
	int confidence[PARA_READ+1];
	{
		double r = 0.01;
		double p = 1+r*r-5*r/3;
		double var = p*(1-p);
		for(int i=0;i<=PARA_READ;++i)
			confidence[i] = i*p-sqrt(var*i)*Cdf_Cl;
	}
	/////////////////////////////////////////
	int CtgNum = Ctgs.CtgNum;

	////////////////////sort KmerNodes
	int NodeNum = NodePool.getNodeNum();
	KmerNode** SortedNodes = new KmerNode*[NodeNum];
	for(int i=0; i<NodeNum; ++i)
		SortedNodes[i] = NodePool.getRef(i);
	sort(SortedNodes, SortedNodes+NodeNum, cmpReads);
	////////////////////merge nodes only with reads
	omp_lock_t type_lock;
	omp_init_lock(&type_lock);
	int ReadLen = Reads.ReadLen;
#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<NodeNum; ++i)
	{
		NodeIDPosi* vect = SortedNodes[i]->myvector;
		int ctgnum = SortedNodes[i]->CtgNum;
		int sub_size = SortedNodes[i]->VSize;
		int* indexes = new int[sub_size];
		int* positions = new int[sub_size];
		bool* isRev = new bool[sub_size];
		for(int j=0;j<sub_size;++j)
		{
			indexes[j] = vect[j].id;
			positions[j] = vect[j].posi >> 2;
			isRev[j] = vect[j].isRev();
		}
		for(int i2=0;i2<ctgnum;++i2)
		{
			int indexes1 = indexes[i2];
			assert(indexes1<CtgNum);
			BaseStr* ctg1 = Ctgs.contigs[indexes1];

	//		ULLN read2;
			for(int j=ctgnum;j<sub_size;++j)
			{
				int indexes2 = indexes[j];
				assert(indexes2<Reads.ReadNum);
				int uidx2 = indexes2 + CtgNum;
				if(uset.find(indexes1)==uset.find(uidx2) || Reads.Score[indexes2]<MINSCORE)
					continue;
	/*			if(positions[j]-positions[i2]>ReadLen-StrLenLowCover)
					break;*/
				/////////////////////////////////////////////////////////////////////////////////////
				if((uset.getCtgLen(indexes1)<LARGELENGTH && uset.getCtgLen(uidx2)<LARGELENGTH )||(uset.getCtgLen(indexes1)<FragLen) || uset.getCtgLen(uidx2)<FragLen)
				{
	/*				if(isRev[j])	read2 = (*Reads.R2[indexes2]);
					else read2 = (*Reads.R1[indexes2]);*/
					int match,mismatch;
					if(isRev[i2])
					{
						if(isRev[j])
							ctg1->checkRead(positions[i2],(*Reads.R1[indexes2]),ReadLen-1-positions[j],ReadLen,match,mismatch);
						else
							ctg1->checkRead(positions[i2],(*Reads.R2[indexes2]),ReadLen-1-positions[j],ReadLen,match,mismatch);
					}
					else
					{
						if(isRev[j])
							ctg1->checkRead(positions[i2],(*Reads.R2[indexes2]),PARA_KMER+positions[j]-1,ReadLen,match,mismatch);
						else
							ctg1->checkRead(positions[i2],(*Reads.R1[indexes2]),positions[j]+PARA_KMER-1,ReadLen,match,mismatch);
					}
/*					if(INTEST && match < 32)
					{
						omp_set_lock(&type_lock);
						cerr << dec << sub_size << '\t' << ctgnum << '\t' << ((vect[i2].posi)&0x3) << '\t' << ((vect[j].posi)&0x3) << endl;
						cerr << SortedNodes[i]->kmer << endl;
						cerr << dec << "match1 < 32:\t" << match << '\t' <<  indexes1 <<'\t' << positions[i2]<<'\t' << indexes2<<'\t'<<positions[j]+PARA_KMER-1 << endl;
						cerr << ctg1->subStr(positions[i2]-positions[j]+1-PARA_KMER, ReadLen) << endl;
						cerr << ctg1->subStr(positions[i2]+1-PARA_KMER, ReadLen) << endl;
						omp_unset_lock(&type_lock);
					}*/
					assert(match >= 32);
					if(match>=MINSCORE && (match-MINSCORE >= confidence[match+mismatch-MINSCORE]))
					{
						omp_set_lock(&type_lock);
						uset.Union(indexes1,uidx2);
						omp_unset_lock(&type_lock);
					}
				}
			}
		}
		delete[]indexes;delete[]positions;delete[]isRev;
	}
	omp_destroy_lock(&type_lock);
	delete[]SortedNodes;
	///////////////////////////////////////////////////////////////////////
}

void MergePE(const ReadsClass& Reads, USet& uset, AccTester& acc_tester)
{
	/////////////////////////////////////////get pair-end graph (PEGraph);
	map<unsigned, map<unsigned, unsigned> >PEGraph;
	SparseMatrix SPMTX;
	int SPMTXSIZE = SPMTX.Size;
	omp_lock_t* hash_lock = new omp_lock_t[SPMTXSIZE];
	for(int i=0;i<SPMTXSIZE;++i)
		omp_init_lock(&hash_lock[i]);

#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<Reads.TotalNum;i+=2)
	{
		int idx1 = Reads.MatchId[i];
		int idx2 = Reads.MatchId[i+1];
		if(idx1!=idx2 && (idx1>=0 && idx2 >= 0))
		{
			unsigned id = SPMTX.hash(idx1,idx2);
			omp_set_lock(&hash_lock[id]);
			SPMTX.insert(idx1,idx2);
			SPMTX.insert(idx2,idx1);
			omp_unset_lock(&hash_lock[id]);
		}
	}
	for(int i=0;i<SPMTXSIZE;++i)
		omp_destroy_lock(&hash_lock[i]);
	delete[] hash_lock;
	SPMTX.toNeighbor(PEGraph);
	////////////////////////////////////////////

	//////////////////////////////////////merge pair-end reads
	map<unsigned,unsigned> weightsum;
	map<unsigned,unsigned> friends;
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		unsigned tsum = 0;
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
			tsum += itr2->second;
		weightsum[itr1->first] = tsum;
		friends[itr1->first] = itr1->second.size();
	}

	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))>=0.25)
				uset.Union(itr1->first, itr2->first);

	if(INTEST)acc_tester.calAcc(uset);
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
		{
			if(uset.getCtgLen(itr1->first)>10000 && uset.getCtgLen(itr2->first)>10000)
				continue;
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))>=0.20)
				uset.Union(itr1->first, itr2->first);
		}
	}
	if(INTEST)acc_tester.calAcc(uset);
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
		{
			if(uset.getCtgLen(itr1->first)>10000 && uset.getCtgLen(itr2->first)>10000)
				continue;
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))>=0.15)
				uset.Union(itr1->first, itr2->first);
		}
	}
	//////////////////////////////////////////////////////////////////////////////
	////////////////
	/*
	cerr << "start to get purity of each ctg." << endl;
	vector<double> purity(acc_tester.CtgNum,0);
	vector<int>majorid(acc_tester.CtgNum,0);
	for(int i=0;i<acc_tester.CtgNum;++i)
	{
		int tmajor = 0, tmax = 0, sum = 0;
		for(int j=0;j<acc_tester.GenomeNum;++j)
		{
			int t = acc_tester.CtgComp[i][j];
			sum += t;
			if(t > tmax)
			{
				tmajor = j;
				tmax = t;
			}
		}
		purity[i] = tmax/(double)sum;
		majorid[i] = tmajor;
	}
	////////////////
	cerr << "pe-graph\t" << endl;
	for(map<unsigned, map<unsigned, unsigned> >::const_iterator itr1 = PEGraph.begin(); itr1 != PEGraph.end(); ++itr1)
	{
		for(map<unsigned, unsigned>::const_iterator itr2 = itr1->second.begin(); itr2 != itr1->second.end(); ++itr2)
		{
			if((itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first])))<0.1)
				continue;

			if(majorid[itr1->first]==majorid[itr2->first])
				cerr << "true\t";
			else
				cerr << "false\t";
			cerr << itr2->second/(double)weightsum[itr1->first] << '\t';
			cerr << itr2->second/(double)weightsum[itr2->first] << '\t';
			cerr << friends[itr1->first] <<'\t'<<friends[itr2->first]<<'\t';
			cerr << purity[itr1->first] << '\t' << purity[itr2->first] << '\t';
			cerr << itr2->second << '\t';
			cerr << weightsum[itr1->first] << '\t';
			cerr << weightsum[itr2->first] << '\t';
			cerr << uset.getCtgLen(itr1->first) << '\t';
			cerr << uset.getCtgLen(itr2->first) << '\t';
			cerr << itr2->second/(double)(max(weightsum[itr1->first],weightsum[itr2->first]))<<'\t';
			cerr << endl;
		}
	}
	cerr << "end of pe-graph\t" << endl;
	*/
}

void MergeAsStep1(const ContigsClass& Ctgs,const ReadsClass& Reads, KmerNodeAloc& NodePool, USet& uset, AccTester& acc_tester, int pair_end_merge_threshold)
{
	///////////////////////testing
	if(INTEST)
	{
		map<int,int>halfmapscore;
		int unmatch = 0;
		int bothunmapped = 0;
		for(int i=0;i<Reads.TotalNum;i+=2)
		{
			int idx1 = Reads.MatchId[i];
			int idx2 = Reads.MatchId[i+1];
			if(idx1!=idx2 && (idx1>=0 && idx2 >= 0))
				++unmatch;
			if(idx1<0 && idx2<0)
				++bothunmapped;
			if(idx1 >=0 && idx2<0)
				++halfmapscore[Reads.Score[i+1]];
			if(idx2 >=0 && idx1<0)
				++halfmapscore[Reads.Score[i]];
		}
		cerr << "Mixed mapped pair-end reads: \t" << unmatch/(double)Reads.TotalNum << '\t' << unmatch << '\t' << Reads.TotalNum << endl;
		cerr << "buth unmaped pair-end reads: \t" << bothunmapped/(double)Reads.TotalNum << '\t' << bothunmapped << '\t' << Reads.TotalNum << endl;
		cerr << "half mapped scores:\t"<< endl;
		for(map<int,int>::const_iterator itr=halfmapscore.begin();itr!=halfmapscore.end();++itr)
			cerr << itr->first << '\t' << itr->second << endl;;
	}
	/////////////////////merge pair-end reads
/*	for(int i=0;i<Reads.TotalNum;i+=2)
	{
		if(Reads.MatchId[i]<0 && Reads.MatchId[i+1]<0)
			uset.Union(Ctgs.CtgNum+Reads.toNewId(i+1),Ctgs.CtgNum+Reads.toNewId(i));
	}*/
	if(INTEST)acc_tester.calAcc(uset);
	///////////////////////////////////////////////////////////////////////
	MergePE(Reads, uset, acc_tester);
	if(INTEST)acc_tester.calAcc(uset);
/*	MergeCtgReads(NodePool, Reads, Ctgs, uset);
	if(INTEST)acc_tester.calAcc(uset);
	MergeReads(NodePool, Reads, uset, Ctgs.CtgNum);
	if(INTEST)acc_tester.calAcc(uset);
	*/
	///////////////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////////////
//process nr database
static const string U = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
int nuc2int(char ch){
	switch(ch){
		case 'A':return 0;
		case 'C':return 1;
		case 'G':return 2;
		case 'T':return 3;
		default: return 1<<10;
	}
	return 1<<10;
}
char Nuc2Prot(string ch){
	assert(U.length()==64);
	int prot = (nuc2int(ch[0])<<4)|(nuc2int(ch[1])<<2)|(nuc2int(ch[2])); 
	if(prot < 0 || prot >= 64)return 'X';
	return U[prot];
}
string toProt(string& contig){
	int len = (contig.length())/3;
	char* buf = new char[len+1];
	for(int i=0;i<len;++i)
		buf[i] = Nuc2Prot(contig.substr(i*3,3));
	buf[len] = 0;
	string ans(buf);delete[]buf;
	return ans;
}

void getProtSearch(const ContigsClass& Ctgs,const string ProtIdxPath,vector<map<int,double> >&ans){
	NRBWTs nrbwts(ProtIdxPath);
	cerr << "nrbwts loaded" << endl;
	int N = Ctgs.CtgNum;
	ans.resize(N);
	omp_lock_t type_lock;
	omp_init_lock(&type_lock);
	//////////////////////////////////////////////////////////////
#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<N;++i){
		omp_set_lock(&type_lock);
		cerr <<"processing contig " << i << " with nr. "<< endl;
		omp_unset_lock(&type_lock);
		////////////////////////////////////////////////////////
		string str[6];
		for(int j=0;j<3;++j){
			string tstr = Ctgs.contigs[i]->str.substr(j);
			str[j] = toProt(tstr);
		}
		string strnuc = Ctgs.contigs[i]->str;
		{
			char* rstrnuc = new char[strnuc.length()+1];
			for(int j=0;j<strnuc.length();++j){
				char ch = strnuc[strnuc.length()-1-j];
				switch(ch){
					case 'A':rstrnuc[j]='T';break;
					case 'C':rstrnuc[j]='G';break;
					case 'G':rstrnuc[j]='C';break;
					case 'T':rstrnuc[j]='A';break;
					default:rstrnuc[j]='N';
				}
			}
			rstrnuc[strnuc.length()] = 0;
			for(int j=0;j<3;++j){
				string tstr(rstrnuc+j);
				str[j+3] = toProt(tstr);
			}
			delete[] rstrnuc;
		}
		////////////////////////////////////////////////////////
		map<int,double> MTX[6];
		for(int j=0;j<6;++j){
			string ProtP = str[j];

			set<string>kmers;
			int strend = (ProtP.length()-(PARA_ANNOTATION_NR_KMER-1));
			for(int l=0;l<strend;++l)
				kmers.insert(ProtP.substr(l,PARA_ANNOTATION_NR_KMER));
			for(set<string>::const_iterator itr = kmers.begin();itr!=kmers.end();++itr){
				set<pair<int,int> > gross_cand_taxid1;
				nrbwts.search(*itr, gross_cand_taxid1);

				set<int> stid_cand;
				for(set<pair<int,int> >::const_iterator itr2=gross_cand_taxid1.begin();itr2!=gross_cand_taxid1.end();++itr2){
					if(!INTEST)
						stid_cand.insert(itr2->second);
				}
				double curScore = 1.0/stid_cand.size();
				for(set<int>::const_iterator itr2=stid_cand.begin();itr2!=stid_cand.end();++itr2)
					MTX[j][*itr2] += curScore;
			}
		}
		map<int,double>maxMTX;
		for(int j=0;j<6;++j){
			for(map<int,double>::iterator itr=MTX[j].begin();itr!=MTX[j].end();++itr){
				maxMTX[itr->first] = max(maxMTX[itr->first],itr->second);
			}
		}
//		omp_set_lock(&type_lock);
		for(map<int,double>::const_iterator itr=maxMTX.begin();itr!=maxMTX.end();++itr)
			ans[i][itr->first] += itr->second;
//		omp_unset_lock(&type_lock);
	}
}

///////////////////////////////////////////////////////////////////////
vector<double> getTaxonConfidence(vector<pair<double,double> >normalPara,vector<int> normalPrior,int level,double val){
	vector<double> ans(8,0);
	if(normalPara.size()==0 || normalPara.size() <= level || fabs(normalPara[level].first) < 1e-9)return ans;
	double totalscore = 0;
	for(int i=0;i<normalPara.size();++i){
		double currentConfi = 0;
		if(fabs(normalPara[i].first)>=1e-9)
			currentConfi = gaussian(normalPara[i].first, normalPara[i].second, val)*normalPrior[i];
		ans[i] = currentConfi;
		if(i>0)ans[i] += ans[i-1];
//		if(currentConfi < 0)cout << "confi->0: "<<currentConfi<<'\t' << normalPara[i].first<<'\t'<<normalPara[i].second << '\t'<<val << endl;
	}
	totalscore = ans[normalPara.size()-1];
	if(totalscore<=1e-100)return ans;
	for(int i=0;i<normalPara.size();++i)
		ans[i] /= totalscore;
	return ans;
}

void getScoreV(int** GenoDBTaxo, const vector<int>&CtgId, const int taxoLevel, int& AnsLevel, double& AnsScore)
{
	map<int,int> TaxoStatis;
	for(vector<int>::const_iterator itr=CtgId.begin();itr!=CtgId.end();++itr)
		++TaxoStatis[GenoDBTaxo[*itr][taxoLevel]];
	double ans = 0;
	for(map<int,int>::const_iterator itr=TaxoStatis.begin();itr!=TaxoStatis.end();++itr)
	{
		double p = itr->second/(double)CtgId.size();
		ans += p*log(p);
	}
	int maxId = 0,maxC = 0;
	for(map<int,int>::const_iterator itr=TaxoStatis.begin();itr!=TaxoStatis.end();++itr)
		if(itr->second > maxC)
		{
			maxId = itr->first;
			maxC = itr->second;
		}
	AnsLevel = maxId;
	AnsScore = -ans;
}

void getTaxo4Clust(int** GenoDBTaxo, const vector<int>&CtgId, double Thresh, int& AnsLevel, double& AnsScore)
{
	for(int i=6;i>=0;--i)
	{
		getScoreV(GenoDBTaxo, CtgId, i, AnsLevel, AnsScore);
		if(AnsScore >= Thresh)
			return;
	}
}

double getScoreEntropy(const map<int,long long> &TaxoComp)
{
	long long sum = 0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
		sum += itr->second;
	if(sum==0) return 1e9;
	double ans = 0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
	{
		double tp = (itr->second)/(double)sum;
		ans += tp*log(tp);
	}
	return ans/(log(2));
}
double getScoreMax(const map<int,long long> &TaxoComp)
{
	long long sum = 0,max=0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
	{
		sum += itr->second;
		if(itr->second > max)
			max = itr->second;
	}
	if(sum==0)return 1e9;
	return max/(double)sum;
}
double getScoreMax2(const map<int,long long> &TaxoComp)
{
	long long max1=0,max2=0,sum=0;
	for(map<int,long long>::const_iterator itr = TaxoComp.begin(); itr!=TaxoComp.end(); ++itr)
	{
		sum += itr->second;
		if(itr->second > max1)
		{
			max2 = max1;
			max1 = itr->second;
		}
		else if(itr->second > max2)
			max2 = itr->second;
	}
	if(sum==0)return 1e9;
	return (max1-max2)/(double)sum;
}

void outputCluster(int* MatchId, long long ReadN, USet& uset,int* metabest, MCPara& mcpara,string infile)
{
	map<int,vector<unsigned> > R;
	vector<unsigned> Gid(ReadN,0);
	for(long long i=0;i<ReadN;++i)
	{
		int ctgid = MatchId[i];
		if(ctgid<0)
			continue;
		ctgid = uset.find(ctgid);
		long long metaid;
		if((metaid=mcpara.getNewId(ctgid))<0)
			continue;
		R[metabest[metaid]].push_back(i);
	}
	ifstream ifs(infile.c_str());
	vector<string> readtitle(ReadN,"");
	if(!ifs.fail())
	{
		const int BMax = 10000;
		char Buf[BMax];
		for(long long i=0;i<ReadN;++i)
		{
			if(ifs.eof())break;
			ifs.getline(Buf,BMax);
			readtitle[i] = Buf;
			ifs.getline(Buf,BMax);
		}
	}
	
	ofstream ofs((infile+".cluster").c_str());
	for(map<int,vector<unsigned> >::const_iterator itr1 = R.begin();itr1!=R.end();++itr1)
	{
		ofs << "cluster " << itr1->first << '\t' << itr1->second.size() << endl;
/*		ofs << "Taxo:\t";
 *				if(itr1->first >=0 && itr1->first < taxoofclust.size())
 *							for(int i=0;i<7;++i)
 *											ofs << taxoofclust[itr1->first].taxo[L7[i]].taxid << ':' << taxoofclust[itr1->first].taxo[L7[i]].score << '\t';
 *													ofs << endl;*/
		int i = itr1->first;
	    ofs << "Taxo (unknown level) :1";
		ofs << endl;
		for(vector<unsigned>::const_iterator itr2 = itr1->second.begin();itr2!=itr1->second.end();++itr2)
		{
			if(readtitle[*itr2].length()==0)
				ofs << "< read " << (*itr2) << endl;
			else 
				ofs << readtitle[*itr2] << endl;
		}
	}
	ofs.close();
}

void outputResultReads(vector<double>&ClustCoverage,int* MatchId, long long ReadN, USet& uset,int* metabest, MCPara& mcpara,string infile,vector<ClustTaxoInfoClass>& taxoofclust)
{
	map<int,vector<unsigned> > R;
	vector<unsigned> Gid(ReadN,0);
	for(long long i=0;i<ReadN;++i)
	{
		int ctgid = MatchId[i];
		if(ctgid<0)
			continue;
		ctgid = uset.find(ctgid);
		long long metaid;
		if((metaid=mcpara.getNewId(ctgid))<0)
			continue;
		R[metabest[metaid]].push_back(i);
	}
	ifstream ifs(infile.c_str());
	vector<string> readtitle(ReadN,"");
	if(!ifs.fail())
	{
		const int BMax = 10000;
		char Buf[BMax];
		for(long long i=0;i<ReadN;++i)
		{
			if(ifs.eof())break;
			ifs.getline(Buf,BMax);
			readtitle[i] = Buf;
			ifs.getline(Buf,BMax);
		}
	}
	
 	const int L7[] = {4,8,12,17,21,24,28};
	const string TaxonLevelName[] = {"Species","Genus","Family","Order","Class","Phylum","Kingdom","Unknown"};
	ofstream ofs((infile+".annotation").c_str());
	for(map<int,vector<unsigned> >::const_iterator itr1 = R.begin();itr1!=R.end();++itr1)
	{
		ofs << "cluster " << itr1->first << '\t' << itr1->second.size() << endl;
/*		ofs << "Taxo:\t";
		if(itr1->first >=0 && itr1->first < taxoofclust.size())
			for(int i=0;i<7;++i)
				ofs << taxoofclust[itr1->first].taxo[L7[i]].taxid << ':' << taxoofclust[itr1->first].taxo[L7[i]].score << '\t';
		ofs << endl;*/
		int i = itr1->first;
		if(taxoofclust[i].level < 7 && taxoofclust[i].level >= 0)
			ofs << "Taxo ("<<TaxonLevelName[taxoofclust[i].level]<<":" << taxoofclust[i].taxon <<")\t"<<"(Mean Coverage: "<<ClustCoverage[i]<<"):\t";
		else
			ofs << "Taxo (unknown level "<<taxoofclust[i].level<<":" << taxoofclust[i].taxon <<")\t";
		ofs << endl;

		for(vector<unsigned>::const_iterator itr2 = itr1->second.begin();itr2!=itr1->second.end();++itr2)
		{
			if(readtitle[*itr2].length()==0)
				ofs << "< read " << (*itr2) << endl;
			else 
				ofs << readtitle[*itr2] << endl;
		}
	}
	ofs.close();
}

void outputResultContigs(vector<double>&ClustCoverage,ContigsClass&Ctgs,const int VCtgNum,USet& uset,int* metabest, MCPara& mcpara,const string infile,vector<ClustTaxoInfoClass>& taxoofclust){
 	const int L7[] = {4,8,12,17,21,24,28};
	int ClustNum = 0;
	for(int i=0;i<VCtgNum;++i)
		ClustNum = max(metabest[i],ClustNum);
	++ClustNum;
	cerr<<"o test 1:" << VCtgNum<<'\t'<<ClustNum <<  endl;
	vector<vector<int> >CtgIdInCluster(ClustNum+1);
	for(int i=0;i<Ctgs.CtgNum;++i){
		int cid = mcpara.getNewId(uset.find(i));
		if(cid < 0 || metabest[cid]<0)
			CtgIdInCluster[ClustNum].push_back(i);
		else CtgIdInCluster[metabest[cid]].push_back(i);
	}
	cerr<<"o test 2" << endl;
	const string TaxonLevelName[] = {"Species","Genus","Family","Order","Class","Phylum","Kingdom","Unknown"};
	ofstream ofs((infile+".annotation").c_str());
	for(int i=0;i<ClustNum;++i){
		ofs << "cluster " << i << '\t' << CtgIdInCluster[i].size() << endl;
		if(taxoofclust[i].level < 7 && taxoofclust[i].level >= 0)
			ofs << "Taxo ("<<TaxonLevelName[taxoofclust[i].level]<<":" << taxoofclust[i].taxon <<")\t"<<"(Mean Coverage: "<<ClustCoverage[i]<<"):\t";
		else
			ofs << "Taxo (unknown level "<<taxoofclust[i].level<<":" << taxoofclust[i].taxon <<")\t";
		/*
		for(int j=0;j<7;++j)
			ofs << taxoofclust[i].taxo[L7[j]].taxid << ':' << taxoofclust[i].taxo[L7[j]].score << '\t';
		*/
		ofs << endl;
		for(vector<int>::const_iterator itr=CtgIdInCluster[i].begin();itr!=CtgIdInCluster[i].end();++itr){
			ofs << Ctgs.info[*itr] << endl;
		}
	}
	cerr<<"o test 3" << endl;
	ofs << "Unannotated contigs: " << endl;
	for(vector<int>::const_iterator itr=CtgIdInCluster[ClustNum].begin();itr!=CtgIdInCluster[ClustNum].end();++itr)
		ofs << Ctgs.info[*itr] << endl;
	ofs.close();
	cerr<<"o test 4" << endl;
}

void outputReads(const ReadsClass& Reads,string filename){
     ofstream fout(filename.c_str());
     for(int i=0;i<Reads.TotalNum;i++)
	     fout<<i<<" "<<Reads.MatchId[i]<<" "<<endl;
     fout.close();
}

void outputUset(const ContigsClass& Ctgs, USet& uset, string filename){

     ofstream fout(filename.c_str());
     for(int i=0;i<Ctgs.CtgNum;i++){
             fout<<i<<" "<<uset.find(i)<<endl;
     }
     fout.close();
}

void createbin(USet& uset,map<int,int>& bin,string filename){
     ifstream fin(filename.c_str());
     int ctgid,ctgnum;
     while(!fin.eof()){
             fin>>ctgid>>ctgnum;
             if(fin.eof())break;
             bin[uset.find(ctgid)]++;

     }
     fin.close();
    
}
