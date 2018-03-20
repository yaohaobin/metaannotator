#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<stack>
#include<algorithm>
#include<sstream>
#include<stdlib.h>

#include <sdsl/suffix_trees.hpp>

using namespace std;
using namespace sdsl;



void gencommon(vector<string>&seqs){
    cst_sct3<csa_bitcompressed<>> cst;
	string seqtext = "";
	int_vector<8>belong;
    for(unsigned i = 0;i<seqs.size();i++){
        
		//gentext(seqtext,seq);
		seqtext+=seqs[i];
		unsigned long oldsize= belong.size();
		belong.resize(seqtext.length()+1);
		for(unsigned long j = oldsize;j<belong.size();j++)
            belong[j] = i;		
		//seqtext.resize(seqtext.size()+1);
		seqtext+='Z';		
		//seqtext[seqtext.size()-1] = 66;		
	
	}    
	cout<<seqtext.length()<<" "<<seqs.size()<<endl;
	construct_im(cst,seqtext,1);
	cout<<size_in_mega_bytes(cst)<<" "<<cst.csa.sigma<<endl;
	vector<unsigned int>common(seqtext.length(),0);
	
	map<unsigned long,unsigned int>coverpos;
	set<unsigned long>coverseq;
	
	unsigned int cur_depth = 0;
	bool prune = false;
	
	vector<pair<unsigned long,unsigned long> > intervals;
	
	unsigned long lb=1;
	
	for(size_t i = 2;i<cst.lcp.size();i++){
             

             unsigned long endpos = cst.csa[i-1]+ cst.lcp[i] -1;
            
             if(common[endpos] >= cst.lcp[i]){
                    lb = i;
                    cur_depth = 0;
                    continue;
             }
             else
                   common[endpos] = cst.lcp[i];
             
	     if(cur_depth == 0){					 
		     cur_depth = cst.lcp[i];
			 if(cur_depth !=0)
			     lb = i-1;
			 else
			     lb = i;
		 }
		 else{
		     if(cur_depth != cst.lcp[i]){
			     if( i-lb>=seqs.size() )
			         intervals.push_back(make_pair(lb,i-1));
				 cur_depth = cst.lcp[i];
				 lb = i-1;
			 }
		}
                
               
			     
	}
	if(cur_depth>0 && cst.lcp.size()-lb>=seqs.size())
	    intervals.push_back(make_pair(lb,cst.lcp.size()-1));
	cout<<"interval: "<<intervals.size()<<endl;	
	vector<pair<unsigned long,unsigned long> >commonseed;
	for(size_t i=0;i<intervals.size();i++){
	    coverseq.clear();
		for(size_t j=intervals[i].first;j<=intervals[i].second;j++)
		    coverseq.insert(belong[cst.csa[j]]);
		if(coverseq.size() == seqs.size()){
	            commonseed.push_back(intervals[i]);
                    //cout<<cst.csa[intervals[i].first]<<" "<<cst.csa[intervals[i].second]<<endl;
                }
	}
	cout<<"commonseed: "<<commonseed.size()<<endl;
	intervals.clear();
	bool repeat = false;

        
        	
	for(size_t i=0;i<commonseed.size();i++){
	    prune = false;
		cur_depth = cst.lcp[commonseed[i].first + 1];
	    for(size_t j=commonseed[i].first;j<=commonseed[i].second;j++){
		    unsigned long endpos = cst.csa[j] + cur_depth - 1;
			if(cur_depth < common[endpos] ) 
			    prune = true;
			else{
			    common[endpos] = cur_depth;
			}
		}
		if(!prune){
		    intervals.push_back(commonseed[i]);
                    //cout<<commonseed[i].first<<" "<<commonseed[i].second<<" "<<cur_depth<<endl;
                }
	}


		
	cout<<"final inter: "<<intervals.size()<<endl;	
    	
	
      


}

void extract(vector<string>&dirs,string prefix){
    
	//string prefix(argv[2]);
	
	vector<string>seqs;
	for(unsigned i = 0;i<dirs.size();i++){
	        string filename = dirs[i];
		ifstream seqfile(filename.c_str());
		string line;
		getline(seqfile,line);
		string seq = "";
		while(seqfile>>line)
		    seq+=line;
		seqfile.close();
		
		seqs.push_back(seq);
		cout<<seq.length()<<endl;
	}
	
	
	gencommon(seqs);
	
	

	


}


