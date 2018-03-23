#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<stack>
#include<algorithm>
#include<sstream>
#include<stdlib.h>
#include<string>
#include<omp.h>

#include <sdsl/suffix_trees.hpp>
#define k 12

using namespace std;
using namespace sdsl;

template<class t_csa, class t_pat=typename t_csa::string_type>
uint64_t search(t_csa& csa,  t_pat pat){
   uint64_t lb=0, rb=csa.size()-1;
   backward_search(csa, lb, rb, pat.begin(), pat.end(), lb, rb);
   return rb+1-lb;

}



void genquery(string& querytext,string& query){
   
    //unsigned long insertpos = seqtext.size();
  


    uint8_t text=0,revcomp=0,temp=0;
    int idx = 0;


    for(size_t i=0;i<query.length();i++){



        switch(query[i]){
            //case 'a': text[i] = 1;break;
            case 'A': temp = 0;break;
            case 'C': temp = 1;break;
            case 'T': temp = 2;break;
            case 'G': temp = 3;break;
            default: temp = 4;

        }

                if(temp<4){
                        text= (text<<2) | temp;
                        revcomp = revcomp | temp<<(2*idx);

                        idx++;
                        if(idx>2){

                                querytext+= text<revcomp?text:revcomp + 1;
                                
                                text = 0;
                                revcomp = 0;
                                idx = 0;
                        }
                }



    }

        if(idx !=0){
                querytext+= text<revcomp?text:revcomp + 1;
            
        }
   

}

template<class t_csa>
int queryseq(string& seq,t_csa& csa){
  int occur = 0;
  set<string>kmers;
  for(unsigned int i=0;i<seq.length()-k;i++){
    
    string querytext="";
    string kmer = seq.substr(i,k);
    genquery(querytext,kmer);
    kmers.insert(querytext);
 
  }
  cout<<kmers.size()<<endl;
  for(set<string>::iterator itr = kmers.begin();itr!=kmers.end();itr++){
    uint64_t hits= locate(csa,*itr).size();
    if(hits>10) occur++;
    
  }
  return occur;
}

int main(int argc,char* argv[]){
  if(argc < 2)return 0;
  csa_bitcompressed<>csa;
  load_from_file(csa,argv[1]);
  ifstream seqfile(argv[2]);
  string line;
  vector<string>seqs;
  while(!seqfile.eof() ){
    
    getline(seqfile,line);
    if(line.length() == 0)break;
    getline(seqfile,line);
    
    seqs.push_back(line);
    //cout<<line.length()<<endl;
    //cout<<queryseq(line,csa)<<endl;
  }
  cout<<"seq num: "<<seqs.size()<<endl;
  #pragma omp parallel{
      #pragma omp for
      for(unsigned int i=0;i<seqs.size();i++)
          cout<<queryseq(seqs[i],csa);
  }        
  seqfile.close();
  
  return 0; 
}
  
  
