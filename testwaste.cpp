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
#define exk 12
using namespace std;
using namespace sdsl;

template<class t_csa, class t_pat=typename t_csa::string_type>
uint64_t search(t_csa& csa,  t_pat pat,t_pat ext){
   uint64_t lb=0, rb=csa.size()-1;
   uint64_t last_lb = lb,uint64_t last_rb = rb;   
   if (backward_search(csa, lb, rb, pat.begin(), pat.end(), lb, rb)>0){
        for (auto it=ext.end(); it != ext.begin() and lb <= rb;) {
            --it;
            if (backward_search(cst.csa, lb, rb, (typename t_cst::char_type)*it, last_lb, last_rb) > 0) {
                lb = last_lb;
                rb = last_rb;
            }
        }
   }
   
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
  map<string,string>kmers;
  for(unsigned int i=0;i<seq.length()-k-exk;i++){
    
    string querytext="",queryext="";
    string ext = seq.substr(i,exk);
    string kmer = seq.substr(i+exk,k);
    genquery(querytext,kmer);
    genquery(queryext,ext);
    kmers[querytext] = extension;
 
  }
  //cout<<kmers.size()<<endl;
  for(map<string,string>::iterator itr = kmers.begin();itr!=kmers.end();itr++){
    uint64_t hits= search(csa,itr->first,itr->second);
    if(hits>5) occur++;
    
  }
  
  return occur;
}

int main(int argc,char* argv[]){
  if(argc < 2)return 0;
  csa_wt<>csa;
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

  
  int_vector<16> result(seqs.size(),0);
  
  #pragma omp parallel for
  
      
      for(unsigned int i=0;i<seqs.size();i++){
           result[i] = queryseq(seqs[i],csa);
           //cout<<result[i]<<endl;
     }
           
           
  unsigned long classify = 0;
  for(unsigned int i=0;i<result.size();i++)
      if(result[i] != 0)
          classify++;
  cout<<"classify: "<<classify<<endl; 
  return 0;
  int_vector<16>common(csa.size(),0);
  cout<<size_in_mega_bytes(common);
  for(unsigned long int i = 0;i<csa.size();i++){
     common[i] = (uint16_t)csa[i];
  }       
  seqfile.close();
  
  return 0; 
}
  
  
