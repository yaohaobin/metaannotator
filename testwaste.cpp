#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<stack>
#include<algorithm>
#include<sstream>
#include<stdlib.h>


#include <sdsl/suffix_trees.hpp>
#define k 9

using namespace std;
using namespace sdsl;


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


int queryseq(string& seq,csa_wt<>& csa){
  int occur = 0;
  set<string>kmers;
  for(unsigned int i=0;i<seq.length()-k;i++){
    
    string querytext="";
    genquery(querytext,seq.substr(i,k));
    kmers.insert(string);
 
  }
  
  for(set<string>::itearator itr = kmers.begin();itr!=kmers.end();itr++)
    if(locate(csa,querytext).size() > 0)
        occur++;
  return occur;
}

int main(int argc,char* argv[]){
  if(argc < 2)continue;
  csa_wt<>csa;
  load_from_file(csa,argv[1]);
  ifstream seqfile(argv[2]);
  string line;
  for(unsigned int i=0;i<100;i++){
    seqfile>>line;
    seqfile>>line;
    cout<<queryseq(line,csa)<<endl;
  }
  
  seqfile.close();
  
  return 0; 
}
  
  
