#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<stack>
#include<algorithm>
#include<sstream>
#include<stdlib.h>
#include"gennode.h"

#include <sdsl/suffix_trees.hpp>
#define N 2

using namespace std;
using namespace sdsl;

int read_old(int_vector<8>& text,string filename){
    
    ifstream fin(filename.c_str());
    string line;
    getline(fin,line);
    string seq="";
    string segment;
    int num = 0;
    while(fin>>segment)
        seq+=segment;
    fin.close();
    long start=text.size();
    text.resize(text.size()+seq.length()+1);
    for(size_t i=0;i<seq.length();i++){
     
      
        
        switch(seq[i]){
            //case 'a': text[i] = 1;break;
            case 'A': text[start+i]=1;break;
            case 'C': text[start+i]=2;break;
            case 'G': text[start+i]=3;break;
            case 'T': text[start+i]=4;break;
            default: text[start+i] = 5;num++;
            
        }
      

    }
    text[text.size()-1] = 6;
    
    return num;

}


int read(vector<unsigned char>& seqtext,string filename){
    
    ifstream fin(filename.c_str());
    string line;
    getline(fin,line);
    string seq=""; 
    string segment; 
    int num = 0;
    while(fin>>segment)
        seq+=segment;
    fin.close();
  
    
	unsigned char text=0,revcomp=0,temp=0;
	int idx = 0;
	
	
    for(size_t i=0;i<seq.length();i++){
         
         
        
        switch(seq[i]){
            //case 'a': text[i] = 1;break;
            case 'A': temp = 0;break;
            case 'C': temp = 1;break;
            case 'T': temp = 2;break;
            case 'G': temp = 3;break;
            default: temp = 4;num++;

        }
		
		if(temp<4){
			text= (text<<2) | temp;
			revcomp = revcomp | temp<<(2*idx);
			
			idx++;
			if(idx>2){
				
				seqtext.push_back( (text<revcomp?text:revcomp) + 1);
				
				text = 0;
				revcomp = 0;
				idx = 0;
			}
		}
		
      

    }
	
	if(idx !=0)
		seqtext.push_back( (text<revcomp?text:revcomp) + 1);
	
    seqtext.push_back(65);

    return num;

}


int gentext(int_vector<8>& seqtext,string& seq){
    int num = 0;
    unsigned long insertpos = seqtext.size();
    seqtext.resize(seqtext.size()+seq.size()+1);


    uint8_t text=0,revcomp=0,temp=0;
    int idx = 0;


    for(size_t i=0;i<seq.length();i++){



        switch(seq[i]){
            //case 'a': text[i] = 1;break;
            case 'A': temp = 0;break;
            case 'C': temp = 1;break;
            case 'T': temp = 2;break;
            case 'G': temp = 3;break;
            default: temp = 4;num++;

        }

                if(temp<4){
                        text= (text<<2) | temp;
                        revcomp = revcomp | temp<<(2*idx);

                        idx++;
                        if(idx>2){

                                seqtext[insertpos] =text<revcomp?text:revcomp + 1;
                                insertpos++;
                                text = 0;
                                revcomp = 0;
                                idx = 0;
                        }
                }



    }

        if(idx !=0){
                seqtext[insertpos] = text<revcomp?text:revcomp + 1;
            insertpos++;
        }
    seqtext[insertpos] = 65;
    if(seqtext.size()>(insertpos+1) )
        seqtext.resize(insertpos+1);


   return num;
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

int read_once(int_vector<8>& seqtext,string filename){
    
    ifstream fin(filename.c_str());
    string line;
    getline(fin,line);
    string seq=""; 
    string segment; 
    int num = 0;
    while(fin>>segment)
        seq+=segment;
    fin.close();
    return gentext(seqtext,seq);

    

}


void split(string& str, char delim,vector<string>& result)
{
    std::size_t current, previous = 0;
    current = str.find(delim);
    while (current != std::string::npos) {
        result.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    if (previous != str.length())
        result.push_back(str.substr(previous, current - previous));
}
void uniqtax(vector<string>& tax ){
    
    for(unsigned i=1;i<tax.size();i++){
        stringstream ss;
        ss<<tax[i]<<i;
        tax[i]=ss.str();
     
    }
    
}


void gencommon(string& commonstr,string lenfilename,string commonfilename){
    ifstream lenfile(lenfilename.c_str());
    string gi;
    int len;
    map<string,int>gi_len;
    int totallen = 0;
    while(lenfile>>gi>>len){
        gi_len[gi] = len+1;
        totallen += len+1;
    }
    lenfile.close();
    int seqnum = gi_len.size();
    vector<int> common(totallen+1,0);
    ifstream commonfile(commonfilename.c_str());
    float mean = 0;
    int num = 0;    
    int depth,lb,rb,lstart,rstart;
    while(commonfile>>depth>>lb>>rb>>lstart>>rstart){
        if (depth > common[lstart+depth-1])
            common[lstart+depth-1] = depth;
    }

    commonfile.close();
    string output(commonfilename.c_str());
    output += ".clean";

    ofstream commonout(output.c_str());
    for(int i=0;i<common.size();i++){
        if (common[i]){
            commonout<<i<<" "<<common[i]<<endl;
            mean = (mean*num + common[i])/(num+1);

            num++;
        }
    }
    cout<<mean<<" "<<num<<endl;
    commonout<<mean<<" "<<num<<endl;
    commonout.close();
  
    

}

int loadtree(string dbfile,vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir,string prefix){
    string taxcode;
    taxtree.resize(6);
    ifstream db(dbfile.c_str());
    string line;
    int nodenum = 0;
    int normalseq = 0;

    long totalsize = 0;
    
    
    int_vector<8>seq;
    while(!db.eof()){
        getline(db,line);
        if(line.length() == 0)break;
        vector<string> info;
        split(line,'\t',info);
        int complete = stoull(info[3]);
        if (complete <0)continue; 

        string filename = prefix+'/'+info[5]+'/'+info[0]+".fa";
        gbkdir[info[1]] = filename;
        

        vector<string> tax;
        //cout<<info.size()<<endl;
        split(info[2],';',tax);
        if(tax.size() == 0)continue;
        if(tax[0] !="Bacteria") continue;
        uniqtax(tax);
        tax.push_back(info[4]);
        
        tax.push_back(info[1]);
      
            

        int unnormal = read_once(seq,filename);
        //if (unnormal>0)cout<<filename<<" "<<unnormal<<endl;
            
        cout<<seq.size()<<endl;
        
      

        int i=0; 
        if (tax.size()-1 > taxtree.size() ) taxtree.resize(tax.size()-1);
        for(vector<string>::iterator itr=tax.begin();(itr+1)!=tax.end();itr++){
            if ( (*itr).length() == 0)
                cout<< line <<endl;    
            if (taxtree[i].find(*itr) != taxtree[i].end())
                taxtree[i][*itr].insert(*(itr+1));
            else{

                set<string> newparent;
                newparent.insert(*(itr+1));
                taxtree[i][*itr] = newparent;
            }
            
            i++;
        }
   }
   
   
   int_vector<8>seqforsa;
   seqforsa.resize(seq.size());
   for(unsigned long i=0;i<seq.size();i++)
       seqforsa[i] = seq[i];
   
   cout<<(int)seq[2345]<<endl;
  
      
   csa_bitcompressed<> csa;
   for(unsigned long i = 0;i<seq.size();i++)
   {
        if( (int)seq[i] == 0)cout<<i<<endl;
        break;
   }
   construct_im(csa,seq);
   cout<<size_in_mega_bytes(csa)<<endl;
   store_to_file(csa,"bitcomplete.csa");  
   return 0;
   string query;
   cin>>query;
   string querytext="";
   genquery(querytext,query);
   cout<<locate(csa,querytext ).size()<<endl;
   
}



void genTree(vector<pair<int,vector<int> > >& commontree,vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir){
    int nodeidx = 0;
    for(int i=taxtree.size()-1;i>=0;i--){
        for(map<string,set<string> >::iterator itr = taxtree[i].begin();itr!=taxtree[i].end();itr++){
            if(itr->second.size() < 2)continue;
                        
        }
    }

}

int main(int argc,char* argv[]){
    //string commonstr="";
    //string lenfilename(argv[1]);
    //string commonfilename(argv[2]);
    vector<map<string,set<string> > > taxtree;
    map<string,string>gbkdir;
    taxtree.resize(6);
    loadtree(argv[1],taxtree,gbkdir,argv[2]);
    return 0;
    for(int i=0;i<taxtree.size();i++)
        cout<<taxtree[i].size()<<endl;

    int level = taxtree.size()-1;   
    cout<<"last"<<endl;
    for(map<string,set<string> >::iterator itr = taxtree[level].begin();itr !=taxtree[level].end();itr++){
        cout<<itr->first<<" "<<itr->second.size()<<endl;
        
    }   
    
    cout<<"node"<<endl;
    Subphytree indextree;
    indextree.genTree(taxtree,gbkdir);
    string prefix(argv[2]);
    indextree.common(prefix);         
    return 0;



}


        







