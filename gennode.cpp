#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<sstream>
#include"gennode.h"
#include <sdsl/suffix_trees.hpp>
using namespace std;
using namespace sdsl;



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

int loadtree(string dbfile,vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir){
    string taxcode;
    taxtree.resize(6);
    ifstream db(dbfile.c_str());
    string line;
    int nodenum = 0;
    while(!db.eof()){
        getline(db,line);
        if(line.length() == 0)break;
        vector<string> info;
        split(line,'\t',info);
        int complete = stoull(info[3]);
        if (complete <0)continue; 
        gbkdir[info[1]] = info[5];      
        vector<string> tax;
        //cout<<info.size()<<endl;
        split(info[2],';',tax);
        if(tax.size() == 0)continue;
        if(tax[0] !="Bacteria") continue;
        uniqtax(tax);
        tax.push_back(info[4]);
        
        tax.push_back(info[1]);
      
         
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
    loadtree(argv[1],taxtree,gbkdir);
     
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
    return 0;



}


        







