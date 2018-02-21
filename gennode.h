#include<vector>
#include<map>
#include<algorithm>
string findcommon_seq(vector<string> dirs){
	
	cout<<dir<<endl;
	return ""+(char)dirs.size();
}
string findcommon_node(vector<int>& children,vector<string >&commonstring){
	
	cout<<"common: "<<children.size()<<endl;
	
	return ""+(char)children.size();
}

string get_seq(string dir){
	
	return dir;
}

bool compnode(node* &n1,node* &n2){
	return n1->leafnum > n2->leafnum;
}


class commonnode{
public:
   	node* parent;
	
	vector<commonnode*> children;
	
	int id;
	
	bool isLeaf;
	
	int leafnum;
	
	string dir;
	
}

class pathnode{
public:
    vector<int>path;
	

	
	bool isLeaf = false;
	
	int id;
	
};

class pathtree{
public:
    vector<pathnode*> paths;
	pathnode* root;
	
	
};

class Subphytree{
public:
	vector<node*> commontree; 
	vector<string> commonstring;
	pathtree heavy;
	map<int,string> id_name;

   
	
	commonnode* root;
	int rootid;
	
	void genTree(vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir);
	
    void heavyPath();
	
	void preorder();
	
	void sortchild();
	
	
	~Subphytree(){
		for(int i=0;i<commontree.size();i++)
			delete commontree[i];
		for(int i=0;i<heavy.paths.size();i++)
			delete heavy.paths[i];
	}
	
private:
   void label(commonnode* node,int idx);
   void constructHeavypath(commonnode* node,pathnode* newpath,int idx);
};

//construct tree where each node presents commonstrings among subtrees,return the number of nodes (root id = number of nodes - 1)
void Subphytree::genTree(vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir){

	
	
	
	
	
	
	map<string,commonnode*> tempmap;
	
    for(int i=taxtree.size()-1;i>=0;i--){
        for(map<string,set<string> >::iterator itr = taxtree[i].begin();itr!=taxtree[i].end();itr++){
          
			
			commonnode* parentnode = new commonnode;
			parentnode->isLeaf = false;
			parentnode->leafnum = 0;
			tempmap[itr->first] = parentnode;
			
			for(int j=0;j<itr->second.size();j++){
				if(gbkdir.find(itr->second[j]) != gbkdir.end() ){
					 commonnode* node = new commonnode;
					 node->isLeaf = true;
					 node->leafnum = 1;
					 node->dir = gbkdir[itr->second[j]];
					 tempmap[itr->second[j]] = node; 
					 
					 parentnode->children.push_back(node);
					 node->parent = parentnode;
					 parentnode->leafnum++；
					 
					 
					 
				}
				else{
					 commonnode* node = tempmap[itr->second[j]];
					 parent->children.push_back(node);
					 parentnode->leafnum += node->leafnum;
					 node->parent = parentnode;
				}
			} 
			
			
            
        }
    }
    
	root = tempmap[taxtree[0].begin()->first]; 
	
	preorder();
	
	for(map<string,commonnode*>::iterator itr=tempmap.begin();itr!=tempmap.end();itr++){
		id_name[itr->second->id] = itr->first;
	}
}



/*
int Subphytree::genTree(int num_nodes,vector<map<string,set<string> > >& taxtree,map<string,string>& gbkdir){
    int nodeidx = 0;
	
	
	map<string,int> name_id;
	
	isLeaf = bit_vector(num_nodes,0);
	
	
    for(int i=taxtree.size()-1;i>=0;i--){
        for(map<string,set<string> >::iterator itr = taxtree[i].begin();itr!=taxtree[i].end();itr++){
            id_name[nodeidx] = itr->first;
			name_id[itr->first] = nodeidx;
			if( ){
				if(itr->second.size() < 2){
					commonstring.push_back( get_seq(gbkdir[itr->first]) );
					
				}
				else{
					vector<string> temp;
					for(int j = 0;j<itr->second.size();j++)
						temp.push_back(gbkdir(itr->second[i]));
					commonstring.push_back(findcommon_seq( temp));
				}
				
				
				isLeaf[nodeidx] = 1;
				commontree.push_back(empty);
				
				commonnode newnode;
				newnode.id = nodeidx;
				commontree.push_back(newnode);
			}
			else{
				commonnode internal;
				internal.id = nodeidx;
				if(itr->second.size() < 2){
																														
					internal.children.push_back(name_id[itr->second[0] ]);
					
					childnode = commontree[name_id[itr->second[0] ] ];
					childnode.parent = nodeidx;
					commonstring.push_back( commonstring[ name_id[ itr->second[0] ] ]);
				}
				else{
					
															
					vector<int> children;
					for(int j = 0;j<itr->second.size();j++){
						int childid = name_id[itr->second[i]];
						children.push_back(childid);
					    childnode = commontree[childid];
						childnode.parent = nodeidx;
						internal.children.push_back(childid);
					}
					commonstring.push_back(findcommon_node( children,commonstring));
					
				}
				commontree.push_back(internal);
				
			}
			nodeidx++;
            
        }
    }
    
	return nodeidx; 
}

*/
void Subphytree::sortchild(){
	for(int i=0;i<commontree.size();i++)
		sort(commontree[i]->children.begin(),commontree[i]->children.end(),compnode);
}

void Subphytree::preorder(){
	
	
	label(root,0);
	
}

int Subphytree::label(commonnode* node,int idx){
	
	for(int i=0;i<node->children.size();i++){
		commonnode* child = node->children[i];
	    if(child.isLeaf){
		    commontree.push_back(child);
			child->id = idx;
			idx++;
		}
		else{
			idx = label(child,idx);
			
			idx++;
		}							
	}
	node->id = idx;
	commontree.push_back(node);
	
	return idx;
	
}


void Subphytree::constructHeavypath(commonnode* node,pathnode* newpath,int idx){
	
	
	
	newpath->path.push_back(commonnode->id);
	if(node->isLeaf) return;
	
	bool heavied = false
	for(int i=0;i<node->children.size();i++){
		
		if(!heavied){
			heavied = true;
			constructHeavypath(node->children[i],newpath);
		}
		else{
			idx++;
			
			pathnode* anotherpath = new pathnode;
			anotherpath->id = idx;
			heavy.paths.push_back(anotherpath);
			constructHeavypath(node->children[i],anotherpath);
		}
	}
}
void Subphytree::heavyPath(){
	
	pathnode* newpath = new pathnode;
	newpath->id = 0;
	heavy.paths.push_back(newpath);
	
	constructHeavypath(root,newpath,0);
	
	
		
}







