#include "ExoPieElement/TreeMaker/interface/baseTree.h"

baseTree::baseTree(std::string identifier, TTree* tree){
  identifier_  = identifier; 
  tree_ = tree;
}


void baseTree::AddBranch(bool* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/O").c_str());
}


void baseTree::AddBranch(int* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}


void baseTree::AddBranch(unsigned int* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/i").c_str());
}


void baseTree::AddBranch(unsigned long* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/l").c_str());
}


void baseTree::AddBranch(float* x, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),x,(brName+"/F").c_str());
}


void baseTree::AddBranch(std::vector<bool>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(), vec);
}


void baseTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),vec);
}


void baseTree::AddBranch(std::vector<float>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),vec);
}


void baseTree::AddBranch(std::vector<std::string>* vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(), vec);
}

void baseTree::AddBranch(TClonesArray** vec, std::string name){
  std::string brName = identifier_+name;
  tree_->Branch(brName.c_str(),"TClonesArray", vec,32000,0);
}

void baseTree::AddBranch(std::vector<std::vector<int> >* vec, std::string name){
   std::string brName = identifier_+name;
   tree_->Branch(brName.c_str(), vec);
}


void baseTree::AddBranch(std::vector<std::vector<float> >* vec, std::string name){
    std::string brName = identifier_+name;
   tree_->Branch(brName.c_str(), vec);
}

