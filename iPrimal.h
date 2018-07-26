
///////////////////////////////////////////////////////////////////////////////
// File name: iPrimal.h
// This file defines the main algorithm of IMPROVED primal assignment Method.
// The time complexity is improved to O(n^3 lg n), via merging multiple steps
// into single stages of swap loop searching. 
// Min-heaps and trees are used to maintain searching results.
// Lantao Liu, Jan 17, 2011
///////////////////////////////////////////////////////////////////////////////

#ifndef PRIMAL_H
#define PRIMAL_H

#include "Define.h"
#include "Assignment.h"
#include <math.h>
#include <set>

template<typename T> class less_comp; //defined later


class Primal{
public:
  Primal(){ cerr<<"Not allowed."<<endl; exit(0); }
  Primal(mat& m){ 
	row_size = m.size();
	col_size = m[0].size();
	assert(row_size <= col_size);  //currently assume they are equal 
	orig_matrix = m;
	oprt_matrix = m;
  }
  ~Primal(){}

  //some accessors
  inline void ClearExploited(void){ 
	exploited_rows.clear(); exploited_cols.clear(); }
  inline void ClearSwapLoop(void){
	swap_loop.clear(); }
  inline uint GetAssignedCol(uint _row_id){ return assignment[_row_id]; }
  inline uint GetAssignedRow(uint _col_id); 
  vector<uint> GetAsgnVec(void){ return assignment; }
  uint GetSwapLoopLength(void){ return swap_loop.size(); }
  void InitAsgnPtrs(void){ asgn_ptrs.clear(); asgn_ptrs.resize(row_size);
	for(uint i=0; i<row_size; i++) asgn_ptrs[i] = NULL; }
  void InitSibPtrs(void){ sib_ptrs.clear(); sib_ptrs.resize(row_size);
	for(uint i=0; i<row_size; i++) sib_ptrs[i] = NULL; }

  // initiations
  void RandSolution(void);
  void InitSolution(const vector<uint>&);
  void InitStage(void);
  
  //approximations, remove these two does not hurt the algo(safely ignore them)
  void GreedyApprox1(vector<uint>&); //local greedy
  void GreedyApprox2(vector<uint>&); //global greedy

  // pre-process and create the heaps
  void Preprocess(void);

  // get the next starting cell, if returns (-1, -1), meaning all are >=0
  cell NextStartEntry() const;
  cell NextStartEntry(uint) const;
  //void UpdateDualVars(const vector<double>& _deltas);
  void UpdateReducedCostMatrix(const vector<double>& _dual_rows, const vector<double>& _dual_cols);

  // search a swap loop using BFS
  bool SearchSwapLoop(cell& _start);

  // dual updates, return: -1 if simply update; 0 if no loop but feasible; 1 if a swap loop is found during the update, ie, root is changed & already covered.
  int DualUpdates(cell&, queue<tree_node*>&);

  // given a swap loop, tasks are swapped, matrix & assign-vec are updated
  void SwapTasks(void);

  //primal algo containing all compoments to get the ultimate solution
  //must init a solution before calling it
  void PrimalAlgo(void);

  // get current sum of costs based on current assignemnt solution 
  double ComputeCostSum(const mat& _m, const vector<uint>& _as) const;

  // some display functions
  void DisplayMatrix(const mat&) const;
  void DisplayMatrix(const mat&, const vector<uint>&) const;
  void DisplayAssignment(void) const;
  void DisplaySet(const set<uint>&) const;
  void DisplayTree(tree_node* _root);
  template<typename T>
  void DisplayVec(const vector<T>& _vec);

  // destroy the whole tree, return number of nodes;
  uint DeleteTree(tree_node* _root);

  // for testing & debugging, not used in final algorithm
  void PrimalTest(void);

  // double check if all values in reduced cost matrix have become feasible
  void DoubleCheck(void); 


private:
  //basic data members
  uint row_size;
  uint col_size;

  vector<double> dual_row_vars;	// dual row variables, ie, row labelling values 
  vector<double> dual_col_vars;	// dual col variables, ie, col labelling values 
  vector<double> deltas;	// vector to store accumulated updates

  set<uint> exploited_rows;	// rows traversed during searching loop
  set<uint> exploited_cols;	// cols traversed during searching loop

  mat orig_matrix;   		// a copy of original matrix
  mat oprt_matrix;   		// the so-called reduced cost matrix

  vector<
	priority_queue<
        	pair<uint, double>,
        	vector<pair<uint,double> >,
        	less_comp<uint> >
	> min_heaps;		// heaps to maintain searching status
 
  vector<uint> assignment;	// indices of vector <-> values in the vector
  list<cell> swap_loop;

  vector<tree_node*> asgn_ptrs;	// store the pointers only to assigned cells
  vector<tree_node*> sib_ptrs;  // store the pointers of last siblings each row 

  //set<uint> greedy_cols;	// store greedily selected starting columns

};


// for priority queue elements' comparison
template<typename T>
class less_comp{
public:
  less_comp(bool __switch = false){ _switch = __switch; };
  bool operator() (const pair<T, double>& a,
                const pair<T, double>& b) const {
  if(!_switch) // output in increasing order
    return a.second > b.second ? true : false ;
  else       // output in decreasing order
    return a.second < b.second ? true : false ;
  }

private:
  bool _switch;

};


#endif


