
///////////////////////////////////////////////////////////////////////////////
// File name: Define.h
// This file defines some commonly used basic data types.
// Lantao Liu, Dec 3, 2011 
///////////////////////////////////////////////////////////////////////////////

#ifndef DEFINE_H
#define DEFINE_H

#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <assert.h>

//#define DEBUG		// to toggle the debug mode (verbose info.)

#define USE_HEAP	// to toggle if use the min-heaps or not

#define USE_GREEDY	// whether to use greedy approx initial solution


#ifdef DEBUG
  #define _cout(expr) std::cout<<expr
#else
  #define _cout(expr)
#endif


#define ROW_SIZE 5
#define COL_SIZE 5
#define DOUBLE_EPSILON 1e-7 //For double type comparison: <= as valid, > invalid

#define SEED 0
#define OUTPUT_FILE "out.txt"
#define DEBUG_FILE "debug.txt"
#define VERBOSE_LEVEL 1
#define DISPLAY_WIDTH 10

#define MAX_RANDOM 1000

#define POS_INF 10e8
#define NEG_INF -10e8
#define EXCEPTION_BROKEN -0x01
#define EXCEPTION_WRONG  -0x02

/*
#ifndef min 
  #define min(x, y) (((x) > (y)) ? (y) : (x))
#endif

#ifndef max
  #define max(x, y) (((x) > (y)) ? (x) : (y))
#endif
*/

using namespace std;

typedef unsigned int uint;

typedef vector<vector<double> > mat;

typedef struct _cell {
  int row_idx;
  int col_idx;
} cell;

typedef struct _tree_node {

  cell c;
  int asgn;		// -1: root, 0: new found asgn, 1: old asgn

  // n-ary tree is constructed with siblings instead of direct children
  _tree_node* child; 
  _tree_node* parent;
  _tree_node* next_sibling; 

  // constructor is a must!
  _tree_node(){ child=NULL; parent=NULL; next_sibling=NULL; }
  ~_tree_node(){}

} tree_node;

/*
typedef struct _cell_node {
  cell c;
  bool is_assigned;
  list<cell> history;
} cell_node;
*/

#endif


