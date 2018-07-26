
#include "iPrimal.h"
#include <limits>
#include "Utils.h"

uint
Primal::GetAssignedRow(uint _col_id){

  uint row = row_size;
  for(uint i=0; i<row_size; i++)
    if(assignment[i] == _col_id) 
      row = i;

  assert(row<row_size);
  return row; 
}


void
Primal::RandSolution(void){

  // assign the diagnol as the initial feasible solution
  assignment.resize(row_size);
  for(uint i=0; i<row_size; i++)
    assignment[i]=i;

  // dual row vars are 0s
  dual_row_vars.resize(row_size, 0);

  // dual col vars are the values on the diagnol
  dual_col_vars.resize(row_size);
  for(uint i=0; i<col_size; i++)
    dual_col_vars[i]=orig_matrix[i][i];

  // update the oprt matrix
  for(uint i=0; i<row_size; i++){
    for(uint j=0; j<col_size; j++){
      oprt_matrix[i][j] -= dual_col_vars[j];
    }
  }

  _cout("NOTE: the index of row/col is from 0, NOT 1 !"<<endl);
  _cout("Update the matrix: "<<endl);
  DisplayMatrix(oprt_matrix);

}


void
Primal::InitSolution(const vector<uint>& _asgn){

  assert(_asgn.size()==row_size);
  // assign the diagnol as the initial feasible solution
  assignment = _asgn;
  //DisplayVec<uint>(assignment);

  // dual row vars are 0s
  dual_row_vars.resize(row_size, 0);

  // dual col vars are the values on the diagnol
  dual_col_vars.resize(col_size);
  for(uint i=0; i<col_size; i++)
    dual_col_vars[i]=orig_matrix[GetAssignedRow(i)][i];

  // update the oprt matrix
  for(uint i=0; i<row_size; i++){
    for(uint j=0; j<col_size; j++){
      oprt_matrix[i][j] -= dual_col_vars[j];
    }
  }

  _cout("NOTE: the index of row/col is from 0, NOT 1 !"<<endl);
  _cout("Update the matrix: "<<endl);
  DisplayMatrix(oprt_matrix);

}



void
Primal::GreedyApprox1(vector<uint>& _asgn){

  _asgn.clear();
  _asgn.resize(row_size);
  set<uint> mem_rows;
  //column-wise searching
  for(uint j=0; j<col_size; j++){
    double min_val = std::numeric_limits<double>::infinity(); 
    uint min_row = row_size;
    for(uint i=0; i<row_size; i++){
      if(orig_matrix[i][j] < min_val &&
      		mem_rows.find(i) == mem_rows.end()){
	min_val = orig_matrix[i][j];
	min_row = i;	
      }
    }
    _asgn[min_row] = j;
    mem_rows.insert(min_row);
  }

  //double check if the assignment is meaningful
  set<uint> test(_asgn.begin(), _asgn.end());
  assert(test.size()==row_size);
  _cout("Approximate assignment: "<<endl);
  DisplayVec<uint>(_asgn);

}


void
Primal::GreedyApprox2(vector<uint>& _asgn){

  _asgn.clear();
  _asgn.resize(row_size);
  set<uint> mem_rows;
  set<uint> mem_cols;
  priority_queue<
                pair<cell, double>,
                vector<pair<cell,double> >,
                less_comp<cell>
		> col_mins;

  //column-wise searching, push the smallest of each col in priority queue
  for(uint j=0; j<col_size; j++){
    double min_val = std::numeric_limits<double>::infinity(); 
    uint min_row = row_size;
    for(uint i=0; i<row_size; i++){
      if(orig_matrix[i][j] < min_val){
	min_val = orig_matrix[i][j];
	min_row = i;	
      }
    }
    cell c;
    c.row_idx = min_row;
    c.col_idx = j;
    col_mins.push(pair<cell, double>(c, min_val));
  }

  //assign those global minimals (conflicted) first
  while(!col_mins.empty()){
    pair<cell, double> p = col_mins.top();
    col_mins.pop();
    if(mem_rows.find(p.first.row_idx) == mem_rows.end()){
      mem_rows.insert(p.first.row_idx);
      mem_cols.insert(p.first.col_idx);
      _asgn[p.first.row_idx] = p.first.col_idx;
    }
  }

  //then assign other un-assigned
  //column-wise searching
  for(uint j=0; j<col_size; j++){
    if(mem_cols.find(j) == mem_cols.end()){
      double min_val = std::numeric_limits<double>::infinity(); 
      uint min_row = row_size;
      for(uint i=0; i<row_size; i++){
        if(orig_matrix[i][j] < min_val &&
      		mem_rows.find(i) == mem_rows.end()){
	  min_val = orig_matrix[i][j];
	  min_row = i;	
        }
      }
      _asgn[min_row] = j;
      mem_rows.insert(min_row);
    }//if
  }//for

  //double check if the assignment is meaningful
  set<uint> test(_asgn.begin(), _asgn.end());
  assert(test.size()==row_size);
  _cout("Approximate assignment: "<<endl);
  DisplayVec<uint>(_asgn);

}



void
Primal::InitStage(void){

  deltas.clear();
  deltas.resize(col_size, 0);
  ClearExploited();
  ClearSwapLoop();
  InitAsgnPtrs();
  InitSibPtrs();
  //init min-heaps;
  min_heaps.clear();  
  min_heaps.resize(row_size);

}


/* preprocess is not used in this code, since the min_heaps are created on the fly, and they are created only when they are required, i.e., immediately after a new row is colored */
void
Primal::Preprocess(void){

  //update oprt_matrix, already done in SwapTask of previous stage
  //UpdateReducedCostMatrix(dual_row_vars, dual_col_vars);  

  min_heaps.clear();
  min_heaps.resize(row_size);

  for(uint i=0; i<row_size; i++)
    for(uint j=0; j<col_size; j++){
      if(oprt_matrix[i][j]>=0 && j!=assignment[i]){
	pair<uint, double> p(j, oprt_matrix[i][j]);
        min_heaps[i].push(p);
      }
    }

}



cell
Primal::NextStartEntry(void) const{
  
  double min_v = 0;
  int min_i, min_j;
  min_i = min_j = -1;
  for(uint i=0; i<row_size; i++)
    for(uint j=0; j<col_size; j++)
      if(oprt_matrix[i][j] < min_v){
	min_v = oprt_matrix[i][j];
	min_i = i;
	min_j = j;
      }

  cell _c;
  _c.row_idx = min_i;
  _c.col_idx = min_j;

  return _c; 
}


cell
Primal::NextStartEntry(uint _col_id) const {
  
  double min_v = 0;
  int min_i, min_j;
  min_i = min_j = -1;
  for(uint i=0; i<row_size; i++)
    if(oprt_matrix[i][_col_id] < min_v){
	min_v = oprt_matrix[i][_col_id];
	min_i = i;
	min_j = _col_id;
      }

  cell _c;
  _c.row_idx = min_i;
  _c.col_idx = min_j;
  _cout("Column ["<<_col_id<<"]: starting from ("<<_c.row_idx<<", "<<_c.col_idx<<") with reduced cost: "<<min_v<<endl);

  return _c; 
}


void
Primal::UpdateReducedCostMatrix(const vector<double>& _dual_rows, const vector<double>& _dual_cols){
 
  for(uint i=0; i<row_size; i++)
    for(uint j=0; j<col_size; j++)
      oprt_matrix[i][j] = orig_matrix[i][j] - _dual_rows[i] - _dual_cols[j];

}


bool
Primal::SearchSwapLoop(cell& _start) {

  //travs.resize(100, 0);

  queue<tree_node*> bfs;  	//used for searching shortest loop
  tree_node* tree_root = new tree_node;  
  tree_root->c = _start;
  tree_root->asgn = -1;

  //vector<double> deltas(col_size, 0);

  //init the current tree_node
  tree_node* p_cur_tree_node = tree_root;
  bfs.push(p_cur_tree_node);
  
  bool found_status = false;
  int dual_update_stat = -1;	//3 cases, can be -1, 0, 1

  uint loop_cnt = 0;
  //continue to search until find a loop or meet a dead end
  while(bfs.size() && !found_status){

    if(loop_cnt++ > col_size*col_size){
      //cout<<" row_size: "<<row_size<<" col_size: "<<col_size<<" loop cnt: "<<loop_cnt<<endl;
      throw EXCEPTION_BROKEN;
    }

    //get the front cell_node and pop it off queue
    p_cur_tree_node = bfs.front();
    bfs.pop();

    // if we are now at an assigned entry, we need to search other feasible
    // entries in the same row
    if(p_cur_tree_node->asgn == 1){
      // record this row as exploited
      exploited_rows.insert(p_cur_tree_node->c.row_idx);

      // create the min_heap on the fly when needed
      for(uint j=0; j<col_size; j++){
        if(oprt_matrix[p_cur_tree_node->c.row_idx][j]>=0 
		&& j!=assignment[p_cur_tree_node->c.row_idx]){
          pair<uint, double> p(j, oprt_matrix[p_cur_tree_node->c.row_idx][j]);
          min_heaps[p_cur_tree_node->c.row_idx].push(p);
        }
      }

      for(uint col_index=0; col_index<col_size; col_index++) {
        // can find multiple satisfying new cell_nodes in a for loop
        // if == 0
        if(fabs(oprt_matrix[p_cur_tree_node->c.row_idx][col_index] + deltas[col_index] - deltas[GetAssignedCol(p_cur_tree_node->c.row_idx)]) < DOUBLE_EPSILON
		&& exploited_cols.find(col_index) == exploited_cols.end()
		&& col_index != (uint)p_cur_tree_node->c.col_idx){
          //current new cell
          cell cur_c;
          cur_c.row_idx = p_cur_tree_node->c.row_idx;
          cur_c.col_idx = col_index;
          // below is pretty useful for tracking or debugging the process
          _cout("  .found a new 0: ("
	  	<<cur_c.row_idx<<", "<<cur_c.col_idx<<")"<<endl);

	  //create a new tree_node
          tree_node* new_tree_node = new tree_node;
	  new_tree_node->c = cur_c;
	  new_tree_node->asgn = 0;
	  new_tree_node->parent = p_cur_tree_node;
          // check if it already has siblings
	  if(p_cur_tree_node->child == NULL){
            p_cur_tree_node->child = new_tree_node;
	    //memorize such unique node into sib_ptrs, for quick localization 
	    sib_ptrs[cur_c.row_idx] = new_tree_node;
          }
          else{
	    assert(sib_ptrs[cur_c.row_idx] != NULL);
	    sib_ptrs[cur_c.row_idx]->next_sibling = new_tree_node;
	    sib_ptrs[cur_c.row_idx] = new_tree_node;
	  }

	  //push the new cell_node in the BFS queue
	  bfs.push(new_tree_node);
        }//if
      }//for i
    }
    // otherwise we must be at an un-assigned entry,
    // and we need search the col and push the one assigned entry out
    else{
      //in iprimal-v2, below is moved to bottom, in the same col, only one asgned cell need be pushed in Q 
      // record this col as exploited 
      //exploited_cols.insert(p_cur_tree_node->c.col_idx);

      for(uint row_index=0; row_index<row_size; row_index++){
        // can find only one satisfying new cell_node in a for loop
        if(assignment[row_index] == (uint)p_cur_tree_node->c.col_idx
		&& exploited_rows.find(row_index) == exploited_rows.end() 
      		//added in iprimal-v2, in the same col, only one asgned cell need be pushed in Q 
		&& exploited_cols.find((uint)p_cur_tree_node->c.col_idx) == exploited_cols.end() 
		){
          //current new cell
          cell cur_c;
          cur_c.row_idx = row_index;
          cur_c.col_idx = p_cur_tree_node->c.col_idx;
          // below is pretty useful for tracking or debugging the process
          _cout("  .found a new *: ("
	  	<<cur_c.row_idx<<", "<<cur_c.col_idx<<")"<<endl);

	  //create a new tree_node
	  tree_node* new_tree_node = new tree_node;
          new_tree_node->c = cur_c;
	  new_tree_node->asgn = 1;
	  new_tree_node->parent = p_cur_tree_node;
          p_cur_tree_node->child = new_tree_node;

	  //memorize such unique node into asgn_ptrs, for quick localization 
          asgn_ptrs[row_index] = new_tree_node;

    	  //do a check if the traversal forms a loop
          if(row_index == (uint)_start.row_idx){
	    // get swap loop
	    tree_node* tmp_ptr = new_tree_node;
            swap_loop.clear();
	    while(tmp_ptr->asgn != -1){
              swap_loop.push_back(tmp_ptr->c);
	      tmp_ptr=tmp_ptr->parent;
	    }
            swap_loop.push_back(tmp_ptr->c); // the root

	    // record row as exploited then exit loop, actually unnecessary
            exploited_rows.insert(row_index);

	    found_status = true;
          }
          else{
	    //push the new cell_node in the BFS queue
	    bfs.push(new_tree_node);
          }

          // record this col as exploited
          exploited_cols.insert(p_cur_tree_node->c.col_idx);
          break;	//break this for loop
        }//if
      }//for i
    }//else
    
    // embed the dual update here, to combine iterations into one stage
    if(bfs.empty() && !found_status){
      _cout("BFS queue is empty, need dual update."<<endl);
      dual_update_stat = DualUpdates(_start, bfs);
      // swap loop already found inside the dual updates
      if(dual_update_stat >=0) 
	break;
    }
  }//while

#ifdef DEBUG
  _cout("Exploited rows: ");
  DisplaySet(exploited_rows);

  _cout("          cols: ");
  DisplaySet(exploited_cols);

  DisplayTree(tree_root);

  // the loop should be an even number
  assert( swap_loop.size()%2 == 0 );

  //if(swap_loop.size() && dual_update_stat != 0){
  if(swap_loop.size()){
     _cout("Found shortest loop: ");
     for(list<cell>::reverse_iterator itr=swap_loop.rbegin(); 
		itr!=swap_loop.rend(); itr++){
       _cout("("<<itr->row_idx<<", "<<itr->col_idx<<")->");
     }
     _cout("head"<<endl<<endl);
  }
#endif

  //destroy the tree here
  //travs[_start.col_idx] += DeleteTree(tree_root);
  return found_status;

}


int
Primal::DualUpdates(cell& _start_cell, queue<tree_node*>& _bfs){

  // search the epsilon (minimal value) for varibles adjustments
  // it must satisfy 3 conditions: 
  // 1. >=0 in oprt_matrix; 2. in exploited rows; 3. in un-exploited cols

  // define and init epsilon as a largest positive value
  double epsilon = std::numeric_limits<double>::infinity();
  bool updated_epsilon = false;

#ifdef USE_HEAP

  /**** lines below are version using min-heaps, relatively faster ****/
  // search for min epsilon
  for(uint row_index=0; row_index<row_size; row_index++){
    if(exploited_rows.find(row_index) != exploited_rows.end()){
      //expose the first unvisited entry
      while(!min_heaps[row_index].empty()){ 
	if(exploited_cols.find(min_heaps[row_index].top().first) != exploited_cols.end())
	  min_heaps[row_index].pop();
        else
	  break;
      }
      if(!min_heaps[row_index].empty()){
        double cell_val = min_heaps[row_index].top().second + deltas[min_heaps[row_index].top().first] - deltas[GetAssignedCol(row_index)];
        if(cell_val < epsilon){
	  epsilon = cell_val;
          updated_epsilon = true;
        }
      }//if !min
    }//if exploited row
  }//for row_index

  // one more round traversal to find all qualified entries
  queue<cell> tmp_q;
  cell qual_c;
  if(updated_epsilon){
    for(uint row_index=0; row_index<row_size; row_index++)
      if(exploited_rows.find(row_index) != exploited_rows.end()){
        while(!min_heaps[row_index].empty()){ 
	  if(exploited_cols.find(min_heaps[row_index].top().first) != exploited_cols.end()){
            min_heaps[row_index].pop();
          }
          else{
            double cell_val = min_heaps[row_index].top().second + deltas[min_heaps[row_index].top().first] - deltas[GetAssignedCol(row_index)];
            if(fabs(cell_val - epsilon)<DOUBLE_EPSILON){
              qual_c.row_idx = row_index;
	      qual_c.col_idx = min_heaps[row_index].top().first;
              tmp_q.push(qual_c);
	      //this entry will be colored as visited later, so remove it here
	      // and this also expose next available entry for judging
	      min_heaps[row_index].pop();
            }
	    else if( cell_val - epsilon > DOUBLE_EPSILON){
              // stop at encountering the first uncolored entry that > epsilon
	      break;
	    }
            else {
	      cerr<<"Very weired thing happens."<<endl;
	      assert(0);
            }
          }//else
        }//while
      }//if exploited rows
  }//if updated epsilon

#else

  /**** lines below are version NOT using min-heaps, relatively slower ****/
  for(uint row_index=0; row_index<row_size; row_index++)
    if(exploited_rows.find(row_index) != exploited_rows.end())
      for(uint col_index=0; col_index<col_size; col_index++)
        if(exploited_cols.find(col_index) == exploited_cols.end()){
	  //double cell_val = oprt_matrix[row_index][col_index] - deltas[col_index];
	  double cell_val = oprt_matrix[row_index][col_index] + deltas[col_index] - deltas[GetAssignedCol(row_index)];
          if(cell_val > 0 && cell_val < epsilon){
            epsilon = cell_val;
            updated_epsilon = true;
          }
        }

  // one more round traversal to find all qualified entries
  queue<cell> tmp_q;
  cell qual_c;
  if(updated_epsilon){
    for(uint row_index=0; row_index<row_size; row_index++)
      if(exploited_rows.find(row_index) != exploited_rows.end())
        for(uint col_index=0; col_index<col_size; col_index++)
          if(exploited_cols.find(col_index) == exploited_cols.end()){
	    double cell_val = oprt_matrix[row_index][col_index] + deltas[col_index] - deltas[GetAssignedCol(row_index)];
            if(fabs(cell_val - epsilon)<DOUBLE_EPSILON){
              qual_c.row_idx = row_index;
	      qual_c.col_idx = col_index;
              tmp_q.push(qual_c);
            }
          }
  }//if

#endif

  // if found nothing, i.e., all available entries are <=0, leaving a void set
  if(!updated_epsilon)
    epsilon = - (oprt_matrix[_start_cell.row_idx][_start_cell.col_idx] 
		+ deltas[_start_cell.col_idx] - deltas[GetAssignedCol(_start_cell.row_idx)]);
  assert(epsilon > 0 && epsilon < std::numeric_limits<double>::infinity());

  _cout("# epsilon: "<<epsilon<<endl);

  // delta update
  for(uint col_index=0; col_index<col_size; col_index++)
    if(exploited_cols.find(col_index) != exploited_cols.end())
      deltas[col_index] += epsilon;  

  _cout("Updated deltas: "<<endl);
  DisplayVec<double>(deltas);

  vector<double> tmp_d(row_size);
  
  for(uint row_index=0; row_index<row_size; row_index++){
    tmp_d[row_index] = oprt_matrix[row_index][_start_cell.col_idx] + deltas[_start_cell.col_idx] - deltas[GetAssignedCol(row_index)];
  }

  _cout("Updated column ["<<_start_cell.col_idx<<"]:"<<endl);
  DisplayVec<double>(tmp_d);

  double min_d = std::numeric_limits<double>::infinity();
  uint new_start_row = row_size;
  for(uint i=0; i<row_size; i++)
    if(tmp_d[i] < min_d){ 
      min_d = tmp_d[i];
      new_start_row = i;
    }
  assert(new_start_row != row_size);
  _cout("Updated root: ("<<new_start_row<<", "<<_start_cell.col_idx<<")"<<endl);

  //terminate this stage if conditions below occur
  if(min_d >= 0){
    _cout("No loop, but the whole column is feasible!"<<endl);
    swap_loop.clear();

    return 0;
  }
  else if(new_start_row !=(uint)_start_cell.row_idx)
    if(exploited_rows.find(new_start_row) != exploited_rows.end()){
      _cout("Loop already exists!"<<endl);
      //update the starting cell
      _start_cell.row_idx = new_start_row;
      //now search the loop or notify the upper funciton to search it
      assert(asgn_ptrs[new_start_row]!=NULL);
      tree_node* tmp_ptr = asgn_ptrs[new_start_row];
      swap_loop.clear();
      while(tmp_ptr->asgn != -1){
        swap_loop.push_back(tmp_ptr->c);
        tmp_ptr=tmp_ptr->parent;
      }
      //finally append the updated root 
      swap_loop.push_back(_start_cell); 

      return 1;
    }
 
  while(!tmp_q.empty()){
    cell cur_c = tmp_q.front();
    tmp_q.pop();
    tree_node* p_cur_tree_node = asgn_ptrs[cur_c.row_idx];
    _cout("Current dead-end entry is: ("<<p_cur_tree_node->c.row_idx<<", "<<p_cur_tree_node->c.col_idx<<")"<<endl);
/*
    for(uint col_index=0; col_index<col_size; col_index++){
      //if(exploited_cols.find(col_index) == exploited_cols.end() && oprt_matrix[epsilon_row_idx][col_index] == oprt_matrix[epsilon_row_idx][epsilon_col_idx] ){
      if(fabs(oprt_matrix[epsilon_row_idx][col_index] + deltas[col_index] - deltas[GetAssignedCol(epsilon_row_idx)]) < DOUBLE_EPSILON
                && exploited_cols.find(col_index) == exploited_cols.end()
                && col_index != (uint)p_cur_tree_node->c.col_idx){
          //current new cell
          cell cur_c;
          cur_c.row_idx = epsilon_row_idx;
          cur_c.col_idx = col_index;
*/
          _cout("  found a new 0: ("
	  	<<cur_c.row_idx<<", "<<cur_c.col_idx<<")"<<endl);

	  //create a new tree_node
          tree_node* new_tree_node = new tree_node;
	  new_tree_node->c = cur_c;
	  new_tree_node->asgn = 0;
	  new_tree_node->parent = p_cur_tree_node;
          // check if it already has siblings
	  if(p_cur_tree_node->child == NULL){
            p_cur_tree_node->child = new_tree_node;
	    sib_ptrs[cur_c.row_idx] = new_tree_node;
          }
          else{
	    assert(sib_ptrs[cur_c.row_idx] != NULL);
	    sib_ptrs[cur_c.row_idx]->next_sibling = new_tree_node;
	    sib_ptrs[cur_c.row_idx] = new_tree_node;
	  }

	  //push the new cell_node in the BFS queue
	  _bfs.push(new_tree_node);
      //}//if
    //}//for
  }//while tmp_q

  return -1;  // did not form loop, countinue searching outside
}


void
Primal::SwapTasks(void){

  double augment_delta = 0;
  cell c_first = swap_loop.back();
  if(!swap_loop.empty()){
    // the first cell in the loop is special
    // only the column corresponding the first cell need to be adjusted

    // get the first cell, which was stored reversely
    augment_delta = - (oprt_matrix[c_first.row_idx][c_first.col_idx] + deltas[c_first.col_idx] - deltas[GetAssignedCol(c_first.row_idx)] );
    assert(augment_delta>=0);
    _cout("Augment delta: "<<augment_delta<<endl);
  }

  _cout("Before swaps:"<<endl);
  DisplayMatrix(oprt_matrix);

  // update dual variables
  for(uint i=0; i<col_size; i++){
    dual_col_vars[i] -= deltas[i];
    dual_row_vars[i] += deltas[GetAssignedCol(i)];
  }
  if(!swap_loop.empty()){
    // also specially process the starting col, to augment the sum
    dual_col_vars[c_first.col_idx] -= augment_delta;
  }

  _cout("Update dual cols: ");
  DisplayVec<double>(dual_col_vars);
  _cout("Update dual rows: ");
  DisplayVec<double>(dual_row_vars);

  UpdateReducedCostMatrix(dual_row_vars, dual_col_vars);

  // swap the tasks in the assignment vector
  // only need the odd-numbered cell from the loop
  for(list<cell>::reverse_iterator itr=swap_loop.rbegin(); itr!=swap_loop.rend(); itr++){
    cell c_odd = *itr++;
    assert(oprt_matrix[c_odd.row_idx][c_odd.col_idx] < DOUBLE_EPSILON);
    assignment[c_odd.row_idx] = c_odd.col_idx;
  } 

  _cout("After swaps:"<<endl);
  DisplayMatrix(oprt_matrix);
  _cout(endl);

}



void
Primal::PrimalAlgo(void){

  //initiation is moved outside
  //RandInitAlgo();

  Utils utils;
  cell next;

  for(uint i=0; i<col_size; i++){
    
      InitStage();
      next = NextStartEntry(i);

      //_cout("Root cell_node: ("<<next.row_idx<<", "<<next.col_idx<<")"<<endl);
      if(next.row_idx != -1){

/*
#ifdef USE_HEAP
        Preprocess();
#endif
*/
        SearchSwapLoop(next);
        SwapTasks();

      }//if

  }//for

  // double check the reduced cost to see if all become feasible
  DoubleCheck();

//#ifdef DEBUG
  cout<<endl<<"***********************************************"<<endl;
  cout<<"Final solution:"<<endl;
  DisplayMatrix(orig_matrix, assignment);
  cout<<"Assigned row-col pairs: "<<endl;
  DisplayAssignment();
  cout<<"Minimization result: "<<ComputeCostSum(orig_matrix, assignment)<<endl;
  cout<<"***********************************************"<<endl<<endl;
//#else
//  cout<<"Minimization result: "<<ComputeCostSum(orig_matrix, assignment)<<endl;
//#endif

}


/* just for testing, not used in the algo */
void 
Primal::PrimalTest(void){

  // applies to the demo in Balinsk's paper
  RandSolution();
  InitStage();
  cell next = NextStartEntry(0);
  cout<<"Starting cell_node: ("<<next.row_idx<<", "<<next.col_idx<<")"<<endl;
  next = NextStartEntry(1);
  cout<<"Starting cell_node: ("<<next.row_idx<<", "<<next.col_idx<<")"<<endl;

  SearchSwapLoop(next);
  SwapTasks();

  next = NextStartEntry(2);
  cout<<"Starting cell_node: ("<<next.row_idx<<", "<<next.col_idx<<")"<<endl;
  InitStage();
  SearchSwapLoop(next);
  SwapTasks();

  next = NextStartEntry(3);
  cout<<"Starting cell_node: ("<<next.row_idx<<", "<<next.col_idx<<")"<<endl;
  InitStage();
  SearchSwapLoop(next);
  SwapTasks();

}


void
Primal::DoubleCheck(void){

  for(uint i=0; i<row_size; i++)
    for(uint j=0; j<col_size; j++)
      if(oprt_matrix[i][j] < 0 && 
		fabs(oprt_matrix[i][j]) > DOUBLE_EPSILON){
        throw EXCEPTION_WRONG;
      }

}



double
Primal::ComputeCostSum(const mat& _m, const vector<uint>& _as) const{

  double sum = 0;
  for(uint i=0; i<_m.size(); i++){
    for(uint j=0; j<_m[0].size(); j++)
      if(_as[i] == j)
	sum += _m[i][j];
  }
 
  return sum;

}


void
Primal::DisplayMatrix(const mat& _m) const{

  if(_m[0].size() > DISPLAY_WIDTH){
    _cout("Matrix is big, not displaying."<<endl);
    return;
  }
  
  for(uint i=0; i<_m.size(); i++){
    for(uint j=0; j<_m[0].size(); j++)
      _cout(" "<<_m[i][j]<<"\t");
    _cout(endl);
  }

}

void
Primal::DisplayMatrix(const mat& _m, const vector<uint>& _as) const{

  if(_m[0].size() > DISPLAY_WIDTH){
    cout<<"Matrix is big, not displaying."<<endl;
    return;
  }

  for(uint i=0; i<_m.size(); i++){
    for(uint j=0; j<_m[0].size(); j++){
      if(_as[i] == j)
        cout<<" "<<_m[i][j]<<"*\t";
      else
        cout<<" "<<_m[i][j]<<"\t";
    }
    cout<<endl;
  }

}


void
Primal::DisplayAssignment(void) const{

  for(uint i=0; i<assignment.size(); i++)
    cout<<"("<<i<<","<<assignment[i]<<") ";
  cout<<endl;

}


void
Primal::DisplaySet(const set<uint>& _s) const{

  for(set<uint>::iterator itr=_s.begin(); itr!=_s.end(); itr++)
    _cout(*itr<<" "); 
  _cout(endl);

}


//templated function not defined in header, so can only used in this .cpp
template<typename T>
void
Primal::DisplayVec(const vector<T>& _vec){

  for(typename vector<T>::const_iterator itr=_vec.begin(); itr!=_vec.end(); itr++)
    _cout(*itr<<" ");
  _cout(endl);

}


void
Primal::DisplayTree(tree_node* _root){

  queue<tree_node*> bfs;
  bfs.push(_root);
  while(!bfs.empty()){
    tree_node* p_n = bfs.front();
    bfs.pop(); 
    if(p_n->asgn==-1)
      _cout("  -the root is ^: ("
		<<p_n->c.row_idx<<", "<<p_n->c.col_idx<<")"<<endl);
    else if (p_n->asgn == 0)
      _cout("  -found a new 0: ("
		<<p_n->c.row_idx<<", "<<p_n->c.col_idx<<")"<<endl);
    else
      _cout("  -found a new *: ("
		<<p_n->c.row_idx<<", "<<p_n->c.col_idx<<")"<<endl);

    tree_node* p_sibling;
    if(p_n->child != NULL){
      bfs.push(p_n->child); 
      p_sibling = p_n->child;
      while(p_sibling->next_sibling != NULL){
        p_sibling = p_sibling->next_sibling;
        bfs.push(p_sibling);
      }
    }
  }

}


uint
Primal::DeleteTree(tree_node* _root){

  uint num = 0;
  queue<tree_node*> bfs;
  bfs.push(_root);
  while(!bfs.empty()){
    tree_node* p_n = bfs.front();
    bfs.pop(); 
    num++;
    tree_node* p_sibling;
    if(p_n->child != NULL){
      bfs.push(p_n->child); 
      p_sibling = p_n->child;
      while(p_sibling->next_sibling != NULL){
        p_sibling = p_sibling->next_sibling;
        bfs.push(p_sibling);
      }
    }
    delete p_n;
  }

  return num;

}




