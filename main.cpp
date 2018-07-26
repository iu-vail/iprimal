
///////////////////////////////////////////////////////////////////////////////
// File name: main.cpp
// First obtain assignment, then call the algorithm.
// Lantao Liu, Dec 5, 2011
/////////////////////////////////////////////////////////////////////////////// 

#include "Define.h"
#include "CmdParser.h"
#include "Assignment.h"
#include "iPrimal.h"
#include "Utils.h"


int main(int argc, char** argv)
{

  //parse command line and generate an assignment
  CmdParser parser(argc, argv);
  parser.DisplayInput();

  Assignment as;
  if(parser.GetInputFileName().empty()){
    if(parser.GetSeed())
      as.SetSeed(parser.GetSeed());
    else
      as.SetSeed(time(NULL));
    _cout(endl<<"  *Seed for random generator: "<<as.GetSeed()<<endl);
    as.RandomGenerate( parser.GetAssignmentSize(), 
			 parser.GetAssignmentSize(), 
			 MAX_RANDOM, 
			 parser.GetSeed() );
  }
  else{
    ifstream myfile(parser.GetInputFileName().c_str());
    as.ImportMatrix(myfile);
  }
  // negate matrix: uncomment below to convert from Minimization to Maximization
  //as.NegateMatrix();

  as.DisplayMatrix(as.GetMatrix());

  if(as.GetRowSize() != as.GetColSize()){
     cerr<<endl<<"  The input problem must be a SQUARE matrix. You can add dummy rows/columns of 0s to convert a rectangular matrix to a square matrix."<<endl<<endl;
     exit(0);
  }

//  vector<uint> asgn_vec;
//  ifstream fasgn((char*)"asgn.txt");
//  as.ImportVec(fasgn, asgn_vec);

  Utils utils;

  try{

    Primal primal(as.GetMatrix());

    double time_start = utils.GetCurTime();
    double time = 0;
  
    //initiation
#ifdef USE_GREEDY
    vector<uint> asgn_vec;
    primal.GreedyApprox2(asgn_vec);
    primal.InitSolution(asgn_vec);
#else
    primal.RandSolution();
#endif

    primal.PrimalAlgo();

    time += utils.GetCurTime() - time_start;
    cout<<"Time used: "<<time<<endl;

    //utils.WriteVec<uint>(primal.GetAsgnVec(), (char*)"asgn.txt");

    //for testing
    //primal.PrimalTest();
  }
  catch(int e){

    if(e==EXCEPTION_BROKEN){
      cerr<<"I suspect there is something wrong..."<<endl;
      cerr<<"Current seed: "<<as.GetSeed()<<endl;
      cerr<<"The problem has been written in to \"debug.txt\". Send this file to code maintainer?  "<<endl;
      utils.WriteMatrix(as);
    }
    else if(e==EXCEPTION_WRONG){
      cerr<<"!!Wrong solution, double check it. "<<endl;
      cerr<<"Current seed: "<<as.GetSeed()<<endl;
      cerr<<"The problem has been written in to \"debug.txt\"."<<endl;
      utils.WriteMatrix(as);
    }
    else
      cerr<<"Unknown exception caught. "<<endl;
  }

  return 0;

}


