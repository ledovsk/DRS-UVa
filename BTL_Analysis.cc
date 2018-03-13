#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "pulse.hh"

std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
	{
	  if ( tmp.find( "=" )  != std::string::npos ) return tmp.substr( tmp.find_last_of("=") + 1 );
	  if ( tmp.find( "--" ) != std::string::npos ) return "yes";
	}
    }
  
  return "";
};


int main( int argc, char** argv)
{

  std::string inputRootFile = ParseCommandLine( argc, argv, "--inputRootFile" );
  if ( inputRootFile == "" )
    {
      std::cerr << "[ERROR]: Please provide a valid ROOT file to run on!! Please use --inputRootFile=<your_file_name>\n";
      return -1;
    }

  TFile* fin = new TFile(inputRootFile.c_str(), "READ");
  TTree* myTree = (TTree*)fin->Get("pulse");
  pulse* myPulse = new pulse( myTree );


  // Check performance of one DRS channel. Few simple distributoins
  myPulse->SimpleCheck(10);


  // Calibrate Ratios for one channel
  //  myPulse->CalibrateOneChannel(10);
  
  return 0;
}
