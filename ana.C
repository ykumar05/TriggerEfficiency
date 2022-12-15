#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

void ana(int sample=1)
{
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events");
  //Declare an instance of our code class
  VLLAna_Trigger m_selec;
  
  if(sample==1){
    //Add one file to chain. This is the input file.
    chain->Add("SingleMuon_2018A_1.root");
    //Set names of output files.
    hstfilename = "hst_SingleMuon.root";
    sumfilename = "sum_SingleMuon.txt";
    m_selec.SetLep(1);
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
  if(sample==2){
    //Add one file to chain. This is the input file.
    chain->Add("VLL_DYJetsToLL_M50_97.root");
    //Set names of output files.
    hstfilename = "hst_DY.root";
    sumfilename = "sum_DY.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
    m_selec.SetLep(1);
  }
  if(sample==3){
    //Add one file to chain. This is the input file.
    chain->Add("VLL_DYJetsToLL_M50_97.root");
    //Set names of output files.
    hstfilename = "hst_eleDY.root";
    sumfilename = "sum_eleDY.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
    m_selec.SetLep(0);
  }
  if(sample==4){
    //Add one file to chain. This is the input file.
    chain->Add("EGamma_2018A_10.root");
    //Set names of output files.
    hstfilename = "hst_EGamma2018A.root";
    sumfilename = "sum_EGamma2018A.txt";
    m_selec.SetLep(0);
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2018);
  }
    
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);
}
