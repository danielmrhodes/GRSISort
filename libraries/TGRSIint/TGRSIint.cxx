#include "TGRSIint.h"

#include "GRootGuiFactory.h"
#include "GVersion.h"
#include "Getline.h"
#include "Globals.h"
#include "TDataParser.h"
#include "TGRSIOptions.h"
#include "TGRSIUtilities.h"
#include "GValue.h"
#include "TROOT.h"
#include "GCanvas.h"

#include "StoppableThread.h"
#include "TAnalysisHistLoop.h"
#include "TAnalysisWriteLoop.h"
#include "TDataLoop.h"
#include "TFilterLoop.h"
#include "TDetBuildingLoop.h"
#include "TEventBuildingLoop.h"
#include "TFragHistLoop.h"
#include "TFragWriteLoop.h"
#include "TFragmentChainLoop.h"
#include "TAnalysisChainLoop.h"
#include "TTerminalLoop.h"
#include "TUnpackingLoop.h"
#include "TPPG.h"
#include "TSortingDiagnostics.h"
#include "TParserLibrary.h"

#include "GRootCommands.h"
#include "TRunInfo.h"

#include "TInterpreter.h"
#include "TGHtmlBrowser.h"
//#include <pstream.h>
#include "TPython.h"

#include <thread>
#include <utility>

#include <pwd.h>

/// \cond CLASSIMP
ClassImp(TGRSIint)
/// \endcond

extern void PopupLogo(bool);
extern void WaitLogo();

TGRSIint* TGRSIint::fTGRSIint = nullptr;

TEnv* TGRSIint::fGRSIEnv = nullptr;

void ReadTheNews();

TGRSIint* TGRSIint::instance(int argc, char** argv, void* options, int numOptions, bool, const char* appClassName)
{
  /// Singleton constructor instance
  if(fTGRSIint == nullptr) {
    fTGRSIint = new TGRSIint(argc, argv, options, numOptions, true, appClassName);
    fTGRSIint->ApplyOptions();
  }
  return fTGRSIint;
}

TGRSIint::TGRSIint(int argc, char** argv, void* options, int numOptions, bool noLogo, const char* appClassName)
  : TRint(appClassName, &argc, argv, options, numOptions, noLogo), fKeepAliveTimer(nullptr),
    main_thread_id(std::this_thread::get_id()), fIsTabComplete(false), fAllowedToTerminate(true), fRootFilesOpened(0),
    fRawFilesOpened(0)
{
  /// Singleton constructor
  fGRSIEnv = gEnv;

  GetSignalHandler()->Remove();
  auto* ih = new TGRSIInterruptHandler();
  ih->Add();

  try {
    TGRSIOptions* opt = TGRSIOptions::Get(argc, argv);
    if(opt->ShouldExit()) {
      exit(0);
    }
  } catch(ParseError& e) {
    exit(1);
  }

  PrintLogo(TGRSIOptions::Get()->ShowLogo());
  SetPrompt("GRSI [%d] ");
  std::string grsipath = getenv("GRSISYS");
  gInterpreter->AddIncludePath(Form("%s/include", grsipath.c_str()));
}

void TGRSIint::ApplyOptions()
{
  /// Applies options from TGRSIOptions. This include things such as batch sorting,
  /// reading material, and logo. Also includes the setup of what to do with mid
  /// and root files that are input.
  TGRSIOptions* opt = TGRSIOptions::Get();

  if(opt->Batch()) {
    MakeBatch();
  }

  bool missing_raw_file = !all_files_exist(opt->InputFiles());

  if(!false) { // this will be change to something like, if(!ClassicRoot)
    LoadGROOTGraphics();
  }

  if(opt->ReadingMaterial()) {
    std::thread fnews = std::thread(ReadTheNews);
    fnews.detach();
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  //////           Start of redoing section        ///////
  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  // load parser library if provided
  if(!opt->ParserLibrary().empty()) {
    try {
      TParserLibrary::Get()->Load();
    } catch(std::runtime_error& e) {
      // if we failed to load the library, try to continue w/o it
      std::cerr<<DRED<<e.what()<<RESET_COLOR<<std::endl;
    }
  }

  TRunInfo::ClearVersion();
  TRunInfo::SetVersion(GRSI_RELEASE);

  TRunInfo::ClearFullVersion();
  TRunInfo::SetFullVersion(GRSI_GIT_COMMIT);

  TRunInfo::ClearDate();
  TRunInfo::SetDate(GRSI_GIT_COMMIT_TIME);

  for(auto& rawFile : opt->InputFiles()) {
    OpenRawFile(rawFile);
  }

  for(const auto& filename : opt->RootInputFiles()) {
    // this will populate gChain if able.
    //   TChannels from the root file will be loaded as file is opened.
    //   GValues from the root file will be loaded as file is opened.
    OpenRootFile(filename);
  }

  SetupPipeline();

  if(opt->StartGui()) {
    StartGUI();
  }

  for(auto& filename : opt->MacroInputFiles()) {
    RunMacroFile(filename);
  }

  std::cout<<StoppableThread::AllThreadHeader()<<std::endl;
  LoopUntilDone();
  if(opt->CloseAfterSort()) {
    int exit_status = missing_raw_file ? 1 : 0;
    Terminate(exit_status);
  }
}

void TGRSIint::LoopUntilDone()
{
  /// Outputs the thread status until all of the threads are complete.
  int iter = 0;
  while(StoppableThread::AnyThreadRunning()) {
    std::this_thread::sleep_for(std::chrono::seconds(1));

    // We need to process events in case a different thread is asking for a file to be opened.
    // However, if there is no stdin, ProcessEvents() will call Terminate().
    // This prevents the terminate from taking effect while in this context.
    fAllowedToTerminate = false;
    gSystem->ProcessEvents();
    fAllowedToTerminate = true;

    ++iter;
    if(TGRSIOptions::Get()->StatusInterval() > 0 && iter % TGRSIOptions::Get()->StatusInterval() == 0) {
      std::cout<<"\r"<<StoppableThread::AllThreadStatus()<<std::endl;
    }
    std::cout<<"\r"<<StoppableThread::AllThreadProgress()<<std::flush;
  }
  std::cout<<std::endl;
}

TGRSIint::~TGRSIint()
{
  /// Default dtor.
}

bool TGRSIint::HandleTermInput()
{
  /// Handles terminal input via TRint
  return TRint::HandleTermInput();
}

void TGRSIint::Terminate(Int_t status)
{
  /// Kills all of the threads if the process is allowed to terminate. This
  /// sends an error to TSortingDiagnostics if an analysis tree is being created
  if(!fAllowedToTerminate) {
    std::cout<<"Not allowed to terminate, sorry!"<<std::endl;
    return;
  }

  StoppableThread::SendStop();
  LoopUntilDone();
  StoppableThread::StopAll();

  if(TGRSIOptions::Get()->MakeAnalysisTree()) {
    TSortingDiagnostics::Get()->Print("error");
  }

  if((clock() % 60) == 0) {
    printf("DING!");
    fflush(stdout);
    gSystem->Sleep(500);
    printf("\r              \r");
    fflush(stdout);
  }

  TSeqCollection* canvases = gROOT->GetListOfCanvases();
  while(canvases->GetEntries() > 0) {
    static_cast<GCanvas*>(canvases->At(0))->Close();
  }

  // TChannel::DeleteAllChannels();
  TRint::Terminate(status);
}

Int_t TGRSIint::TabCompletionHook(char* buf, int* pLoc, std::ostream& out)
{
  /// Tries to do a tab completion. Returns false if unsuccsessful
  fIsTabComplete = true;
  auto result    = TRint::TabCompletionHook(buf, pLoc, out);
  fIsTabComplete = false;
  return result;
}

Long_t TGRSIint::ProcessLine(const char* line, Bool_t sync, Int_t* error)
{
  /// This takes over the native root command line. There are two main reasons for this
  /// 1. To keep the command line thread-safe.
  /// 2. To block TCanvas from opening, and to instead use our GCanvas.

  // If you print while fIsTabComplete is true, you will break tab complete.
  // Any diagnostic print statements should be done after this if statement.
  if(fIsTabComplete) {
    long res = TRint::ProcessLine(line, sync, error);
    return res;
  }

  const char* canvas = strstr(line, "TCanvas");
  if(canvas != nullptr) {
    const_cast<char*>(canvas)[0] = 'G';
  }

  if(std::this_thread::get_id() != main_thread_id) {
    return DelayedProcessLine(line);
  }

  return TRint::ProcessLine(line, sync, error);
}

void ReadTheNews()
{
  /// Opens a random wikipedia page for your enjoyment
#ifdef __APPLE__
  gROOT->ProcessLine(".! open http://en.wikipedia.org/wiki/Special:Random > /dev/null 2>&1;");
#else
  gROOT->ProcessLine(".! xdg-open http://en.wikipedia.org/wiki/Special:Random > /dev/null 2>&1;");
#endif
}

void TGRSIint::PrintLogo(bool print)
{
  /// Prints the GRSISort logo to terminal
  if(print) {
#ifdef LINUX
    const std::string& ref       = ProgramName();
    const unsigned int reflength = ref.length() - 78;
#else
    const std::string& ref       = "Sorting Program for Online and Offline Nuclear Data";
    const unsigned int reflength = 53;
#endif

    const unsigned int width = reflength + (reflength % 2);
    printf("\t*%s*\n", std::string(width, '*').c_str());
    printf("\t*%*s%*s*\n", width / 2 + 4, "GRSI Sort", width / 2 - 4, "");
    printf("\t*%*s%*s*\n", width / 2 + 12, "a remake of GRSI SPOON", width / 2 - 12, "");
    printf("\t*%*s%*s*\n", width / 2 + reflength / 2, ref.c_str(), width / 2 - reflength / 2, "");
    printf("\t*%*s%*s*\n", width / 2 + 14, "A lean, mean sorting machine", width / 2 - 14, "");
    printf("\t*%*s%*s*\n", width / 2 + 9, "version " GRSI_RELEASE, width / 2 - 9, "");
    printf("\t*%s*\n", std::string(width, '*').c_str());

    std::thread drawlogo(&TGRSIint::DrawLogo, this);
    drawlogo.detach();
  } else {
    std::cout<<"\tgrsisort version "<<GRSI_RELEASE<<std::endl;
  }
}

TFile* TGRSIint::OpenRootFile(const std::string& filename, Option_t* opt)
{
  /// Opens root files provided on the command line. Also tells you where these files
  /// are stored (ie _file0). If these files are analysis or fragment trees, they are
  /// automatically chained into chains called gFragment and gAnalysis. Once this is
  /// complete, the TChannels, GValues and RunInfo are also read in.
  TString sopt(opt);
  sopt.ToLower();

  TFile* file = nullptr;
  if(sopt.Contains("recreate") || sopt.Contains("new")) {
    // We are being asked to make a new file.
    file = new TFile(filename.c_str(), "RECREATE");
    if(file != nullptr && file->IsOpen()) {
      // Give access to the file inside the interpreter.
      const char* command = Form("TFile* _file%i = (TFile*)%luL;", fRootFilesOpened, (unsigned long)file);
      TRint::ProcessLine(command);
      fRootFilesOpened++;
    } else {
      std::cout<<"Could not create "<<filename<<std::endl;
    }
  } else {
    // Open an already existing file.
    file = new TFile(filename.c_str(), opt);
    if(file != nullptr && file->IsOpen()) {
      // Give access to the file inside the interpreter.
      const char* command = Form("TFile* _file%i = (TFile*)%luL;", fRootFilesOpened, (unsigned long)file);
      TRint::ProcessLine(command);
      std::cout<<"\tfile "<<BLUE<<file->GetName()<<RESET_COLOR<<" opened as "<<BLUE<<"_file"
	       <<fRootFilesOpened<<RESET_COLOR<<std::endl;

      // If FragmentTree exists, add the file to the chain.
      if(file->FindObjectAny("FragmentTree") != nullptr) {
	if(gFragment == nullptr) {
	  // TODO: Once we have a notifier set up
	  gFragment = new TChain("FragmentChain");
	  // gFragment->SetNotify(GrutNotifier::Get());
	}
	printf("file %s added to gFragment.\n", file->GetName());
#if MAJOR_ROOT_VERSION < 6
	gFragment->AddFile(file->GetName(), TChain::kBigNumber, "FragmentTree");
#else
	gFragment->AddFile(file->GetName(), TTree::kMaxEntries, "FragmentTree");
#endif
      }

      // If AnalysisTree exists, add the file to the chain.
      if(file->FindObjectAny("AnalysisTree") != nullptr) {
	if(gAnalysis == nullptr) {
	  gAnalysis = new TChain("AnalysisChain");
	  // TODO: Once we have a notifier set up
	  // gAnalysis->SetNotify(GrutNotifier::Get());
	}
	printf("file %s added to gAnalysis.\n", file->GetName());
#if MAJOR_ROOT_VERSION < 6
	gAnalysis->AddFile(file->GetName(), TChain::kBigNumber, "AnalysisTree");
#else
	gAnalysis->AddFile(file->GetName(), TTree::kMaxEntries, "AnalysisTree");
#endif
      }

      if(file->FindObjectAny("Channel") != nullptr) {
	file->Get("Channel"); // this calls TChannel::Streamer
      }
      if(file->FindObjectAny("Values") != nullptr) {
	file->Get("Values");
      }
      fRootFilesOpened++;
    } else {
      std::cout<<"Could not open "<<filename<<std::endl;
    }
  }

  AddFileToGUI(file);
  TRunInfo::ReadInfoFromFile();

  return file;
}

TRawFile* TGRSIint::OpenRawFile(const std::string& filename)
{
  /// Opens Raw input file and stores them in _raw if successfuly opened.
  if(!file_exists(filename.c_str())) {
    std::cerr<<R"(File ")"<<filename<<R"(" does not exist)"<<std::endl;
					 return nullptr;
					 }

  // create new raw file
  auto* file = TParserLibrary::Get()->CreateRawFile(filename);
  fRawFiles.push_back(file);

  if(file != nullptr) {
    const char* command = Form("TRawFile* _raw%i = (TRawFile*)%luL;", fRawFilesOpened, (unsigned long)file);
    ProcessLine(command);

    std::cout<<"\tfile "<<BLUE<<filename<<RESET_COLOR<<" opened as "
	     <<BLUE<<"_raw"<<fRawFilesOpened<<RESET_COLOR<<std::endl;
  }
  fRawFilesOpened++;
  return file;
}

void TGRSIint::SetupPipeline()
{
  /// Finds all of the files input as well as flags provided and makes all
  /// of the decisions about what to sort and what order to open everything up
  /// in. This also creates the output files. Starts the threads and gets the
  /// sorting going. This is really the brains of the command line sorting routine.
  TGRSIOptions* opt = TGRSIOptions::Get();

  // Determining which parts of the pipeline need to be set up.

  bool missing_raw_file = false;
  for(auto& filename : opt->InputFiles()) {
    if(!file_exists(filename.c_str())) {
      missing_raw_file = true;
      std::cerr<<"File not found: "<<filename<<std::endl;
    }
  }

  // Which input files do we have
  bool has_raw_file = !opt->InputFiles().empty() && opt->SortRaw() && !missing_raw_file;
  bool filter_data = opt->CompiledFilterFile().length();
	
  bool has_input_fragment_tree = gFragment != nullptr; // && opt->SortRoot();
  bool has_input_analysis_tree = gAnalysis != nullptr; // && opt->SortRoot();

  // Which output files are could possibly be made
  bool able_to_write_fragment_histograms =
    ((has_raw_file || has_input_fragment_tree) && opt->FragmentHistogramLib().length() > 0);
  bool able_to_write_fragment_tree       = (has_raw_file && !missing_raw_file);
  bool able_to_write_analysis_histograms = ((has_raw_file || has_input_fragment_tree || has_input_analysis_tree) &&
					    opt->AnalysisHistogramLib().length() > 0);
  bool able_to_write_analysis_tree = (able_to_write_fragment_tree || has_input_fragment_tree);

  // Which output files will we make
  bool write_fragment_histograms = (able_to_write_fragment_histograms && opt->MakeHistos());
  bool write_fragment_tree       = able_to_write_fragment_tree && opt->WriteFragmentTree();
  bool write_analysis_histograms =
    (able_to_write_analysis_histograms  && opt->MakeHistos()); //TODO: make it so we aren't always trying to generate both frag and analysis histograms
  bool write_analysis_tree = (able_to_write_analysis_tree && opt->MakeAnalysisTree());

  // Which steps need to be performed to get from the inputs to the outputs
  bool self_stopping = opt->CloseAfterSort();

  bool read_from_raw = (has_raw_file && (write_fragment_histograms || write_fragment_tree ||
					 write_analysis_histograms || write_analysis_tree));

  bool read_from_fragment_tree =
    (has_input_fragment_tree && (write_fragment_histograms || write_analysis_histograms || write_analysis_tree));

  bool generate_analysis_data =
    ((read_from_raw || read_from_fragment_tree) && (write_analysis_histograms || write_analysis_tree));

  bool read_from_analysis_tree =
    (has_input_analysis_tree && (write_analysis_histograms || write_analysis_tree) && !generate_analysis_data);

  // Extract the run number and sub run number from whatever we were given
  int run_number     = 0;
  int sub_run_number = 0;
  if(read_from_raw) {
    run_number     = fRawFiles[0]->GetRunNumber();
    sub_run_number = fRawFiles[0]->GetSubRunNumber();
  } else if(read_from_fragment_tree) {
    auto run_title = gFragment->GetListOfFiles()->At(0)->GetTitle();
    run_number     = GetRunNumber(run_title);
    sub_run_number = GetSubRunNumber(run_title);
  } else if(read_from_analysis_tree) {
    auto run_title = gAnalysis->GetListOfFiles()->At(0)->GetTitle();
    run_number     = GetRunNumber(run_title);
    sub_run_number = GetSubRunNumber(run_title);
  }

  // Choose output file names for the 4 possible output files
  std::string output_fragment_tree_filename = opt->OutputFragmentFile();
  if(output_fragment_tree_filename.length() == 0) {
    if(sub_run_number == -1) {
      output_fragment_tree_filename = Form("fragment%05i.root", run_number);
    } else {
      output_fragment_tree_filename = Form("fragment%05i_%03i.root", run_number, sub_run_number);
    }
  }

  std::string output_fragment_hist_filename = opt->OutputFragmentHistogramFile();
  if(output_fragment_hist_filename.length() == 0) {
    if(sub_run_number == -1) {
      output_fragment_hist_filename = Form("hist_fragment%05i.root", run_number);
    } else {
      output_fragment_hist_filename = Form("hist_fragment%05i_%03i.root", run_number, sub_run_number);
    }
  }

  std::string output_analysis_tree_filename = opt->OutputAnalysisFile();
  if(output_analysis_tree_filename.length() == 0) {
    if(sub_run_number == -1) {
      output_analysis_tree_filename = Form("analysis%05i.root", run_number);
    } else {
      output_analysis_tree_filename = Form("analysis%05i_%03i.root", run_number, sub_run_number);
    }
  }

  std::string output_analysis_hist_filename = opt->OutputAnalysisHistogramFile();
  if(output_analysis_hist_filename.length() == 0) {
    if(sub_run_number == -1) {
      output_analysis_hist_filename = Form("hist_analysis%05i.root", run_number);
    } else {
      output_analysis_hist_filename = Form("hist_analysis%05i_%03i.root", run_number, sub_run_number);
    }
  }

  std::vector<TFile*> cut_files;
  for(auto filename : opt->InputCutFiles()) {
    TFile* tfile = OpenRootFile(filename);
    cut_files.push_back(tfile);
    std::cout << "loading cuts\n" << std::endl;
    //AddFileToGUI(tfile);
    /*
    if(tfile && GUIIsRunning()){
      std::cout << "HERE" << std::endl;
      TPython::Bind(tfile,"tdir");
      ProcessLine("TPython::Exec(\"window.LoadCutFile(tdir)\");");
    }
    */
  }

  // No need to set up all the loops if we are just opening the interpreter or viewing histograms
  if(!read_from_raw && !read_from_fragment_tree && !read_from_analysis_tree) {
    return;
  }

  ////////////////////////////////////////////////////
  ////////////  Setting up the loops  ////////////////
  ////////////////////////////////////////////////////

  if(!write_fragment_histograms && !write_fragment_tree && !write_analysis_histograms && !write_analysis_tree) {
    // We still might want to read the calibration, values, or run info files
    for(const auto& cal_filename : opt->CalInputFiles()) {
      TChannel::ReadCalFile(cal_filename.c_str());
    }
    for(const auto& val_filename : opt->ValInputFiles()) {
      GValue::ReadValFile(val_filename.c_str());
    }
    for(const auto& info_filename : opt->ExternalRunInfo()) {
      TRunInfo::Get()->ReadInfoFile(info_filename.c_str());
    }
    return;
  }
  // Set the width of the status and each column
  StoppableThread::ColumnWidth(TGRSIOptions::Get()->ColumnWidth());
  StoppableThread::StatusWidth(TGRSIOptions::Get()->StatusWidth());

  // Different queues that can show up
  std::vector<std::shared_ptr<ThreadsafeQueue<std::shared_ptr<const TFragment>>>>    fragmentQueues;
  std::vector<std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TEpicsFrag>>>>         scalerQueues;
  std::vector<std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>>     analysisQueues;

  // The different loops that can run
  TDataLoop*          dataLoop          = nullptr;
  TFilterLoop*        FilterLoop        = nullptr;
  TUnpackingLoop*     unpackLoop        = nullptr;
  TFragmentChainLoop* fragmentChainLoop = nullptr;
  TEventBuildingLoop* eventBuildingLoop = nullptr;
  TDetBuildingLoop*   detBuildingLoop   = nullptr;
  TAnalysisChainLoop* analysisChainLoop = nullptr;

  int num_loops = 1;

  // If needed, read from the raw file
  if(read_from_raw) {
    if(fRawFiles.size() > 1) {
      std::cerr<<"I'm going to ignore all but first .mid"<<std::endl;
    }

    //dataLoop = TDataLoop::Get("1_input_loop", fRawFiles[0]);
    dataLoop = TDataLoop::Get(Form("%d_input_loop",num_loops), fRawFiles[0]);
    dataLoop->SetSelfStopping(self_stopping);
    num_loops++;

    //unpackLoop               = TUnpackingLoop::Get("2_unpack_loop");
    unpackLoop = TUnpackingLoop::Get(Form("%d_unpack_loop",num_loops));
    unpackLoop->InputQueue() = dataLoop->OutputQueue();
    num_loops++;
  }

  // If needed, read from the fragment tree
  else if(read_from_fragment_tree) {
    //fragmentChainLoop = TFragmentChainLoop::Get("1_chain_loop",gFragment);
    fragmentChainLoop = TFragmentChainLoop::Get(Form("%d_chain_loop",num_loops),gFragment);
    fragmentChainLoop->SetSelfStopping(self_stopping);
    num_loops++;
  }

  else if(read_from_analysis_tree) {
    //analysisChainLoop = TAnalysisChainLoop::Get("1_analysis_loop",gAnalysis);
    analysisChainLoop = TAnalysisChainLoop::Get(Form("%d_analysis_loop",num_loops),gAnalysis);
    analysisChainLoop->SetSelfStopping(self_stopping);
    num_loops++;
  }

  // if I am passed any calibrations, lets load those, this
  // will overwrite any with the same address previously read in.
  for(const auto& cal_filename : opt->CalInputFiles()) {
    TChannel::ReadCalFile(cal_filename.c_str());
  }
  // Set the run number and sub-run number
  if(!fRawFiles.empty()) {
    TRunInfo::Get()->SetRunInfo(fRawFiles[0]->GetRunNumber(), fRawFiles[0]->GetSubRunNumber());
  } else {
    TRunInfo::Get()->SetRunInfo(0, -1);
  }

  for(const auto& val_filename : opt->ValInputFiles()) {
    GValue::ReadValFile(val_filename.c_str());
  }
  for(const auto& info_filename : opt->ExternalRunInfo()) {
    TRunInfo::Get()->ReadInfoFile(info_filename.c_str());
  }

  TEventBuildingLoop::EBuildMode event_build_mode = TEventBuildingLoop::EBuildMode::kDefault;
  if(TRunInfo::Get()->GetDetectorInformation() != nullptr) {
    // call DetectorInformation::Set again, in case the calibration files we've loaded now changed things
    TRunInfo::Get()->GetDetectorInformation()->Set();
    event_build_mode = TRunInfo::Get()->GetDetectorInformation()->BuildMode();
  } else {
    std::cout<<"no detector information, can't set build mode"<<std::endl;
  }

  // If requested, write the fragment histograms
  if(write_fragment_histograms) {
    //TFragHistLoop* loop = TFragHistLoop::Get("3_frag_hist_loop");
    TFragHistLoop* loop = TFragHistLoop::Get(Form("%d_frag_hist_loop",num_loops));
    loop->SetOutputFilename(output_fragment_hist_filename);
    num_loops++;

    if(unpackLoop != nullptr) {
      loop->InputQueue() = unpackLoop->AddGoodOutputQueue();
    }
    if(fragmentChainLoop != nullptr) {
      loop->InputQueue() = fragmentChainLoop->AddOutputQueue();
    }
    fragmentQueues.push_back(loop->InputQueue());
  }

  // If requested, write the fragment tree
  if(write_fragment_tree) {
    
    //TFragWriteLoop* loop = TFragWriteLoop::Get("4_frag_write_loop", output_fragment_tree_filename);
    TFragWriteLoop* loop = TFragWriteLoop::Get(Form("%d_frag_write_loop",num_loops),
					       output_fragment_tree_filename);
    num_loops++;

    fNewFragmentFile     = output_fragment_tree_filename;
    if(unpackLoop != nullptr) {
      loop->InputQueue() = unpackLoop->AddGoodOutputQueue(TGRSIOptions::Get()->FragmentWriteQueueSize());
      loop->BadInputQueue() = unpackLoop->BadOutputQueue();
      loop->ScalerInputQueue() = unpackLoop->ScalerOutputQueue();
      scalerQueues.push_back(loop->ScalerInputQueue());
    }
    if(fragmentChainLoop != nullptr) {
      loop->InputQueue() = fragmentChainLoop->AddOutputQueue();
    }
    fragmentQueues.push_back(loop->InputQueue());
  }

  // If needed, generate the individual detectors from the TFragments
  if(generate_analysis_data) {
    
    TGRSIOptions::AnalysisOptions()->Print();
    //eventBuildingLoop = TEventBuildingLoop::Get("5_event_build_loop", event_build_mode, opt->AnalysisOptions()->BuildWindow());
    eventBuildingLoop = TEventBuildingLoop::Get(Form("%d_event_build_loop",num_loops),event_build_mode,
						opt->AnalysisOptions()->BuildWindow());
    num_loops++;
    
    eventBuildingLoop->SetSortDepth(opt->SortDepth());
    if(unpackLoop != nullptr) {
      eventBuildingLoop->InputQueue() = unpackLoop->AddGoodOutputQueue();
    }
    if(fragmentChainLoop != nullptr) {
      eventBuildingLoop->InputQueue() = fragmentChainLoop->AddOutputQueue();
    }
    fragmentQueues.push_back(eventBuildingLoop->InputQueue());

    //detBuildingLoop = TDetBuildingLoop::Get("6_det_build_loop");
    detBuildingLoop = TDetBuildingLoop::Get(Form("%d_det_build_loop",num_loops));
    detBuildingLoop->InputQueue() = eventBuildingLoop->OutputQueue();
    num_loops++;
  }

  // If requested, filter the built data
  if(filter_data) {

    FilterLoop = TFilterLoop::Get(Form("%d_filter_loop",num_loops));
    num_loops++;
    
    //if(raw_filtered_output) {
    //FilterLoop->OpenRawOutputFile(opt->OutputFilteredFile());
    //}

    for (auto cut_file: cut_files) {
      FilterLoop->AddCutFile(cut_file);
    }
	  
    if(detBuildingLoop != nullptr) {
	    
      FilterLoop->InputQueue() = detBuildingLoop->AddOutputQueue();;
    }
    else {
      std::cerr << DRED << "Must currently have a det building to loop to filter data" 
		<< RESET_COLOR <<std::endl;
      exit(1);
    }

    //FilterLoop->SetSelfStopping(self_stopping);
  }
  
  // If requested, write the analysis histograms
  if(write_analysis_histograms) {
    //TAnalysisHistLoop* loop = TAnalysisHistLoop::Get("7_analysis_hist_loop");
    TAnalysisHistLoop* loop = TAnalysisHistLoop::Get(Form("%d_analysis_hist_loop",num_loops));
    num_loops++;
		
    for(auto cut_file : cut_files) {
      loop->AddCutFile(cut_file);
    }

    loop->SetOutputFilename(output_analysis_hist_filename);
    if(detBuildingLoop != nullptr) {
      if(FilterLoop != nullptr) {
	loop->InputQueue() = FilterLoop->AddOutputQueue();
      }
      else {
	loop->InputQueue() = detBuildingLoop->AddOutputQueue();
      }
    }
    else {
      loop->InputQueue() = analysisChainLoop->OutputQueue();
    }

    analysisQueues.push_back(loop->InputQueue());
  }

  // If requested, write the analysis tree
  if(write_analysis_tree) {
    
    //TAnalysisWriteLoop* loop = TAnalysisWriteLoop::Get("8_analysis_write_loop",output_analysis_tree_filename);
    TAnalysisWriteLoop* loop = TAnalysisWriteLoop::Get(Form("%d_analysis_write_loop",num_loops),
						       output_analysis_tree_filename);
    num_loops++;

    if(FilterLoop != nullptr) {
      loop->InputQueue() = FilterLoop->AddOutputQueue(TGRSIOptions::Get()->AnalysisWriteQueueSize());
    }
    else {
      loop->InputQueue() = detBuildingLoop->AddOutputQueue(TGRSIOptions::Get()->AnalysisWriteQueueSize());
    }

    if(TGRSIOptions::Get()->SeparateOutOfOrder()) {
      loop->OutOfOrderQueue() = eventBuildingLoop->OutOfOrderQueue();
    }
    analysisQueues.push_back(loop->InputQueue());
  }
  
  StoppableThread::ResumeAll();
}

void TGRSIint::RunMacroFile(const std::string& filename)
{
  /// Runs a macro file. This happens when a .C file is provided on the command line
  if(file_exists(filename.c_str())) {
    const char* command = Form(".x %s", filename.c_str());
    ProcessLine(command);
  } else {
    // check if commandline arguments were supplied
    size_t beginning_pos = filename.find_first_of('(');
    if(beginning_pos != std::string::npos && filename.back() == ')') {
      std::string trueFilename = filename.substr(0,beginning_pos);
      std::string arguments = filename.substr(beginning_pos, std::string::npos);
      if(file_exists(trueFilename.c_str())) {
	const char* command = Form(".L %s", trueFilename.c_str());
	ProcessLine(command);
	command = Form("%s%s", trueFilename.substr(0,filename.find_first_of('.')).c_str(), arguments.c_str());
	ProcessLine(command);
      } else {
	std::cerr<<R"(File ")"<<trueFilename<<R"(" does not exist)"<<std::endl;
						 }
    } else {
      std::cerr<<R"(File ")"<<filename<<R"(" does not exist)"<<std::endl;
					   }
  }
}

void TGRSIint::DrawLogo()
{
  /// Draws the logo. Can be suppressed with -l
  PopupLogo(false);
  WaitLogo();
}

void TGRSIint::LoadGROOTGraphics()
{
  /// Loads root graphics in unless -b is used for batch mode.
  if(gROOT->IsBatch()) {
    return;
  }
  // force Canvas to load, this ensures global GUI Factory ptr exists.
  gROOT->LoadClass("TCanvas", "Gpad");
  gGuiFactory = new GRootGuiFactory();
}

void TGRSIint::LoadGCutG(GCutG* cutg) {
  if(GUIIsRunning()) {
    TPython::Bind(cutg, "cutg");
    ProcessLine("TPython::Exec(\"window.LoadCutG(cutg)\");");
  }
}

void TGRSIint::PrintHelp(bool print)
{
  /// Prints the help. Not sure this is used anymore.
  if(print) {
    printf(DRED BG_WHITE "     Sending Help!!     " RESET_COLOR "\n");
    new TGHtmlBrowser(gSystem->ExpandPathName("${GRSISYS}/README.html"));
  }
  return;
}

bool TGRSIInterruptHandler::Notify()
{
  /// When ctrl-c is pressed, this takes over. This can be used in the future
  /// for safe cleanup.
  if(!StoppableThread::AnyThreadRunning()) {
    std::cout<<std::endl<<DRED<<BG_WHITE<<"   Control-c was pressed in interactive mode.   "<<RESET_COLOR<<std::endl;
    exit(1);
  }
  static int timesPressed = 0;
  timesPressed++;
  switch(timesPressed) {
  case 1:
    std::cout<<std::endl<<DRED<<BG_WHITE<<"   Control-c was pressed, terminating input loop.   "<<RESET_COLOR<<std::endl;
    TGRSIint::instance()->Terminate();
    break;
  case 2:
    std::cout<<std::endl<<DRED<<BG_WHITE<<"   Control-c was pressed, stopping all queues.   "<<RESET_COLOR<<std::endl;
    StoppableThread::ClearAllQueues();
    break;
  default:
    std::cout<<std::endl<<DRED<<BG_WHITE<<"   No you shutup!   "<<RESET_COLOR<<std::endl;
    exit(1);
  }
  return true;
}

//   These variables are to be accessed only from DelayedProcessLine
// and DelayedProcessLine_Action, nowhere else.
//   These are used to pass information between the two functions.
// DelayedProcessLine_Action() should ONLY be called from a TTimer running
// on the main thread.  DelayedProcessLine may ONLY be called
// from inside ProcessLine().
namespace {
  std::mutex              g__CommandListMutex;
  std::mutex              g__ResultListMutex;
  std::mutex              g__CommandWaitingMutex;
  std::condition_variable g__NewResult;

  std::string g__LineToProcess;
  bool        g__ProcessingNeeded;

  Long_t g__CommandResult;
  bool   g__CommandFinished;
}  // namespace

Long_t TGRSIint::DelayedProcessLine(std::string command)
{
  std::lock_guard<std::mutex> any_command_lock(g__CommandWaitingMutex);

  g__LineToProcess    = std::move(command);
  g__CommandFinished  = false;
  g__ProcessingNeeded = true;
  TTimer::SingleShot(0, "TGRSIint", this, "DelayedProcessLine_Action()");

  std::unique_lock<std::mutex> lock(g__ResultListMutex);
  while(!g__CommandFinished) {
    g__NewResult.wait(lock);
  }

  return g__CommandResult;
}

void TGRSIint::DelayedProcessLine_Action()
{
  std::string message;
  {
    std::lock_guard<std::mutex> lock(g__CommandListMutex);
    if(!g__ProcessingNeeded) {
      return;
    }
    message = g__LineToProcess;
  }

  Long_t result = ProcessLine(message.c_str());
  Getlinem(EGetLineMode::kInit, (static_cast<TRint*>(gApplication))->GetPrompt());

  {
    std::lock_guard<std::mutex> lock(g__ResultListMutex);
    g__CommandResult    = result;
    g__CommandFinished  = true;
    g__ProcessingNeeded = false;
  }

  g__NewResult.notify_one();
}
