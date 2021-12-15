#include "TAnalysisChainLoop.h"
#include "GRootCommands.h"

#include <chrono>
#include <thread>

#include "TClass.h"
#include "TFile.h"
#include "TThread.h"

#include "TDetector.h"
#include "TGRSIint.h"
#include "TUnpackedEvent.h"

TAnalysisChainLoop* TAnalysisChainLoop::Get(std::string name, TChain *chain){
  if(name.length()==0){
    name = "chain_loop";
  }

  StoppableThread* thread = StoppableThread::Get(name);
  if(!thread) {
    if(!chain && !gAnalysis){
      return 0;
    } else if(!chain) {
      chain = gAnalysis;
    }
    thread = new TAnalysisChainLoop(name,chain);
  }
  return dynamic_cast<TAnalysisChainLoop*>(thread);
}

TAnalysisChainLoop::TAnalysisChainLoop(std::string name, TChain *chain)
  : StoppableThread(name),
    fEntriesRead(0), fEntriesTotal(chain->GetEntries()),
    input_chain(chain),
    output_queue(std::make_shared<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>()),
    fSelfStopping(true) {
  SetupChain();
}

TAnalysisChainLoop::~TAnalysisChainLoop() { }

void TAnalysisChainLoop::ClearQueue() {
  while(output_queue->Size()){
    std::shared_ptr<TUnpackedEvent> event = NULL;
    output_queue->Pop(event);
    //if(event){
    //delete event;
    //}
  }
}

int TAnalysisChainLoop::SetupChain() {
  if(!input_chain)
    return 0;

  TObjArray *array = input_chain->GetListOfBranches();
  for(int x=0;x<array->GetSize();x++) {
    TBranch *b = (TBranch*)array->At(x);
    if(b) {
      TClass *c = TClass::GetClass(b->GetName());
      if(c) {
        printf("Found  %s!\n",b->GetName());
        TDetector** det = new TDetector*;
	*det = NULL;
        det_map[c] = det;
        input_chain->SetBranchAddress(b->GetName(),det_map[c]);
      }
    }
  }
  return 0;
}

std::string TAnalysisChainLoop::Status() {
  return Form("Event: %ld / %ld", long(fEntriesRead), fEntriesTotal);
}


void TAnalysisChainLoop::Restart() {
  fEntriesRead = 0;
}

void TAnalysisChainLoop::OnEnd() {
  output_queue->SetFinished();
}

bool TAnalysisChainLoop::Iteration() {
  
  if(fEntriesRead >= fEntriesTotal) {
    if(fSelfStopping) {
      return false;
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
      return true;
    }
  }

  for(auto& elem : det_map){
    *elem.second = (TDetector*)elem.first->New();
  }
  input_chain->GetEntry(fEntriesRead++);

  std::shared_ptr<TUnpackedEvent> event = std::make_shared<TUnpackedEvent>();
  for(auto& elem : det_map) {
    TDetector* det = *elem.second;
    /*
    //if(!det->TestBit(TDetector::kUnbuilt)){
    if(det->IsBuilt()) {
      //event->AddDetector(det);
      std::cout << "HERE" << std::endl;
      //event->AddDetector(std::shared_ptr<TDetector>(det));
    } 
    else {
      
      //if(det->Timestamp()!=-1 && det->Size()!=0) {
      //std::cout << det->IsA()->GetName() << " was not present in this event (TS="
      //	  << det->Timestamp() << ")"
      //	  << std::endl;
      //}
      
      //std::cout << "I don't know what this condition means..." << std::endl;
      delete det;
    }
    */
    event->AddDetector(std::shared_ptr<TDetector>(det));
  }
    
  output_queue->Push(event);
  fInputSize = fEntriesTotal - fEntriesRead;
    
  return true;
}
