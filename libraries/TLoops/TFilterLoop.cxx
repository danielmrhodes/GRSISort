#include "TFilterLoop.h"

//#include "MakeUnique.h"
#include "TGRSIOptions.h"
//#include "TRawFileOut.h"

TFilterLoop* TFilterLoop::Get(std::string name) {
  if(name.length()==0)
    name = "filter_loop";
  TFilterLoop* loop = dynamic_cast<TFilterLoop*>(StoppableThread::Get(name));
  if(!loop)
    loop = new TFilterLoop(name);
  return loop;
}

TFilterLoop::TFilterLoop(std::string name)
  : StoppableThread(name),
    input_queue(std::make_shared<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>())/*,
    output_queue(std::make_shared<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>()),
    filtered_output(nullptr)*/{

  LoadLib(TGRSIOptions::Get()->CompiledFilterFile());
}

TFilterLoop::~TFilterLoop() { }

void TFilterLoop::ClearQueue() {
  while(input_queue->Size()){
    std::shared_ptr<TUnpackedEvent> event;
    input_queue->Pop(event);
    //if(event){
    //delete event;
    //}
  }

  for(const auto& outQueue : fOutputQueues) {
      while(outQueue->Size() != 0u) {
         std::shared_ptr<TUnpackedEvent> event;
         outQueue->Pop(event);
      }
   }

  /*
  while(output_queue->Size()){
    std::shared_ptr<TUnpackedEvent> event;
    output_queue->Pop(event);
    //if(event){
    //delete event;
    //}
  }
  */
  
  return;
}

bool TFilterLoop::Iteration() {
  
  std::shared_ptr<TUnpackedEvent> event = std::make_shared<TUnpackedEvent>();
  fInputSize = input_queue->Pop(event);
  
  if(fInputSize < 0) {
    fInputSize = 0;

    if(input_queue->IsFinished()) {
      //output_queue->SetFinished();
      for(const auto& outQueue : fOutputQueues) {
	outQueue->SetFinished();
      }
      return false;
    }
    else {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      return true;
    }
  }
  else {
    ++fItemsPopped;
  }

  if(event != nullptr) {
    HandleEvent(event);
  } 
  
  return true;
}

void TFilterLoop::HandleEvent(std::shared_ptr<TUnpackedEvent> event) {

  if(compiled_filter.MatchesCondition(event)) {
    /*
    if(filtered_output) {
      filtered_output->Write(*event);
    }
    */

    // No need to keep the raw data around after writing it to disk,
    // can save on RAM at this point.
    event->ClearRawData();

    //output_queue->Push(event);
    for(const auto& outQueue : fOutputQueues) {
      outQueue->Push(event);
   }
    
  } 
  //else {
    //delete event;
  //}

  return;
}

/*
void TFilterLoop::OpenRawOutputFile(const std::string& output_filename) {
  filtered_output = make_unique<TRawFileOut>(output_filename);
}
*/

void TFilterLoop::AddCutFile(TFile* cut_file) {
  compiled_filter.AddCutFile(cut_file);
}

void TFilterLoop::LoadLibrary(std::string library) {
  compiled_filter.Load(library);
}

std::string TFilterLoop::GetLibraryName() const {
  return compiled_filter.GetLibraryName();
}
