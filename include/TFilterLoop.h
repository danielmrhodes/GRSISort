#ifndef _TFILTERLOOP_H_
#define _TFILTERLOOP_H_

/** \addtogroup Loops
 *  @{
 */

#ifndef __CINT__
#include <memory>
#endif
#include <string>

#include "StoppableThread.h"
#include "TCompiledFilter.h"
#include "ThreadsafeQueue.h"

class TFile;
//class TRawFileOut;

class TFilterLoop : public StoppableThread {
public:
  static TFilterLoop* Get(std::string name="");

  ~TFilterLoop() override;

#ifndef __CINT__
  std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>& InputQueue() { return input_queue; }
  std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>& OutputQueue(unsigned int i) {
     //if(i >= fOutputQueues.size())
     //return std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>(NULL);
     return fOutputQueues.at(i);
   }
  std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>& AddOutputQueue(size_t maxSize = 50000)
   {
      std::stringstream name;
      name << "filter_queue_" << fOutputQueues.size();
      fOutputQueues.push_back(std::make_shared<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>(name.str(), maxSize));
      return fOutputQueues.back();
   }
#endif

  void LoadLibrary(std::string library);
  std::string GetLibraryName() const;

  //void OpenRawOutputFile(const std::string& output_filename);

  void AddCutFile(TFile* cut_file);

  void ClearQueue() override;

  size_t GetItemsPushed() override
   {
      if(fOutputQueues.size() > 0) {
         return fOutputQueues.back()->ItemsPushed();
      }
      return std::numeric_limits<size_t>::max();
   } // this should work fine as all loops are always filled at the same time

   size_t GetItemsPopped() override { return 0; }  // fOutputQueue->ItemsPopped(); }
   size_t GetItemsCurrent() override { return 0; } // fOutputQueue->Size();        }
   size_t GetRate() override { return 0; }

protected:
  bool Iteration() override;
  void LoadLib(std::string libname) { compiled_filter.Load(libname); }

private:
  TFilterLoop(std::string name);

  void HandleEvent(std::shared_ptr<TUnpackedEvent> event);

  TCompiledFilter compiled_filter;

  //std::string output_filename;

#ifndef __CINT__
  std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>> input_queue;
  //std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>> output_queue;
  std::vector<std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>>  fOutputQueues;
  //std::unique_ptr<TRawFileOut> filtered_output;
#endif

  ClassDefOverride(TFilterLoop,0);
};

/*! @} */
#endif /* _TFILTERLOOP_H_ */
