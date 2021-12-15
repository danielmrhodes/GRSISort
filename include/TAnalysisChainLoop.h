#ifndef _TANALYSISCHAINLOOP_H_
#define _TANALYSISCHAINLOOP_H_

/** \addtogroup Loops
 *  @{
 */

#ifndef __CINT__
#include <atomic>
#endif

#include <map>

#include "TChain.h"
#include "TClass.h"

#include "TUnpackingLoop.h"
#include "StoppableThread.h"
#include "ThreadsafeQueue.h"
#include "TUnpackedEvent.h"

class TUnpackedEvent;

class TAnalysisChainLoop : public StoppableThread {
public:
  static TAnalysisChainLoop* Get(std::string name="",TChain *chain=nullptr);
  ~TAnalysisChainLoop() override;

#ifndef __CINT__
  std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>>& OutputQueue() { return output_queue; }
#endif

  size_t GetItemsPushed()  override { return fEntriesRead;}
  size_t GetItemsPopped()  override { return 0;}
  size_t GetItemsCurrent() override { return fEntriesTotal;}
  size_t GetRate()         override { return 0;}

  std::string Status() override;
  void ClearQueue() override;

  void OnEnd() override;

  void SetSelfStopping(bool self_stopping) { fSelfStopping = self_stopping; }
  bool GetSelfStopping() const { return fSelfStopping; }
  void Restart();

  bool Iteration() override;

private:
  TAnalysisChainLoop(std::string name, TChain *chain);

#ifndef __CINT__
  std::atomic_long fEntriesRead;
#endif
  long fEntriesTotal;

  TChain *input_chain;
#ifndef __CINT__
  std::shared_ptr<ThreadsafeQueue<std::shared_ptr<TUnpackedEvent>>> output_queue;
#endif

  bool fSelfStopping;

  int SetupChain();
  std::map<TClass*, TDetector**> det_map;
  //std::map<TClass*,std::shared_ptr<TDetector*>> det_map;
  ClassDefOverride(TAnalysisChainLoop, 0);
};

/*! @} */
#endif /* _TANALYSISCHAINLOOP_H_ */
