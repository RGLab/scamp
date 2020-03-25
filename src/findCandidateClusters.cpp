#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "scmp.h"
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <random>

void ccSearchThread(const std::vector<std::vector<double>>& dmRef,
		    const std::vector<std::vector<int>>& rvRef,
		    const double& dipVal,
		    const int& cLowerBound,
		    const bool& repAllow,
		    const int& maxSearchD,
		    const long& mclNum,
		    const long& maxNumGat,
		    std::mutex& mutRef,
		    searchResults& ccdsRef,
		    const bool& searchRandom,
		    const bool& searchRestricted,
		    long& bsCounter,
		    const bool& annotationForest,
		    const long& pathologyLimit,
		    std::vector<bool>& sigilVector,
		    const int indexSigil,//important to copy this index from main thread -- else iteration points past end of sigilVector
		    unsigned long long rSeed,
		    const bool& searchParEx,
		    int parExStartRoot,//also need to copy in, otherwise roots repeat.
		    std::condition_variable& condRef,
		    int& threadCounter,
		    unsigned long subSampleThreshold,
		    unsigned long subSampleSize,
		    unsigned long subSampleIterations)
{
  //randomly sample from the space of clustering trees, subject to the restrictions in rvRef
  searchResults parResult =  candidateClusterSearch(dmRef,
						    rvRef,
						    dipVal,
						    cLowerBound,
						    repAllow,
						    maxSearchD,
						    mclNum,
						    maxNumGat,
						    searchRandom,
						    searchRestricted,
						    annotationForest,
						    rSeed,
						    searchParEx,
						    parExStartRoot,
						    subSampleThreshold,
						    subSampleSize,
						    subSampleIterations);
  
  //having completed a random candidate cluster search, lock for updating
  std::lock_guard<std::mutex> guard(mutRef); 
  if ((parResult.abortIteration) && (((dmRef[0]).size()) > (100*pathologyLimit))) {
    //this condition is only true if the search didn't find a candidate cluster
    //and the subset size exceed a multiple of the pathology limit. The size check is b/c
    //we want to allow residual clusters to be lumped together.
    ++bsCounter;
  }
  else {
    std::vector<std::vector<long>> parCands = parResult.candidates;
    for (auto cc : parCands) {
      ccdsRef.candidates.push_back(cc);
    }
    std::vector<gateInfo> parGateLocs = parResult.gateLocations;
    //std::cout << "gate info found in thread " << indexSigil << std::endl;
    for (auto gL : parGateLocs) {
      ccdsRef.gateLocations.push_back(gL);
      //std::cout << "colNumber: " << gL.colNumber << "; numGates: " << gL.numGates << "; gateDepth: " << gL.gateDepth << std::endl;
    }
    ccdsRef.abortIteration = false;
  }
  //update the sigil vector to notify main thread this thread can be relaunched.
  sigilVector[indexSigil] = true;
  //decrement the thread count, and notify scheduling thread.
  --threadCounter;
  condRef.notify_one();
  return;
}


searchResults findCandidateClusters(const std::vector<std::vector<double>>& fDataVals,
				    const std::vector<std::vector<int>>& fRestrictedVals,
				    const double& fDipT,
				    const int& fClusterLB,
				    const bool& fRepeatsAllowed,
				    const int& fMaxSearchDepth,
				    const long& fMaxClusterNum,
				    const long& fMaxNumberOfGates,
				    const bool& fUseRestrictedValue,
				    const bool& fRandomSearch,
				    const double& maxAllowedTime,
				    const bool& printDebugInfo,
				    const bool& annotationForestRun,
				    unsigned long numThreadsTotal,
				    unsigned long long& randomSeed,
				    unsigned long subSampleThreshold,
				    unsigned long subSampleSize,
				    unsigned long subSampleIterations)
{
  searchResults fCandClusters;
  bool underMaxTime = true;
  bool noPathology = true;
  auto startSearch = std::chrono::steady_clock::now();
  auto checkPoint = std::chrono::steady_clock::now();
  std::chrono::duration<double> currentDuration;
  double elapsedTime;
  long badSearchCounter = 0;
  const auto numInSubset = ((fDataVals[0]).size());
  const long pathologyLimit = (fClusterLB * 10); //only check for vacuous clusterings if the subset is large enough.
  searchResults parResult;
  bool noAnnForestInterrupt = true;
  //auto duration_s = std::chrono::duration<double>(dcast);
  if (printDebugInfo) {
    std::cout << "Max allowed time for search: " << maxAllowedTime << std::endl;
  }
  std::random_device rdCCS;
  unsigned long long currentSeed = randomSeed;
  bool ccSearchParex = false; // flag for exhaustive parallel search.
  int ccSearchParexRoot = 0; //used to signal starting root for exhaustive parallel search.
  if ((fRandomSearch) && (numThreadsTotal > 2)) {
    std::mutex ccMutex;
    std::unique_lock<std::mutex> guardCandidateClusters(ccMutex);
    //std::unique_lock<std::mutex> initLock(ccMutex);
    std::condition_variable searchCompleteCV;
    int launchedThreads = 0, activeThreads = 0;
    //std::mutex sigMutex;
    std::vector<bool> sigVec((numThreadsTotal-1),false);
    std::vector<std::thread> threadsVec;
    bool newThreadLaunched = false;
    //launch the initial worker threads
    for (int i = 0; i < (numThreadsTotal-1); ++i) {
      // if the user has set a random seed, progressively increment for each thread.
      if (randomSeed > 0) {
	currentSeed += 1;
      }
      //otherwise, set seed for thread from random device.
      else {
	currentSeed = rdCCS();
      }
      ++activeThreads;
      ++launchedThreads;
      threadsVec.push_back(std::thread(ccSearchThread,
				       std::ref(fDataVals),
				       std::ref(fRestrictedVals),
				       std::ref(fDipT),
				       std::ref(fClusterLB),
				       std::ref(fRepeatsAllowed),
				       std::ref(fMaxSearchDepth),
				       std::ref(fMaxClusterNum),
				       std::ref(fMaxNumberOfGates),
				       std::ref(ccMutex),
				       std::ref(fCandClusters),
				       std::ref(fRandomSearch),
				       std::ref(fUseRestrictedValue),
				       std::ref(badSearchCounter),
				       std::ref(annotationForestRun),
				       std::ref(pathologyLimit),
				       std::ref(sigVec),
				       i,
				       currentSeed,
				       std::ref(ccSearchParex),
				       ccSearchParexRoot,
				       std::ref(searchCompleteCV),
				       std::ref(activeThreads),
				       subSampleThreshold,
				       subSampleSize,
				       subSampleIterations));
      if (printDebugInfo) {
	std::cout << "launched thread " << i << std::endl;
	std::cout << "Used seed " << currentSeed << std::endl;
      }
    }
    //continue to spawn threads until we find enough candidate clusters.
    while (launchedThreads && underMaxTime && noPathology && noAnnForestInterrupt) {
      //each thread terminates by decrementing activeThreads.
      //by starting a loop iteration with a wait, and then setting launchedThreads to activeThreads
      //loop will exit once fMaxClusterNum activeThreads have completed running.
      searchCompleteCV.wait(guardCandidateClusters,
			    [&](){return !(activeThreads == launchedThreads);});
      launchedThreads = activeThreads; //activeThreads is decremented by terminated thread.
      newThreadLaunched = false;
      for (int i = 0; i < (numThreadsTotal-1); ++i) {
	if ((sigVec[i]) && (fCandClusters.candidates.size() < fMaxClusterNum)) {
	  //thread i has signaled completion, and we have not found enough candidate clusters.
	  //reset sigil and launch new thread.
	  newThreadLaunched = true;
	  sigVec[i] = false;
	  //check if the thread is joinable. if so, join it.
	  if ((threadsVec[i]).joinable()) {
	    (threadsVec[i]).join();
	  }
	  // if the user has set a random seed, progressively increment for next thread.
	  if (randomSeed > 0) {
	    currentSeed += 1;
	  }
	  //otherwise, set seed for next thread from random device.
	  else {
	    currentSeed = rdCCS();
	  }
	  ++activeThreads;
	  ++launchedThreads;
	  //launch new thread.
	  threadsVec[i] = std::thread(ccSearchThread,
				      std::ref(fDataVals),
				      std::ref(fRestrictedVals),
				      std::ref(fDipT),
				      std::ref(fClusterLB),
				      std::ref(fRepeatsAllowed),
				      std::ref(fMaxSearchDepth),
				      std::ref(fMaxClusterNum),
				      std::ref(fMaxNumberOfGates),
				      std::ref(ccMutex),
				      std::ref(fCandClusters),
				      std::ref(fRandomSearch),
				      std::ref(fUseRestrictedValue),
				      std::ref(badSearchCounter),
				      std::ref(annotationForestRun),
				      std::ref(pathologyLimit),
				      std::ref(sigVec),
				      i,
				      currentSeed,
				      std::ref(ccSearchParex),
				      ccSearchParexRoot,
				      std::ref(searchCompleteCV),
				      std::ref(activeThreads),
				      subSampleThreshold,
				      subSampleSize,
				      subSampleIterations);
	  if (printDebugInfo) {
	    std::cout << "Launched a new thread in slot: " << i << std::endl;
	    std::cout << "Used seed: " << currentSeed << std::endl;
	  }
	}
      }
      //check if we've had too many bad searches.
      //this check targets residual searches: numInSubset small.
      if (badSearchCounter > (9*(fMaxClusterNum/10))) {
	if (printDebugInfo) {
	  std::cout << "Search pathology: vacuous results. Abort and re-noise." << std::endl;
	}
	noPathology = false;
      }
      if (newThreadLaunched) {
	checkPoint = std::chrono::steady_clock::now();
	currentDuration = checkPoint-startSearch;
	elapsedTime = currentDuration.count();
	if (printDebugInfo) {
	  std::cout << "Elapsed search time: " << elapsedTime << std::endl;
	  if (annotationForestRun) {
	    std::cout << "Gates found: " << fCandClusters.gateLocations.size() << std::endl;
	  }
	  else{
	    std::cout << "Candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
	  }
	}
	if (elapsedTime > maxAllowedTime) {
	  underMaxTime = false;
	}
	if (annotationForestRun) {
	  if (fCandClusters.gateLocations.size() > fMaxClusterNum) {
	    noAnnForestInterrupt = false;
	  }
	}
      }
    }
    //join any dangling threads
    for (int i = 0; i < (numThreadsTotal-1); ++i) {
      if ((threadsVec[i]).joinable()) {
	(threadsVec[i]).join();
      }
    }
    //get final stats on search.
    if (printDebugInfo) {
      checkPoint = std::chrono::steady_clock::now();
      currentDuration = checkPoint-startSearch;
      elapsedTime = currentDuration.count();
      std::cout << "Final search time: " << elapsedTime << std::endl;
      if (annotationForestRun) {
	std::cout << "Total gates found: " << fCandClusters.gateLocations.size() << std::endl;
      }
      else{
	std::cout << "Total candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
      }
    }
    //finally, void existing results if a search pathology or time violation is signaled.
    if ((!noPathology) || (!underMaxTime)) {
      fCandClusters.candidates.clear();
      fCandClusters.gateLocations.clear();
      fCandClusters.abortIteration = true;
    }
  }
  else if (numThreadsTotal > 2) {
    std::mutex ccMutex;
    std::unique_lock<std::mutex> guardCandidateClusters(ccMutex);
    std::condition_variable searchCompleteCV;
    int launchedThreads = 0, activeThreads = 0;
    std::vector<bool> sigVec((numThreadsTotal-1),false);
    std::vector<std::thread> threadsVec;
    bool newThreadLaunched = false;
    //build up vector of starting roots.
    std::vector<bool> vRoots = determineAnnotationColumns(fDataVals,fDataVals.size(),
							  fDipT,fRestrictedVals,
							  fUseRestrictedValue);
    std::vector<int> vRootsIndex;
    for (int vri=0; vri !=vRoots.size();++vri)
      if (vRoots[vri])
	vRootsIndex.push_back(vri);
    if (printDebugInfo) {
      std::cout << "Searching across " << vRootsIndex.size() << " roots in forest." << std::endl;
    }
    int allRootNum = vRootsIndex.size();
    if (allRootNum <= 1) {
      guardCandidateClusters.unlock();
      //in the case of no viable or a single viable root, pass to exhaustiveSearch.
      fCandClusters =  candidateClusterSearch(fDataVals,
					      fRestrictedVals,
					      fDipT,
					      fClusterLB,
					      fRepeatsAllowed,
					      fMaxSearchDepth,
					      fMaxClusterNum,
					      fMaxNumberOfGates,
					      fRandomSearch,
					      fUseRestrictedValue,
					      annotationForestRun,
					      (randomSeed+1),
					      ccSearchParex,
					      ccSearchParexRoot,
					      subSampleThreshold,
					      subSampleSize,
					      subSampleIterations);
    }
    else {
      ccSearchParex = true;
      //check if there are fewer roots in the forest than threads available
      int usableThreadNum = std::min((allRootNum+1),static_cast<int>(numThreadsTotal));
      for (int i = 0; i < (usableThreadNum-1); ++i) {
	ccSearchParexRoot = vRootsIndex.back();
	vRootsIndex.pop_back();
	++activeThreads;
	++launchedThreads;
	if (printDebugInfo) {
	  std::cout << "fRandomSearch at launch: " << fRandomSearch << std::endl;
	  std::cout << "ccSearchParex at launch: " << ccSearchParex << std::endl;
	  std::cout << "ccSearchParexRoot at launch: " << ccSearchParexRoot << "." << std::endl;
	  std::cout << "launchedThreads at launch: " << launchedThreads << "." << std::endl;
	  std::cout << "activeThreads at launch: " << activeThreads << "." << std::endl;
	}
	threadsVec.push_back(std::thread(ccSearchThread,
					 std::ref(fDataVals),
					 std::ref(fRestrictedVals),
					 std::ref(fDipT),
					 std::ref(fClusterLB),
					 std::ref(fRepeatsAllowed),
					 std::ref(fMaxSearchDepth),
					 std::ref(fMaxClusterNum),
					 std::ref(fMaxNumberOfGates),
					 std::ref(ccMutex),
					 std::ref(fCandClusters),
					 std::ref(fRandomSearch),
					 std::ref(fUseRestrictedValue),
					 std::ref(badSearchCounter),
					 std::ref(annotationForestRun),
					 std::ref(pathologyLimit),
					 std::ref(sigVec),
					 i,
					 (randomSeed+1),
					 std::ref(ccSearchParex),
					 ccSearchParexRoot,
					 std::ref(searchCompleteCV),
					 std::ref(activeThreads),
					 subSampleThreshold,
					 subSampleSize,
					 subSampleIterations));
	if (printDebugInfo) {
	  std::cout << "launched thread " << i << std::endl;
	  std::cout << "Used seed " << (randomSeed+1) << std::endl;
	}
      }
      if (vRootsIndex.size() == 0) {
	if (printDebugInfo) {
	  std::cout << "More theads than roots. Spin off." << std::endl;
	  std::cout << "launchedThreads at start: " << launchedThreads << "." << std::endl;
	  std::cout << "activeThreads at start: " << activeThreads << "." << std::endl;
	}
	//we have launched searches along the entire space. wait for them to complete.
	while (launchedThreads) {
	  searchCompleteCV.wait(guardCandidateClusters,
				[&](){return !(activeThreads == launchedThreads);});
	  if (printDebugInfo) {
	    std::cout << " after alert, launchedThreads: " << launchedThreads << std::endl;
	    std::cout << " after alert, activeThreads: " << activeThreads << std::endl;
	  }
	  launchedThreads = activeThreads;
	  if (printDebugInfo) {
	    std::cout << " after reconcile, launchedThreads: " << launchedThreads << std::endl;
	    std::cout << " after reconcile, activeThreads: " << activeThreads << std::endl;
	  }
	}
	if (printDebugInfo)
	  std::cout << "search complete: " << launchedThreads << " remain." << std::endl;
      }
      else {
	//monitor thread pool, launching new threads along new roots as each search completes.
	while (launchedThreads && underMaxTime && noPathology && noAnnForestInterrupt) {
	  //as in the previous case, each thread terminates by decrementing activeThreads.
	  //here, by starting a loop iteration with a wait, and then setting launchedThreads to activeThreads
	  //loop will exit once activeThreads have explored the remainder of vRootsIndex.
	  searchCompleteCV.wait(guardCandidateClusters,
				[&](){return !(activeThreads == launchedThreads);});
	  launchedThreads = activeThreads;
	  newThreadLaunched = false;
	  for (int i = 0; i < (usableThreadNum-1); ++i) {
	    if ((sigVec[i]) && (vRootsIndex.size())) {
	      //thread i has signaled completion, and there are still unexplored
	      //roots in the forest. launch another thread
	      ccSearchParexRoot = vRootsIndex.back();
	      vRootsIndex.pop_back();
	      ++activeThreads;
	      ++launchedThreads;
	      //reset sigil and launch new thread.
	      newThreadLaunched = true;
	      sigVec[i] = false;
	      //check if the thread is joinable. if so, join it.
	      if ((threadsVec[i]).joinable()) {
		(threadsVec[i]).join();
	      }
	      //launch new thread.
	      threadsVec[i] = std::thread(ccSearchThread,
					  std::ref(fDataVals),
					  std::ref(fRestrictedVals),
					  std::ref(fDipT),
					  std::ref(fClusterLB),
					  std::ref(fRepeatsAllowed),
					  std::ref(fMaxSearchDepth),
					  std::ref(fMaxClusterNum),
					  std::ref(fMaxNumberOfGates),
					  std::ref(ccMutex),
					  std::ref(fCandClusters),
					  std::ref(fRandomSearch),
					  std::ref(fUseRestrictedValue),
					  std::ref(badSearchCounter),
					  std::ref(annotationForestRun),
					  std::ref(pathologyLimit),
					  std::ref(sigVec),
					  i,
					  (randomSeed+1),
					  std::ref(ccSearchParex),
					  ccSearchParexRoot,
					  std::ref(searchCompleteCV),
					  std::ref(activeThreads),
					  subSampleThreshold,
					  subSampleSize,
					  subSampleIterations);
	      if (printDebugInfo) {
		std::cout << "Launched a new thread in slot: " << i << std::endl;
		std::cout << "Used seed: " << (randomSeed+1) << std::endl;
	      }
	    }
	  }
	  //check if we've had too many bad searches.
	  //this check targets residual searches: numInSubset small.
	  if (badSearchCounter > (9*(fMaxClusterNum/10))) {
	    if (printDebugInfo) {
	      std::cout << "Search pathology: vacuous results. Abort and re-noise." << std::endl;
	    }
	    noPathology = false;
	  }
	  if (newThreadLaunched) {
	    checkPoint = std::chrono::steady_clock::now();
	    currentDuration = checkPoint-startSearch;
	    elapsedTime = currentDuration.count();
	    if (printDebugInfo) {
	      std::cout << "Elapsed search time: " << elapsedTime << std::endl;
	      if (annotationForestRun) {
		std::cout << "Gates found: " << fCandClusters.gateLocations.size() << std::endl;
	      }
	      else{
		std::cout << "Candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
	      }
	    }
	    if (elapsedTime > maxAllowedTime) {
	      underMaxTime = false;
	    }
	    if (annotationForestRun) {
	      if (fCandClusters.gateLocations.size() > fMaxClusterNum) {
		noAnnForestInterrupt = false;
	      }
	    }
	  }
	}
      }
      //join any dangling threads
      for (int i = 0; i < (usableThreadNum-1); ++i) {
	if ((threadsVec[i]).joinable()) {
	  (threadsVec[i]).join();
	}
      }
      //get final stats
      if (printDebugInfo) {
	checkPoint = std::chrono::steady_clock::now();
	currentDuration = checkPoint-startSearch;
	elapsedTime = currentDuration.count();
	std::cout << "Final search time: " << elapsedTime << std::endl;
	if (annotationForestRun) {
	  std::cout << "Total gates found: " << fCandClusters.gateLocations.size() << std::endl;
	}
	else{
	  std::cout << "Total candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
	}
      }
      //finally, void existing results if a search pathology or time violation is signaled.
      if ((!noPathology) || (!underMaxTime)) {
	fCandClusters.candidates.clear();
	fCandClusters.gateLocations.clear();
	fCandClusters.abortIteration = true;
      }
    }
  }
  else if (fRandomSearch) { 
    while ((fCandClusters.candidates.size() < fMaxClusterNum) && underMaxTime && noPathology && noAnnForestInterrupt) {
      // if the user has set a random seed, progressively increment for next search.
      if (randomSeed > 0) {
	currentSeed += 1;
      }
      //otherwise, set seed for next search from random device.
      else {
	currentSeed = rdCCS();
      }
      parResult = candidateClusterSearch(fDataVals,
					 fRestrictedVals,
					 fDipT,
					 fClusterLB,
					 fRepeatsAllowed,
					 fMaxSearchDepth,
					 fMaxClusterNum,
					 fMaxNumberOfGates,
					 fRandomSearch,
					 fUseRestrictedValue,
					 annotationForestRun,
					 currentSeed,
					 ccSearchParex,
					 ccSearchParexRoot,
					 subSampleThreshold,
					 subSampleSize,
					 subSampleIterations);
      //check for search pathology
      if ((parResult.abortIteration) && (numInSubset > pathologyLimit)) {
	//this condition is only true if the search didn't find a candidate cluster.
	++badSearchCounter;
	if (badSearchCounter > (9*(fMaxClusterNum/10))) {
	  //if too many searches fail, give up and re-noise
	  if (printDebugInfo) {
	    std::cout << "Search pathology: vacuous results. Abort and re-noise." << std::endl;
	  }
	  fCandClusters.candidates.clear();
	  fCandClusters.gateLocations.clear();
	  fCandClusters.abortIteration = true;
	  noPathology = false;
	}
      }
      //check time elapsed for search.
      checkPoint = std::chrono::steady_clock::now();
      currentDuration = checkPoint-startSearch;
      elapsedTime = currentDuration.count();
      if (elapsedTime > maxAllowedTime) {
	//if have run too long, terminate search
	fCandClusters.candidates.clear();
	fCandClusters.gateLocations.clear();
	fCandClusters.abortIteration = true;
	underMaxTime = false;
      }
      else {
	std::vector<std::vector<long>> parCands = parResult.candidates;
	std::vector<gateInfo> parGateLocs = parResult.gateLocations;
	for (auto cc : parCands) {
	  fCandClusters.candidates.push_back(cc);
	}
	for (auto gL : parGateLocs) {
	  fCandClusters.gateLocations.push_back(gL);
	}
	if (noPathology) {
	  fCandClusters.abortIteration = false;
	}
      }
      if (printDebugInfo) {
	std::cout << "Elapsed restricted search time: " << elapsedTime << std::endl;
	std::cout << "Current seed: " << currentSeed << std::endl;
	if (annotationForestRun) {
	  std::cout << "Gates found: " << fCandClusters.gateLocations.size() << std::endl;
	}
	else{
	  std::cout << "Candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
	}
      }
      if (annotationForestRun) {
	if (fCandClusters.gateLocations.size() > fMaxClusterNum) {
	  noAnnForestInterrupt = false;
	}
      }
    }
  }
  else {
    //pass in the random seed directly in this event -- it won't be used.
    fCandClusters =  candidateClusterSearch(fDataVals,
					    fRestrictedVals,
					    fDipT,
					    fClusterLB,
					    fRepeatsAllowed,
					    fMaxSearchDepth,
					    fMaxClusterNum,
					    fMaxNumberOfGates,
					    fRandomSearch,
					    fUseRestrictedValue,
					    annotationForestRun,
					    (randomSeed+1),
					    ccSearchParex,
					    ccSearchParexRoot,
					    subSampleThreshold,
					    subSampleSize,
					    subSampleIterations);
  }
  //whatever the search parameter, update the random seed with the current seed.
  if (randomSeed > 0) {
    randomSeed = currentSeed;
  }
  return fCandClusters;
}
