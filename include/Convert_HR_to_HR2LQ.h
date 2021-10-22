#ifndef CONVERT_HR_TO_HR2LQ
#define CONVERT_HR_TO_HR2LQ

#include "MatchingAlgorithm.h"
#include "PreferenceList.h"
#include "VertexBookkeeping.h"
#include <vector>
#include <queue>
#include <ostream>

// Converts input HR instance to HR2LQ and prints the HR2LQ instance in output file
class Convert_HR_to_HR2LQ : public MatchingAlgorithm {
public:
    explicit Convert_HR_to_HR2LQ(std::shared_ptr<BipartiteGraph> G, bool A_proposing = false);
    ~Convert_HR_to_HR2LQ() override = default;
    
    std::shared_ptr<MatchedPairListType> compute_matching() override;
};

#endif
