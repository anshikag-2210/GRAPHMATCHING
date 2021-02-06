#ifndef RELAXED_STABLE_H
#define RELAXED_STABLE_H

#include "MatchingAlgorithm.h"
#include "PreferenceList.h"
#include "VertexBookkeeping.h"
#include <stack>
#include <ostream>

// to compute 3/2 approximation of MAXRSM
class RelaxedStable : public MatchingAlgorithm {
private:
    // typedef list for unmatched vertices
    typedef std::stack<VertexPtr> FreeListType;

public:
    explicit RelaxedStable(std::shared_ptr<BipartiteGraph> G, bool A_proposing);
    ~RelaxedStable() override = default;

    std::unique_ptr<BipartiteGraph> get_modified_graph();
    std::shared_ptr<MatchedPairListType> compute_matching() override;
    bool is_relaxed_stable(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchedPairListType> M);
};

#endif

