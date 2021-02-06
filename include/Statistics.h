#ifndef STATISTICS_H
#define STATISTICS_H

#include "BipartiteGraph.h"
#include "MatchingAlgorithm.h"
#include "PartnerList.h"
#include <vector>

class Statistics {
public:
    //size of matching
    unsigned long int S;
    //number of blocking pairs
    unsigned long int BPC;
    //number of blocking residents
    unsigned long int BR;
    //number of residents matched to their rank 1 hospitals
    unsigned long int R1;
    //Deficiency of graph
    unsigned long int def;

public:
    Statistics();
    virtual ~Statistics();
    void get_statistics(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M);
    void get_smfq_statistics(std::shared_ptr<BipartiteGraph> G,
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M, 
        std::shared_ptr<MatchingAlgorithm::MatchedPairListType> Ms,
        char *alg, std::map<VertexPtr, unsigned int> &cost,
        std::vector<std::vector<int>>& additional_output);
};

#endif
