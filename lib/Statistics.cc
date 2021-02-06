#include "Statistics.h"
#include "StableMarriage.h"
#include "Vertex.h"
#include "PartnerList.h"
#include "Utils.h"
#include <map>
#include <vector>
#include <iostream>

Statistics::Statistics()
{
    //initially set all statistics to zero
    S = 0;
    BPC = 0;
    BR = 0;
    R1 = 0;
    def = 0;
}

Statistics::~Statistics()
{}

//To get the statistics of matching with respect to graph
void Statistics::get_statistics(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M) {
    
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    auto A_partition = G->get_A_partition();

    //for each resident
    //set in_free_list to false for every resident
    // we use in_free_list as in_blocking_pair
    for (auto& it : A_partition) {
        auto v = it.second;
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
        bookkeep_data[v].in_free_list = false;
    }
    //for each resident check if it can be in a blocking pair
    for (auto& it : A_partition) {
        //resident
        auto u = it.second;
        auto& u_pref_list = u->get_preference_list();

        // while u hasn't exhausted its preference list
        while (not bookkeep_data[u].is_exhausted()) {
            //hospital to be checked if (r,h) is a blocking pair 
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();
            //incrementing propose pointer of u
            bookkeep_data[u].begin += 1;

            //if this hospital is matched to this resident then no need to check for later hospitals
            auto M_u = M->find(u);
            if (M_u != M->end()) {
                auto& partners = M_u->second;
                // v is found in u's matchings
                // no need to check for other hospitals in this pref list
                if (partners.find(v) != partners.cend()) {
                    // v's rank on u's preference list
                    auto v_rank = compute_rank(v, u_pref_list);
                    //if v is ranked 1 in u's preference list
                    if (v_rank == 1) {
                        R1++;
                    }
                    break;
                }
            }

            //if v is fullysubscribed
            if (number_of_partners(M, v) == v->get_upper_quota()) {
                // v's worst partner
                auto v_worst_partner = M->at(v).get_least_preferred();
                auto v_worst_partner_rank = compute_rank(v_worst_partner.vertex, v_pref_list);

                auto u_rank = compute_rank(u, v_pref_list);

                //if v prefers its worst partner to u, check for next hospital
                if (v_worst_partner_rank < u_rank) {
                    continue;
                }
            }
            //if v is not fully subscribed or
            //if v prefers u to its worst partner then (u,v) is a blocking pair
            bookkeep_data[u].in_free_list = true;
            //increment number of blocking pairs
            BPC++;
        }
        // if u is matched, increment size of matching
        if (number_of_partners(M, u) > 0) {
            S++;
        }
        // if this resident is in some blocking pair
        if (bookkeep_data[u].in_free_list == true) {
            BR++;
        }
    }

    StableMarriage alg(G, true);
    auto M1 = alg.compute_matching();

    //for each hospital
    auto B_partition = G->get_B_partition();

    unsigned long int sum_of_lower_quota = 0;
    for (auto& B1 : B_partition) {
        auto v = B1.second;
        auto v_lower_quota = v->get_lower_quota();
        sum_of_lower_quota += v_lower_quota;
        auto M_v = M1->find(v);
        //if v is matched to some vertices
        if (M_v != M1->end()) {
            //partner list of v
            auto& partners = M_v->second;
            if (partners.size() < v_lower_quota) {
                def += v_lower_quota - partners.size();
            }
        }
        //if v is not matched to any vertex
        else {
            def += v_lower_quota;
        }
    }
    std::cout << S << "," << BPC << "," << BR << "," << R1 << "," << def << ","<<sum_of_lower_quota<<"\n";
}

//To get the statistics of smfq matching with respect to graph
void Statistics::get_smfq_statistics(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> Ms,
    char* alg, std::map<VertexPtr, unsigned int>& cost,
    std::vector<std::vector<int>>& additional_output) {

    auto A_partition = G->get_A_partition();
    auto B_partition = G->get_B_partition();

    unsigned long int size_of_matching = 0;
    unsigned long int cost_of_matching = 0;
    unsigned long int rank1 = 0;
    unsigned long int rank2 = 0;
    unsigned long int max_dev = 0;
    unsigned long int avg_dev = 0;

    //for each resident r
    for (const auto& it : A_partition) {
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        auto M_r = M->find(r);

        // if r is matched, increment size of matching
        if (number_of_partners(M, r) > 0) {
            size_of_matching++;
        }

        if (M_r != M->end()) {
            auto& partners = M_r->second;

            for (const auto& i : partners) {
                //hospital  h matched to r
                auto h = i.vertex;
                cost_of_matching += cost[h];
                // h's rank on r's preference list
                auto h_rank = compute_rank(h, r_pref_list);
                //if h is ranked 1 in r's preference list
                if (h_rank == 1) {
                    rank1++;
                }
                //if h is ranked 2 in r's preference list
                if (h_rank == 2) {
                    rank2++;
                }
            }
        }
    }

    // to store number of allocations of each hospital
    std::vector<int> alloc_vector;
    //for each hospital h
    for (const auto& it : B_partition) {
        auto h = it.second;
        alloc_vector.push_back(number_of_partners(M, h));
        int dev_of_h = number_of_partners(M, h) - number_of_partners(Ms, h);
        if (dev_of_h < 0)dev_of_h = dev_of_h * -1;

        avg_dev += dev_of_h;
        if (dev_of_h > max_dev) {
            max_dev = dev_of_h;
        }
    }

    //push alloc_vector into additional output
    additional_output.push_back(alloc_vector);

    avg_dev = avg_dev / B_partition.size();
    std::cout << alg << "," << size_of_matching << "," << cost_of_matching << "," 
        << rank1 << "," << rank2 << "," << max_dev << "," << avg_dev << "\n";
}