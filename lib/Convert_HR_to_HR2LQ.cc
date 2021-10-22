#include "Convert_HR_to_HR2LQ.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include "StableMarriage.h"
#include "Popular.h"
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <climits>
#include <iostream>

Convert_HR_to_HR2LQ::Convert_HR_to_HR2LQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> Convert_HR_to_HR2LQ::compute_matching() {

    // Input graph
    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();

    // compute stable matching
    StableMarriage alg(G, false);
    auto Ms = alg.compute_matching();

    // compute maximum matching
    PopularAmongMaxCard alg1(G, false);
    auto M_max = alg1.compute_matching();

    // For each resident 
    std::cout<<"@PartitionA\n";
    auto size = G->get_A_partition().size();
    for (auto& it : G->get_A_partition()) {
        auto r = it.second;

        //If resident is unmatched in M_s
        // and is matched in M_max, print it as a lq resident
        if (number_of_partners(Ms, r) == 0 && number_of_partners(M_max, r) > 0) {
            std::cout<<r->get_id()<<" ("<<number_of_partners(M_max, r)<<", "<<r->get_upper_quota()<<")";
        }
        // else print it as non lq resident
        else{
            std::cout<<r->get_id()<<" ("<<r->get_lower_quota()<<", "<<r->get_upper_quota()<<")";
        }
        size--;
        if(size){
            std::cout<<", ";
        }
        else{
            std::cout<<";\n@End\n\n";
        }
    }

    // For each hospital
    std::cout<<"@PartitionB\n";
    size = G->get_B_partition().size();
    for (auto& it : G->get_B_partition()) {
        auto h = it.second;
        
        // if hospital has less matched partners in M_s
        // than in M_max, print its lq value in between those 2 values
        if (number_of_partners(Ms, h) < number_of_partners(M_max, h)) {
            std::cout<<h->get_id()<<" ("<<number_of_partners(M_max, h)<<", "<<h->get_upper_quota()<<")";
        }
        else{
            std::cout<<h->get_id()<<" ("<<h->get_lower_quota()<<", "<<h->get_upper_quota()<<")";
        }
        size--;
        if(size){
            std::cout<<", ";
        }
        else{
            std::cout<<";\n@End\n\n";
        }
    }

    // print preference lists of residents
    std::cout<<"@PreferenceListsA\n";
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        std::cout<<r->get_id()<<": ";

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            std::cout << h->get_id();
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;

            if(not bookkeep_data[r].is_exhausted()){
                std::cout<<", ";
            }
            else{
                std::cout<<";\n";
            }
        }
    }
    std::cout<<"@End\n\n";

    // print preference lists of hospitals
    std::cout<<"@PreferenceListsB\n";
    for (auto& it : G->get_B_partition()) {
        //hospital
        auto h = it.second;
        auto& h_pref_list = h->get_preference_list();
        bookkeep_data[h] = VertexBookkeeping(0, h->get_preference_list().size());

        std::cout<<h->get_id()<<": ";

        // while h hasn't exhausted its preference list
        while (not bookkeep_data[h].is_exhausted()) {
            //resident in h's preflist
            auto r = h_pref_list.at(bookkeep_data[h].begin).vertex;
            std::cout << r->get_id();
            //incrementing propose pointer of r
            bookkeep_data[h].begin += 1;

            if(not bookkeep_data[h].is_exhausted()){
                std::cout<<", ";
            }
            else{
                std::cout<<";\n";
            }
        }
    }
    std::cout<<"@End";

    return M;
}
