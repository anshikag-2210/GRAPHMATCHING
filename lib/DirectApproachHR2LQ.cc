#include "DirectApproachHR2LQ.h"
#include "BipartiteGraph.h"
#include "StableMarriage.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include <set>
#include <iostream>

DirectApproachHR2LQ::DirectApproachHR2LQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

bool DirectApproachHR2LQ::level_proposing(std::shared_ptr<MatchingAlgorithm::MatchedPairListType>& M,
    const BipartiteGraph::ContainerType& proposing_partition) {

    // queue to maintain proposing vertices
    FreeListType free_list;

    // VertexBookkeeping for maintaining propose pointer and level of proposing vertices
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    // input graph
    std::shared_ptr<BipartiteGraph> G = get_graph();

    // to maintain max level possible of proposing vertices
    unsigned int max_level = 0;

    // set the level of every proposing vertex to 0
    // mark all deficient proposing vetices free (by pushing into the free_list)
    // set max_level to sum of lq's of all proposing vetices
    for (auto& it : proposing_partition) {
        auto v = it.second;
        max_level += v->get_lower_quota();
        // sets level 0
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
        // if v is deficient and not present in free list(queue)
        if (number_of_partners(M, v) < v->get_lower_quota()) {
            if (bookkeep_data[v].in_free_list == false) {
                // add it to free list
                bookkeep_data[v].in_free_list = true;
                free_list.push(v);
            }
        }
    }

    // while there is at least one proposing vertex in the free list
    while (not free_list.empty()) {
        // proposing vertex at the front in free list
        auto u = free_list.front();
        //remove u from freelist
        free_list.pop();
        bookkeep_data[u].in_free_list = false;

        auto& u_pref_list = u->get_preference_list();

        // if u hasn't exhausted its preference list
        if (not bookkeep_data[u].is_exhausted()) {
            // highest ranked resident/hospital to whom u has not yet proposed in this level
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();

            //if v is already fully subscribed
            if (number_of_partners(M, v) == v->get_upper_quota()) {
                // v's least preferred partner
                auto v_partner = M->at(v).get_least_preferred();
                auto possible_partner = Partner(u, compute_rank(u, v_pref_list), bookkeep_data[u].level);

                // least preferred partner is compared to present proposing vertex
                // Note : < operator is overloaded
                // if v prefers u over least preferred in its matchings
                if (v_partner < possible_partner) {
                    // remove v_partner from v's matchings and vice versa
                    M->at(v).remove_least_preferred();
                    M->at(v_partner.vertex).remove(v);

                    // add v and u to the matching
                    add_partner(M, u, v, compute_rank(v, u_pref_list), 0);
                    add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);

                    // add v_partner to free_list if it is not already present
                    if (bookkeep_data[v_partner.vertex].in_free_list == false) {
                        // add it to free list
                        bookkeep_data[v_partner.vertex].in_free_list = true;
                        free_list.push(v_partner.vertex);
                    }
                }
                // if v doesnt prefer u
                else {
                    // add u to free list
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            //if v has vacancies
            else {
                // add v and u to the matching
                add_partner(M, u, v, compute_rank(v, u_pref_list), 0);
                add_partner(M, v, u, compute_rank(u, v_pref_list), bookkeep_data[u].level);
            }
            // if u has vacancies and its level is 0
            if (bookkeep_data[u].level == 0 && number_of_partners(M, u) < u->get_upper_quota()) {
                if (bookkeep_data[u].in_free_list == false) {
                    // add it to free list
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            // if u is deficient
            else if (u->get_lower_quota() > 0 && number_of_partners(M, u) < u->get_lower_quota()) {
                if (bookkeep_data[u].in_free_list == false) {
                    // add it to free list
                    bookkeep_data[u].in_free_list = true;
                    free_list.push(u);
                }
            }
            //increment u's propose pointer
            bookkeep_data[u].begin++;
        }
        // if u is a lq hospital and exhausted its pref list and
        // its level is less than max level
        else if (u->get_lower_quota() > 0 && bookkeep_data[u].level < max_level) {
            // increment h's level
            bookkeep_data[u].level += 1;
            // reset proposal index
            bookkeep_data[u].begin = 0;
            // add it to free list
            bookkeep_data[u].in_free_list = true;
            free_list.push(u);
        }

    }

    // if any proposing vertex is deficient
    // then feasible matching not possible
    for (auto& it : proposing_partition) {
        auto v = it.second;
        if (number_of_partners(M, v) < v->get_lower_quota()) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> DirectApproachHR2LQ::compute_matching() {
    
    std::shared_ptr<BipartiteGraph> G = get_graph();
    
    //find stable matching
    StableMarriage alg(G, false);
    auto M = alg.compute_matching();

    // hospital proposing phase
    if (!level_proposing(M, G->get_B_partition())) {
        std::cout << "Feasible matching not possible\n";
    }

    // resident proposing phase
    if (!level_proposing(M, G->get_A_partition())) {
        std::cout << "Feasible matching not possible\n";
    }
    
    return M;
}
