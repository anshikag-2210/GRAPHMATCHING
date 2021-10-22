#include "SEAPopularHRLQ.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include <set>
#include <iostream>

SEAPopularHRLQ::SEAPopularHRLQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

std::shared_ptr<MatchingAlgorithm::MatchedPairListType> SEAPopularHRLQ::compute_matching() {
    // queue
    FreeListType free_list;
    // VertexBookkeeping for maintaining propose pointer and level of hospitals
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();
    max_level = G->get_A_partition().size();

    // set the level of every hospital to 0
    // mark all hospitals free (by pushing into the free_list)
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        free_list.push(v);
        // sets level 0
        bookkeep_data[v] = VertexBookkeeping(0, v->get_preference_list().size());
        bookkeep_data[v].in_free_list = true;
    }

    // while there is at least one hospital in the free list
    while (not free_list.empty()) {
        // hospital at the front in free list
        auto h = free_list.front();
        auto h_id = h->get_id();
        //remove h from freelist
        free_list.pop();
        bookkeep_data[h].in_free_list = false;

        auto& h_pref_list = h->get_preference_list();

        // if h hasn't exhausted its preference list
        if (not bookkeep_data[h].is_exhausted()) {
            // highest ranked resident to whom h has not yet proposed
            auto r = h_pref_list.at(bookkeep_data[h].begin).vertex;
            auto r_pref_list = r->get_preference_list();
            
            //if resident r is already matched
            if (number_of_partners(M, r) == r->get_upper_quota()) {
                auto r_partner = M->at(r).get_least_preferred();
                auto possible_partner = Partner(h, compute_rank(h, r_pref_list), bookkeep_data[h].level);

                // if level[resident] < level[hospital]
                // or if both levels same and resident prefers hospital h
                // to its partner
                // Note : < operator is overloaded
                if (r_partner < possible_partner) {
                    //////invariant
                    //if (possible_partner.level - r_partner.level > 2) {
                      //  std::cout << r->get_id() <<" jumped from level "
                       //     <<r_partner.level<<" to "<<possible_partner.level<<"\n";
                    //}
                    /////////////////////////////////////////////////

                    // remove r_partner from r and vice versa
                    M->at(r).remove_least_preferred();
                    M->at(r_partner.vertex).remove(r);

                    // add r and h to the matching
                    add_partner(M, h, r, compute_rank(r, h_pref_list), 0);
                    add_partner(M, r, h, compute_rank(h, r_pref_list), bookkeep_data[h].level);

                    // add r_partner to free_list if it is not already present
                    if (bookkeep_data[r_partner.vertex].in_free_list == false) {
                        // add it to free list
                        bookkeep_data[r_partner.vertex].in_free_list = true;
                        free_list.push(r_partner.vertex);
                    }
                }
                // if r doesnt prefer h
                else {
                    // add it to free list
                    bookkeep_data[h].in_free_list = true;
                    free_list.push(h);
                }
                //if resident r is not matched
            }
            else {
                // add r and h to the matching
                add_partner(M, h, r, compute_rank(r, h_pref_list), 0);
                add_partner(M, r, h, compute_rank(h, r_pref_list), bookkeep_data[h].level);
            }
            // if h has vacancies and its level is 0
            if (bookkeep_data[h].level == 0 && number_of_partners(M, h) < h->get_upper_quota()) {
                if (bookkeep_data[h].in_free_list == false) {
                    // add it to free list
                    bookkeep_data[h].in_free_list = true;
                    free_list.push(h);
                }
            }
            // if h is lower quota hospital and is deficient
            else if (h->get_lower_quota() > 0 && number_of_partners(M, h) < h->get_lower_quota()) {
                if (bookkeep_data[h].in_free_list == false) {
                    // add it to free list
                    bookkeep_data[h].in_free_list = true;
                    free_list.push(h);
                }
            }
            //increment h's propose pointer
            bookkeep_data[h].begin++;
        }
        // if h is a lq hospital and exhausted its pref list and
        // its level is less than number of residents
        else if (h->get_lower_quota() > 0 && bookkeep_data[h].level < max_level) {
            // increment h's level
            bookkeep_data[h].level += 1;

            // reset proposal index
            bookkeep_data[h].begin = 0;
            // add it to free list
            bookkeep_data[h].in_free_list = true;
            free_list.push(h);
        }

        //////////////// invariant
        //auto m_h = m->find(h);
        //std::set<int> s;
        ////if h is matched to some residents
        //if (m_h != m->end()) {
        //    //partner list of h
        //    auto& partners = m_h->second;
        //    for (const auto& i : partners) {
        //        auto res_in_h = i.vertex;
        //        auto r_partner = m->at(res_in_h).get_least_preferred();
        //        s.insert(r_partner.level);
        //        if (s.size() == 2) {
        //            auto first = s.begin(); 
        //            std::advance(first, 1);
        //            if (abs(*(s.begin()) - *(first)) > 1) {
        //                std::cout << "error\n";
        //                break;
        //            }
        //        }
        //    }
        //    if (s.size() > 2) {
        //        std::cout << "error\n";
        //    }
        //}
        /////////////////////////////////////
    }

    return M;
}
