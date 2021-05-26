#include "LpApproxSMFQ.h"
#include "HApproxSMFQ.h"
#include "Exact_Exponential_SMFQ.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include "StableMarriage.h"
#include "Statistics.h"
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <climits>
#include <iostream>

LpApproxSMFQ::LpApproxSMFQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

//check if output matching is relaxed stable
bool LpApproxSMFQ::is_envy_free(std::shared_ptr<BipartiteGraph> G,
    std::shared_ptr<MatchingAlgorithm::MatchedPairListType> M) {

    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;
    auto A_partition = G->get_A_partition();

    //for each resident check if it can be in a blocking pair
    for (auto& it : A_partition) {
        //resident
        auto u = it.second;
        auto& u_pref_list = u->get_preference_list();

        bookkeep_data[u] = VertexBookkeeping(0, u->get_preference_list().size());

        // while u hasn't exhausted its preference list
        while (not bookkeep_data[u].is_exhausted()) {
            //hospital to be checked if (r,h) is a blocking pair 
            auto v = u_pref_list.at(bookkeep_data[u].begin).vertex;
            auto v_pref_list = v->get_preference_list();

            //if this hospital is matched to this resident then no need to check for later hospitals
            auto M_u = M->find(u);
            if (M_u != M->end()) {
                auto& partners = M_u->second;
                // v is found in u's matchings
                if (partners.find(v) != partners.cend()) {
                    break;
                }
            }
         
            // if v has a matching compare u with v's worst partner
            if (M->find(v) != M->end()) {
                // v's worst partner
                auto v_worst_partner = M->at(v).get_least_preferred();
                auto v_worst_partner_rank = compute_rank(v_worst_partner.vertex, v_pref_list);
                auto u_rank = compute_rank(u, v_pref_list);
                //if v prefers its worst partner to u check for next hospital
                if (v_worst_partner_rank > u_rank) {
                    return false;
                }
            }
            
            //incrementing propose pointer of u
            bookkeep_data[u].begin += 1;
        }
    }
    return true;
}

int LpApproxSMFQ::get_lower_bound1(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& cost, std::map<VertexPtr, unsigned int> &index,
    std::vector<std::vector<bool>> &edges) {

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    int lb1 = 0;
   
    // for each resident find least cost hospital in its pref list
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto index_of_r = index[r];
        auto& r_pref_list = r->get_preference_list();

        int least_cost = INT_MAX;

        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto index_of_h = index[h];

            if (edges[index_of_r][index_of_h] == true) {
                if (cost[h] < least_cost) {
                    least_cost = cost[h];
                }
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }

        if(least_cost != INT_MAX)lb1 += least_cost;
    }
    return lb1;
}

int LpApproxSMFQ::prune_graph(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& cost, std::map<VertexPtr, unsigned int>& index,
    std::vector<std::vector<bool>> edges, VertexPtr r, VertexPtr h) {

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    auto index_of_r = index[r];
    auto& r_pref_list = r->get_preference_list();

    auto index_of_h = index[h];
    auto& h_pref_list = h->get_preference_list();
    
    bookkeep_data[h] = VertexBookkeeping(0, h->get_preference_list().size());

    // go through h's preference list
    while (not bookkeep_data[h].is_exhausted()) {
        //resident r1 in h's preflist
        auto r1 = h_pref_list.at(bookkeep_data[h].begin).vertex;
        auto index_of_r1 = index[r1];
        auto& r1_pref_list = r1->get_preference_list();

        // if r1 is r, stop the process
        if (index_of_r1 == index_of_r) {
            break;
        }
        // h prefers r1 over r
        else {
            bookkeep_data[r1] = VertexBookkeeping(0, r1->get_preference_list().size());
            bool start_deleting = false;

            // go through r1's preference list
            while (not bookkeep_data[r1].is_exhausted()) {
                //hospital h1 in r1's preflist
                auto h1 = r1_pref_list.at(bookkeep_data[r1].begin).vertex;
                auto index_of_h1 = index[h1];

                // h1 is less preferred than h by r1
                if (start_deleting) {
                    // delete the  edge (r1, h1)
                    edges[index_of_r1][index_of_h1] = false;
                }

                // if h is occurred in r1's preference list
                // then all next hospitals are less preferred than h
                if (index_of_h1 == index_of_h) {
                    start_deleting = true;
                }

                //incrementing propose pointer of r
                bookkeep_data[r1].begin += 1;
            }
        }

        //incrementing propose pointer of h
        bookkeep_data[h].begin += 1;
    }

    bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

    // go through r's preference list
    while (not bookkeep_data[r].is_exhausted()) {
        //hospital h1 in r's preflist
        auto h1 = r_pref_list.at(bookkeep_data[r].begin).vertex;
        auto index_of_h1 = index[h1];
        auto& h1_pref_list = h1->get_preference_list();

        // if h1 is h, stop the process
        if (index_of_h1 == index_of_h) {
            break;
        }
        // r prefers h1 over h
        else {
            bookkeep_data[h1] = VertexBookkeeping(0, h1->get_preference_list().size());
            bool start_deleting = false;

            // go through h1's preference list
            while (not bookkeep_data[h1].is_exhausted()) {
                //resident r1 in h1's preflist
                auto r1 = h1_pref_list.at(bookkeep_data[h1].begin).vertex;
                auto index_of_r1 = index[r1];

                // r1 is less preferred than r by h1
                if (start_deleting) {
                    // delete the  edge (r1, h1)
                    edges[index_of_r1][index_of_h1] = false;
                }

                // if r is occurred in h1's preference list
                // then all next residents are less preferred than r
                if (index_of_r1 == index_of_r) {
                    start_deleting = true;
                }

                //incrementing propose pointer of r
                bookkeep_data[h1].begin += 1;
            }
        }

        //incrementing propose pointer of r
        bookkeep_data[r].begin += 1;
    }

    // remove all other edges for r except h
    bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

    // go through r's preference list
    while (not bookkeep_data[r].is_exhausted()) {
        //hospital h1 in r's preflist
        auto h1 = r_pref_list.at(bookkeep_data[r].begin).vertex;
        auto index_of_h1 = index[h1];
        
        // if h1 is h, stop the process
        if (index_of_h1 == index_of_h) {
            edges[index_of_r][index_of_h1] = true;
        }
        else {
            edges[index_of_r][index_of_h1] = false;
        }

        //incrementing propose pointer of r
        bookkeep_data[r].begin += 1;
    }

    // pruning is completed
    // return lb1 of this pruned graph
    return get_lower_bound1(G, cost, index, edges);
}

void LpApproxSMFQ::print_lower_bounds(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& cost) {

    unsigned int R, H;
    R = G->get_A_partition().size();
    H = G->get_B_partition().size();

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    //to maintain indices of residents and hospitals
    std::map<VertexPtr, unsigned int> index;

    // to maintain edges of pruned graph for every tuple
    std::vector<std::vector<bool>> edges(R, std::vector<bool>(H, false));

    int index_count = 0;

    // find indices of residents
    for (auto& A1 : G->get_A_partition()) {
        auto u = A1.second;
        index[u] = index_count;
        index_count++;
    }

    index_count = 0;

    // find indices of hospitals
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        index[v] = index_count;
        index_count++;
    }

    // fill edges adjacency matrix
    // for each resident
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto index_of_r = index[r];
        auto& r_pref_list = r->get_preference_list();

        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto index_of_h = index[h];

            edges[index_of_r][index_of_h] = true;
               
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }

    // finding lb1 on initial graph
    int lb1 = get_lower_bound1(G, cost, index, edges);
    std::cout << " ,LB1 = " << lb1;

    // finding lb2 on initial graph
    int lb2 = 0;

    // for each resident
    for (auto& it : G->get_A_partition()) {
        //resident
        auto r = it.second;
        auto index_of_r = index[r];
        auto& r_pref_list = r->get_preference_list();

        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        int lb_of_r = INT_MAX;

        // for each hospital in r's pref list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto index_of_h = index[h];

            // prune graph for edge (r,h) 
            // and find lb1 on pruned graph
            int lb1_pruned_graph = prune_graph(G, cost, index, edges, r, h);

            if (lb1_pruned_graph < lb_of_r) {
                lb_of_r = lb1_pruned_graph;
            }

            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }

        if (lb_of_r > lb2) {
            lb2 = lb_of_r;
        }
    }

    std::cout << " ,LB2 = " << lb2;
}

// GCD of 'a' and 'b' 
void LpApproxSMFQ::print_additional_output(std::shared_ptr<BipartiteGraph> G
    , std::vector<std::vector<int>> &additional_output
    , std::vector<std::string> &additional_output_names)
{
    std::cout << "\nAdditional output\n";
    std::cout << "Hospital";
    //print every hospital name
    for (auto& B1 : G->get_B_partition()) {
        auto h = B1.second;
        std::cout << ", " << h->get_id();
    }
    std::cout << "\n";

    // print name and then additional vector
    for (int i = 0; i < additional_output.size(); i++) {
        std::cout << additional_output_names[i];
        for (int j = 0; j < additional_output[i].size(); j++) {
            std::cout << ", " << additional_output[i][j];
        }
        std::cout << "\n";
    }
}

// GCD of 'a' and 'b' 
int LpApproxSMFQ::gcd(int a, int b)
{
    if (b == 0)return a;
    return gcd(b, a % b);
}

/*
// to find costs of hospitals
// and to print input parameters
int LpApproxSMFQ::find_costs(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int> &cost, 
    std::vector<std::vector<int>> &additional_output) {

    unsigned int R, H, LR = 0, LH = 0, UC;

    R = G->get_A_partition().size();
    H = G->get_B_partition().size();

    // find max length of resident pref list
    for (auto& A1 : G->get_A_partition()) {
        auto u = A1.second;
        // if u's pref list is empty
        if (u->get_preference_list().size() == 0) {
            std::cout << "Resident " << u->get_id() << " has no preference list\n";
            return -1;
        }
        if (u->get_preference_list().size() > LR) {
            LR = u->get_preference_list().size();
        }
    }

    unsigned long int lcm = 0;
    // finding lcm's of upper quota's of all hospitals
    // find max length of hospital pref list
    for (auto& B1 : G->get_B_partition()) {
        auto v = B1.second;
        // if v's pref list is empty
        if (v->get_preference_list().size() == 0) {
            std::cout << "Hospital " << v->get_id() << " has no preference list\n";
            return -1;
        }
        if (v->get_preference_list().size() > LH) {
            LH = v->get_preference_list().size();
        }

        auto v_upper_quota = v->get_upper_quota();
        //for first hospital
        if (lcm == 0) {
            lcm = v_upper_quota;
        }
        else {
            lcm = (((v_upper_quota * lcm)) / (gcd(v_upper_quota, lcm)));
        }
    }

    // set the cost of every hospital
    // find number of distinct costs
    // collect the additional data about capacity and cost vectors 
    std::set<int, std::greater<int> > s1;
    std::vector<int> capacity_vector;
    std::vector<int> cost_vector;
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        cost[v] = (lcm / (v->get_upper_quota()));
        s1.insert(cost[v]);
        capacity_vector.push_back(v->get_upper_quota());
        cost_vector.push_back(cost[v]);
    }
    // number of distinct costs will be size of the set
    UC = s1.size();
    // insert capcity and cost vectors in additional output
    additional_output.push_back(capacity_vector);
    additional_output.push_back(cost_vector);

    //print in desired format
    std::cout << "Input Parameters\n";
    std::cout << "R = "<< R;
    std::cout << " ,H = " << H;
    std::cout << " ,LR = " << LR;
    std::cout << " ,LH = " << LH;
    std::cout << " ,UC = " << UC;
    print_lower_bounds(G, cost);
    std::cout << "\n\n";
    std::cout << "Output\n";
    std::cout << "Algo, Size, Cost, Rank1, Rank2, Max-Dev, Avg-Dev\n";
    return 1;
}*/

// to find costs of hospitals
// and to print input parameters
int LpApproxSMFQ::find_costs(std::shared_ptr<BipartiteGraph> G,
    std::map<VertexPtr, unsigned int>& cost,
    std::vector<std::vector<int>>& additional_output) {

    unsigned int R, H, LR = 0, LH = 0, UC;

    R = G->get_A_partition().size();
    H = G->get_B_partition().size();

    // find max length of resident pref list
    for (auto& A1 : G->get_A_partition()) {
        auto u = A1.second;
        // if u's pref list is empty
        if (u->get_preference_list().size() == 0) {
            std::cout << "Resident " << u->get_id() << " has no preference list\n";
            return -1;
        }
        if (u->get_preference_list().size() > LR) {
            LR = u->get_preference_list().size();
        }
    }

    std::vector<int> capacities;

    // finding maximum of upper quota's of all hospitals
    // find max length of hospital pref list
    for (auto& B1 : G->get_B_partition()) {
        auto v = B1.second;
        // if v's pref list is empty
        if (v->get_preference_list().size() == 0) {
            std::cout << "Hospital " << v->get_id() << " has no preference list\n";
            return -1;
        }
        if (v->get_preference_list().size() > LH) {
            LH = v->get_preference_list().size();
        }

        int v_upper_quota = v->get_upper_quota();
        if (std::find(capacities.begin(), capacities.end(), v_upper_quota) == capacities.end()) {
            capacities.push_back(v_upper_quota);
        }
    }

    // sort the distinct_costs_of_r vector
    sort(capacities.begin(), capacities.end(), std::greater<int>());
    
    // set the cost of every hospital
    // find number of distinct costs
    // collect the additional data about capacity and cost vectors 
    std::set<int, std::greater<int> > s1;
    std::vector<int> capacity_vector;
    std::vector<int> cost_vector;
    for (auto& it : G->get_B_partition()) {
        auto v = it.second;
        int v_upper_quota = v->get_upper_quota();
        cost[v] = std::find(capacities.begin(), capacities.end(), v_upper_quota) - capacities.begin() + 1;
        s1.insert(cost[v]);
        capacity_vector.push_back(v->get_upper_quota());
        cost_vector.push_back(cost[v]);
    }
    // number of distinct costs will be size of the set
    UC = s1.size();
    // insert capcity and cost vectors in additional output
    additional_output.push_back(capacity_vector);
    additional_output.push_back(cost_vector);

    //print in desired format
    std::cout << "Input Parameters\n";
    std::cout << "R = " << R;
    std::cout << " ,H = " << H;
    std::cout << " ,LR = " << LR;
    std::cout << " ,LH = " << LH;
    std::cout << " ,UC = " << UC;
    print_lower_bounds(G, cost);
    return 1;
}


std::shared_ptr<MatchingAlgorithm::MatchedPairListType> LpApproxSMFQ::compute_matching() {
    
    std::shared_ptr<BipartiteGraph> G = get_graph();
    auto M = std::make_shared<MatchingAlgorithm::MatchedPairListType>();

    //to maintain costs of hospitals
    std::map<VertexPtr, unsigned int> cost;

    //to maintain if a hospital is min_cost hospital for some resident
    std::map<VertexPtr, bool> min_cost_hospital;

    // to maintain capacity vector, cost vector for each algorithm
    std::vector<std::vector<int>> additional_output;
    
    // to maintain row names for additional output
    std::vector<std::string> additional_output_names;

    //find cost of hospitals
    //if hospital or resident has no pref list, find_costs method returns -1
    if (!find_costs(G, cost, additional_output)) {
        return M;
    }

    unsigned int min_max_lb = 0;
    //find H approx matching
    HApproxSMFQ alg1(G, true);
    auto Mh = alg1.compute_matching1(min_max_lb);
    std::cout << " ,MINMAX_LB = " << min_max_lb;

    // printing output heading
    std::cout << "\n\n";
    std::cout << "Output\n";
    std::cout << "Algo, Size, Cost, Rank1, Rank2, Max-Dev, Avg-Dev\n";

    additional_output_names.push_back("Capacity");
    additional_output_names.push_back("Cost");

    Statistics s;

    //find stable matching
    StableMarriage alg(G, true);
    auto Ms = alg.compute_matching();
    // print stable matching statistics
    s.get_smfq_statistics(G, Ms, Ms, "SM", cost, additional_output);
    additional_output_names.push_back("SM");

    // VertexBookkeeping for maintaining propose pointer
    // to know current proposing position in pref list 
    std::map<VertexPtr, VertexBookkeeping> bookkeep_data;

    auto A_partition = G->get_A_partition();
 
    // mark all hospitals as not min cost for any resident
    for (auto& it : G->get_B_partition()) {
        auto h = it.second;
        min_cost_hospital[h] = false;
    }

    // for each resident find min cost hospital
    for (auto& it : A_partition) {
        //resident r
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        //min cost of any hospital in this resident's pref list
        unsigned long int min_cost = 0;

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital h in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            
            // if h has less than min cost 
            // or it is the first hospital
            // set min_cost as cost of that hospital
            if (cost[h] < min_cost || min_cost == 0) {
                min_cost = cost[h];
            }
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
        //resetting proposing pointer of resident r
        bookkeep_data[r].begin = 0;

        // now mark the min cost hospital
        // as appearing first in r's pref list with min_cost
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital h in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            //if h has min cost
            if (cost[h] == min_cost) {
                min_cost_hospital[h] = true;
                break;
            }
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }

    // match each resident to most preferred hospital h
    // which has min_cost_hospital[h] true
    for (auto& it : A_partition) {
        //resident
        auto r = it.second;
        auto& r_pref_list = r->get_preference_list();
        bookkeep_data[r] = VertexBookkeeping(0, r->get_preference_list().size());

        // while r hasn't exhausted its preference list
        while (not bookkeep_data[r].is_exhausted()) {
            //hospital in r's preflist
            auto h = r_pref_list.at(bookkeep_data[r].begin).vertex;
            auto h_pref_list = h->get_preference_list();

            //if h has min cost
            if (min_cost_hospital[h] == true) {
                // add r and h to the matching
                add_partner(M, h, r, compute_rank(r, h_pref_list), 0);
                add_partner(M, r, h, compute_rank(h, r_pref_list), 0);
                break;
            }
            //incrementing propose pointer of r
            bookkeep_data[r].begin += 1;
        }
    }
    s.get_smfq_statistics(G, M, Ms, "Lp-approx", cost, additional_output);
    additional_output_names.push_back("Lp-approx");

    /*
    //find exact exponential matching
    Exact_Exponential_SMFQ alg1(G, true);
    auto Me = alg1.compute_matching();
    // print exact exponential matching statistics
    s.get_smfq_statistics(G, Me, Ms, "Exact Exp", cost, additional_output);
    additional_output_names.push_back("Exact Exp");
    */

    // print H Appx matching statistics
    s.get_smfq_statistics(G, Mh, Ms, "H Approx", cost, additional_output);
    additional_output_names.push_back("H Approx");

    //printing additional outputs
    print_additional_output(G, additional_output, additional_output_names);

    // printing matchings
    std::cout << "\nStable Matching\n";
    if (!is_envy_free(G, Ms)) {
        std::cout << "Not Envy Free\n";
    }
    print_matching(G, Ms, std::cout);

    std::cout << "\nLp_Approx Matching\n";
    if (!is_envy_free(G, M)) {
        std::cout << "Not Envy Free\n";
    }
    print_matching(G, M, std::cout);

    /*
    std::cout << "\nExact Exp Matching\n";
    if (!is_envy_free(G, Me)) {
        std::cout << "Not Envy Free\n";
    }
    print_matching(G, Me, std::cout);
    */

    std::cout << "\nH_Approx Matching\n";
    if (!is_envy_free(G, Mh)) {
        std::cout << "Not Envy Free\n";
    }
    print_matching(G, Mh, std::cout);

    return M;
}
