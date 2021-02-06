#include "LpApproxSMFQ.h"
#include "Vertex.h"
#include "Partner.h"
#include "Utils.h"
#include "StableMarriage.h"
#include "Statistics.h"
#include <set>
#include <vector>
#include <string>
#include <iterator>
#include <iostream>

LpApproxSMFQ::LpApproxSMFQ(std::shared_ptr<BipartiteGraph> G,
    bool A_proposing)
    : MatchingAlgorithm(std::move(G), A_proposing)
{}

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
    std::cout << " ,UC = " << UC<<"\n\n";
    std::cout << "Output\n";
    std::cout << "Algo, Size, Cost, Rank1, Rank2, Max-Dev, Avg-Dev\n";
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

    //printing additional outputs
    print_additional_output(G, additional_output, additional_output_names);

    // printing matchings
    std::cout << "\nStable Matching\n";
    print_matching(G, Ms, std::cout);
    std::cout << "\nLp_Approx Matching\n";
    print_matching(G, M, std::cout);

    return M;
}
