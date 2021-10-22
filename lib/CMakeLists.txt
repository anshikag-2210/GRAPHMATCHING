cmake_minimum_required(VERSION 3.0)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(graph
                MatchingAlgorithm.cc NProposingMatching.cc
                YokoiEnvyfreeHRLQ.cc MaximalEnvyfreeHRLQ.cc
                HHeuristicHRLQ.cc RHeuristicHRLQ.cc
                BipartiteGraph.cc Partner.cc PartnerList.cc PreferenceList.cc Vertex.cc
                GraphReader.cc Utils.cc 
                RelaxedStable.cc Statistics.cc SEAPopularHRLQ.cc 
                LpApproxSMFQ.cc Exact_Exponential_SMFQ.cc DirectApproachHR2LQ.cc
                Node.cc Edge.cc Graph.cc HApproxSMFQ.cc Convert_HR_to_HR2LQ.cc)
