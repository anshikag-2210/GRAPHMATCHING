#ifndef T_DEFS_H
#define T_DEFS_H

#include <map>
#include <memory>

// forward declaration
class Vertex;
class Node;

/// Id type for a vertex
typedef std::string IdType;

// Rank representation for a vertex in the preference/partner list
// this must be unique for all the vertices in a list
typedef int RankType;

// Pointer type for vertices
typedef std::shared_ptr<Vertex> VertexPtr;

/// Id type for a node in Graph
typedef int NodeType;

// Pointer type for nodes
typedef std::shared_ptr<Node> NodePtr;

#endif
