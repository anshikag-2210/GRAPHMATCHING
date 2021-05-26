#include "Node.h"

Node::Node()
{}

Node::Node(IdType node_name, NodeType node_id) 
	: node_name_(node_name), node_id_(node_id)
{}

Node::~Node()
{}

const IdType& Node::get_name() const {
	return node_name_;
}

const NodeType& Node::get_id() const {
	return node_id_;
}
