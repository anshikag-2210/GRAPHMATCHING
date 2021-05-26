#include "Edge.h"

Edge::Edge()
{}

Edge::Edge(NodePtr U, NodePtr V, int weight /* = 0 */, bool matched /* = false */)
	: U_(U), V_(V), weight_(weight), matched_(matched)
{}

Edge::~Edge()
{}

const NodePtr& Edge::get_U() const {
	return U_;
}

const NodePtr& Edge::get_V() const {
	return V_;
}

const int& Edge::get_weight() const {
	return weight_;
}

void Edge::set_weight(int weight) {
	weight_ = weight;
}

bool Edge::is_matched() {
	return matched_;
}

void Edge::set_matched(bool matched) {
	matched_ = matched;
}