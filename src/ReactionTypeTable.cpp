#include <ReactionTypeTable.hpp>

namespace motility
{

ReactionTypeTable::ReactionTypeTable() {}

ReactionTypeTable::ReactionType::ReactionType()
{
	forward_const = 0;
	backward_const = 0;
}

void ReactionTypeTable::ReactionType::clear()
{
	reactants.clear();
	products.clear();
	forward_const = 0;
	backward_const = 0;
}

ReactionTypeTable::Table& ReactionTypeTable::instance()
{
	static ReactionTypeTable rtt;
	return rtt.table;
}

}
