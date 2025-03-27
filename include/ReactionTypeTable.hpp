#ifndef REACTIONTYPETABLE_HPP_
#define REACTIONTYPETABLE_HPP_

#include <string>
#include <list>
#include <map>

namespace motility
{

class ReactionTypeTable
{

  public:

	struct ReactionType
	{
		std::list<std::string> reactants, products;
		double forward_const, backward_const;

		ReactionType();
		void clear();
	};

	typedef std::map<std::string, ReactionType> Table;

  private:

	Table table;

	ReactionTypeTable();

  public:

	static Table& instance();

};

}

#endif /*REACTIONTYPETABLE_HPP_*/
