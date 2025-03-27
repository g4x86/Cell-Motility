#ifndef PARAMETERTABLE_HPP_
#define PARAMETERTABLE_HPP_

#include <string>
#include <map>

namespace motility
{

class ParameterTable
{

	std::map<std::string, std::string> table;

	ParameterTable();

	virtual ~ParameterTable();

  public:

	typedef std::map<std::string, std::string> Table;

	static std::map<std::string, std::string>& instance();

};

}

#endif /*PARAMETERTABLE_HPP_*/
