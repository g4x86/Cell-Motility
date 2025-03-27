#include <ParameterTable.hpp>

namespace motility
{

ParameterTable::ParameterTable() {}

ParameterTable::~ParameterTable() {}

std::map<std::string, std::string>& ParameterTable::instance()
{
	static ParameterTable pt;
	return pt.table;
}

}
