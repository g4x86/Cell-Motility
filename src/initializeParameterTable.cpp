#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ParameterTable.hpp>
#include <algorithms.hpp>
#include <TokenIterator.hpp>
#include <initializeParameterTable.hpp>

namespace motility
{

void initializeParameterTable(std::ifstream& input)
{
	ParameterTable::Table& param_table = ParameterTable::instance();
	while (!input.eof())
	{
		std::string buf;
		getline(input, buf);
		if(buf.size() > 0)
		{
			if(buf[0] != '[' && buf[0] != ';')
			{
				Delimiters delim("=");
				TokenIterator<std::string::const_iterator, Delimiters> wordIter(buf.begin(), buf.end(), delim);
				std::string keyword = *wordIter++;
				std::string keyvalue = *wordIter++;
				param_table[keyword] = keyvalue;
			}
		}
	}
}

}
