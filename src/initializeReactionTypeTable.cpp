#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ReactionTypeTable.hpp>
#include <algorithms.hpp>
#include <TokenIterator.hpp>
#include <initializeReactionTypeTable.hpp>

namespace motility
{

void getWord(const std::string& buf, std::string& word, const char* delims, size_t pos)
{
	Delimiters delimiters(delims);
	TokenIterator<std::string::const_iterator, Delimiters> wordIter(buf.begin(), buf.end(), delimiters);
	for(size_t i = 0; i < pos; ++i) *wordIter++;
	word = *wordIter++;
}

void getWords(const std::string& buf, std::list<std::string>& words, const char* delims, size_t pos)
{
	Delimiters delimiters(delims);
	TokenIterator<std::string::const_iterator, Delimiters> wordIter(buf.begin(), buf.end(), delimiters), end;
	for(size_t i = 0; i < pos; ++i) *wordIter++;
	while(wordIter != end) words.push_back(*wordIter++);
}

void initializeReactionTypeTable(std::ifstream& input)
{
	ReactionTypeTable::Table& reactionTable = ReactionTypeTable::instance();
	ReactionTypeTable::ReactionType rt;
	size_t line_cnt = 0;
	std::string name, rc;
	while (!input.eof())
	{
		std::string buf;
		getline(input, buf);
		if(buf.size() > 0 && buf[0] != '[' && buf[0] != ';')
		{
			switch(line_cnt)
			{
				case 0 :
				{
					getWord(buf, name, "=", 1);
					++line_cnt;
					break;
				}
				case 1 :
				{
					getWords(buf, rt.reactants, "=,", 1);
					++line_cnt;
					break;
				}
				case 2 :
				{
					getWords(buf, rt.products, "=,", 1);
					++line_cnt;
					break;
				}
				case 3 :
				{
					getWord(buf, rc, "=", 1);
					rt.forward_const = strtod(rc);
					rc.clear();
					++line_cnt;
					break;
				}
				case 4 :
				{
					getWord(buf, rc, "=", 1);
					rt.backward_const = strtod(rc);
					rc.clear();
					++line_cnt;
					break;
				}
			}
			if(line_cnt > 4)
			{
				reactionTable[name] = rt;
				name.clear();
				rt.clear();
				line_cnt = 0;
			}
		}
	}
}

}
