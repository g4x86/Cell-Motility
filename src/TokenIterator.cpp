#include <string>
#include <TokenIterator.hpp>


namespace motility
{

bool Isalpha::operator()(char c)
{ 
	return std::isalpha(c); 
}

Delimiters::Delimiters() {}

Delimiters::Delimiters(const std::string& excl)
{
	exclude = excl;
}

bool Delimiters::operator()(char c)
{
	return exclude.find(c) == std::string::npos;
}

}
