#ifndef TOKENITERATOR_HPP_
#define TOKENITERATOR_HPP_

#include <string>
#include <iterator>
#include <functional>
#include <algorithm>

namespace motility
{

// This algorithm is provided by Thinking in C++ (Vol2)

struct Isalpha
{ 
	bool operator()(char c);
};

class Delimiters
{
	std::string exclude;

  public:

	Delimiters();

	Delimiters(const std::string& excl);

	bool operator()(char c);
};

template <typename InputIter, typename Pred = Isalpha>
class TokenIterator
{
	InputIter first;
	InputIter last;
	std::string word;
	Pred predicate;

  public:

    using iterator_category = std::input_iterator_tag;
    using value_type = std::string;
    using difference_type = std::ptrdiff_t;
    using pointer = std::string*;
    using reference = std::string&;

	TokenIterator(InputIter begin, InputIter end, Pred pred = Pred()) 
	{
		first = begin;
		last = end;
		predicate = pred;
		operator++();
	}

	TokenIterator() {}

	TokenIterator& operator++()
	{
		word.resize(0);
		first = find_if(first, last, predicate);
		while (first != last && predicate(*first)) word += *first++;
		return *this;
	}

	class Proxy
	{
		std::string word;

	  public:

		Proxy(const std::string& w)
		{
			word = w;
		}

		std::string operator*()
		{
			return word;
		}
	};

	Proxy operator++(int)
	{
		Proxy d(word);
		operator++();
		return d; 
	}

	std::string operator*() const
	{
		return word;
	}

	std::string* operator->() const
	{
		return &(operator*()); 
	}

	bool operator==(const TokenIterator&)
	{
		return (word.size() == 0 && first == last);
	}

	bool operator!=(const TokenIterator& rv)
	{
		return !(*this == rv);
	}
};

}

#endif /*TOKENITERATOR_HPP_*/
