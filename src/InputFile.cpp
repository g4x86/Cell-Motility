#include <algorithms.hpp>
#include <InputFile.hpp>

namespace motility
{

InputFile::InputFile(const std::string& fn)
{
	name = fn;
	stream.open(name.c_str());
	if(!stream.is_open())
	{
		std::stringstream msg;
		msg << "failed to read from " << name;
		handleErrorEvent(msg);
	}
}

InputFile::~InputFile()
{
	if(stream.is_open()) stream.close();
}

std::string InputFile::getName() const
{
	return name;
}

std::ifstream& InputFile::getStream()
{
	return stream;
}

}
