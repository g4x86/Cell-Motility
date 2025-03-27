#include <algorithms.hpp>
#include <OutputFile.hpp>

namespace motility
{

OutputFile::OutputFile(const std::string& fn)
{
	name = fn;
	stream.open(name.c_str());
	if(!stream.is_open())
	{
		std::stringstream msg;
        msg << "failed to write to " << name;
		handleErrorEvent(msg);
	}
}

OutputFile::~OutputFile()
{
	if(stream.is_open()) stream.close();
}

std::string OutputFile::getName() const
{
	return name;
}

std::ofstream& OutputFile::getStream()
{
	return stream;
}

}
