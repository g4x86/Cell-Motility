#ifndef OUTPUTFILE_HPP_
#define OUTPUTFILE_HPP_

#include <string>
#include <fstream>

namespace motility
{

class OutputFile
{
  private:

	std::string name;

	std::ofstream stream;

  public:

	OutputFile(const std::string& fn);

	virtual ~OutputFile();

	std::string getName() const;

	std::ofstream& getStream();
};

}

#endif /*OUTPUTFILE_HPP_*/
