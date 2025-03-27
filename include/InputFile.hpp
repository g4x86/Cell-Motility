#ifndef INPUTFILE_HPP_
#define INPUTFILE_HPP_

#include <string>
#include <fstream>

namespace motility
{

class InputFile
{
  private:

	std::string name;

	std::ifstream stream;

  public:

	InputFile(const std::string& fn);

	virtual ~InputFile();

	std::string getName() const;

	std::ifstream& getStream();
};

}

#endif /*INPUTFILE_HPP_*/
