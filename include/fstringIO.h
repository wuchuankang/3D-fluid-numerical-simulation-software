
#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include <string>

using std::ifstream;
using std::ofstream;
using std::runtime_error;
using std::string;

//read a text line from file, srips any comment, leading and trailing spaces
//assume that the comment starts with the harsh(#) sign

inline string ReadLine(ifstream &stream) 
{
	string str;
	string::size_type ic;	//iterator (int)
	if (stream.good())
	{
		getline(stream, str);

		//get rid of comment
		ic = str.find_first_of('#');
		if (ic != string::npos) 
			str.erase( str.begin()+ic, str.end() );
		
		//get rid of leading spaces
		ic = str.find_first_not_of(' ');
		if ((int)ic > 0)
			str.erase( str.begin(), str.begin()+ic );

		//get of rid of trailing spaces
		ic = str.find_last_not_of(' ');
		if ((int)ic != string::npos )
			str.erase( str.begin()+ic+1, str.end());
	}
	else
		throw runtime_error("could not read line from file");

	return str;
}
