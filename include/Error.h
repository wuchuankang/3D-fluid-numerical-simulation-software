
#ifndef ERROR_H
#define ERROR_H
class Error
{
public:
  /// Writes given error message to the standard output and exits the program.
  ///
  /// @param message  text of the error message
  ///
  static void Message( const char *message )
  {
    using namespace std;
    cout << "Error: " << message << endl << endl;
    exit(EXIT_FAILURE);
  }
};

#endif