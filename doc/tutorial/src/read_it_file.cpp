#include "itcomm.h"

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declare a vector:
  vec a;	

  //Declare the it_file class:
  it_file ff;

  //Open the file it_file_test.it for reading:
  ff.open("it_file_test.it");

  //Read the variable a from the file. Put result in vector a:
  ff >> Name("a") >> a;

  //Print the result:
  cout << "a = " << a << endl;

  //Exit the program:
  return 0;

}
