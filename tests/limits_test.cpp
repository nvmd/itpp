#include <limits>
#include <iostream>

using std::cout;
using std::endl;
using std::numeric_limits;

int main()
{
 
 cout << "================================" << endl;
 cout << "    Test of Numerical Limits   " << endl;
 cout << "================================" << endl;
 
 cout<< "char = " << numeric_limits<char>::digits << "bits" << endl;
 cout<< "unsigned char = " << numeric_limits<unsigned char>::digits << "bits" << endl;
 
 cout<< "short = " << numeric_limits<short>::digits << "bits" << '\n';
 cout<< "unsigned short = " << numeric_limits<unsigned short>::digits << "bits" << endl;
 
 cout<< "int = " << numeric_limits<int>::digits << "bits" << '\n';
 cout<< "unsigned int = " << numeric_limits<unsigned int>::digits << "bits" << endl;
 
 cout<< "long int = " << numeric_limits<long int>::digits << "bits" << endl;
 
 cout<< "long = " << numeric_limits<long>::digits << "bits" << endl;
 cout<< "unsigned long = " << numeric_limits<unsigned long>::digits << "bits" << endl;
 
 cout<< "float = " << numeric_limits<float>::digits << "bits" << endl;
 cout<< "double = " << numeric_limits<double>::digits << "bits" << endl;
 //cout<< "long double = " << numeric_limits<long double>::digits << "bits" << endl;
 
 return 0;

}
