#include <itpp/itcomm.h>


using std::cout;
using std::endl;
using namespace itpp;


bvec set_errors(const bvec &input, const ivec errpos) {
 bvec output=input ;
for (int i=0;i<errpos.length();i++)
 output(errpos(i))^=1;
return output ;
}


int main(void) {
 BCH bch(31, 21, 2, "3 5 5 1");

 cout << "==================================" << endl;
 cout << "      Test of BCH encoder/decoder" << endl;
 cout << "==================================" << endl;

 cout << "A two error case (should be corrected)" << endl;
 bvec input=randb(21);
 bvec encoded=bch.encode(input);
 bvec err=set_errors(encoded, (ivec) "1 2") ; // error positions
 bvec decoded=bch.decode(err);

 cout << "input=   " <<input << endl ;
 cout << "encoded= " << encoded << endl ;
 cout << "err=     " << err <<  endl ;
 cout << "decoded= " << decoded << endl ;
 
 input=randb(21);
 encoded=bch.encode(input);
 err=set_errors(encoded, (ivec) "1 2 27") ; // error positions
 decoded=bch.decode(err);

 cout << "A three error case (will cause decoding errors)" << endl;
 cout << "input=   " <<input << endl ;
 cout << "encoded= " << encoded << endl ;
 cout << "err=     " << err <<  endl ;
 cout << "decoded= " << decoded << endl ;
 
 return 0;
}

