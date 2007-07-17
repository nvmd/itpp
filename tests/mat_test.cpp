#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main(void)
{
 mat A = randn(3,4);
 mat B = randn(3,4);
 mat C = randn(4,3);
 vec v = randn(4);
 vec u = randn(3);

 mat E = randn(1,3);

 mat M;
 imat iM;
 smat sM;
 cmat cM;
 bmat bM;
 mat D;
 vec s;

 cout << "====================================================" << endl;
 cout << "mat_test: Test of matrix operations in mat.h/mat.cpp" << endl;
 cout << "====================================================" << endl;
 cout << "A = " << A << endl;
 cout << "B = " << B << endl;
 cout << "C = " << C << endl;
 cout << "E = " << E << endl << endl;

 cout << "Testing indexing" << endl;
 cout << "================" << endl;
 cout << "A(1,2) = " << A(1,2) << endl;
 cout << "A(2,3) = " << A(2,3) << endl;
 cout << "A(6) = " << A(6) << endl;
 cout << "A(0,2,1,3) = " << A(0,2,1,3) << endl;
 cout << "A.get_row(1) = " << A.get_row(1) << endl;
 cout << "A.get_rows(1,2) = " << A.get_rows(1,2) << endl;
 cout << "A.get_col(2) = " << A.get_col(2) << endl;
 cout << "A.get_cols(2,3) = " << A.get_cols(2,3) << endl << endl;

 cout << "Testing initialisation with string" << endl;
 cout << "==================================" << endl;
 cout << "sM = \"0xFF, -021 ,   100; 0,-0x01; 0xA, 10 012;  \"" << endl;
 sM = "0xFF, -021 ,   100; 0,-0x01; 0xA, 10 012;  ";
 cout << "sM = " << sM << endl;
 cout << "cM = \" 1+3i, (.33,1) ;  (333,-1) 2-0.2E6i\"" << endl;
 cM = " 1+3i, (.33,1) ;  (333,-1) 2-0.2E6i";
 cout << "cM = " << cM << endl;
 cout << "bM = \" 1 1 0; 0 1; 1 1 ,1,   1; ; 0 1\"" << endl;
 bM = " 1 1 0; 0 1; 1 1 ,1,   1; ; 0 1";
 cout << "bM = " << bM << endl << endl;

 cout << "Testing setting/copying/swapping" << endl;
 cout << "================================" << endl;
 cout << "v = " << v << endl;
 cout << "u = " << u << endl;
 D = A;
 cout << "D = " << D << endl;
 D.set_row(1, v);
 cout << "D.set_row(1, v):" << endl;
 cout << "D = " << D << endl;
 D.set_col(2, u);
 cout << "D.set_col(2, u):" << endl;
 cout << "D = " << D << endl;
 D.copy_row(1, 2);
 cout << "D.copy_row(1, 2):" << endl;
 cout << "D = " << D << endl;
 D.copy_col(2, 3);
 cout << "D.copy_col(2, 3):" << endl;
 cout << "D = " << D << endl;
 D.swap_rows(0, 2);
 cout << "D.swap_rows(0, 2):" << endl;
 cout << "D = " << D << endl;
 D.swap_cols(0, 3);
 cout << "D.swap_cols(0, 3):" << endl;
 cout << "D = " << D << endl;

 D.set_submatrix(1,2,2,3, A(0,1,0,1));
 cout << "D.set_submatrix(1,2,2,3, A(0,1,0,1): " << endl;
 cout << "D = " << D << endl;
 D.set_submatrix(1,2,2,3, 5.3);
 cout << "D.set_submatrix(1,2,2,3, 5.3): " << endl;
 cout << "D = " << D << endl;

 cout << "A = " << A << endl;
 cout << "B = " << B << endl;
 D = concat_horizontal(A,B);
 cout << "D = concat_horizontal(A,B):" << endl;
 cout << "D = " << D << endl;
 D = concat_vertical(A,B);
 cout << "D = concat_vertical(A,B):" << endl;
 cout << "D = " << D << endl << endl;

 cout << "Testing operators (=, +, - *, /)" << endl;
 cout << "================================" << endl;
 D = A;
 cout << "D = A = " << D << endl;
 D = 5.3;
 cout << "D = 5.3 = " << D << endl;

 D = A+B;
 cout << "D = A+B = " << D << endl;
 D += A;
 cout << "D += A = " << D << endl;
 D = A+5.3;
 cout << "D = A+5.3 = " << D << endl;
 D += 5.3;
 cout << "D += 5.3 = " << D << endl;

 D = A-B;
 cout << "D = A-B = " << D << endl;
 D -= A;
 cout << "D -= A = " << D << endl;
 D = A-5.3;
 cout << "D = A-5.3 = " << D << endl;
 D -= 5.3;
 cout << "D -= 5.3 = " << D << endl;

 D = A*C;
 cout << "D = A*C = " << D << endl;
 D *= A;
 cout << "D *= A = " << D << endl;
 D = A*5.3;
 cout << "D = A*5.3 = " << D << endl;
 D = 3.9*A;
 cout << "D = 3.9*A = " << D << endl;
 D *= 5.3;
 cout << "D *= 5.3 = " << D << endl;
 D = elem_mult(A,B);
 cout << "D = elem_mult(A,B) = " << D << endl;

 s = A*v;
 cout << "s = A*v = " << s << endl;
 D = u*E;
 cout << "D = u*E = " << D << endl;

 D = A/5.3;
 cout << "D = A/5.3 = " << D << endl;
 D /= 5.3;
 cout << "D /= 5.3 = " << D << endl;
 D = elem_div(A,B);
 cout << "D = elem_div(A,B) = " << D << endl;
 D /= A;
 cout << "D /= A = " << D << endl;

 return 0;

}
