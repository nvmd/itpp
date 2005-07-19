#include "itpp/itbase.h"

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
    CPU_Timer t1;
    Real_Timer t2;

    tic();
    t1.start();
    t2.start();

    while (t2.get_time() < 4.0)
			;

    t1.stop();
    t2.stop();
    double t3 = toc();

    if (t2.get_time() > 4.0)
      cout << "Real_Timer is OK" << endl;

    if (t1.get_time() < t2.get_time())
      cout << "CPU_Timer is OK" << endl;

    if (t3 >= t2.get_time())
      cout << "tic() and toc() is OK" << endl;

    return 0;
}
