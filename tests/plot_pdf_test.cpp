#include <itpp/itgraphics.h>

using std::cout;
using std::endl;
using namespace itpp;


#ifdef HAVE_HARU

int main()
{
  PlotPDF plot("plot");

  vec x=linspace(0,1.0, 100);
  vec y1= sin(2*pi*x);
  vec y2= 0.7*sin(4*pi*x);

  plot.set_plot(x,y1);
  plot.set_plot(x,y2);
  
  // line style dashed
  plot.set_line_style(1, DASHED);

  // set axis range
  plot.set_axis(0, 1.0, -1.0, 1.0);

  // grid on
  plot.grid_on();

  plot.set_xlabel("x");
  plot.set_ylabel("sin(2*pi*x) and 0.7sin(4*pi*x)");
  plot.set_title("sin(2*pi*x) and 0.7sin(4*pi*x)");

  plot.save_pdf();
}

#else
int main() { cout << "libharu is needed for this test program" << endl; }

#endif
