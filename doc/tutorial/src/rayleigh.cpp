#include "itcomm.h"

//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declare scalars and vectors:
  double noise_var, norm_dopp;
  int Nsamples;
  Array<cvec> chan;

  //Open the file rayleigh_test.it
  it_file ff("rayleigh_test.it");

  //Declare the variable my_channel as an instance of the class Rayleigh_Channel:
  TDL_Channel my_channel;

  //The normalized Doppler frequency is set to 0.1. 
  norm_dopp = 0.1;
  my_channel.set_norm_doppler(norm_dopp);

  //Generate Nsamples of the fading process:
  Nsamples = 10000;
  my_channel.generate(Nsamples, chan);

  //Save the results to file:
  ff << Name("chan") << chan;
  ff.close();

  //Exit program:
  return 0;

}
