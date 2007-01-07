#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

extern int main(int argc, char **argv)
{ 
  int Nbits=1000*1000*5000;          // maximum number of bits simulated for any SNR point
  int Nbers=2000;                   // target number of bit errors per SNR point
  double BERmin=1e-6;               // BER at which to terminate simulation
  vec EbN0db = "0.6:0.2:5";

  LDPC_Code C(argv[1]);
  bool single_snr_mode=false;
  if (argc==3) {
    double x;
    sscanf(argv[2],"%lf",&x);
    EbN0db.set_size(1);
    EbN0db(0)=x;
    single_snr_mode=true;
  }

  cout << "Running with Eb/N0: " << EbN0db << endl;

  // High performance: 2500 iterations, high resolution LLR algebra
  C.setup_decoder("bp","2500 1 0");  

  // Alt. setting -- High speed: 50 iterations, logmax approximation
  // C.setup_decoder("bp","50 1 0",LLR_calc_unit(12,0,7));  

  cout << C << endl;

  int N = C.get_nvar();             // number of bits per codeword
  BPSK Mod;
  bvec bitsin = zeros_b(N);
  vec s = Mod.modulate_bits(bitsin);
  
  RNG_randomize();
  for (int j=0; j<length(EbN0db); j++) {
    // noise variance is N0/2 per dimension
    double N0 = pow(10.0,-EbN0db(j)/10.0) / C.get_rate();
    AWGN_Channel chan(N0/2);
    BERC berc; // counters for coded and uncoded BER
    BLERC ferc; // counter for coded FER
    ferc.set_blocksize(N);
    for (long int i=0; i<Nbits; i+=C.get_nvar())   {
      // Received data
      vec x = chan(s);
      
      // Demodulate
      vec softbits=Mod.demodulate_soft_bits(x,N0);

      //Decode the received bits
      bvec bitsout=C.decode(softbits);

      //Count the number of errors
      berc.count(bitsin,bitsout);
      ferc.count(bitsin,bitsout);      

      if (single_snr_mode) {
	cout << "Eb/N0=" << EbN0db(j) << "  Simulated " 
	     << ferc.get_total_blocks() << " frames and " 
	     << berc.get_total_bits() << " bits. " 
	     << "Obtained " << berc.get_errors() << " bit errors. " 
	     << " BER: " << berc.get_errorrate() 
	     << " FER: " << ferc.get_errorrate() << endl;
	cout.flush();
      }      else 	{
	  if (berc.get_errors()>Nbers) { break;}	
      }
    }

    cout << "Eb/N0=" << EbN0db(j) << "  Simulated " 
	 << ferc.get_total_blocks() << " frames and " 
	 << berc.get_total_bits() << " bits. " 
	 << "Obtained " << berc.get_errors() << " bit errors. " 
	 << " BER: " << berc.get_errorrate() 
	 << " FER: " << ferc.get_errorrate() << endl;
    cout.flush();
    if (berc.get_errorrate()<BERmin)  {  break;  }
  } 
  return 0;
}
