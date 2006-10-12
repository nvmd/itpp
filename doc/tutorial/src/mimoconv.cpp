
#include <itpp/itcomm.h>

using std::cout;
using std::endl;
using namespace itpp;
using namespace std;

// ---- Zero-forcing detector with soft information --------

void ZF_demod(ND_UPAM &channel, ivec &LLR_apr,  ivec &LLR_apost, double sigma2, mat &H, vec &y)
{
  it_assert(H.rows()>=H.cols(),"ZF_demod() - underdetermined systems not tolerated");
  vec shat=ls_solve_od(H,y);
  vec Sigma2=diag(inv(H.transpose()*H))*sigma2;
  vec h(length(shat));
  for (int i=0; i<length(shat); i++) {
    shat(i) = shat(i)/sqrt(Sigma2(i));
    h(i) = 1.0/sqrt(Sigma2(i));
  }
  channel.map_demod(LLR_apr,LLR_apost,1.0,h,shat);
}


extern int main(int argc, char **argv)
{ 
  // -- modulation and channel parameters --
  int nC;                    // type of constellation  (1=QPSK, 2=16-QAM, 3=64-QAM)
  int nRx;                   // number of receive antennas
  int nTx;                   // number of transmit antennas 
  int Tc;                    // coherence time (number of channel vectors with same H)

  if (argc!=5) {
    cout << "Usage: cm nTx nRx nC Tc" << endl << "Example: cm 2 2 1 100000 (2x2 QPSK MIMO on slow fading channel)" << endl;
    exit(1); 
  } else {
    sscanf(argv[1],"%i",&nTx);
    sscanf(argv[2],"%i",&nRx);
    sscanf(argv[3],"%i",&nC);
    sscanf(argv[4],"%i",&Tc);
  }

  cout << "Initializing.. " << nTx << " TX antennas, " << nRx << " RX antennas, " 
       << (1<<nC) << "-PAM per dimension, coherence time " << Tc << endl;

  // -- simulation control parameters --
  const vec EbN0db = "-5:2:50";
  const int Nmethods =2;               // number of demodulators to try
  const long int Nbitsmax=500*1000*1000;     // maximum number of bits to ever simulate per SNR point

  int Nbers, Nfers;              // target number of bit/frame errors per SNR point
  double BERmin, FERmin;         // BER/FER at which to terminate simulation
  if (Tc==1) {           // Fast fading channel, BER is of primary interest here
    BERmin = 0.001;     // stop simulating a given method if BER<this value
    FERmin = 1.0e-10;    // stop simulating a given method if FER<this value
    Nbers = 1000;        // move to next SNR point after counting 1000 bit errors
    Nfers = 100000;      // move to next SNR point after counting 100000 frame errors (in practice never)
  } else {               // Slow fading channel, FER is of primary interest here
    BERmin = 1.0e-15;    // stop simulating a given method if BER<this value
    FERmin = 0.01;      // stop simulating a given method if FER<this value
    Nbers = -1;          // do not stop on this condition
    Nfers = 200;         // move to next SNR point after counting 200 frame errors
  }

  // -- Channel code parameters --
  Convolutional_Code code;
  ivec generator(3);
  double rate=1.0/3.0;
  generator(0)=0133;
  generator(1)=0165;
  generator(2)=0171;
  code.set_generator_polynomials(generator, 7);
  const int Nu = 100;  // block length before coding
  bvec dummy;
  code.encode_tail(randb(Nu),dummy);
  const int Nc = length(dummy);

  // ============= Initialize ====================================

  const int Nctx = (int) (2*nC*nTx*ceil(double(Nc)/double(2*nC*nTx)));   // Total number of bits to transmit   
  const int Nvec = Nctx/(2*nC*nTx);                         // Number of TX vectors
  const int Nbitspvec = 2*nC*nTx;                            // Number of bits per TX vector

  // initialize MIMO channel
  ND_UPAM chan;
  chan.set_Gray_PAM(2*nTx,1<<nC);  
  cout << chan << endl;

  //  RNG_randomize();

  Array<vec> Y(Nvec);
  Array<mat> H(Nvec/Tc+1);

  ivec Contflag = ones_i(Nmethods);
  if (pow(2.0,nC*2.0*nTx)>256) {   // ML decoder too complex..
    Contflag(1)=0;  
  }
  
  // ================== Run simulation =======================
  for (int nsnr=0; nsnr<length(EbN0db); nsnr++) {
    double N0 = pow(10.0,-EbN0db(nsnr)/10.0) / rate;
    Array<BERC> berc(Nmethods);  // counter for coded BER
    Array<BERC> bercu(Nmethods); // counter for uncoded BER
    Array<BLERC> ferc(Nmethods); // counter for coded FER

    for (int i=0; i<Nmethods; i++) {
      ferc(i).set_blocksize(Nu);
    }

    long int nbits=0;
    while (nbits<Nbitsmax) {
      nbits += Nu;

      // generate and encode random data
      bvec inputbits = randb(Nu);                
      bvec txbits;
      code.encode_tail(inputbits, txbits);
      txbits=concat(txbits,zeros_b(Nctx-Nc));

      // -- generate channel and data ----
      for (int k=0; k<Nvec; k++) {
	if (k%Tc==0) {       // generate a new channel realization
	  mat Hr = 1/sqrt(2.0)*randn(nRx,nTx);
	  mat Hi = 1/sqrt(2.0)*randn(nRx,nTx);
	  H(k/Tc).set_size(2*nRx,2*nTx);
	  H(k/Tc).set_submatrix(0,nRx-1,0,nTx-1,Hr);
	  H(k/Tc).set_submatrix(nRx,2*nRx-1,0,nTx-1,Hi);
	  H(k/Tc).set_submatrix(0,nRx-1,nTx,2*nTx-1,-Hi);
	  H(k/Tc).set_submatrix(nRx,2*nRx-1,nTx,2*nTx-1,Hr);
	}
	
	// modulate and transmit bits
	bvec bitstmp = txbits(k*2*nTx*nC,(k+1)*2*nTx*nC-1);
	vec x=chan.modulate_bits(bitstmp);
	vec e=sqrt(N0/2.0)*randn(2*nRx);
	Y(k) = H(k/Tc)*x+e;
      }

      // -- run decoder loop ---     

      Array<QLLRvec> LLRin(Nmethods);     
      for (int i=0; i<Nmethods; i++) {
	LLRin(i) = zeros_i(Nctx);
      }

      for (int k=0; k<Nvec; k++) {           // RX data vectors
	QLLRvec llr_apr =zeros_i(nC*2*nTx);
	QLLRvec llr_apost =zeros_i(nC*2*nTx);
	     
	// zero forcing demodulation
	if (Contflag(0)) {
	  ZF_demod(chan,llr_apr,llr_apost,N0/2.0,H(k/Tc),Y(k));
	  LLRin(0).set_subvector(k*Nbitspvec,(k+1)*Nbitspvec-1,llr_apost);
	}
	  
	// ML demodulation
	if (Contflag(1)) { 
	  chan.map_demod(llr_apr, llr_apost, N0/2.0, H(k/Tc), Y(k));
	  LLRin(1).set_subvector(k*Nbitspvec,(k+1)*Nbitspvec-1,llr_apost);
	}	  
      }
            
      // Run decoder and count errors
      for (int i=0; i<Nmethods; i++) {
	bvec decoded_bits;
	if (Contflag(i)) {
	  // QLLR values must be converted to real numbers since decoder wants this
	  vec llr=chan.get_llrcalc().to_double(LLRin(i).left(Nc));
	  bercu(i).count(txbits(0,Nc-1),llr(0,Nc-1)>0);
	  code.decode_tail(-llr,decoded_bits);
	  berc(i).count(inputbits(0,Nu-1),decoded_bits(0,Nu-1));
	  ferc(i).count(inputbits(0,Nu-1),decoded_bits(0,Nu-1));
	}
      }	
      
      // check whether it is time to terminate the simulation
      // (terminate when all methods still running have counted at
      // least Nbers or Nfers bit resp. frame errors)
      int minber=1000000;
      int minfer=1000000;
      for (int i=0; i<Nmethods; i++) {
	if (Contflag(i)) {
	  minber=min(minber,round_i(berc(i).get_errors()));   
	  minfer=min(minfer,round_i(ferc(i).get_errors()));  
	}
      }
      if (Nbers>0 && minber>Nbers) { break;}
      if (Nfers>0 && minfer>Nfers) { break;}

    }
    
    cout << "----------------------" << endl;
    cout << "Eb/N0: " << EbN0db(nsnr) << " simulated " << nbits << " bits." << endl;
    cout << " Uncoded BER: " << bercu(0).get_errorrate() << " (ZF);     " << bercu(1).get_errorrate() << " (ML)" << endl;
    cout << " Coded BER: " << berc(0).get_errorrate() << " (ZF);     " << berc(1).get_errorrate() << " (ML)" << endl;
    cout << " FER: " << ferc(0).get_errorrate() << " (ZF);     " << ferc(1).get_errorrate() << " (ML)" << endl;
    cout.flush();

    // check wheter it is time to stop simulating (stop when all
    // methods have reached the min BER/FER of interest)
    int contflag=0;
    for (int i=0; i<Nmethods; i++) {
      if (Contflag(i)) {
	if (berc(i).get_errorrate()>BERmin)  {  contflag=1;  } else { Contflag(i)= 0; } 
	if (ferc(i).get_errorrate()>FERmin)  {  contflag=1;  } else { Contflag(i)= 0; } 
      }
    }
    if (contflag) { continue; } else {break; }

  }
  
  return 0;
}
