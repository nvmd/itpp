#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

int main(int argc, char **argv)
{
	//receiver parameters
	string demapper_method[] = {"Hassibi_maxlogMAP", "GA", "sGA", "mmsePIC", "zfPIC"};
	string code_name[] = {"V-BLAST_MxN", "imp_V-BLAST_MxN", "Golden_2x2", "Damen_2x2"};
	int const_size = 4;
	int rec_antennas = 2;
	int perm_len = pow2i(10);
	int em_antennas = 2;
	int channel_uses = 2;
	double sigma2 = 1e-3;//the variance on each dimension

	//QAM modulator class
	QAM mod(const_size);

	//Space-Time code parameters
	STC st_block_code(code_name[0], const_size, em_antennas, channel_uses);
	em_antennas = st_block_code.get_nb_emission_antenna();
	channel_uses = st_block_code.get_channel_uses();
	int nb_symb = perm_len/mod.bits_per_symbol();
	int nb_subblocks = nb_symb/st_block_code.get_nb_symbols_per_block();
	int tx_duration = channel_uses*nb_subblocks;

	cout << "Number of symbols " << nb_symb << endl;
	cout << "Number of subblocks " << nb_subblocks << endl;
	cout << "Channel uses [symbols] " << channel_uses << endl;
	cout << "Tx duration [symbol durations] " << tx_duration << endl;

	//MIMO channel
	if (em_antennas != rec_antennas)
	{
		cout << "The number of emission antenna must be equal to the number of "
				"reception antenna" << endl;
		return EXIT_FAILURE;
	}
	cmat ch_matrix = eye_c(em_antennas);
	cmat ch_attenuations(em_antennas*rec_antennas, nb_subblocks);
	for (int ns=0; ns < nb_subblocks; ++ns)
	{
		ch_attenuations.set_col(ns, cvectorize(ch_matrix));
	}

	//SISO block
	SISO siso;
	siso.set_constellation(mod.bits_per_symbol(), mod.get_symbols(),
			mod.get_bits2symbols());
	siso.set_noise(sigma2);

	//bits generation
	bvec bits = randb(perm_len);

	//modulation
	cvec em = mod.modulate_bits(bits)/sqrt(em_antennas);

	cmat rec(tx_duration,rec_antennas);
	vec demapper_apriori_data(perm_len);
	demapper_apriori_data.zeros();
	vec demapper_extrinsic_data(perm_len);
	BERC berc;
	for (unsigned int c = 0; c < sizeof(code_name)/sizeof(code_name[0]); ++c)
	{
		//ST code
		st_block_code.setup(code_name[c], const_size, em_antennas, channel_uses);
		cmat S = st_block_code.encode(em);
		cout << "\n" << code_name[c] << endl;

		//MIMO channel
		for (int ns=0; ns < nb_subblocks; ++ns)
		{
			rec.set_submatrix(ns*channel_uses, 0,
					S(ns*channel_uses, (ns+1)*channel_uses-1, 0,
							em_antennas-1)*ch_matrix);
		}
		rec += sqrt(2*sigma2)*randn_c(tx_duration, rec_antennas);

		//SISO demapper
		siso.set_impulse_response(ch_attenuations);
		siso.set_st_block_code(st_block_code.get_nb_symbols_per_block(),
				st_block_code.get_1st_gen_matrix(),
				st_block_code.get_2nd_gen_matrix(),
				rec_antennas);
		for (unsigned int d = 0; d < sizeof(demapper_method)/sizeof(demapper_method[0]); ++d)
		{
			siso.set_demapper_method(demapper_method[d]);
			siso.demapper(demapper_extrinsic_data, rec, demapper_apriori_data);
			//BER
			berc.clear();
			berc.count(bits, demapper_extrinsic_data>0);
			cout << demapper_method[d] <<", BER = " << berc.get_errorrate() << endl;
		}
	}

	return EXIT_SUCCESS;
}
