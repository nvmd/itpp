/*!
 * \file
 * \brief Definitions for Space Time Codes (STC) class
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2011  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef STC_H
#define STC_H

#include <itpp/itbase.h> //IT++ base module

namespace itpp
{

/*!
  \ingroup misccommfunc
  \brief Space Time block Codes (STC) class

  Implements space time block codes using Hassibi's model. The following codes
  are available:
  - V-BLAST_MxN
  - imp_V-BLAST_MxN
  - Alamouti_2xN
  - Switched_Alamouti_4xN
  - Double_Alamouti_4xN
  - Jafarkhani_4xN
  - Golden_2x2
  - Damen_2x2
  - 34ortho_3xN
  - 36LD_3xN
  - 37LD_3xN
  - 39LD_3xN

  The code name and the constellation size are always needed to generate the
  requested code. The number of emission antenna and the channel uses are required
  by some codes to be provided by the user as input parameters, for other codes
  they have predefined values and don't need to be specified by the user.
  The number of symbols per block is always set internally. Therefore, it is
  recommended that after setting up the generator matrices (either through the
  constructor or through the setup() method) to call getters in order to obtain
  the number of emission antenna, the channel uses and the number of symbols per
  block.

  Usage example:
  \code
  STC stc(code_name, const_size);
  nb_em_antenna = stc.get_nb_emission_antenna();
  channel_uses = stc.get_channel_uses();
  symb_block = stc.get_nb_symbols_per_block();
  //symbol generation
  enc_symb = stc.encode(symb);
  \endcode

  Reference: B. Hassibi and B. M. Hochwald, ''High-rate codes that are linear in space and time,``
  IEEE Transactions on Information Theory, vol. 48, pp. 1804-1824, July 2002
 */
class STC
{
public:
	//! Space Time Code constructor (sets up the generator matrices using Hassibi's method)
	inline STC(const std::string &in_code_name, //!< code name (see available codes)
    		int in_const_size, //!< constellation size (should be at least two)
    		int in_em_antenna = 0, //!< number of emission antenna (for some codes it is set internally and should be obtained with get_nb_emission_antenna())
    		int in_channel_uses = 0 //!< number of channel uses (for some codes it is set internally and should be obtained with get_channel_uses())
    		)
	{
		setup(in_code_name, in_const_size, in_em_antenna, in_channel_uses);
	}
    //! Sets up the generator matrices using Hassibi's method (can be used to obtain new generator matrices, e.g. for a different code)
    inline void setup(const std::string &in_code_name, //!< code name (see available codes)
    		int in_const_size, //!< constellation size (should be at least two)
    		int in_em_antenna = 0, //!< number of emission antenna (for some codes it is set internally and should be obtained with get_nb_emission_antenna())
    		int in_channel_uses = 0 //!< number of channel uses (for some codes it is set internally and should be obtained with get_channel_uses())
    		)
    {
    	code_name = in_code_name;
    	it_assert(in_const_size >= 2, "Constellation size should be at least two");
    	const_size = in_const_size;
        em_antenna = in_em_antenna;
        channel_uses = in_channel_uses;
        Hassibi_block_code();
    }
    //! Gets the number of emission antenna (for some codes this is a predefined parameter)
    inline int get_nb_emission_antenna(void) const
    {
    	return em_antenna;
    }
    //! Gets the channel uses  (for some codes this is a predefined parameter)
    inline int get_channel_uses(void) const
    {
    	return channel_uses;
    }
    //! Gets the number of symbols per block (for all codes this is an output parameter)
    inline int get_nb_symbols_per_block(void) const
    {
    	return symb_block;
    }
    //! Gets the first generator matrix of the ST code following Hassibi's approach
    inline itpp::cmat get_1st_gen_matrix(void) const
    {
        return A;
    }
    //! Gets the second generator matrix of the ST code following Hassibi's approach
    inline itpp::cmat get_2nd_gen_matrix(void) const
    {
        return B;
    }
    //! Encodes input symbols according to the specified ST code
    itpp::cmat encode(const itpp::cvec &symb);
private:
    STC(const STC&);//not used
    STC& operator=(const STC&);//not used
    void Hassibi_block_code(void);
    itpp::cmat diag_pow(const itpp::cmat &in_mat, double in_exp);
    itpp::mat mat_pow(const itpp::mat &in_mat, int in_exp);
    std::string code_name;
    int const_size;
    int em_antenna;
    int channel_uses;
    int symb_block;
    itpp::cmat A;
    itpp::cmat B;
};

}
#endif /* STC_H_ */
