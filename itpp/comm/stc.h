/*!
 * \file
 * \brief Definitions for Space Time Codes (STC) class
 * \author Bogdan Cristea
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

  Implements Space Time block Codes using Hassibi's model

  Reference: B. Hassibi and B. M. Hochwald, ''High-rate codes that are linear in space and time,``
  IEEE Transactions on Information Theory, vol. 48, pp. 1804-1824, July 2002
 */
class STC
{
public:
	//! Space Time Code constructor
	inline STC(int &in_out_em_antennas, int &in_out_channel_uses,
			int &out_symb_block,
    		const std::string &in_code_name, int in_const_size) :
    		em_antennas(in_out_em_antennas),
    		channel_uses(in_out_channel_uses),
    		symb_block(out_symb_block)
	{
		setup(in_out_em_antennas, in_out_channel_uses, out_symb_block,
				in_code_name, in_const_size);
	}
    //! Setup ST block codes (Hassibi's method is used)
    inline void setup(int &in_out_em_antennas, int &in_out_channel_uses,
    		int &out_symb_block,
    		const std::string &in_code_name, int in_const_size)
    {
        em_antennas = in_out_em_antennas;
        channel_uses = in_out_channel_uses;
        symb_block = out_symb_block;
        code_name = in_code_name;
        const_size = in_const_size;
        Hassibi_block_code();
    };
    //! Encodes input symbols according to specified ST code
    itpp::cmat encode(const itpp::cvec &symb)
    {
        return Hassibi_encode(symb);
    };
    //! Gets the first generator matrix of the ST code following Hassibi's approach
    const itpp::cmat& get_1st_gen_matrix(void) const
    {
        return A;
    };
    //! Gets the second generator matrix of the ST code following Hassibi's approach
    const itpp::cmat& get_2nd_gen_matrix(void) const
    {
        return B;
    };
private:
    STC(const STC&);
    STC& operator=(const STC&);
    void Hassibi_block_code(void);
    itpp::cmat Hassibi_encode(const itpp::cvec &symb);
    itpp::cmat diag_pow(const itpp::cmat &in_mat, double in_exp);
    itpp::mat mat_pow(const itpp::mat &in_mat, int in_exp);
    int const_size;
    itpp::cmat A;
    itpp::cmat B;
    int &em_antennas;
    int &channel_uses;
    int &symb_block;
    std::string code_name;
};

}
#endif /* STC_H_ */
