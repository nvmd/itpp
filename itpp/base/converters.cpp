/*!
 * \file 
 * \brief Implementation of converters between different vector and matrix types
 * \author Tony Ottosson, Tobias Ringstrom and Pal Frenger
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/base/converters.h>
#include <itpp/base/matfunc.h>


namespace itpp {

  //--------------- Converters for vectors ------------------

  template <class T>
  bvec to_bvec(const Vec<T> &v)
  {
    bvec temp(v.length());
    for (int i=0;i<v.length();i++){
      //temp(i)=static_cast<bin>(v(i));
      temp(i)=bin(v(i));
    }
    return temp;
  }

  template <class T>
  svec to_svec(const Vec<T> &v)
  {
    svec temp(v.length());
    for (int i=0;i<v.length();i++){
      //temp(i)=static_cast<short>(v(i));
      temp(i)=short(v(i));
    }
    return temp;
  }

  template <class T>
  ivec to_ivec(const Vec<T> &v)
  {
    ivec temp(v.length());
    for (int i=0;i<v.length();i++){
      //temp(i)=static_cast<int>(v(i));
      temp(i)=int(v(i));
    }
    return temp;
  }

  template <class T>
  vec to_vec(const Vec<T> &v)
  {
    vec temp(v.length());
    for (int i=0;i<v.length();i++){
      //temp(i)=static_cast<double>(v(i));
      temp(i)=double(v(i));
    }
    return temp;
  }

  template <class T>
  cvec to_cvec(const Vec<T> &v)
  {
    cvec temp(v.length());
    for (int i=0;i<v.length();i++){
      //temp(i)=static_cast<complex<double> >(v(i));
      temp(i)=std::complex<double>(v(i));
    }
    return temp;
  }

  cvec to_cvec(const bvec &v)
  {
    cvec temp(v.length());
    for (int i=0;i<v.length();i++){
      if (v(i)==bin(0)) { 
	temp(i) = std::complex<double>(0.0,0.0);
      } else {
	temp(i) = std::complex<double>(1.0,0.0);
      }
    }
    return temp;
  }

  ivec to_ivec(int s) {ivec out(1); out(0) = s; return out;}

  vec to_vec(double s) {vec out(1); out(0) = s; return out;}

  cvec to_cvec(double real, double imag) {cvec out(1); out(0) = std::complex<double>(real,imag); return out;}

  //------------- end of converters for vectors -------------

  //--------------- Converters for matricies ------------------

  template <class T>
  bmat to_bmat(const Mat<T> &m)
  {
    bmat temp(m.rows(),m.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	//temp(i,j)=static_cast<bin>(m(i,j));
	temp(i,j)=bin(m(i,j));
      }
    }

    return temp;
  }

  template <class T>
  smat to_smat(const Mat<T> &m)
  {
    smat temp(m.rows(),m.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	//temp(i,j)=static_cast<short>(m(i,j));
	temp(i,j)=short(m(i,j));
      }
    }

    return temp;
  }

  template <class T>
  imat to_imat(const Mat<T> &m)
  {
    imat temp(m.rows(),m.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	//temp(i,j)=static_cast<int>(m(i,j));
	temp(i,j)=int(m(i,j));
      }
    }

    return temp;
  }

  template <class T>
  mat to_mat(const Mat<T> &m)
  {
    mat temp(m.rows(),m.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	//temp(i,j)=static_cast<double>(m(i,j));
	temp(i,j)=double(m(i,j));
      }
    }

    return temp;
  }

  template <class T>
  cmat to_cmat(const Mat<T> &m)
  {
    cmat temp(m.rows(),m.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	//temp(i,j)=static_cast<std::complex<double> >(m(i,j));
	temp(i,j)=std::complex<double>(m(i,j));
      }
    }

    return temp;
  }

  cmat to_cmat(const bmat &m)
  {
    cmat temp(m.rows(),m.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	if (m(i,j)==bin(0)) { 
	  temp(i,j) = std::complex<double>(0.0,0.0);
	} else {
	  temp(i,j) = std::complex<double>(1.0,0.0);
	}
      }
    }

    return temp;
  }

  //------------- end of converters for matricies -------------

  template <class T>
  cvec to_cvec(const Vec<T> &real, const Vec<T> &imag)
  {
    assert(real.length()==imag.length());
    cvec temp(real.length());
    for(int i=0;i<real.length();i++){
      temp(i)=std::complex<double>(real(i),imag(i));
    }
		
    return temp;
  }

  template <class T>
  cmat to_cmat(const Mat<T> &real, const Mat<T> &imag)
  {
    it_assert1(real.rows()==imag.rows() && real.cols()==imag.cols(), "to_cmat::real and imag part sizes does not match");
    cmat temp(real.rows(),real.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	temp(i,j)=std::complex<double>(real(i,j), imag(i,j));
      }
    }

    return temp;
  }

  bvec dec2bin(int length, int index)
  {
    int i, bintemp = index;
    bvec temp(length);

    for (i=length-1; i>=0; i--) {
      temp(i) = bin(bintemp & 1);
      bintemp = (bintemp >> 1);
    }
    return temp;
  }

  bvec dec2bin(int index, bool msb_first)
  {
    int length = needed_bits(index);
    int i, bintemp = index;
    bvec temp(length);

    for (i=length-1; i>=0; i--) {
      temp(i) = bin(bintemp & 1);
      bintemp = (bintemp >> 1);
    }
    if (msb_first){
      return temp;
    } else{
      return reverse(temp);
    }
  }

  void dec2bin(int index, bvec &v)
  {
    int i, bintemp = index;
    v.set_size(needed_bits(index), false);

    for (i=v.size()-1; i>=0; i--) {
      v(i) = bin(bintemp & 1);
      bintemp = (bintemp >> 1);
    }
  }

  int bin2dec(const bvec &inbvec, bool msb_first)
  {
    int i, temp=0;
    int sizebvec=inbvec.length();
    if (msb_first) {
      for (i=0; i<sizebvec; i++) {
	temp+=pow2i(sizebvec-i-1)*int(inbvec(i));
      }
    } else {
      for (i=0; i<sizebvec; i++) {
	temp+=pow2i(i)*int(inbvec(i));
      }
    }
    return temp;
  }

  bvec oct2bin(const ivec &octalindex, short keepzeros)
  {
    int length = octalindex.length(), i;
    bvec out(3*length);
    for (i=0; i<length; i++) {
      out.replace_mid(3*i,dec2bin(3,octalindex(i)));
    }
    //remove zeros if keepzeros = 0
    if (keepzeros == 0) {
      for (i=0; i<out.length(); i++) {
	if ( (short)out(i) != 0) { 
	  return out.right(out.length()-i);
	  break;
	}
      }
      return bvec("0");
    } else {
      return out;
    }
  }

  ivec bin2oct(const bvec &inbits)
  {
    int start, Itterations = (int)ceil(float(inbits.length())/3.0);
    ivec out(Itterations);
    for (int i=Itterations-1; i>0; i--) {
      start = 3*i - ( 3*Itterations - inbits.length() );
      out(i) = bin2dec(inbits.mid(start,3));
    }
    out(0) = bin2dec( inbits.left(inbits.length()-((Itterations-1)*3)) );
    return out;
  }

  ivec bin2pol(const bvec &inbvec)
  {
    return 1-2*to_ivec(inbvec);
  }

  bvec pol2bin(const ivec &inpol)
  {
    return to_bvec((1-inpol)/2);
  }

  /*
  template <>
  string to_str(const double &i)
  {
    std::ostringstream ss;
    ss.precision(8);
    ss.setf(ostringstream::scientific,ostringstream::floatfield);
    ss << i;
    return ss.str();
  }
*/  
  std::string to_str(const double &i, const int precision)
  {
    std::ostringstream ss;
    ss.precision(precision);
    ss.setf(std::ostringstream::scientific, std::ostringstream::floatfield);
    ss << i;
    return ss.str();
  }
  
  // ---------------------- Instantiations -----------------------------------------

  template bvec to_bvec(const svec &v);
  //#ifndef _MSC_VER
  template bvec to_bvec(const Vec<int> &v);
  //#endif
  template svec to_svec(const bvec &v);
  template svec to_svec(const ivec &v);
  template svec to_svec(const svec &v);
  template svec to_svec(const vec &v);

  //#ifndef _MSC_VER
  template ivec to_ivec(const bvec &v);
  //#endif
  template ivec to_ivec(const svec &v);
  template ivec to_ivec(const ivec &v);
  template ivec to_ivec(const vec &v);

  template vec to_vec(const bvec &v);
  template vec to_vec(const svec &v);
  template vec to_vec(const ivec &v);
  template vec to_vec(const vec &v);

  // Template instantiation of to_cvec
  //template cvec to_cvec(const bvec &v); //Specialization created above

  template cvec to_cvec(const svec &v);
  template cvec to_cvec(const ivec &v);
  template cvec to_cvec(const vec &v);
  template cvec to_cvec(const cvec &v);
  template cvec to_cvec(const bvec &real, const bvec &imag);
  template cvec to_cvec(const svec &real, const svec &imag);
  template cvec to_cvec(const ivec &real, const ivec &imag);
  template cvec to_cvec(const vec &real, const vec &imag);

  template bmat to_bmat(const smat &m);
  template bmat to_bmat(const imat &m);
  template smat to_smat(const bmat &m);
  template smat to_smat(const imat &m);

  template imat to_imat(const bmat &m);
  template imat to_imat(const smat &m);
  template imat to_imat(const imat &m);
  // Template instantiation of to_imat
  template imat to_imat(const mat &m);

  template mat to_mat(const bmat &m);
  template mat to_mat(const smat &m);
  template mat to_mat(const imat &m);
  template mat to_mat(const mat &m);

  // Template instantiation of to_cmat
  // template cmat to_cmat(const bmat &m); //Specialization created above

  template cmat to_cmat(const smat &m);
  template cmat to_cmat(const imat &m);
  template cmat to_cmat(const mat &m);
  template cmat to_cmat(const cmat &m);
  template cmat to_cmat(const bmat &real, const bmat &imag);
  template cmat to_cmat(const smat &real, const smat &imag);
  template cmat to_cmat(const imat &real, const imat &imag);
  template cmat to_cmat(const mat &real, const mat &imag);

} // namespace itpp
