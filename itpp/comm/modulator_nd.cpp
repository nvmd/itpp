/*!
 * \file
 * \brief Implementation of vector (MIMO) modulator classes
 * \author Erik G. Larsson and Adam Piatyszek
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

#include <itpp/comm/modulator_nd.h>
#include <itpp/comm/commfunc.h>
#include <itpp/base/algebra/cholesky.h>
#include <itpp/base/algebra/inv.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/converters.h>
#include <itpp/base/itcompat.h>

namespace itpp
{

// ----------------------------------------------------------------------
// Modulator_ND
// ----------------------------------------------------------------------

QLLRvec Modulator_ND::probabilities(QLLR l)
{
  QLLRvec result(2);

  if (l < 0) {        // this can be done more efficiently
    result(1) = -llrcalc.jaclog(0, -l);
    result(0) = result(1) - l;
  }
  else {
    result(0) = -llrcalc.jaclog(0, l);
    result(1) = result(0) + l;
  }
  return result;
}

Array<QLLRvec> Modulator_ND::probabilities(const QLLRvec &l)
{
  Array<QLLRvec> result(length(l));
  for (int i = 0; i < length(l); i++) {
    result(i) = probabilities(l(i));
  }
  return result;
}

void Modulator_ND::update_LLR(const Array<QLLRvec> &logP_apriori, int s,
                              QLLR scaled_norm, int j, QLLRvec &p1,
                              QLLRvec &p0)
{
  QLLR log_apriori_prob_const_point = 0;
  int b = 0;
  for (int i = 0; i < k(j); i++) {
    log_apriori_prob_const_point +=
      ((bitmap(j)(s, i) == 0) ? logP_apriori(b)(1) : logP_apriori(b)(0));
    b++;
  }

  b = 0;
  for (int i = 0; i < k(j); i++) {
    if (bitmap(j)(s, i) == 0) {
      p1(b) = llrcalc.jaclog(p1(b), scaled_norm
                             + log_apriori_prob_const_point);
    }
    else {
      p0(b) = llrcalc.jaclog(p0(b), scaled_norm
                             + log_apriori_prob_const_point);
    }
    b++;
  }
}

void Modulator_ND::update_LLR(const Array<QLLRvec> &logP_apriori,
                              const ivec &s, QLLR scaled_norm,
                              QLLRvec &p1, QLLRvec &p0)
{
  QLLR log_apriori_prob_const_point = 0;
  int b = 0;
  for (int j = 0; j < nt; j++) {
    for (int i = 0; i < k(j); i++) {
      log_apriori_prob_const_point +=
        ((bitmap(j)(s[j], i) == 0) ? logP_apriori(b)(1) : logP_apriori(b)(0));
      b++;
    }
  }

  b = 0;
  for (int j = 0; j < nt; j++) {
    for (int i = 0; i < k(j); i++) {
      if (bitmap(j)(s[j], i) == 0) {
        p1(b) = llrcalc.jaclog(p1(b), scaled_norm
                               + log_apriori_prob_const_point);
      }
      else {
        p0(b) = llrcalc.jaclog(p0(b), scaled_norm
                               + log_apriori_prob_const_point);
      }
      b++;
    }
  }
}

// ----------------------------------------------------------------------
// Modulator_NRD
// ----------------------------------------------------------------------

void Modulator_NRD::modulate_bits(const bvec &bits, vec &out_symbols) const
{
  it_assert(length(bits) == sum(k), "Modulator_NRD::modulate_bits(): "
            "The number of input bits does not match.");

  out_symbols.set_size(nt);

  int b = 0;
  for (int i = 0; i < nt; ++i) {
    int symb = bin2dec(bits.mid(b, k(i)));
    out_symbols(i) = symbols(i)(bits2symbols(i)(symb));
    b += k(i);
  }
}

vec Modulator_NRD::modulate_bits(const bvec &bits) const
{
  vec result(nt);
  modulate_bits(bits, result);
  return result;
}


void Modulator_NRD::demodulate_soft_bits(const vec &y, const mat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori,
    Soft_Demod_Method method)
{
  switch (method) {
  case FULL_ENUM_LOGMAP:
    demodulate_soft_bits(y, H, sigma2, LLR_apriori, LLR_aposteriori);
    break;
  case ZF_LOGMAP: {
    it_assert(H.rows() >= H.cols(), "Modulator_NRD::demodulate_soft_bits():"
              " ZF demodulation impossible for undetermined systems");
    // Set up the ZF detector
    mat Ht = H.T();
    mat inv_HtH = inv(Ht * H);
    vec shat = inv_HtH * Ht * y;
    vec h = ones(shat.size());
    for (int i = 0; i < shat.size(); ++i) {
      // noise covariance of shat
      double sigma_zf = std::sqrt(inv_HtH(i, i) * sigma2);
      shat(i) /= sigma_zf;
      h(i) /= sigma_zf;
    }
    demodulate_soft_bits(shat, h, 1.0, zeros_i(sum(k)), LLR_aposteriori);
  }
  break;
  default:
    it_error("Modulator_NRD::demodulate_soft_bits(): Improper soft "
             "demodulation method");
  }
}

QLLRvec Modulator_NRD::demodulate_soft_bits(const vec &y, const mat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    Soft_Demod_Method method)
{
  QLLRvec result;
  demodulate_soft_bits(y, H, sigma2, LLR_apriori, result, method);
  return result;
}

void Modulator_NRD::demodulate_soft_bits(const vec &y, const vec &h,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori)
{
  it_assert(length(LLR_apriori) == sum(k),
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");
  it_assert((length(h) == length(y)) && (length(h) == nt),
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");

  // set size of the output vector
  LLR_aposteriori.set_size(LLR_apriori.size());

  // normalisation constant "minus one over two sigma^2"
  double moo2s2 = -1.0 / (2.0 * sigma2);

  int b = 0;
  for (int i = 0; i < nt; ++i) {
    QLLRvec bnum = -QLLR_MAX * ones_i(k(i));
    QLLRvec bdenom = bnum;
    Array<QLLRvec> logP_apriori = probabilities(LLR_apriori(b, b + k(i) - 1));
    for (int j = 0; j < M(i); ++j) {
      double norm2 = moo2s2 * sqr(y(i) - h(i) * symbols(i)(j));
      QLLR scaled_norm = llrcalc.to_qllr(norm2);
      update_LLR(logP_apriori, j, scaled_norm, i, bnum, bdenom);
    }
    LLR_aposteriori.set_subvector(b, bnum - bdenom);
    b += k(i);
  }
}

void Modulator_NRD::demodulate_soft_bits(const vec &y, const mat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori)
{
  int np = sum(k); // number of bits in total
  int nr = H.rows();
  it_assert(length(LLR_apriori) == np,
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");
  it_assert((H.rows() == length(y)) && (H.cols() == nt),
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");

  // set size of the output vector
  LLR_aposteriori.set_size(LLR_apriori.size());

  // normalisation constant "minus one over two sigma^2"
  double moo2s2 = -1.0 / (2.0 * sigma2);

  bool diff_update = false;
  for (int i = 0; i < length(M); ++i) {
    // differential update only pays off for larger dimensions
    if (nt * M(i) > 4) {
      diff_update = true;
      break;
    }
  }

  Array<QLLRvec> logP_apriori = probabilities(LLR_apriori);

  mat Ht = H.T();
  mat HtH = Ht * H;
  vec ytH = Ht * y;

  QLLRvec bnum = -QLLR_MAX * ones_i(np);
  QLLRvec bdenom = bnum;
  ivec s = zeros_i(nt);
  double norm = 0.0;

  // Go over all constellation points  (r=dimension, s=vector of symbols)
  int r = nt - 1;
  while (true) {

    if (diff_update) {
      update_norm(norm, r, s(r), 0, ytH, HtH, s);
    }
    s(r) = 0;

    while (true) {
      if (s(r) > M(r) - 1) {
        if (r == nt - 1) {
          goto exit_point;
        }
        r++;
      }
      else {
        if (r == 0) {
          if (!diff_update) {
            norm = 0.0;
            for (int p = 0; p < nr; ++p) {
              double d = y(p);
              for (int i = 0; i < nt; ++i) {
                d -= H(p, i) * symbols(i)[s[i]];
              }
              norm += sqr(d);
            }
          }
          QLLR scaled_norm = llrcalc.to_qllr(norm * moo2s2);
          update_LLR(logP_apriori, s, scaled_norm, bnum, bdenom);
        }
        else {
          r--;
          break;
        }
      }
      if (diff_update) {
        update_norm(norm, r, s(r), s(r) + 1, ytH, HtH, s);
      }
      s(r)++;
    }
  }

exit_point:
  LLR_aposteriori = bnum - bdenom;

}


void Modulator_NRD::update_norm(double &norm, int k, int sold, int snew,
                                const vec &ytH, const mat &HtH,
                                const ivec &s)
{

  int m = length(s);
  double cdiff = symbols(k)[snew] - symbols(k)[sold];

  norm += sqr(cdiff) * HtH(k, k);
  cdiff *= 2.0;
  norm -= cdiff * ytH[k];
  for (int i = 0; i < m; ++i) {
    norm += cdiff * HtH(i, k) * symbols(k)[s[i]];
  }
}


std::ostream &operator<<(std::ostream &os, const Modulator_NRD &mod)
{
  os << "--- REAL MIMO (NRD) CHANNEL ---------" << std::endl;
  os << "Dimension (nt):           " << mod.nt << std::endl;
  os << "Bits per dimension (k):   " << mod.k << std::endl;
  os << "Symbols per dimension (M):" << mod.M << std::endl;
  for (int i = 0; i < mod.nt; i++) {
    os << "Bitmap for dimension " << i << ": " << mod.bitmap(i) << std::endl;
    // skip printing the trailing zero
    os << "Symbol coordinates for dimension " << i << ": " << mod.symbols(i).left(mod.M(i)) << std::endl;
  }
  os << mod.get_llrcalc() << std::endl;
  return os;
}

// ----------------------------------------------------------------------
// Modulator_NCD
// ----------------------------------------------------------------------

void Modulator_NCD::modulate_bits(const bvec &bits, cvec &out_symbols) const
{
  it_assert(length(bits) == sum(k), "Modulator_NCD::modulate_bits(): "
            "The number of input bits does not match.");

  out_symbols.set_size(nt);

  int b = 0;
  for (int i = 0; i < nt; ++i) {
    int symb = bin2dec(bits.mid(b, k(i)));
    out_symbols(i) = symbols(i)(bits2symbols(i)(symb));
    b += k(i);
  }
}

cvec Modulator_NCD::modulate_bits(const bvec &bits) const
{
  cvec result(nt);
  modulate_bits(bits, result);
  return result;
}


void Modulator_NCD::demodulate_soft_bits(const cvec &y, const cmat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori,
    Soft_Demod_Method method)
{
  switch (method) {
  case FULL_ENUM_LOGMAP:
    demodulate_soft_bits(y, H, sigma2, LLR_apriori, LLR_aposteriori);
    break;
  case ZF_LOGMAP: {
    it_assert(H.rows() >= H.cols(), "Modulator_NCD::demodulate_soft_bits():"
              " ZF demodulation impossible for undetermined systems");
    // Set up the ZF detector
    cmat Hht = H.H();
    cmat inv_HhtH = inv(Hht * H);
    cvec shat = inv_HhtH * Hht * y;
    cvec h = ones_c(shat.size());
    for (int i = 0; i < shat.size(); ++i) {
      double sigma_zf = std::sqrt(real(inv_HhtH(i, i)) * sigma2);
      shat(i) /= sigma_zf;
      h(i) /= sigma_zf;
    }
    demodulate_soft_bits(shat, h, 1.0, zeros_i(sum(k)), LLR_aposteriori);
  }
  break;
  default:
    it_error("Modulator_NCD::demodulate_soft_bits(): Improper soft "
             "demodulation method");
  }
}

QLLRvec Modulator_NCD::demodulate_soft_bits(const cvec &y, const cmat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    Soft_Demod_Method method)
{
  QLLRvec result;
  demodulate_soft_bits(y, H, sigma2, LLR_apriori, result, method);
  return result;
}


void Modulator_NCD::demodulate_soft_bits(const cvec &y, const cvec &h,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori)
{
  it_assert(length(LLR_apriori) == sum(k),
            "Modulator_NCD::demodulate_soft_bits(): Wrong sizes");
  it_assert((length(h) == length(y)) && (length(h) == nt),
            "Modulator_NCD::demodulate_soft_bits(): Wrong sizes");

  // set size of the output vector
  LLR_aposteriori.set_size(LLR_apriori.size());

  // normalisation constant "minus one over sigma^2"
  double moos2 = -1.0 / sigma2;

  int b = 0;
  for (int i = 0; i < nt; ++i) {
    QLLRvec bnum = -QLLR_MAX * ones_i(k(i));
    QLLRvec bdenom = -QLLR_MAX * ones_i(k(i));
    Array<QLLRvec> logP_apriori = probabilities(LLR_apriori(b, b + k(i) - 1));
    for (int j = 0; j < M(i); ++j) {
      double norm2 = moos2 * sqr(y(i) - h(i) * symbols(i)(j));
      QLLR scaled_norm = llrcalc.to_qllr(norm2);
      update_LLR(logP_apriori, j, scaled_norm, i, bnum, bdenom);
    }
    LLR_aposteriori.set_subvector(b, bnum - bdenom);
    b += k(i);
  }
}

void Modulator_NCD::demodulate_soft_bits(const cvec &y, const cmat &H,
    double sigma2,
    const QLLRvec &LLR_apriori,
    QLLRvec &LLR_aposteriori)
{
  int np = sum(k); // number of bits in total
  int nr = H.rows();
  it_assert(length(LLR_apriori) == np,
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");
  it_assert((H.rows() == length(y)) && (H.cols() == nt),
            "Modulator_NRD::demodulate_soft_bits(): Wrong sizes");

  // set size of the output vector
  LLR_aposteriori.set_size(LLR_apriori.size());

  // normalisation constant "minus one over sigma^2"
  double moos2 = -1.0 / sigma2;

  bool diff_update = false;
  for (int i = 0; i < length(M); ++i) {
    // differential update only pays off for larger dimensions
    if (nt * M(i) > 4) {
      diff_update = true;
    }
  }

  Array<QLLRvec> logP_apriori = probabilities(LLR_apriori);

  cmat HtH = H.hermitian_transpose() * H;
  cvec ytH = conj(H.hermitian_transpose() * y);

  QLLRvec bnum = -QLLR_MAX * ones_i(np);
  QLLRvec bdenom = -QLLR_MAX * ones_i(np);
  ivec s(nt);
  s.clear();
  double norm = 0.0;

  // Go over all constellation points  (r=dimension, s=vector of symbols)
  int r = nt - 1;
  while (true) {

    if (diff_update) {
      update_norm(norm, r, s(r), 0, ytH, HtH, s);
    }
    s(r) = 0;

    while (true) {
      if (s(r) > M(r) - 1)  {
        if (r == nt - 1) {
          goto exit_point;
        }
        r++;
      }
      else {
        if (r == 0) {
          if (!diff_update) {
            norm = 0.0;
            for (int p = 0; p < nr; ++p) {
              std::complex<double> d = y(p);
              for (int i = 0; i < nt; ++i) {
                d -= H(p, i) * symbols(i)[s[i]];
              }
              norm += sqr(d);
            }
          }
          QLLR scaled_norm = llrcalc.to_qllr(norm * moos2);
          update_LLR(logP_apriori, s, scaled_norm, bnum, bdenom);
        }
        else {
          r--;
          break;
        }
      }
      if (diff_update) {
        update_norm(norm, r, s(r), s(r) + 1, ytH, HtH, s);
      }
      s(r)++;
    }
  }

exit_point:
  LLR_aposteriori = bnum - bdenom;
}


void Modulator_NCD::update_norm(double &norm, int k, int sold, int snew,
                                const cvec &ytH, const cmat &HtH,
                                const ivec &s)
{
  int m = length(s);
  std::complex<double> cdiff = symbols(k)[snew] - symbols(k)[sold];

  norm += sqr(cdiff) * (HtH(k, k).real());
  cdiff *= 2.0;
  norm -= (cdiff.real() * ytH[k].real() - cdiff.imag() * ytH[k].imag());
  for (int i = 0; i < m; i++) {
    norm += (cdiff * HtH(i, k) * conj(symbols(k)[s[i]])).real();
  }
}


std::ostream &operator<<(std::ostream &os, const Modulator_NCD &mod)
{
  os << "--- COMPLEX MIMO (NCD) CHANNEL --------" << std::endl;
  os << "Dimension (nt):           " << mod.nt << std::endl;
  os << "Bits per dimension (k):   " << mod.k << std::endl;
  os << "Symbols per dimension (M):" << mod.M << std::endl;
  for (int i = 0; i < mod.nt; i++) {
    os << "Bitmap for dimension " << i << ": "
    << mod.bitmap(i) << std::endl;
    os << "Symbol coordinates for dimension " << i << ": "
    << mod.symbols(i).left(mod.M(i)) << std::endl;
  }
  os << mod.get_llrcalc() << std::endl;
  return os;
}

// ----------------------------------------------------------------------
// ND_UPAM
// ----------------------------------------------------------------------

ND_UPAM::ND_UPAM(int nt, int Mary)
{
  set_M(nt, Mary);
}

void ND_UPAM::set_M(int nt_in, int Mary)
{
  nt = nt_in;
  ivec Mary_temp(nt);
  Mary_temp = Mary;
  set_M(nt, Mary_temp);
}

void ND_UPAM::set_M(int nt_in, ivec Mary)
{
  nt = nt_in;
  it_assert(length(Mary) == nt, "ND_UPAM::set_M(): Mary has wrong length");
  k.set_size(nt);
  M = Mary;
  bitmap.set_size(nt);
  symbols.set_size(nt);
  bits2symbols.set_size(nt);
  spacing.set_size(nt);

  for (int i = 0; i < nt; i++) {
    k(i) = round_i(::log2(static_cast<double>(M(i))));
    it_assert((k(i) > 0) && ((1 << k(i)) == M(i)),
              "ND_UPAM::set_M(): M is not a power of 2.");

    symbols(i).set_size(M(i) + 1);
    bits2symbols(i).set_size(M(i));
    bitmap(i) = graycode(k(i));
    double average_energy = (M(i) * M(i) - 1) / 3.0;
    double scaling_factor = std::sqrt(average_energy);

    for (int j = 0; j < M(i); ++j) {
      symbols(i)(j) = ((M(i) - 1) - j * 2) / scaling_factor;
      bits2symbols(i)(bin2dec(bitmap(i).get_row(j))) = j;
    }

    // the "symbols" vector must end with a zero; only for a trick
    // exploited in update_norm()
    symbols(i)(M(i)) = 0.0;

    spacing(i) = 2.0 / scaling_factor;
  }
}

int ND_UPAM::sphere_search_SE(const vec &y_in, const mat &H,
                              const imat &zrange, double r, ivec &zhat)
{
  // The implementation of this function basically follows the
  // Schnorr-Eucner algorithm described in Agrell et al. (IEEE
  // Trans.  IT, 2002), but taking into account constellation
  // boundaries, see the "accelerated sphere decoder" in Boutros et
  // al. (IEEE Globecom, 2003).  No lattice reduction is performed.
  // Potentially the function can be speeded up by performing
  // lattice reduction, but it seems difficult to keep track of
  // constellation boundaries.

  mat R = chol(H.transpose() * H);
  mat Ri = inv(R);
  mat Q = H * Ri;
  vec y = Q.transpose() * y_in;
  mat Vi = Ri.transpose();

  int n = H.cols();
  vec dist(n);
  dist(n - 1) = 0;
  double bestdist = r * r;
  int status = -1; // search failed

  mat E = zeros(n, n);
  for (int i = 0; i < n; i++) {   // E(k,:) = y*Vi;
    for (int j = 0; j < n; j++) {
      E(i*n + n - 1) += y(j) * Vi(j + n * i);
    }
  }

  ivec z(n);
  zhat.set_size(n);
  z(n - 1) = floor_i(0.5 + E(n * n - 1));
  z(n - 1) = std::max(z(n - 1), zrange(n - 1, 0));
  z(n - 1) = std::min(z(n - 1), zrange(n - 1, 1));
  double p = (E(n * n - 1) - z(n - 1)) / Vi(n * n - 1);
  ivec step(n);
  step(n - 1) = sign_nozero_i(p);

  // Run search loop
  int k = n - 1;  // k uses natural indexing, goes from 0 to n-1

  while (true) {
    double newdist = dist(k) + p * p;

    if ((newdist < bestdist) && (k != 0)) {
      for (int i = 0; i < k; i++) {
        E(k - 1 + i*n) = E(k + i * n) - p * Vi(k + i * n);
      }

      k--;
      dist(k) = newdist;
      z(k) = floor_i(0.5 + E(k + k * n));
      z(k) = std::max(z(k), zrange(k, 0));
      z(k) = std::min(z(k), zrange(k, 1));
      p = (E(k + k * n) - z(k)) / Vi(k + k * n);

      step(k) = sign_nozero_i(p);
    }
    else {
      while (true) {
        if (newdist < bestdist) {
          zhat = z;
          bestdist = newdist;
          status = 0;
        }
        else if (k == n - 1) {
          goto exit_point;
        }
        else {
          k++;
        }

        z(k) += step(k);

        if ((z(k) < zrange(k, 0)) || (z(k) > zrange(k, 1))) {
          step(k) = (-step(k) - sign_nozero_i(step(k)));
          z(k) += step(k);
        }

        if ((z(k) >= zrange(k, 0)) && (z(k) <= zrange(k, 1))) {
          break;
        }
      }

      p = (E(k + k * n) - z(k)) / Vi(k + k * n);
      step(k) = (-step(k) - sign_nozero_i(step(k)));
    }
  }

exit_point:
  return status;
}


int ND_UPAM::sphere_decoding(const vec &y, const mat &H, double rstart,
                             double rmax, double stepup,
                             QLLRvec &detected_bits)
{
  it_assert(H.rows() == length(y),
            "ND_UPAM::sphere_decoding(): dimension mismatch");
  it_assert(H.cols() == nt,
            "ND_UPAM::sphere_decoding(): dimension mismatch");
  it_assert(rstart > 0, "ND_UPAM::sphere_decoding(): radius error");
  it_assert(rmax > rstart, "ND_UPAM::sphere_decoding(): radius error");

  // This function can be improved, e.g., by using an ordered search.

  vec ytemp = y;
  mat Htemp(H.rows(), H.cols());
  for (int i = 0; i < H.cols(); i++) {
    Htemp.set_col(i, H.get_col(i)*spacing(i));
    ytemp += Htemp.get_col(i) * 0.5 * (M(i) - 1.0);
  }

  imat crange(nt, 2);
  for (int i = 0; i < nt; i++) {
    crange(i, 0) = 0;
    crange(i, 1) = M(i) - 1;
  }

  int status = 0;
  double r = rstart;
  ivec s(sum(M));
  while (r <= rmax) {
    status = sphere_search_SE(ytemp, Htemp, crange, r, s);

    if (status == 0) { // search successful
      detected_bits.set_size(sum(k));
      int b = 0;
      for (int j = 0; j < nt; j++) {
        for (int i = 0; i < k(j); i++) {
          if (bitmap(j)((M(j) - 1 - s[j]), i) == 0) {
            detected_bits(b) = 1000;
          }
          else {
            detected_bits(b) = -1000;
          }
          b++;
        }
      }

      return status;
    }
    r = r * stepup;
  }

  return status;
}

// ----------------------------------------------------------------------
// ND_UQAM
// ----------------------------------------------------------------------

// The ND_UQAM (MIMO with uniform QAM) class could alternatively
// have been implemented by using a ND_UPAM class of twice the
// dimension, but this does not fit as elegantly into the class
// structure

ND_UQAM::ND_UQAM(int nt, int Mary)
{
  set_M(nt, Mary);
}

void ND_UQAM::set_M(int nt_in, int Mary)
{
  nt = nt_in;
  ivec Mary_temp(nt);
  Mary_temp = Mary;
  set_M(nt, Mary_temp);
}

void ND_UQAM::set_M(int nt_in, ivec Mary)
{
  nt = nt_in;
  it_assert(length(Mary) == nt, "ND_UQAM::set_M(): Mary has wrong length");
  k.set_size(nt);
  M = Mary;
  L.set_size(nt);
  bitmap.set_size(nt);
  symbols.set_size(nt);
  bits2symbols.set_size(nt);

  for (int i = 0; i < nt; ++i) {
    k(i) = round_i(::log2(static_cast<double>(M(i))));
    it_assert((k(i) > 0) && ((1 << k(i)) == M(i)),
              "ND_UQAM::set_M(): M is not a power of 2");

    L(i) = round_i(std::sqrt(static_cast<double>(M(i))));
    it_assert(L(i)*L(i) == M(i), "ND_UQAM: constellation M must be square");

    symbols(i).set_size(M(i) + 1);
    bitmap(i).set_size(M(i), k(i));
    bits2symbols(i).set_size(M(i));
    double average_energy = (M(i) - 1) * 2.0 / 3.0;
    double scaling_factor = std::sqrt(average_energy);
    bmat gray_code = graycode(levels2bits(L(i)));

    for (int j1 = 0; j1 < L(i); ++j1) {
      for (int j2 = 0; j2 < L(i); ++j2) {
        symbols(i)(j1*L(i) + j2) =
          std::complex<double>(((L(i) - 1) - j2 * 2.0) / scaling_factor,
                               ((L(i) - 1) - j1 * 2.0) / scaling_factor);
        bitmap(i).set_row(j1*L(i) + j2, concat(gray_code.get_row(j1),
                                               gray_code.get_row(j2)));
        bits2symbols(i)(bin2dec(bitmap(i).get_row(j1*L(i) + j2)))
        = j1 * L(i) + j2;
      }
    }

    // must end with a zero; only for a trick exploited in
    // update_norm()
    symbols(i)(M(i)) = 0.0;
  }
}

// ----------------------------------------------------------------------
// ND_UPSK
// ----------------------------------------------------------------------

ND_UPSK::ND_UPSK(int nt, int Mary)
{
  set_M(nt, Mary);
}

void ND_UPSK::set_M(int nt_in, int Mary)
{
  nt = nt_in;
  ivec Mary_temp(nt);
  Mary_temp = Mary;
  set_M(nt, Mary_temp);
}

void ND_UPSK::set_M(int nt_in, ivec Mary)
{
  nt = nt_in;
  it_assert(length(Mary) == nt, "ND_UPSK::set_M() Mary has wrong length");
  k.set_size(nt);
  M = Mary;
  bitmap.set_size(nt);
  symbols.set_size(nt);
  bits2symbols.set_size(nt);

  for (int i = 0; i < nt; ++i) {
    k(i) = round_i(::log2(static_cast<double>(M(i))));
    it_assert((k(i) > 0) && ((1 << k(i)) == M(i)),
              "ND_UPSK::set_M(): M is not a power of 2");

    symbols(i).set_size(M(i) + 1);
    bits2symbols(i).set_size(M(i));
    bitmap(i) = graycode(k(i));

    double delta = 2.0 * pi / M(i);
    double epsilon = delta / 10000.0;

    for (int j = 0; j < M(i); ++j) {
      std::complex<double> symb
      = std::complex<double>(std::polar(1.0, delta * j));

      if (std::abs(std::real(symb)) < epsilon) {
        symbols(i)(j) = std::complex<double>(0.0, std::imag(symb));
      }
      else if (std::abs(std::imag(symb)) < epsilon) {
        symbols(i)(j) = std::complex<double>(std::real(symb), 0.0);
      }
      else {
        symbols(i)(j) = symb;
      }

      bits2symbols(i)(bin2dec(bitmap(i).get_row(j))) = j;
    }

    // must end with a zero; only for a trick exploited in
    // update_norm()
    symbols(i)(M(i)) = 0.0;
  }
}

} // namespace itpp
