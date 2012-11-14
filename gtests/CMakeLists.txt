#
# \file CMakeLists.txt
# \brief cmake configuration file for gtests
# \author Bogdan Cristea
#
# -------------------------------------------------------------------------
#
# Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
#
# This file is part of IT++ - a C++ library of mathematical, signal
# processing, speech processing, and communications classes and functions.
#
# IT++ is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along
# with IT++.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------

file ( GLOB ITPP_TEST_SRCS *.cpp )

include_directories ( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${GTEST_DIR} ${GTEST_DIR}/include )

add_executable ( itpp_gtests ${ITPP_TEST_SRCS} ${GTEST_DIR}/src/gtest-all.cc ${GTEST_DIR}/src/gtest_main.cc )

target_link_libraries ( itpp_gtests ${ITPP_LIBS} ${libitpp_target} pthread )