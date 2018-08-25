//
// Newton's method library to solve simultaneous equations 2018-08-26.02
// https://github.com/trueroad/newton_method
//
// sample-fast.cc:
//   Sample for position calculation like GPS/GNSS receiver
//   (fast functions)
//
// Copyright (C) 2017, 2018 Masamichi Hosoda.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//

#include<iomanip>
#include<iostream>
#include<cmath>

#include"newton_method/newton.hh"

// Satellite position (known, given from received data)
std::vector<std::vector<double>> Si
{
  {   -11327938.990000,     9886884.330000,    21895433.227000 },
  {     4755496.711000,    19362623.328000,    18112665.323000 },
  {    -7506201.243000,    24076860.073000,     7092793.940000 },
  {   -23085789.286000,    12409399.010000,     4602891.246000 },
  {   -21893190.888000,    -2248546.668000,    14796664.928000 },
  {   -24893247.395000,     3827508.606000,    -8794926.751000 },
  {   -12971740.598000,   -10587013.898000,    21061849.442000 },
  {     7069732.127000,    22267387.067000,    12627670.276000 },
};

// Observed distance (known, measurement result)
std::vector<double> Ri
{
     20690632.972000,
     23225588.018000,
     21288081.687000,
     21187099.471000,
     21833271.739000,
     24393427.283000,
     24031767.538000,
     23630886.925000,
};

// Observed weight (known, measurement result)
std::vector<double> W
{
          0.34246575,
          0.41806020,
          0.44662796,
          0.33967391,
          0.33411293,
          0.29682398,
          0.30759766,
          0.31046259,
};

// Receiver position (actually unknown)
std::vector<double> unknown_P
{
   -3947719.26542369,   3364403.97164603,   3699487.31861822
};

// Error epsilon (0.1 mm)
double epsilon {0.0001};

inline double square (double x)
{
  return x * x;
}

inline double distance (const std::vector<double> &p1,
                        const std::vector<double> &p2)
{
  return sqrt ( square ( p1[0] - p2[0])
                + square ( p1[1] - p2[1])
                + square ( p1[2] - p2[2]) );
}

inline double distance (const double *p1, const double *p2)
{
  return sqrt ( square ( p1[0] - p2[0])
                + square ( p1[1] - p2[1])
                + square ( p1[2] - p2[2]) );
}

// Calculate cols of F and J
inline void calc_fj_cols (int i, double *fi, double *ji, const double *x)
{
  auto d {distance (Si[i].data (), x)};

  *fi = d + x[3] - Ri[i];

  for (int n = 0; n < 3; ++n)
    {
      ji[n] = - (Si[i][n] - x[n]) / d;
    }
  ji[3] = 1.0;
}

// Calculate F (to be solved as equations) and J (Jacobian matrix of F)
void calc_fj (double *f, double *j, const double *x)
{
  for (int i = 0; i < Si.size (); ++i)
    calc_fj_cols (i, f + i, j + (i * 4), x);
}

int main (void)
{
  std::cout
    << "Newton's method library to solve simultaneous equations" << std::endl
    << "https://github.com/trueroad/newton_method" << std::endl
    << std::endl
    << "Sample for position calculation like GPS/GNSS receiver" << std::endl
    << "(fast functions)" << std::endl
    << std::endl
    << "Copyright (C) 2018 Masamichi Hosoda. All rights reserved."
    << std::endl << std::endl;

  std::cout << std::fixed;

  std::cout << "Number of satellites: " << Si.size () << std::endl;
  std::cout << "Satellite position: Si (x y z)" << std::endl;
  for (auto i: Si)
    {
      std::cout << " (";
      for (auto j: i)
        std::cout << " " << j;
      std::cout << ")" << std::endl;
    }

  std::cout << "Observed distance: Ri" << std::endl;
  for (auto i: Ri)
    std::cout << " " << i << std::endl;

  std::cout << "Observed weight: W" << std::endl;
  for (auto i: W)
    std::cout << " " << i << std::endl;

  std::cout << std::defaultfloat;

  std::cout
    << "Trying newton's method" << std::endl
    << " for solving the receiver position P"
    " and the clock offset distance delta_S" << std::endl
    << " from the satellite position Si"
    " and observed distance Ri..."
    << std::endl;
  std::vector<double> solution;
  try
    {
      newton_method::newton_method nm;

      nm.set_function_fast (calc_fj);
      nm.set_epsilon_F (epsilon);
      nm.set_epsilon_deltaX (epsilon);
      nm.set_weight (W);
      std::vector<double> initial_value {0.0, 0.0, 0.0, 0.0};
      solution = nm.solve_fast<newton_method::least_square::weighted>
        (initial_value, Si.size ());
    }
  catch ( const std::exception &e )
    {
      std::cerr << e.what () << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << "unknown exception!!" << std::endl;
      throw;
    }

  std::cout << std::fixed;
  std::cout << std::endl << "Solution: P (x y z) and delta_S:" << std::endl;
  for (auto i: solution)
    std::cout << " " << i;
  std::cout << std::endl;

  std::cout << "Error distance:" << std::endl;
  double error_distance = distance (solution, unknown_P);
  std::cout << " " << error_distance << std::endl;

  return 0;
}
