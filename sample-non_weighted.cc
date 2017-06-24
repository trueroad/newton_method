//
// Newton's method library to solve simultaneous equations 2017-06-24.16
// https://github.com/trueroad/newton_method/
//
// Sample: Position calculation like GPS/GNSS receiver
//         (overdetermined system, non weighted)
//
// Copyright (C) 2017 Masamichi Hosoda. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY AUTHOR AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
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

#include"newton.hh"

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

// Receiver position (actually unknown)
std::vector<double> unknown_P
{
   -3947719.36876915,   3364403.46661849,   3699487.64248845
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

// Calculate F (to be solved as equations)
std::vector<double> calc_f (const std::vector<double> &Xarg)
{
  std::vector<double> fv;

  for (int i = 0; i < Si.size (); ++i)
    fv.push_back (distance (Si[i], Xarg) + Xarg[3] - Ri[i]);

  return fv;
}

// Calculate J (Jacobian matrix of F)
std::vector<std::vector<double>> calc_j (const std::vector<double> &Xarg)
{
  std::vector<std::vector<double>> jv;

  for (int i = 0; i < Si.size (); ++i)
    {
      std::vector<double> j_row;
      double d = distance (Si[i], Xarg);

      for (int n = 0; n < 3; ++n)
        j_row.push_back (- (Si[i][n] - Xarg[n]) / d);
      j_row.push_back (1.0);

      jv.push_back (j_row);
    }

  return jv;
}

int main (void)
{
  std::cout
    << "Newton's method library to solve simultaneous equations" << std::endl
    << "https://github.com/trueroad/newton_method/" << std::endl
    << std::endl
    << "Sample: Position calculation like GPS/GNSS receiver" << std::endl
    << "        (overdetermined system, non weighted)" << std::endl
    << std::endl
    << "Copyright (C) 2017 Masamichi Hosoda. All rights reserved."
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

      nm.set_function (calc_f, calc_j);
      nm.set_epsilon_F (epsilon);
      nm.set_epsilon_deltaX (epsilon);
      std::vector<double> initial_value {0.0, 0.0, 0.0, 0.0};
      solution = nm.solve (initial_value);
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
