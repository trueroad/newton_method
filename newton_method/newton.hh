//
// Newton's method library to solve simultaneous equations 2018-08-23.15
// https://github.com/trueroad/newton_method
//
// newton.hh: Library header file
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

#ifndef INCLUDE_GUARD_NEWTON_METHOD_HH
#define INCLUDE_GUARD_NEWTON_METHOD_HH

#include<functional>
#include<memory>
#include<vector>

namespace newton_method
{
  enum class algorithm
    {
      PartialPivLU,
      FullPivLU,
      HouseholderQR,
      ColPivHouseholderQR,
      FullPivHouseholderQR,
      LLT,
      LDLT,
      JacobiSVD
    };

  enum class least_square
    {
      through_pass,
      weighted,
      normal_equation,
      weighted_normal_equation
    };

  enum class completion_status
    {
      none,
      epsilon_F,
      epsilon_deltaX,
      max_iteration_count_exceeded
    };

  class newton_method
  {
  public:
    // Set function
    void
    set_function (std::function<std::vector<double>
                    (const std::vector<double> &)> /* calc_function */,
                  std::function<std::vector<std::vector<double>>
                    (const std::vector<double> &)> /* calc_jacobian_matrix */
                  ) noexcept;

    // Set iteration parameters
    void set_max_iteration (int /* k */) noexcept;
    void set_max_iteration_exception (bool /* bte */ ) noexcept;
    void set_epsilon_F (double /* ef */) noexcept;
    void set_epsilon_deltaX (double /* ex */) noexcept;

    // Set algorithm
    void set_algorithm (algorithm /* a */) noexcept;
    void set_least_square (least_square /* ls */) noexcept;
    void set_weight (const std::vector<double> & /* weight */);

    // Solver
    std::vector<double>
    solve (const std::vector<double> & /* initial_value */);

    // Get status
    completion_status get_completion_status () noexcept;

    // Constructor and destructor
    newton_method ();
    ~newton_method ();

  private:
    // Pimpl
    class impl;
    std::unique_ptr<impl> pimpl_;

    // Deleted constructors
    newton_method (newton_method const &) = delete;
    newton_method (newton_method &&) = delete;
    newton_method & operator = (newton_method const &) = delete;
    newton_method & operator = (newton_method &&) = delete;
  };
}

#endif // INCLUDE_GUARD_NEWTON_METHOD_HH
