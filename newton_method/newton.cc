//
// Newton's method library to solve simultaneous equations 2018-08-26.02
// https://github.com/trueroad/newton_method
//
// newton.cc: Inner common implementation
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

#include"newton-private.hh"

//
// Pimpl constructor and destructor
//
namespace newton_method
{
  newton_method::newton_method ()
    :pimpl_ (new impl ())
  {
  }

  newton_method::~newton_method () = default;
}

//
// Pimpl delegater
//
namespace newton_method
{
  void
  newton_method::set_function (std::function<std::vector<double>
                                 (const std::vector<double> &)>
                                 calc_function,
                               std::function<std::vector<std::vector<double>>
                                 (const std::vector<double> &)>
                                 calc_jacobian_matrix) noexcept
  {
    pimpl_->set_function (calc_function, calc_jacobian_matrix);
  }

  void
  newton_method::set_function_fast
  (std::function<void (double *, double *, const double *)>
   calc_fast_f_and_j) noexcept
  {
    pimpl_->set_function_fast (calc_fast_f_and_j);
  }

  void newton_method::set_max_iteration (int k) noexcept
  {
    pimpl_->set_max_iteration (k);
  }

  void newton_method::set_max_iteration_exception (bool bte) noexcept
  {
    pimpl_->set_max_iteration_exception (bte);
  }

  void newton_method::set_epsilon_F (double ef) noexcept
  {
    pimpl_->set_epsilon_F (ef);
  }

  void newton_method::set_epsilon_deltaX (double ex) noexcept
  {
    pimpl_->set_epsilon_deltaX (ex);
  }

  void newton_method::set_weight (const std::vector<double> &weight)
  {
    pimpl_->set_weight (weight);
  }

  completion_status newton_method::get_completion_status () noexcept
  {
    return pimpl_->get_completion_status ();
  }
}
