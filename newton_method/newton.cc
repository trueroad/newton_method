//
// Newton's method library to solve simultaneous equations 2018-08-23.15
// https://github.com/trueroad/newton_method
//
// newton.cc: Inner implementation
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

#include"newton.hh"

#ifdef DEBUG_NEWTON_METHOD
#include<iostream>
#endif
#include<limits>
#include<stdexcept>

#define EIGEN_MPL2_ONLY
#include<Eigen/Dense>

//
// Inner class
//
namespace newton_method
{
  class newton_method::impl
  {
  public:
    // Set function
    inline void
    set_function (std::function<std::vector<double>
                    (const std::vector<double> &)> /* calc_function */,
                  std::function<std::vector<std::vector<double>>
                    (const std::vector<double> &)> /* calc_jacobian_matrix */
                  ) noexcept;

    // Set iteration parameters
    void set_max_iteration (int k) noexcept
    {
      max_iteration_ = k;
    }
    void set_max_iteration_exception (bool bte) noexcept
    {
      b_throw_exception_max_iteration_ = bte;
    }
    void set_epsilon_F (double ef) noexcept
    {
      epsilon_F_ = ef;
    }
    void set_epsilon_deltaX (double ex) noexcept
    {
      epsilon_deltaX_ = ex;
    }

    // Set algorithm
    void set_algorithm (algorithm a) noexcept
    {
      algorithm_ = a;
    }
    void set_least_square (least_square ls) noexcept
    {
      least_square_ = ls;
    }
    inline void set_weight (const std::vector<double> &weight);

    // Solver
    template <least_square LEAST_SQUARE, algorithm ALGORITHM>
    std::vector<double>
    solve (const std::vector<double> & /* initial_value */);

    // Get completion status
    completion_status get_completion_status () noexcept
    {
      return completion_status_;
    }

    // Constructor and destructor
    impl () = default;
    ~impl () = default;

  private:
    // Inner
    inline void check_args (void);
    template <typename T>
    void start_iteration_k (int /* k */, const T & /* X */);
    inline Eigen::VectorXd calc_F (const std::vector<double> & /* x */);
    inline bool check_F (const Eigen::VectorXd & /* F */);
    inline Eigen::MatrixXd calc_J (const std::vector<double> & /* x */);
    template <least_square LEAST_SQUARE, algorithm ALGORITHM,
              typename T_F, typename T_J>
    Eigen::VectorXd calc_deltaX (const T_F &F, const T_J &J);
    template <typename T_X>
    bool check_deltaX (const T_X & /* X */,
                       const Eigen::VectorXd & /* deltaX */);
    inline void exceed_max_iteration (void);

    // Function
    std::function<std::vector<double> (const std::vector<double> &)> f_;
    std::function<std::vector<std::vector<double>>
                              (const std::vector<double> &)> j_;
    bool bset_ = false;

    // Iteration parameters
    int max_iteration_ = 256;
    bool b_throw_exception_max_iteration_ = true;
    double epsilon_F_ = std::numeric_limits<double>::epsilon () * 16;
    double epsilon_deltaX_ = std::numeric_limits<double>::epsilon () * 16;

    // Algorithm
    algorithm algorithm_ = algorithm::ColPivHouseholderQR;
    least_square least_square_ = least_square::through_pass;
    Eigen::MatrixXd W_;

    // Status
    completion_status completion_status_ = completion_status::none;

    // Deleted constructors
    impl (newton_method::impl const &) = delete;
    impl (newton_method::impl &&) = delete;
    impl & operator = (newton_method::impl const &) = delete;
    impl & operator = (newton_method::impl &&) = delete;
  };
}

//
// Inner implementation
//
namespace newton_method
{
  inline void newton_method::impl::set_function
  (std::function<std::vector<double>
     (const std::vector<double> &)> calc_function,
   std::function<std::vector<std::vector<double>>
     (const std::vector<double> &)> calc_jacobian_matrix) noexcept
  {
    f_ = calc_function;
    j_ = calc_jacobian_matrix;
    bset_ = true;
  }

  inline void
  newton_method::impl::set_weight (const std::vector<double> &weight)
  {
    Eigen::Map<const Eigen::VectorXd> w {weight.data(),
        static_cast<int> (weight.size ())};
    W_ = w.asDiagonal ();

#ifdef DEBUG_NEWTON_METHOD
    std::cout << "W =" << std::endl << W_ << std::endl;
#endif
  }

  inline void
  newton_method::impl::check_args (void)
  {
    if ( !bset_ )
      {
        throw std::invalid_argument ("Function is not set.");
      }
    if ( max_iteration_ <= 0 )
      {
        throw std::invalid_argument
          ("Max iteration must be greater than 0.");
      }
  }

  template <typename T>
  void
  newton_method::impl::start_iteration_k (int k, const T &X)
  {
#ifdef DEBUG_NEWTON_METHOD
    std::cout << std::endl
              << "*** iteration " << k << " ***" << std::endl
              << "X =" << std::endl << X << std::endl;
#endif
  }

  inline Eigen::VectorXd
  newton_method::impl::calc_F (const std::vector<double> &x)
  {
    auto fv {f_ (x)};
    Eigen::Map<Eigen::VectorXd> F_mapped {fv.data (),
        static_cast<int> (fv.size ())};
    Eigen::VectorXd F {F_mapped};

#ifdef DEBUG_NEWTON_METHOD
    std::cout << "F =" << std::endl << F << std::endl;
#endif

    return F;
  }

  inline bool
  newton_method::impl::check_F (const Eigen::VectorXd &F)
  {
#ifdef DEBUG_NEWTON_METHOD
    std::cout
      << "F.norm = " << F.norm ();
#endif

    if (F.norm () < epsilon_F_)
      {
#ifdef DEBUG_NEWTON_METHOD
        std::cout
          << " is smaller than epsilon (" << epsilon_F_ << ")."
          << std::endl
          << "Completed." << std::endl;
#endif

        completion_status_ = completion_status::epsilon_F;

        return true;
      }
#ifdef DEBUG_NEWTON_METHOD
    std::cout
      << " is larger than epsilon (" << epsilon_F_ << ")."
      << std::endl
      << "Continuing." << std::endl;
#endif

    return false;
  }

  inline Eigen::MatrixXd
  newton_method::impl::calc_J (const std::vector<double> &x)
  {
    auto jvv {j_ (x)};
    Eigen::MatrixXd J {jvv.size (), jvv[0].size ()};

    for (int i = 0; i < jvv.size (); ++i)
      {
        Eigen::Map<Eigen::VectorXd> row {jvv[i].data (),
            static_cast<int> (jvv[i].size ())};
        J.row (i) = row;
      }

#ifdef DEBUG_NEWTON_METHOD
    std::cout << "J =" << std::endl << J << std::endl;
#endif

    return J;
  }

  template <least_square LEAST_SQUARE, algorithm ALGORITHM,
            typename T_F, typename T_J>
  Eigen::VectorXd
  newton_method::impl::calc_deltaX (const T_F &F, const T_J &J)
  {
    Eigen::VectorXd deltaX;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;

    switch (LEAST_SQUARE)
      {
      case least_square::through_pass:
        A = J;
        b = - F;
        break;
      case least_square::weighted:
        A = W_ * J;
        b = W_ * -F;
        break;
      case least_square::normal_equation:
        A = J.transpose () * J;
        b = - J.transpose () * F;
        break;
      case least_square::weighted_normal_equation:
        A = J.transpose () * W_ * J;
        b = - J.transpose () * W_ * F;
        break;
      default:
        throw std::invalid_argument
          ("Unknown least_square.");
        break;
      }

    switch (ALGORITHM)
      {
      case algorithm::PartialPivLU:
        deltaX = A.partialPivLu ().solve (b);
        break;
      case algorithm::FullPivLU:
        deltaX = A.fullPivLu ().solve (b);
        break;
      case algorithm::HouseholderQR:
        deltaX = A.householderQr ().solve (b);
        break;
      case algorithm::ColPivHouseholderQR:
        deltaX = A.colPivHouseholderQr ().solve (b);
        break;
      case algorithm::FullPivHouseholderQR:
        deltaX = A.fullPivHouseholderQr ().solve (b);
        break;
      case algorithm::LLT:
        deltaX = A.llt ().solve (b);
        break;
      case algorithm::LDLT:
        deltaX = A.ldlt ().solve (b);
        break;
      case algorithm::JacobiSVD:
        deltaX = A.jacobiSvd
          (Eigen::ComputeThinU | Eigen::ComputeThinV).solve (b);
        break;
      default:
        throw std::invalid_argument
          ("Unknown algorithm.");
        break;
      }

#ifdef DEBUG_NEWTON_METHOD
    std::cout << "deltaX =" << std::endl << deltaX << std::endl;
#endif

    return deltaX;
  }

  template <typename T_X>
  bool
  newton_method::impl::check_deltaX (const T_X &X,
                                     const Eigen::VectorXd &deltaX)
  {
    Eigen::VectorXd dX_X {deltaX.array () / X.array ()};
#ifdef DEBUG_NEWTON_METHOD
    std::cout
      << "deltaX.array / X.array = " << std::endl
      << dX_X << std::endl
      << "( deltaX.array / X.array ) norm = " << dX_X.norm ();
#endif

    if (dX_X.norm () < epsilon_deltaX_)
      {
#ifdef DEBUG_NEWTON_METHOD
        std::cout
          << " is smaller than epsilon (" << epsilon_deltaX_ << ")."
          << std::endl
          << "Completed." << std::endl;
#endif

        completion_status_ = completion_status::epsilon_deltaX;

        return true;
      }

#ifdef DEBUG_NEWTON_METHOD
    std::cout
      << " is larger than epsilon (" << epsilon_deltaX_ << ")."
      << std::endl
      << "Continuing." << std::endl;
#endif

    return false;
  }

  inline void
  newton_method::impl::exceed_max_iteration (void)
  {
#ifdef DEBUG_NEWTON_METHOD
    std::cout
      << "Although the number of iterations has exceeded the maximum ("
      << max_iteration_ << "), the solution does not converge."
      << std::endl;
#endif

    completion_status_ = completion_status::max_iteration_count_exceeded;

    if (b_throw_exception_max_iteration_)
      {
        throw std::runtime_error
          ("Although the number of iterations has exceeded the maximum, "
           "the solution does not converge.");
      }
  }

  template <least_square LEAST_SQUARE, algorithm ALGORITHM>
  std::vector<double>
  newton_method::impl::solve (const std::vector<double> &initial_value)
  {
    check_args ();

    // X and xv is linked. When X is updated, xv is also updated.
    std::vector<double> xv {initial_value};
    Eigen::Map<Eigen::VectorXd> X {xv.data (), static_cast<int> (xv.size ())};

    // Iteration
    for ( int k = 0; k < max_iteration_; ++k )
      {
        start_iteration_k (k, X);

        auto F {calc_F (xv)};
        if (check_F (F))
          return xv;

        auto J {calc_J (xv)};

        auto deltaX {calc_deltaX<LEAST_SQUARE, ALGORITHM> (F, J)};

        X = X + deltaX;  // xv is also updated.
        if (check_deltaX (X, deltaX))
          return xv;
      }

    exceed_max_iteration ();
    return xv;
  }
}

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

  void newton_method::set_algorithm (algorithm a) noexcept
  {
    pimpl_->set_algorithm (a);
  }

  void newton_method::set_least_square (least_square ls) noexcept
  {
    pimpl_->set_least_square (ls);
  }

  void newton_method::set_weight (const std::vector<double> &weight)
  {
    pimpl_->set_weight (weight);
  }

  template <least_square LEAST_SQUARE, algorithm ALGORITHM>
  std::vector<double>
  newton_method::solve (const std::vector<double> &initial_value)
  {
    return pimpl_->solve<LEAST_SQUARE, ALGORITHM> (initial_value);
  }

  completion_status newton_method::get_completion_status () noexcept
  {
    return pimpl_->get_completion_status ();
  }
}

//
// Explicit instantiation
//
namespace newton_method
{
  template
  std::vector<double>
  newton_method::solve<least_square::through_pass,
                       algorithm::ColPivHouseholderQR>
  (const std::vector<double> & /* initial_value */);
}
