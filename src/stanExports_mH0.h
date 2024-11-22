// Generated by rstantools.  Do not edit by hand.

/*
    BMABMDR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BMABMDR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BMABMDR.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_mH0_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 51> locations_array__ =
  {" (found before start of program)",
  " (in 'string', line 32, column 2 to column 16)",
  " (in 'string', line 35, column 2 to column 9)",
  " (in 'string', line 36, column 2 to column 17)",
  " (in 'string', line 40, column 4 to column 28)",
  " (in 'string', line 39, column 44 to line 41, column 3)",
  " (in 'string', line 39, column 8 to line 41, column 3)",
  " (in 'string', line 38, column 4 to column 15)",
  " (in 'string', line 37, column 38 to line 39, column 3)",
  " (in 'string', line 37, column 2 to line 41, column 3)",
  " (in 'string', line 42, column 2 to column 24)",
  " (in 'string', line 45, column 3 to column 60)",
  " (in 'string', line 46, column 3 to column 47)",
  " (in 'string', line 53, column 6 to column 137)",
  " (in 'string', line 52, column 18 to line 54, column 5)",
  " (in 'string', line 52, column 4 to line 54, column 5)",
  " (in 'string', line 51, column 45 to line 55, column 4)",
  " (in 'string', line 51, column 9 to line 55, column 4)",
  " (in 'string', line 49, column 7 to column 126)",
  " (in 'string', line 48, column 18 to line 50, column 5)",
  " (in 'string', line 48, column 4 to line 50, column 5)",
  " (in 'string', line 47, column 39 to line 51, column 4)",
  " (in 'string', line 47, column 3 to line 55, column 4)",
  " (in 'string', line 19, column 2 to column 8)",
  " (in 'string', line 20, column 9 to column 10)",
  " (in 'string', line 20, column 2 to column 14)",
  " (in 'string', line 21, column 9 to column 10)",
  " (in 'string', line 21, column 2 to column 14)",
  " (in 'string', line 22, column 9 to column 10)",
  " (in 'string', line 22, column 2 to column 15)",
  " (in 'string', line 23, column 2 to column 13)",
  " (in 'string', line 24, column 2 to column 20)",
  " (in 'string', line 25, column 2 to column 27)",
  " (in 'string', line 26, column 2 to column 15)",
  " (in 'string', line 27, column 2 to column 15)",
  " (in 'string', line 28, column 2 to column 14)",
  " (in 'string', line 29, column 2 to column 16)",
  " (in 'string', line 3, column 4 to column 12)",
  " (in 'string', line 4, column 4 to column 12)",
  " (in 'string', line 5, column 4 to column 12)",
  " (in 'string', line 6, column 4 to column 12)",
  " (in 'string', line 7, column 4 to column 15)",
  " (in 'string', line 8, column 4 to column 14)",
  " (in 'string', line 9, column 4 to column 43)",
  " (in 'string', line 10, column 4 to column 42)",
  " (in 'string', line 11, column 4 to column 39)",
  " (in 'string', line 12, column 4 to column 38)",
  " (in 'string', line 13, column 4 to column 41)",
  " (in 'string', line 14, column 4 to column 28)",
  " (in 'string', line 15, column 4 to column 31)",
  " (in 'string', line 2, column 71 to line 16, column 3)"};
template <bool propto__, typename T0__, typename T1__, typename T2__,
          typename T3__, typename T4__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>,
                              stan::is_stan_scalar<T3__>,
                              stan::is_stan_scalar<T4__>>* = nullptr>
stan::promote_args_t<T0__, T1__, T2__, T3__, T4__>
pert_dist_lpdf(const T0__& theta, const T1__& lb, const T2__& md, const T3__&
               ub, const T4__& gama, std::ostream* pstream__);
template <bool propto__, typename T0__, typename T1__, typename T2__,
          typename T3__, typename T4__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>,
                              stan::is_stan_scalar<T3__>,
                              stan::is_stan_scalar<T4__>>*>
stan::promote_args_t<T0__, T1__, T2__, T3__, T4__>
pert_dist_lpdf(const T0__& theta, const T1__& lb, const T2__& md, const T3__&
               ub, const T4__& gama, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__, T4__>;
  int current_statement__ = 0;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  // suppress unused var warning
  (void) DUMMY_VAR__;
  try {
    local_scalar_t__ x1 = DUMMY_VAR__;
    local_scalar_t__ x2 = DUMMY_VAR__;
    local_scalar_t__ x3 = DUMMY_VAR__;
    local_scalar_t__ x4 = DUMMY_VAR__;
    local_scalar_t__ alpha = DUMMY_VAR__;
    local_scalar_t__ beta = DUMMY_VAR__;
    current_statement__ = 43;
    alpha = (1 + ((gama * (md - lb)) / (ub - lb)));
    current_statement__ = 44;
    beta = (1 + ((gama * (ub - md)) / (ub - lb)));
    current_statement__ = 45;
    x1 = ((alpha - 1) * stan::math::log((theta - lb)));
    current_statement__ = 46;
    x2 = ((beta - 1) * stan::math::log((ub - theta)));
    current_statement__ = 47;
    x3 = (((alpha + beta) - 1) * stan::math::log((ub - lb)));
    current_statement__ = 48;
    x4 = stan::math::lbeta(alpha, beta);
    current_statement__ = 49;
    return (((x1 + x2) - x3) - x4);
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
}
#include <stan_meta_header.hpp>
class model_mH0 final : public model_base_crtp<model_mH0> {
private:
  int N;
  Eigen::Matrix<double,-1,1> n_data__;
  Eigen::Matrix<double,-1,1> m_data__;
  Eigen::Matrix<double,-1,1> s2_data__;
  double shift;
  Eigen::Matrix<double,-1,1> priormu_data__;
  Eigen::Matrix<double,-1,-1> priorSigma_data__;
  double priorlb;
  double priorub;
  double priorg;
  int data_type;
  Eigen::Map<Eigen::Matrix<double,-1,1>> n{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double,-1,1>> m{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double,-1,1>> s2{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double,-1,1>> priormu{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> priorSigma{nullptr, 0, 0};
public:
  ~model_mH0() {}
  model_mH0(stan::io::var_context& context__, unsigned int random_seed__ = 0,
            std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_mH0_namespace::model_mH0";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 23;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 23;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 24;
      stan::math::validate_non_negative_index("n", "N", N);
      current_statement__ = 25;
      context__.validate_dims("data initialization", "n", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      n_data__ = Eigen::Matrix<double,-1,1>::Constant(N,
                   std::numeric_limits<double>::quiet_NaN());
      new (&n) Eigen::Map<Eigen::Matrix<double,-1,1>>(n_data__.data(), N);
      {
        std::vector<local_scalar_t__> n_flat__;
        current_statement__ = 25;
        n_flat__ = context__.vals_r("n");
        current_statement__ = 25;
        pos__ = 1;
        current_statement__ = 25;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 25;
          stan::model::assign(n, n_flat__[(pos__ - 1)],
            "assigning variable n", stan::model::index_uni(sym1__));
          current_statement__ = 25;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 26;
      stan::math::validate_non_negative_index("m", "N", N);
      current_statement__ = 27;
      context__.validate_dims("data initialization", "m", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      m_data__ = Eigen::Matrix<double,-1,1>::Constant(N,
                   std::numeric_limits<double>::quiet_NaN());
      new (&m) Eigen::Map<Eigen::Matrix<double,-1,1>>(m_data__.data(), N);
      {
        std::vector<local_scalar_t__> m_flat__;
        current_statement__ = 27;
        m_flat__ = context__.vals_r("m");
        current_statement__ = 27;
        pos__ = 1;
        current_statement__ = 27;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 27;
          stan::model::assign(m, m_flat__[(pos__ - 1)],
            "assigning variable m", stan::model::index_uni(sym1__));
          current_statement__ = 27;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 28;
      stan::math::validate_non_negative_index("s2", "N", N);
      current_statement__ = 29;
      context__.validate_dims("data initialization", "s2", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      s2_data__ = Eigen::Matrix<double,-1,1>::Constant(N,
                    std::numeric_limits<double>::quiet_NaN());
      new (&s2) Eigen::Map<Eigen::Matrix<double,-1,1>>(s2_data__.data(), N);
      {
        std::vector<local_scalar_t__> s2_flat__;
        current_statement__ = 29;
        s2_flat__ = context__.vals_r("s2");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 29;
          stan::model::assign(s2, s2_flat__[(pos__ - 1)],
            "assigning variable s2", stan::model::index_uni(sym1__));
          current_statement__ = 29;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 30;
      context__.validate_dims("data initialization", "shift", "double",
        std::vector<size_t>{});
      shift = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 30;
      shift = context__.vals_r("shift")[(1 - 1)];
      current_statement__ = 31;
      context__.validate_dims("data initialization", "priormu", "double",
        std::vector<size_t>{static_cast<size_t>(2)});
      priormu_data__ = Eigen::Matrix<double,-1,1>::Constant(2,
                         std::numeric_limits<double>::quiet_NaN());
      new (&priormu)
        Eigen::Map<Eigen::Matrix<double,-1,1>>(priormu_data__.data(), 2);
      {
        std::vector<local_scalar_t__> priormu_flat__;
        current_statement__ = 31;
        priormu_flat__ = context__.vals_r("priormu");
        current_statement__ = 31;
        pos__ = 1;
        current_statement__ = 31;
        for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
          current_statement__ = 31;
          stan::model::assign(priormu, priormu_flat__[(pos__ - 1)],
            "assigning variable priormu", stan::model::index_uni(sym1__));
          current_statement__ = 31;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 32;
      context__.validate_dims("data initialization", "priorSigma", "double",
        std::vector<size_t>{static_cast<size_t>(2), static_cast<size_t>(2)});
      priorSigma_data__ = Eigen::Matrix<double,-1,-1>::Constant(2, 2,
                            std::numeric_limits<double>::quiet_NaN());
      new (&priorSigma)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(priorSigma_data__.data(), 2,
        2);
      {
        std::vector<local_scalar_t__> priorSigma_flat__;
        current_statement__ = 32;
        priorSigma_flat__ = context__.vals_r("priorSigma");
        current_statement__ = 32;
        pos__ = 1;
        current_statement__ = 32;
        for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
          current_statement__ = 32;
          for (int sym2__ = 1; sym2__ <= 2; ++sym2__) {
            current_statement__ = 32;
            stan::model::assign(priorSigma, priorSigma_flat__[(pos__ - 1)],
              "assigning variable priorSigma",
              stan::model::index_uni(sym2__), stan::model::index_uni(sym1__));
            current_statement__ = 32;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 32;
      stan::math::check_cov_matrix(function__, "priorSigma", priorSigma);
      current_statement__ = 33;
      context__.validate_dims("data initialization", "priorlb", "double",
        std::vector<size_t>{});
      priorlb = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 33;
      priorlb = context__.vals_r("priorlb")[(1 - 1)];
      current_statement__ = 34;
      context__.validate_dims("data initialization", "priorub", "double",
        std::vector<size_t>{});
      priorub = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 34;
      priorub = context__.vals_r("priorub")[(1 - 1)];
      current_statement__ = 35;
      context__.validate_dims("data initialization", "priorg", "double",
        std::vector<size_t>{});
      priorg = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 35;
      priorg = context__.vals_r("priorg")[(1 - 1)];
      current_statement__ = 36;
      context__.validate_dims("data initialization", "data_type", "int",
        std::vector<size_t>{});
      data_type = std::numeric_limits<int>::min();
      current_statement__ = 36;
      data_type = context__.vals_i("data_type")[(1 - 1)];
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 2;
  }
  inline std::string model_name() const final {
    return "model_mH0";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ = "model_mH0_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,1> par =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(2, DUMMY_VAR__);
      current_statement__ = 1;
      par = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(2);
      local_scalar_t__ a = DUMMY_VAR__;
      local_scalar_t__ invsigma2 = DUMMY_VAR__;
      current_statement__ = 9;
      if ((stan::math::primitive_value(stan::math::logical_eq(data_type, 1))
          ||
          stan::math::primitive_value(stan::math::logical_eq(data_type, 3)))) {
        current_statement__ = 7;
        a = stan::model::rvalue(par, "par", stan::model::index_uni(1));
      } else {
        current_statement__ = 6;
        if ((stan::math::primitive_value(stan::math::logical_eq(data_type, 2))
            ||
            stan::math::primitive_value(stan::math::logical_eq(data_type, 4)))) {
          current_statement__ = 4;
          a = (stan::math::log(
                 stan::model::rvalue(par, "par", stan::model::index_uni(1)))
            - shift);
        }
      }
      current_statement__ = 10;
      invsigma2 = stan::math::exp(
                    stan::model::rvalue(par, "par", stan::model::index_uni(2)));
      {
        current_statement__ = 11;
        lp_accum__.add(pert_dist_lpdf<propto__>(
                         stan::model::rvalue(par, "par",
                           stan::model::index_uni(1)), priorlb,
                         stan::model::rvalue(priormu, "priormu",
                           stan::model::index_uni(1)), priorub, priorg,
                         pstream__));
        current_statement__ = 12;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(
                         stan::model::rvalue(par, "par",
                           stan::model::index_uni(2)),
                         stan::model::rvalue(priormu, "priormu",
                           stan::model::index_uni(2)),
                         stan::model::rvalue(priorSigma, "priorSigma",
                           stan::model::index_uni(2),
                           stan::model::index_uni(2))));
        current_statement__ = 22;
        if ((stan::math::primitive_value(stan::math::logical_eq(data_type, 1))
            ||
            stan::math::primitive_value(stan::math::logical_eq(data_type, 3)))) {
          current_statement__ = 20;
          for (int i = 1; i <= N; ++i) {
            current_statement__ = 18;
            lp_accum__.add((((((-0.5 *
              stan::model::rvalue(n, "n", stan::model::index_uni(i))) *
              stan::math::log((2 * stan::math::pi()))) + ((0.5 *
              stan::model::rvalue(n, "n", stan::model::index_uni(i))) *
              stan::math::log(invsigma2))) - (((0.5 *
              (stan::model::rvalue(n, "n", stan::model::index_uni(i)) - 1)) *
              stan::model::rvalue(s2, "s2", stan::model::index_uni(i))) *
              invsigma2)) - (((0.5 *
              stan::model::rvalue(n, "n", stan::model::index_uni(i))) *
              stan::math::square(
                (stan::model::rvalue(m, "m", stan::model::index_uni(i)) - a)))
              * invsigma2)));
          }
        } else {
          current_statement__ = 17;
          if ((stan::math::primitive_value(
                 stan::math::logical_eq(data_type, 2))
              ||
              stan::math::primitive_value(
                stan::math::logical_eq(data_type, 4)))) {
            current_statement__ = 15;
            for (int i = 1; i <= N; ++i) {
              current_statement__ = 13;
              lp_accum__.add(((((((-0.5 *
                stan::model::rvalue(n, "n", stan::model::index_uni(i))) *
                stan::math::log((2 * stan::math::pi()))) + ((0.5 *
                stan::model::rvalue(n, "n", stan::model::index_uni(i))) *
                stan::math::log(invsigma2))) - (((0.5 *
                (stan::model::rvalue(n, "n", stan::model::index_uni(i)) - 1))
                * stan::model::rvalue(s2, "s2", stan::model::index_uni(i))) *
                invsigma2)) - (((0.5 *
                stan::model::rvalue(n, "n", stan::model::index_uni(i))) *
                stan::math::square(
                  (stan::model::rvalue(m, "m", stan::model::index_uni(i)) -
                  a))) * invsigma2)) -
                (stan::model::rvalue(m, "m", stan::model::index_uni(i)) *
                stan::model::rvalue(n, "n", stan::model::index_uni(i)))));
            }
          }
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_mH0_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,1> par =
        Eigen::Matrix<double,-1,1>::Constant(2,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      par = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(2);
      double a = std::numeric_limits<double>::quiet_NaN();
      double invsigma2 = std::numeric_limits<double>::quiet_NaN();
      out__.write(par);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 9;
      if ((stan::math::primitive_value(stan::math::logical_eq(data_type, 1))
          ||
          stan::math::primitive_value(stan::math::logical_eq(data_type, 3)))) {
        current_statement__ = 7;
        a = stan::model::rvalue(par, "par", stan::model::index_uni(1));
      } else {
        current_statement__ = 6;
        if ((stan::math::primitive_value(stan::math::logical_eq(data_type, 2))
            ||
            stan::math::primitive_value(stan::math::logical_eq(data_type, 4)))) {
          current_statement__ = 4;
          a = (stan::math::log(
                 stan::model::rvalue(par, "par", stan::model::index_uni(1)))
            - shift);
        }
      }
      current_statement__ = 10;
      invsigma2 = stan::math::exp(
                    stan::model::rvalue(par, "par", stan::model::index_uni(2)));
      if (emit_transformed_parameters__) {
        out__.write(a);
        out__.write(invsigma2);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,1> par =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(2, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(par,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,1>>(2),
        "assigning variable par");
      out__.write(par);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "par", "double",
        std::vector<size_t>{static_cast<size_t>(2)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,1> par =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(2, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> par_flat__;
        current_statement__ = 1;
        par_flat__ = context__.vals_r("par");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
          current_statement__ = 1;
          stan::model::assign(par, par_flat__[(pos__ - 1)],
            "assigning variable par", stan::model::index_uni(sym1__));
          current_statement__ = 1;
          pos__ = (pos__ + 1);
        }
      }
      out__.write(par);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"par"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"a", "invsigma2"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(2)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{}, std::vector<size_t>{}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
      param_names__.emplace_back(std::string() + "par" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      param_names__.emplace_back(std::string() + "a");
      param_names__.emplace_back(std::string() + "invsigma2");
    }
    if (emit_generated_quantities__) {}
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
      param_names__.emplace_back(std::string() + "par" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      param_names__.emplace_back(std::string() + "a");
      param_names__.emplace_back(std::string() + "invsigma2");
    }
    if (emit_generated_quantities__) {}
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"par\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(2) + "},\"block\":\"parameters\"},{\"name\":\"a\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"invsigma2\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"par\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(2) + "},\"block\":\"parameters\"},{\"name\":\"a\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"invsigma2\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = 2;
    const size_t num_transformed = emit_transformed_parameters * ((1 + 1));
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = 2;
    const size_t num_transformed = emit_transformed_parameters * ((1 + 1));
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_mH0_namespace::model_mH0;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_mH0_namespace::profiles__;
}
#endif
#endif
