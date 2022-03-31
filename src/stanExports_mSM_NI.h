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
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_mSM_NI_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_mSM_NI");
    reader.add_event(26, 24, "end", "model_mSM_NI");
    return reader;
}
#include <stan_meta_header.hpp>
class model_mSM_NI
  : public stan::model::model_base_crtp<model_mSM_NI> {
private:
        int N;
        vector_d n;
        vector_d m;
        vector_d s2;
        vector_d priormu;
        matrix_d priorSigma;
public:
    model_mSM_NI(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_mSM_NI(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_mSM_NI_namespace::model_mSM_NI";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 3;
            validate_non_negative_index("n", "N", N);
            context__.validate_dims("data initialization", "n", "vector_d", context__.to_vec(N));
            n = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("n");
            pos__ = 0;
            size_t n_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < n_j_1_max__; ++j_1__) {
                n(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 4;
            validate_non_negative_index("m", "N", N);
            context__.validate_dims("data initialization", "m", "vector_d", context__.to_vec(N));
            m = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("m");
            pos__ = 0;
            size_t m_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < m_j_1_max__; ++j_1__) {
                m(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("s2", "N", N);
            context__.validate_dims("data initialization", "s2", "vector_d", context__.to_vec(N));
            s2 = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("s2");
            pos__ = 0;
            size_t s2_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < s2_j_1_max__; ++j_1__) {
                s2(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("priormu", "(N + 1)", (N + 1));
            context__.validate_dims("data initialization", "priormu", "vector_d", context__.to_vec((N + 1)));
            priormu = Eigen::Matrix<double, Eigen::Dynamic, 1>((N + 1));
            vals_r__ = context__.vals_r("priormu");
            pos__ = 0;
            size_t priormu_j_1_max__ = (N + 1);
            for (size_t j_1__ = 0; j_1__ < priormu_j_1_max__; ++j_1__) {
                priormu(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("priorSigma", "(N + 1)", (N + 1));
            validate_non_negative_index("priorSigma", "(N + 1)", (N + 1));
            context__.validate_dims("data initialization", "priorSigma", "matrix_d", context__.to_vec((N + 1),(N + 1)));
            priorSigma = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>((N + 1), (N + 1));
            vals_r__ = context__.vals_r("priorSigma");
            pos__ = 0;
            size_t priorSigma_j_2_max__ = (N + 1);
            size_t priorSigma_j_1_max__ = (N + 1);
            for (size_t j_2__ = 0; j_2__ < priorSigma_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < priorSigma_j_1_max__; ++j_1__) {
                    priorSigma(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            stan::math::check_cov_matrix(function__, "priorSigma", priorSigma);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 10;
            validate_non_negative_index("par", "(N + 1)", (N + 1));
            num_params_r__ += (N + 1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_mSM_NI() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 10;
        if (!(context__.contains_r("par")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable par missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("par");
        pos__ = 0U;
        validate_non_negative_index("par", "(N + 1)", (N + 1));
        context__.validate_dims("parameter initialization", "par", "vector_d", context__.to_vec((N + 1)));
        Eigen::Matrix<double, Eigen::Dynamic, 1> par((N + 1));
        size_t par_j_1_max__ = (N + 1);
        for (size_t j_1__ = 0; j_1__ < par_j_1_max__; ++j_1__) {
            par(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(par);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable par: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 10;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> par;
            (void) par;  // dummy to suppress unused var warning
            if (jacobian__)
                par = in__.vector_constrain((N + 1), lp__);
            else
                par = in__.vector_constrain((N + 1));
            // transformed parameters
            current_statement_begin__ = 13;
            validate_non_negative_index("a", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a(N);
            stan::math::initialize(a, DUMMY_VAR__);
            stan::math::fill(a, DUMMY_VAR__);
            current_statement_begin__ = 14;
            local_scalar_t__ invsigma2;
            (void) invsigma2;  // dummy to suppress unused var warning
            stan::math::initialize(invsigma2, DUMMY_VAR__);
            stan::math::fill(invsigma2, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 15;
            stan::math::assign(a, stan::math::exp(stan::model::rvalue(par, stan::model::cons_list(stan::model::index_min_max(1, N), stan::model::nil_index_list()), "par")));
            current_statement_begin__ = 16;
            stan::math::assign(invsigma2, stan::math::exp(get_base1(par, (N + 1), "par", 1)));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 13;
            size_t a_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < a_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(a(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: a" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable a: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 14;
            if (stan::math::is_uninitialized(invsigma2)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: invsigma2";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable invsigma2: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            // model body
            current_statement_begin__ = 19;
            lp_accum__.add(multi_normal_log<propto__>(par, priormu, priorSigma));
            current_statement_begin__ = 21;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 22;
                lp_accum__.add((((((-(0.5) * get_base1(n, i, "n", 1)) * stan::math::log((2 * stan::math::pi()))) + ((0.5 * get_base1(n, i, "n", 1)) * stan::math::log(invsigma2))) - (((0.5 * (get_base1(n, i, "n", 1) - 1)) * get_base1(s2, i, "s2", 1)) * invsigma2)) - (((0.5 * get_base1(n, i, "n", 1)) * square((get_base1(m, i, "m", 1) - get_base1(a, i, "a", 1)))) * invsigma2)));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("par");
        names__.push_back("a");
        names__.push_back("invsigma2");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back((N + 1));
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_mSM_NI_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> par = in__.vector_constrain((N + 1));
        size_t par_j_1_max__ = (N + 1);
        for (size_t j_1__ = 0; j_1__ < par_j_1_max__; ++j_1__) {
            vars__.push_back(par(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 13;
            validate_non_negative_index("a", "N", N);
            Eigen::Matrix<double, Eigen::Dynamic, 1> a(N);
            stan::math::initialize(a, DUMMY_VAR__);
            stan::math::fill(a, DUMMY_VAR__);
            current_statement_begin__ = 14;
            double invsigma2;
            (void) invsigma2;  // dummy to suppress unused var warning
            stan::math::initialize(invsigma2, DUMMY_VAR__);
            stan::math::fill(invsigma2, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 15;
            stan::math::assign(a, stan::math::exp(stan::model::rvalue(par, stan::model::cons_list(stan::model::index_min_max(1, N), stan::model::nil_index_list()), "par")));
            current_statement_begin__ = 16;
            stan::math::assign(invsigma2, stan::math::exp(get_base1(par, (N + 1), "par", 1)));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t a_j_1_max__ = N;
                for (size_t j_1__ = 0; j_1__ < a_j_1_max__; ++j_1__) {
                    vars__.push_back(a(j_1__));
                }
                vars__.push_back(invsigma2);
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_mSM_NI";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t par_j_1_max__ = (N + 1);
        for (size_t j_1__ = 0; j_1__ < par_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "par" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t a_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < a_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "a" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            param_name_stream__.str(std::string());
            param_name_stream__ << "invsigma2";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t par_j_1_max__ = (N + 1);
        for (size_t j_1__ = 0; j_1__ < par_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "par" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t a_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < a_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "a" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            param_name_stream__.str(std::string());
            param_name_stream__ << "invsigma2";
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_mSM_NI_namespace::model_mSM_NI stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
