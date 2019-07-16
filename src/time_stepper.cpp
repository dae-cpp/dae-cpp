/*
 * Simple yet efficient Adaptive Time Stepping algorithm
 */

#include <iostream>
#include <cmath>

#include "solver.h"

namespace daecpp_namespace_name
{

/*
 * Simple yet efficient Adaptive Time Stepping.
 * Returns negative status if the time stepping method failed to converge,
 * positive status if the time step should be restarted from the scratch,
 * or 0 in case of success.
 */
int Solver::m_adaptive_time_stepping(state_type &x,
                                     const state_type_matrix &x_prev, int iter)
{
    m_iterator_state.dt_eval = m_iterator_state.dt[0];

    if(m_opt.time_stepping == 1)  // S-SATS (Stability-based time stepping)
    {
        if(iter < m_opt.dt_increase_threshold)
        {
            m_increase_dt();
        }
        else if(iter >= m_opt.dt_decrease_threshold - 1)
        {
            m_decrease_dt();
            if(m_check_dt())
                return -1;  // Method failed to converge
        }
    }
    else if(m_opt.time_stepping == 2)  // V-SATS (Variability-based)
    {
        double norm1 = 0.0;
        double norm2 = 0.0;

        // Estimate NORM(C(n+1) - C(n)) and NORM(C(n))
        if(m_opt.vsats_norm == 2)  // NORM_2
        {
            for(MKL_INT i = 0; i < m_size; i++)
            {
                norm1 += (x[i] - x_prev[0][i]) * (x[i] - x_prev[0][i]);
                norm2 += x_prev[0][i] * x_prev[0][i];
            }

            norm1 = std::sqrt(norm1);
            norm2 = std::sqrt(norm2);
        }
        else  // NORM_inf (experimental)
        {
            for(MKL_INT i = 0; i < m_size; i++)
            {
                double adiff1 = std::abs(x[i] - x_prev[0][i]);
                double adiff2 = std::abs(x_prev[0][i]);

                if(adiff1 > norm1)
                {
                    norm1 = adiff1;

                    if(adiff2 < m_opt.dt_eps_m)
                        norm2 = 1;
                    else
                        norm2 = adiff2;
                }
            }
        }

        // Monitor function
        double eta = norm1 / (norm2 + m_opt.dt_eps_m);

        if(m_opt.verbosity > 1)
            std::cout << " (eta = " << eta << ")";

        // The time step should be reduced, scrape the current time iteration
        if(eta > m_opt.dt_eta_max)
        {
            if(m_opt.verbosity > 0)
                std::cout << " <- redo: dt_eta = " << eta;
            if(m_reset_ti_state(x, x_prev))
                return -2;  // Method failed to converge
            return 2;       // Re-run the current iteration
        }

        // The time step can be increased
        if(eta < m_opt.dt_eta_min)
        {
            m_increase_dt();
        }
    }
    else  // Internal error
    {
        return -10;
    }

    return 0;
}

/*
 * Scrapes the current time iteration and decreases the time step
 */
int Solver::m_reset_ti_state(state_type &x, const state_type_matrix &x_prev)
{
    m_iterator_state.t -= m_iterator_state.dt[0];
    m_iterator_state.step_counter_local--;
    m_steps--;
    m_iterator_state.final_time_step = false;
    m_iterator_state.dt[0] /= m_opt.dt_decrease_factor;
    m_iterator_state.dt_eval = m_iterator_state.dt[0];
    m_iterator_state.current_scheme = m_reset_ti_scheme();

    x = x_prev[0];

    return m_check_dt();
}

/*
 * Updates time integrator scheme when the time step changes
 */
int Solver::m_reset_ti_scheme()
{
    if(m_iterator_state.step_counter_local && m_opt.bdf_order == 2)
        return 2;  // BDF-2
    else
        return 1;  // BDF-1
}

/*
 * Increases the time step
 */
void Solver::m_increase_dt()
{
    m_iterator_state.dt_eval =
        m_iterator_state.dt[0] * m_opt.dt_increase_factor;

    if(!m_check_dt())
    {
        m_iterator_state.current_scheme = m_reset_ti_scheme();

        if(m_opt.verbosity > 0)
            std::cout << '>';
    }
}

/*
 * Decreases the time step
 */
void Solver::m_decrease_dt()
{
    m_iterator_state.dt_eval =
        m_iterator_state.dt[0] / m_opt.dt_decrease_factor;

    m_iterator_state.current_scheme = m_reset_ti_scheme();

    if(!m_check_dt() && m_opt.verbosity > 0)
        std::cout << '<';
}

/*
 * Checks if dt is within the interval defined in solver_options.h
 */
int Solver::m_check_dt()
{
    if(m_iterator_state.dt_eval < m_opt.dt_min)
    {
        std::cout << "\nERROR: The time step was reduced to "
                  << m_iterator_state.dt_eval
                  << " but the scheme failed to converge\n";
        return -1;
    }
    else if(m_iterator_state.dt_eval > m_opt.dt_max)
    {
        m_iterator_state.dt_eval = m_opt.dt_max;
        return 1;
    }
    else
    {
        return 0;
    }
}

}  // namespace daecpp_namespace_name
