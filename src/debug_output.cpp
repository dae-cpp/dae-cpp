/*
 * A set of helper functions to print on screen or write to files
 * Jacobian matrix, Mass matrix, the RHS for debugging purposes.
 */

#include <iostream>  // std::cout
#include <iomanip>   // std::setw etc.
#include <fstream>   // File output
#include <string>    // std::string, std::to_string

#include "RHS.h"
#include "jacobian.h"

namespace daecpp_namespace_name
{

const char delimiter = '\t';  // Delimiter of columns in output text files

/*
 * Helper function to write the RHS vector to a file
 */
void RHS::dump(const state_type &x, const double t)
{
    const MKL_INT size = x.size();

    state_type f(size);  // the vector to be saved

    this->operator()(x, f, t);  // calls the RHS

    std::ofstream outFile;

    outFile.open("dump_RHS_" + std::to_string(m_dump_file_counter++) + ".txt");
    outFile << "t = " << t << ":\n";
    outFile << "i" << delimiter << "x[i]" << delimiter << "RHS[i]" << '\n';
    for(MKL_INT i = 0; i < size; i++)
        outFile << i << delimiter << x[i] << delimiter << f[i] << '\n';
    outFile.close();
}

/*
 * Helper function to show Jacobian structure
 */
void Jacobian::print(const state_type &x, const double t)
{
    if(x.size() > 1000)
    {
        std::cout << "\nJacobian::print -- too much output. Skipped.\n";
        return;
    }

    sparse_matrix_holder J;

    this->operator()(J, x, t);

    std::cout << std::right;
    std::cout << "\nJacobian matrix at time t = " << t << ':';
    std::cout << "\n-----------------------------------------\n";
    std::cout << std::setw(7) << "i" << std::setw(16) << "J.A |"
              << std::setw(10) << "J.ja |" << std::setw(8) << "J.ia";
    std::cout << "\n-----------------------------------------\n";

    std::size_t size = (J.A.size() > J.ia.size()) ? J.A.size() : J.ia.size();

    for(std::size_t i = 0; i < size; i++)
    {
        std::cout << std::setw(7) << i << ": ";
        std::cout << std::setw(12);

        if(i < J.A.size())
            std::cout << J.A[i];
        else
            std::cout << ' ';

        std::cout << " | " << std::setw(7);

        if(i < J.ja.size())
            std::cout << J.ja[i];
        else
            std::cout << ' ';

        std::cout << " | ";

        if(i < J.ia.size())
            std::cout << std::setw(7) << J.ia[i];

        std::cout << std::endl;
    }
}

}  // namespace daecpp_namespace_name
