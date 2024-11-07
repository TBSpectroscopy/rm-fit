/**
-------------------------------------------------------------------------------------------------

 Thibault BERTIN
 Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES), C.P. 160/09
 Universite Libre de Bruxelles
 50 avenue F. D. Roosevelt, B-1050 Brussels, Belgium
 Phone: +32.2.650.24.18 - E-mail: thibault.bertin@ulb.be - Web: http://www.ulb.ac.be/cpm

-------------------------------------------------------------------------------------------------
**/



/**
 * Calculate the absorption cross section using speed-dependent uncorrelated hard collisions relaxation matrix formula developped by Ciurylo and Pine "Speed-dependent line mixing profiles" - JQSRT 67 (2000) 375-393 eq 2.13
 **/

#include <cmath>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
//#include <iostream>       // Uncomment this and all "std::cout" if progress indicator is desired

using namespace std::complex_literals; // Allow "i" "il" and "if" operators for complex

/** Trapezoid integration **/
Eigen::MatrixXcd integrate(const std::vector<Eigen::MatrixXcd> &arr, double inter = 1.0)
{
    int n_nu0 = arr[0].rows();
    Eigen::MatrixXcd integ = Eigen::MatrixXcd::Zero(n_nu0, n_nu0);
    for(unsigned int i = 0; i < arr.size() - 1; ++i)
    {
        Eigen::MatrixXcd temp = arr[i] + arr[i+1];
        integ += temp;
    }
    integ *= inter/2;
    return integ;
}



/** Calculate the absorption cross section for given lines using matrix inversions **/
pybind11::array_t<double> calc_abs(pybind11::array_t<double> vnu, pybind11::array_t<double> v, pybind11::array_t<double> mu, double vp, Eigen::MatrixXd nu0, Eigen::MatrixXd beta, Eigen::MatrixXd rho, Eigen::VectorXd X, std::vector<Eigen::MatrixXcd> W)
{
    int n_nu = vnu.size();
    int n_nu0 = nu0.rows();
    int n_v = v.size();
    int n_mu = mu.size();

    pybind11::buffer_info bnu = vnu.request();
    pybind11::buffer_info bv = v.request();
    pybind11::buffer_info bmu = mu.request();

    double *anu = (double*) bnu.ptr;
    double *av = (double*) bv.ptr;
    double *amu = (double*) bmu.ptr;

    const double c = 29979245800;
    const double c2_constant = 1.438776877;

    std::vector<Eigen::MatrixXcd> G(n_nu);
    std::vector<double> vsigma(n_nu);
    double inter_mu = (2/((double)n_mu - 1));
    double inter_v = av[n_v - 1]/((double)n_v - 1);
    //std::cout << "Calculating G(nu)" << std::endl;
    //std::cout << "0%" << std::flush;
    for(int i = 0; i < n_nu; ++i)
    {
        Eigen::MatrixXd nu = Eigen::MatrixXd::Identity(n_nu0, n_nu0);
        nu.diagonal() *= anu[i];
        std::vector<Eigen::MatrixXcd> integ_mu(n_mu);
        Eigen::MatrixXcd kv = Eigen::MatrixXcd::Identity(n_nu0, n_nu0); // Executed faster if declared here
        for (int j = 0; j < n_mu; ++j)
        {
            integ_mu[j].resize(n_nu0, n_nu0);
        }
        std::vector<Eigen::MatrixXcd> integ_v(n_v);
        for (int j = 0; j < n_v; ++j)
        {
            integ_v[j].resize(n_nu0, n_nu0);
            for (int k = 0; k < n_mu; ++k)
            {
                kv.diagonal() = nu.diagonal() * av[j] * amu[k] / c;
                integ_mu[k] = W[j] + beta - (1i*(nu - nu0 - kv));
                integ_mu[k] = integ_mu[k].inverse();

            }
            integ_v[j].noalias() = integrate(integ_mu, inter_mu) * (2.0 / (std::sqrt(M_PI) * std::pow(vp, 3))) * std::pow(av[j], 2) * std::exp(-std::pow(av[j]/vp, 2));
        }
        G[i].noalias() = integrate(integ_v, inter_v);
        //if (i % static_cast<int>(static_cast<double>(n_nu)/10) == 0 && i != 0)
        //{
        //    std::cout << "\r" << std::round(static_cast<double>(i)/static_cast<double>(n_nu) * 100) << "%" << std::flush;
        //}
        Eigen::MatrixXcd brack;
        brack.noalias() = Eigen::MatrixXcd::Identity(n_nu0, n_nu0) - (G[i]*beta);
        brack = brack.inverse() * G[i] * rho;
        Eigen::VectorXcd vec;
        vec.noalias() = brack * X;
        Eigen::Matrix<std::complex<double>, 1, 1> sigmac;
        sigmac.noalias() = X.transpose() * vec;
        vsigma[i] = sigmac(0,0).real() * anu[i] * (1.0 - std::exp(-c2_constant * anu[i])) / M_PI;
    }
    //std::cout << std::endl;
    pybind11::array sigma =  pybind11::cast(vsigma);
    return sigma;
}



/** Calculate the speed-dependent absorption cross section for given lines using partial matrix diagonalization **/
pybind11::array_t<double> calc_abs_diag(pybind11::array_t<double> vnu, pybind11::array_t<double> v, pybind11::array_t<double> mu, double vp, Eigen::MatrixXd nu0, Eigen::MatrixXd beta, Eigen::MatrixXd rho, Eigen::VectorXd X, std::vector<Eigen::MatrixXcd> W)
{
    int n_nu = vnu.size();
    int n_nu0 = nu0.rows();
    int n_v = v.size();
    int n_mu = mu.size();

    pybind11::buffer_info bnu = vnu.request();
    pybind11::buffer_info bv = v.request();
    pybind11::buffer_info bmu = mu.request();

    double *anu = (double*) bnu.ptr;
    double *av = (double*) bv.ptr;
    double *amu = (double*) bmu.ptr;

    const double c = 29979245800;
    const double c2_constant = 1.438776877;

    std::vector<Eigen::MatrixXcd> G(n_nu);
    std::vector<Eigen::MatrixXcd> A(n_v);
    std::vector<Eigen::MatrixXcd> Ainv(n_v);
    std::vector<Eigen::VectorXcd> Z(n_v);
    std::vector<double> vsigma(n_nu);
    double inter_mu = (2/((double)n_mu - 1));
    double inter_v = av[n_v - 1]/((double)n_v - 1);
    //std::cout << "Calculating G(nu)" << std::endl;
    //std::cout << "0%" << std::flush;

    for (int j = 0; j < n_v; ++j)
    {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es;
        es.compute(nu0 - (1i * (W[j] + beta)), true); // Compute eigenvalues only of (nu1 - iW -ibeta)
        A[j] = es.eigenvectors();
        Ainv[j] = es.eigenvectors().inverse();
        Z[j] = es.eigenvalues();
    }

    for(int i = 0; i < n_nu; ++i)
    {
        Eigen::MatrixXd nu = Eigen::MatrixXd::Identity(n_nu0, n_nu0);
        nu.diagonal() *= anu[i];
        std::vector<Eigen::MatrixXcd> integ_mu(n_mu);
        Eigen::MatrixXcd kv = Eigen::MatrixXcd::Identity(n_nu0, n_nu0); // Executed faster if declared here
        std::vector<Eigen::MatrixXcd> integ_v(n_v);
        for (int j = 0; j < n_v; ++j)
        {
            integ_v[j].resize(n_nu0, n_nu0);
            for (int k = 0; k < n_mu; ++k)
            {
                integ_mu[k] = Eigen::MatrixXcd::Zero(n_nu0, n_nu0);
                kv.diagonal() = nu.diagonal() * av[j] * amu[k] / c;
                integ_mu[k].diagonal().noalias() = ((1i * Z[j]) - (1i * (nu.diagonal() - kv.diagonal()))).cwiseInverse();
                integ_mu[k] = A[j] * integ_mu[k] * Ainv[j];
            }
            integ_v[j].noalias() = integrate(integ_mu, inter_mu) * (2.0 / (std::sqrt(M_PI) * std::pow(vp, 3))) * std::pow(av[j], 2) * std::exp(-std::pow(av[j]/vp, 2));
        }
        G[i].noalias() = integrate(integ_v, inter_v);
        //if (i % static_cast<int>(static_cast<double>(n_nu)/10) == 0 && i != 0)
        //{
        //    std::cout << "\r" << std::round(static_cast<double>(i)/static_cast<double>(n_nu) * 100) << "%" << std::flush;
        //}
        Eigen::MatrixXcd brack;
        brack.noalias() = Eigen::MatrixXcd::Identity(n_nu0, n_nu0) - (G[i]*beta);
        brack = brack.inverse() * G[i] * rho;
        Eigen::VectorXcd vec;
        vec.noalias() = brack * X;
        Eigen::Matrix<std::complex<double>, 1, 1> sigmac;
        sigmac.noalias() = X.transpose() * vec;
        vsigma[i] = sigmac(0,0).real() * anu[i] * (1.0 - std::exp(-c2_constant * anu[i])) / M_PI;
    }
    //std::cout << std::endl;
    pybind11::array sigma =  pybind11::cast(vsigma);
    return sigma;
}

/** Calculate the absorption cross section for given lines using matrix inversions and incorporating the coupling into the correlation matrix **/
pybind11::array_t<double> calc_abs_corr(pybind11::array_t<double> vnu, pybind11::array_t<double> v, pybind11::array_t<double> mu, double vp, Eigen::MatrixXd nu0, Eigen::MatrixXd beta, Eigen::MatrixXd rho, Eigen::VectorXd X, Eigen::MatrixXcd C, std::vector<Eigen::MatrixXcd> W)
{
    int n_nu = vnu.size();
    int n_nu0 = nu0.rows();
    int n_v = v.size();
    int n_mu = mu.size();

    pybind11::buffer_info bnu = vnu.request();
    pybind11::buffer_info bv = v.request();
    pybind11::buffer_info bmu = mu.request();

    double *anu = (double*) bnu.ptr;
    double *av = (double*) bv.ptr;
    double *amu = (double*) bmu.ptr;

    const double c = 29979245800;
    const double c2_constant = 1.438776877;

    std::vector<Eigen::MatrixXcd> G(n_nu);
    std::vector<double> vsigma(n_nu);
    double inter_mu = (2/((double)n_mu - 1));
    double inter_v = av[n_v - 1]/((double)n_v - 1);
    //std::cout << "Calculating G(nu)" << std::endl;
    //std::cout << "0%" << std::flush;
    for(int i = 0; i < n_nu; ++i)
    {
        Eigen::MatrixXd nu = Eigen::MatrixXd::Identity(n_nu0, n_nu0);
        nu.diagonal() *= anu[i];
        std::vector<Eigen::MatrixXcd> integ_mu(n_mu);
        Eigen::MatrixXcd kv = Eigen::MatrixXcd::Identity(n_nu0, n_nu0); // Executed faster if declared here
        for (int j = 0; j < n_mu; ++j)
        {
            integ_mu[j].resize(n_nu0, n_nu0);
        }
        std::vector<Eigen::MatrixXcd> integ_v(n_v);
        for (int j = 0; j < n_v; ++j)
        {
            integ_v[j].resize(n_nu0, n_nu0);
            for (int k = 0; k < n_mu; ++k)
            {
                kv.diagonal() = nu.diagonal() * av[j] * amu[k] / c;
                integ_mu[k] = W[j] + beta - (1i*(nu - nu0 - kv));
                integ_mu[k].diagonal() = integ_mu[k].diagonal().cwiseInverse();

            }
            integ_v[j].noalias() = integrate(integ_mu, inter_mu) * (2.0 / (std::sqrt(M_PI) * std::pow(vp, 3))) * std::pow(av[j], 2) * std::exp(-std::pow(av[j]/vp, 2));
        }
        G[i] = integrate(integ_v, inter_v);
        //if (i % static_cast<int>(static_cast<double>(n_nu)/10) == 0 && i != 0)
        //{
        //    std::cout << "\r" << std::round(static_cast<double>(i)/static_cast<double>(n_nu) * 100) << "%" << std::flush;
        //}
        Eigen::MatrixXcd brack;
        Eigen::MatrixXcd B;
        B.noalias() = beta - C;
        brack.noalias() = Eigen::MatrixXcd::Identity(n_nu0, n_nu0) - (G[i]*B);
        brack = brack.inverse() * G[i] * rho;
        Eigen::VectorXcd vec;
        vec.noalias() = brack * X;
        Eigen::Matrix<std::complex<double>, 1, 1> sigmac;
        sigmac.noalias() = X.transpose() * vec;
        vsigma[i] = sigmac(0,0).real() * anu[i] * (1.0 - std::exp(-c2_constant * anu[i])) / M_PI;
    }
    //std::cout << std::endl;
    pybind11::array sigma =  pybind11::cast(vsigma);
    return sigma;
}


PYBIND11_MODULE(calc_g, handle)
{
    handle.doc() = "C++ function to calculate G and alpha";
    handle.def("calc_abs", &calc_abs);
    handle.def("calc_abs_diag", &calc_abs_diag);
    handle.def("calc_abs_corr", &calc_abs_corr);
}

