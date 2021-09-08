#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "scenario.h"
#include "models.h"
#include "solvers.h"
#include <iostream>
#include <pybind11/eigen.h>
#include <exception>      // std::exception

#include <Eigen/LU>

namespace py = pybind11;

struct Parameters {
    Parameters() {};
    
    py::array_t<double>  Lambda;
    // Taxa de mortalidade equilibrio
    py::array_t<double>  mu_eq;
    // Quanto a infectividade do assimptomatico é menor do que a do simptomatico(entre 0 e 1)
    py::array_t<double>  alpha;
    // Fator de correção para infecções pré-sintomáticas
    py::array_t<double>  ksi;
    // Probabilidade de um exposto se tornar sintomatico
    py::array_t<double>  rho;
    // Porcentagem de infectados que precisarao de hospitalizacao
    py::array_t<double>  phi;
    // Porcentagem de infectados que precisarao de leitos de UTI
    py::array_t<double>  eta;
    // Taxa de conversao de individuos infectados
    py::array_t<double>  a;
    // Taxa de mortalidade covid
    py::array_t<double>  mu_cov;
    // Coeficiente de taxa de recuperacao de assintomaticos (modelo SEAIR)
    py::array_t<double>  theta;
    // relação entre mortalidade de infectados e mortalidade de hospitalizados
    py::array_t<double>  gama_A;
    // Coeficiente de taxa de hospitalização de um infectado (modelo SEAHIR-Qia)
    py::array_t<double>  gama_H;
    // Coeficiente de taxa de recuperação de hospitalizado (modelo SEAHIR-Qia)
    py::array_t<double>  gama_HR;
    // Coeficiente de taxa de recuperação do infectado sintomático não quarentenado (modelo SEAHIR-Qia)
    py::array_t<double>  gama_RI;
    // Coeficiente de taxa de recuperação do infectado assintomático não quarentenado (modelo SEAHIR-Qia)
    py::array_t<double>  gama_RA;
    // Coeficiente de taxa de recuperação do quarentenado sintomático (modelo SEAHIR-Qia)
    py::array_t<double>  gama_RQI;
    // Coeficiente de taxa de recuperação do quarentenado assintomático (modelo SEAHIR-Qia)
    py::array_t<double>  gama_RQA;
    // Coeficiente de taxa de hospitalização do quarentenado sintomático (modelo SEAHIR-Qia)
    py::array_t<double>  gama_HQI;
    // Fração média de óbitos entre infectados
    py::array_t<double>  Tc;
    // Fração média de óbitos entre hospitalizados
    py::array_t<double>  Tlc;

    // ---- Parâmetros que variam no tempo -----

    // Taxa de contaminacao
    // Eigen::Matrix<double, MAX_DAYS, NEA> beta[MAX_DAYS][NEA][NEA];
    // Taxa de recuperacao
    Eigen::Matrix<double, MAX_DAYS, NEA> gama;
    // Coeficiente de taxa de movimento do infectado sintomático para a quarentena (modelo SEAHIR-Qia)
    Eigen::Matrix<double, MAX_DAYS, NEA> gama_QI;
    // Coeficiente de taxa de movimento do infectado assintomático para a quarentena (modelo SEAHIR-Qia)
    Eigen::Matrix<double, MAX_DAYS, NEA> gama_QA;
    // Probabilidade de testar um infectado sintomatico (modelo SEAHIR-Qia)
    Eigen::Matrix<double, MAX_DAYS, NEA> xI;
    // Probabilidade de testar um infectado assintomatico (modelo SEAHIR-Qia)
    Eigen::Matrix<double, MAX_DAYS, NEA> xA;
};

void read_array(double arr[NEA], py::array_t<double> pyArr){
    py::buffer_info info = pyArr.request();
    
    auto ptr = static_cast<double *>(info.ptr);

    int n = 1;
    for (auto r: info.shape) {
      n *= r;
    }

    if (n != NEA) { 
        throw std::exception();
    }

    for (int i = 0; i < NEA; i++) {
        arr[i] = *ptr++;
    }
}

void read_array_per_day(double array[MAX_DAYS][NEA], Eigen::Matrix<double, MAX_DAYS, NEA> matrix) {
    for (int i=0; i<MAX_DAYS; i++)
        for (int j=0; i<NEA; i++)
            array[i][j] = matrix(i, j);
}

void run_model(int model, Eigen::MatrixXd y0, Parameters *p) {
    std::cout << "Model chosen: "<<   model << std::endl;
    
    ScenarioParameters params;

    read_array(params.Lambda, p->Lambda);
    read_array(params.mu_eq, p->mu_eq);
    read_array(params.alpha, p->alpha);
    read_array(params.ksi, p->ksi);
    read_array(params.rho, p->rho);
    read_array(params.phi, p->phi);
    read_array(params.eta, p->eta);
    read_array(params.a, p->a);
    read_array(params.mu_cov, p->mu_cov);
    read_array(params.theta, p->theta);
    read_array(params.gama_A, p->gama_A);
    read_array(params.gama_H, p->gama_H);
    read_array(params.gama_HR, p->gama_HR);
    read_array(params.gama_RI, p->gama_RI);
    read_array(params.gama_RA, p->gama_RA);
    read_array(params.gama_RQI, p->gama_RQI);
    read_array(params.gama_RQA, p->gama_RQA);
    read_array(params.gama_HQI, p->gama_HQI);
    read_array(params.Tc, p->Tc);
    read_array(params.Tlc, p->Tlc);

    read_array_per_day(params.gama, p->gama);
    read_array_per_day(params.gama_QI, p->gama_QI);
    read_array_per_day(params.gama_QA, p->gama_QA);
    read_array_per_day(params.xI, p->xI);
    read_array_per_day(params.xA, p->xA);;

    for (int i = 0; i < NEA; i++) {
        std::cout << params.theta[i];
    }


    // DerivFunc derivs = get_model(Model::Enum(model));

        
    // driver2D_simple(derivs, y0, &params, &output);
    
	// if (write_output(output_filename, &output) == 1) {
	// 	return 1;
	// }
	// else {
	// 	printf("Exececucao bem sucedida. Verificar consistencia dos resultados.\n");
	// 	return 0;
	// }
}


PYBIND11_MODULE(cmodels, m) {
    m.doc() = "pybind11 cmodels plugin"; // optional module docstring

    py::class_<Parameters>(m, "Parameters")
        .def(py::init<>())
        .def_readwrite("Lambda", &Parameters::Lambda)
        .def_readwrite("mu_eq", &Parameters::mu_eq)
        .def_readwrite("alpha", &Parameters::alpha)
        .def_readwrite("ksi", &Parameters::ksi)
        .def_readwrite("rho", &Parameters::rho)
        .def_readwrite("phi", &Parameters::phi)
        .def_readwrite("eta", &Parameters::eta)
        .def_readwrite("a", &Parameters::a)
        .def_readwrite("mu_cov", &Parameters::mu_cov)
        .def_readwrite("theta", &Parameters::theta)
        .def_readwrite("gama_A", &Parameters::gama_A)
        .def_readwrite("gama_H", &Parameters::gama_H)
        .def_readwrite("gama_HR", &Parameters::gama_HR)
        .def_readwrite("gama_RI", &Parameters::gama_RI)
        .def_readwrite("gama_RA", &Parameters::gama_RA)
        .def_readwrite("gama_RQI", &Parameters::gama_RQI)
        .def_readwrite("gama_RQA", &Parameters::gama_RQA)
        .def_readwrite("gama_HQI", &Parameters::gama_HQI)
        .def_readwrite("Tc", &Parameters::Tc)
        .def_readwrite("Tlc", &Parameters::Tlc)
        .def_readwrite("gama", &Parameters::gama)
        .def_readwrite("gama_QI", &Parameters::gama_QI)
        .def_readwrite("gama_QA", &Parameters::gama_QA)
        .def_readwrite("xI", &Parameters::xI)
        .def_readwrite("xA", &Parameters::xA);

    m.def("model", &run_model, "A function which adds two numbers");

}