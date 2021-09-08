#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "scenario.h"
#include "models.h"
#include "solvers.h"
#include <iostream>
#include <pybind11/eigen.h>
#include <exception>      // std::exception

#include <Eigen/LU>
#include <pybind11/stl.h>

namespace py = pybind11;

struct Output {
    std::vector<Eigen::Matrix<double, NA, 1>> Y_sum;
    std::vector<Eigen::Matrix<double, NEA, NA>> YOUT;

    Output(double Y_sum[NA][MAX_DAYS], double YOUT[NEA][NA][MAX_DAYS]) {
        for (int i=0; i<MAX_DAYS; i++){
            Eigen::Matrix<double, NA, 1> Y_sum_daily;
            Eigen::Matrix<double, NEA, NA> YOUT_daily;
            for (int j=0; j<NA; j++){
                Y_sum_daily(j, 0) = Y_sum[j][i];
                for (int k=0; k<NEA; k++){
                    YOUT_daily(k, j) = YOUT[k][j][i];
                }
            }
            this->Y_sum.emplace_back(Y_sum_daily);
            this->YOUT.emplace_back(YOUT_daily);
        }
    };
    Output(){};
};


struct PhaseParameters {
    py::array_t<double> gama;
    py::array_t<double> gama_Qi;
    py::array_t<double> gama_Qa;
    py::array_t<double> xI;
    py::array_t<double> xA;
    Eigen::Matrix<double, NEA, NEA> beta;
    PhaseParameters() {};
};

struct DynamicParameters {
    public:
        DynamicParameters(PhaseParameters *phase0Parameters){ 
            this->vectorPhaseParameters.push_back(phase0Parameters); 
            this->days.push_back(0);
        };
        void addPhase(int startDay, PhaseParameters *phaseParameters){ 
            this->vectorPhaseParameters.push_back(phaseParameters); 
            this->days.push_back(startDay);
        };

        std::vector<PhaseParameters*> vectorPhaseParameters;
        std::vector<int> days;
};

struct StaticParameters {
    StaticParameters() {};
    
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

void read_array_per_day(double array[MAX_DAYS][NEA], py::array_t<double> pyArr, int day) {
    py::buffer_info info = pyArr.request();
    
    auto ptr = static_cast<double *>(info.ptr);

    int n = 1;
    for (auto r: info.shape) {
      n *= r;
    }

    if (n != NEA) { 
        throw std::exception();
    }

    for (int j = 0; j < NEA; j++) 
        for (int k=day; k < MAX_DAYS; k++)
            array[k][j] =  *ptr++;
}

void read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], Eigen::Matrix<double, NEA, NEA>  pyMatrix, int day) {
    if (pyMatrix.rows() != NEA ||pyMatrix.rows() != NEA)
        throw std::exception();
    for (int i = 0; i < NEA; i++) 
        for (int j = 0; j < NEA; j++) 
            for (int k=day; k < MAX_DAYS; k++)
                matrix[k][i][j] =  pyMatrix(i, j);
}

Output run_model(int model, Eigen::Matrix<double, NA, NEA>  y0, StaticParameters *p, DynamicParameters *dp) {    
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

    for (int i =0; i < dp->days.size(); i++){
        int day = dp->days.at(i);
        read_array_per_day(params.gama, dp->vectorPhaseParameters.at(i)->gama, day);
        read_array_per_day(params.gama_QI, dp->vectorPhaseParameters.at(i)->gama_Qi, day);
        read_array_per_day(params.gama_QA, dp->vectorPhaseParameters.at(i)->gama_Qa, day);
        read_array_per_day(params.xI, dp->vectorPhaseParameters.at(i)->xI, day);
        read_array_per_day(params.xA, dp->vectorPhaseParameters.at(i)->xA, day);
        read_matrix_per_day(params.beta, dp->vectorPhaseParameters.at(i)->beta, day);
    }

    for (int i = 0; i < NEA; i++) {
        std::cout << params.theta[i];
    }

    ScenarioOutput output;

    auto model_selector = static_cast<Model::Enum>(model);
    DerivFunc derivs = get_model(model_selector);
    driver2D_eigen(derivs, y0, &params, &output);
    Output out(output.Y_sum, output.YOUT);
    return out;
}


PYBIND11_MODULE(cmodels, m) {
    m.doc() = "pybind11 cmodels plugin"; // optional module docstring

    py::class_<Output>(m, "Output")
        .def(py::init<>())
        .def_readwrite("Y_sum", &Output::Y_sum)
        .def_readwrite("YOUT", &Output::YOUT);

    py::class_<PhaseParameters>(m, "PhaseParameters")
        .def(py::init<>())
        .def_readwrite("gama", &PhaseParameters::gama)
        .def_readwrite("gama_Qi", &PhaseParameters::gama_Qi)
        .def_readwrite("gama_Qa", &PhaseParameters::gama_Qa)        
        .def_readwrite("xI", &PhaseParameters::xI)
        .def_readwrite("beta", &PhaseParameters::beta)
        .def_readwrite("xA", &PhaseParameters::xA);
    
    py::class_<DynamicParameters>(m, "DynamicParameters")
        .def(py::init<PhaseParameters*>())
        .def("addPhase", &DynamicParameters::addPhase);

    py::class_<StaticParameters>(m, "StaticParameters")
        .def(py::init<>())
        .def_readwrite("Lambda", &StaticParameters::Lambda)
        .def_readwrite("mu_eq", &StaticParameters::mu_eq)
        .def_readwrite("alpha", &StaticParameters::alpha)
        .def_readwrite("ksi", &StaticParameters::ksi)
        .def_readwrite("rho", &StaticParameters::rho)
        .def_readwrite("phi", &StaticParameters::phi)
        .def_readwrite("eta", &StaticParameters::eta)
        .def_readwrite("a", &StaticParameters::a)
        .def_readwrite("mu_cov", &StaticParameters::mu_cov)
        .def_readwrite("theta", &StaticParameters::theta)
        .def_readwrite("gama_A", &StaticParameters::gama_A)
        .def_readwrite("gama_H", &StaticParameters::gama_H)
        .def_readwrite("gama_HR", &StaticParameters::gama_HR)
        .def_readwrite("gama_RI", &StaticParameters::gama_RI)
        .def_readwrite("gama_RA", &StaticParameters::gama_RA)
        .def_readwrite("gama_RQI", &StaticParameters::gama_RQI)
        .def_readwrite("gama_RQA", &StaticParameters::gama_RQA)
        .def_readwrite("gama_HQI", &StaticParameters::gama_HQI)
        .def_readwrite("Tc", &StaticParameters::Tc)
        .def_readwrite("Tlc", &StaticParameters::Tlc);

    m.def("model", &run_model, "A function which adds two numbers");

}