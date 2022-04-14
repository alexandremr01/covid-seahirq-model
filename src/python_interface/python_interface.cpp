#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../seahirq_lib/scenario.h"
#include "../seahirq_lib/models.h"
#include "../seahirq_lib/solvers.h"
#include <iostream>
#include <string>
#include <pybind11/eigen.h>
#include <exception>      // std::exception

#include <Eigen/LU>
#include <Eigen/Eigenvalues> 

#include <pybind11/stl.h>

namespace py = pybind11;

struct Output {
    std::vector<Eigen::Matrix<double, NA, 1>> Y_sum;  // for each day: # individuals per disease stage
    std::vector<Eigen::Matrix<double, NEA, NA>> YOUT; // for each day: # invidiuals per age strata and disease stage

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
    Eigen::Matrix<double, NEA, 1>  xI;
    Eigen::Matrix<double, NEA, 1>  xA;
    Eigen::Matrix<double, NEA, NEA> beta;
    PhaseParameters() {};
};

struct DynamicParameters {
    public:
        DynamicParameters(
            Eigen::Matrix<double, NEA, 1>  xI0,
            Eigen::Matrix<double, NEA, 1>  xA0,
            Eigen::Matrix<double, NEA, NEA> beta0
        ){ 
            // std::cout << " Phase 0 " << phase0Parameters.xI;
            PhaseParameters *phase0Parameters = new PhaseParameters();
            phase0Parameters->xI = xI0;
            phase0Parameters->xA = xA0;
            phase0Parameters->beta = beta0;
            this->vectorPhaseParameters.push_back(phase0Parameters); 
            this->days.push_back(0);
        };
        void addPhase(
            int startDay, 
            Eigen::Matrix<double, NEA, 1>  xI,
            Eigen::Matrix<double, NEA, 1>  xA,
            Eigen::Matrix<double, NEA, NEA> beta
        ){ 
            // std::cout << " Phase 0 " << phase0Parameters.xI;
            PhaseParameters *phaseParameters = new PhaseParameters();
            phaseParameters->xI = xI;
            phaseParameters->xA = xA;
            phaseParameters->beta = beta;
            // std::cout << " Phase: " << phaseParameters->xI;
            this->vectorPhaseParameters.push_back(phaseParameters); 
            this->days.push_back(startDay);
        };

        std::vector<PhaseParameters*> vectorPhaseParameters;
        std::vector<int> days;
};

struct StaticParameters {
    StaticParameters() {};
    
    Eigen::Matrix<double, NEA, 1> Lambda;
    // Taxa de mortalidade equilibrio
    Eigen::Matrix<double, NEA, 1> mu_eq;
    // Quanto a infectividade do assimptomatico é menor do que a do simptomatico(entre 0 e 1)
    Eigen::Matrix<double, NEA, 1> alpha;
    // Fator de correção para infecções pré-sintomáticas
    Eigen::Matrix<double, NEA, 1> ksi;
    // Probabilidade de um exposto se tornar sintomatico
    Eigen::Matrix<double, NEA, 1> rho;
    // Porcentagem de infectados que precisarao de hospitalizacao
    Eigen::Matrix<double, NEA, 1> phi;
    // Porcentagem de infectados que precisarao de leitos de UTI
    Eigen::Matrix<double, NEA, 1> eta;
    // Taxa de conversao de individuos infectados
    Eigen::Matrix<double, NEA, 1> a;
    // Taxa de mortalidade covid
    Eigen::Matrix<double, NEA, 1> mu_cov;
    // Coeficiente de taxa de recuperacao de assintomaticos (modelo SEAIR)
    Eigen::Matrix<double, NEA, 1> theta;
    // relação entre mortalidade de infectados e mortalidade de hospitalizados
    Eigen::Matrix<double, NEA, 1> gama_A;
    // Coeficiente de taxa de hospitalização de um infectado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_H;
    // Coeficiente de taxa de recuperação de hospitalizado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_HR;
    // Coeficiente de taxa de recuperação do infectado sintomático não quarentenado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RI;
    // Coeficiente de taxa de recuperação do infectado assintomático não quarentenado (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RA;
    // Coeficiente de taxa de recuperação do quarentenado sintomático (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RQI;
    // Coeficiente de taxa de recuperação do quarentenado assintomático (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_RQA;
    // Coeficiente de taxa de hospitalização do quarentenado sintomático (modelo SEAHIR-Qia)
    Eigen::Matrix<double, NEA, 1> gama_HQI;
    // Fração média de óbitos entre infectados
    Eigen::Matrix<double, NEA, 1> Tc;
    // Fração média de óbitos entre hospitalizados
    Eigen::Matrix<double, NEA, 1> Tlc;


    Eigen::Matrix<double, NEA, 1> gama;

    Eigen::Matrix<double, NEA, 1> rel_pop;
};

void read_array(double arr[NEA], Eigen::Matrix<double, NEA, 1> pyArr){
    for (int i = 0; i < NEA; i++) {
        arr[i] = pyArr(i, 0);
    }
}

void read_array_per_day(double array[MAX_DAYS][NEA], Eigen::Matrix<double, NEA, 1>  pyMatrix, int day) {
    for (int j = 0; j < NEA; j++) 
        for (int k=day; k < MAX_DAYS; k++){ 
            array[k][j] =  pyMatrix(j, 0);
                // std::cout << "Preenchendo dia " << k << " faixa " << j << " com " << pyMatrix(j, 0);
        }
}

void read_matrix_per_day(double matrix[MAX_DAYS][NEA][NEA], Eigen::Matrix<double, NEA, NEA>  pyMatrix, int day) {
    if (pyMatrix.rows() != NEA ||pyMatrix.rows() != NEA)
        throw std::exception();
    for (int i = 0; i < NEA; i++) 
        for (int j = 0; j < NEA; j++) 
            for (int k=day; k < MAX_DAYS; k++)
                matrix[k][i][j] =  pyMatrix(i, j);
}

void view_array(double arr[NEA], std::string name){
    std::cout<< name << ": ";
    for (int i = 0; i < NEA; i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

double spectral_radius(Eigen::Matrix<double, 6*NEA, 6*NEA> matrix){
    double largest = 0;
    Eigen::EigenSolver<Eigen::Matrix<double, 6*NEA, 6*NEA>> solver(6*NEA);
    solver.compute(matrix, false);

    Eigen::VectorXcd eivals = solver.eigenvalues();
    int r = eivals.rows();
    for (int i=0; i<r; i++){
        std::complex<double> z = eivals(i);
        if (std::abs(z) > std::abs(largest))
            largest = std::abs(z);
    }
    return largest;
}
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
#include <fstream>

double calculateR0(StaticParameters *p, DynamicParameters *dp) {
    // # A unica entrada de novas infecções é em E, vindo de A, I ou E
    Eigen::Matrix<double, 6*NEA, 6*NEA>  F;
    F.setZero(6*NEA, 6*NEA);

    struct Mat {
        enum Enum
        {
            E, A, I, Qi, Qa, H
        };
    }; 
    for (int i=0; i < NEA; i++){
        for (int j=0; j < NEA; j++){
            double val = dp->vectorPhaseParameters.at(0)->beta(i, j) * p->rel_pop(i);
            F(Mat::E*NEA + i, j) = p->ksi(j) * val;
            F(Mat::E*NEA +i, NEA + j) = p->alpha(j) * val; // A
            F(Mat::E*NEA +i, 2*NEA + j) = val;
        }
    }
    // # Transições entre compartimentos de infectados
    Eigen::Matrix<double, 6*NEA, 6*NEA>  V;
    V.setZero(6*NEA, 6*NEA);

    double gama_QI[NEA];
    double gama_QA[NEA];

    for (int k = 0; k < NEA; k++)
	{
        gama_QI[k] = dp->vectorPhaseParameters.at(0)->xI[k] * (p->gama_H[k] + p->gama_RI[k]) / (1 - dp->vectorPhaseParameters.at(0)->xI[k]);
        gama_QA[k] = dp->vectorPhaseParameters.at(0)->xA[k] * p->gama_RA[k] / (1 - dp->vectorPhaseParameters.at(0)->xA[k]);
	}

    for (int i = 0; i < NEA; i++) {
        int j;
    // # Exposto -> Exposto
        j = Mat::E + i;
        V(j, j) =  p->mu_eq(i) +  p->a(i);
    // # Assintomático
        j = Mat::A*NEA + i;
        V(j, j) = p->mu_eq(i) +  p->gama_RA(i) +  gama_QA[i];
        V(j, Mat::E*NEA+i) = - p->a(i) * (1- p->rho(i));
    // # Infectado sintomático
        j = Mat::I*NEA + i;
        V(j, j) = p->mu_eq(i) + p->theta(i)*p->mu_cov(i) + p->gama_H(i) + p->gama_RI(i) + gama_QI[i];
        V(j, Mat::E*NEA+i) = -p->a(i)*p->rho(i);
    // # Quarentenado sintomático
        j = Mat::Qi*NEA + i;
        V(j, j) = p->mu_eq(i) + p->gama_HQI(i) + p->gama_RQI(i);
        V(j, Mat::I*NEA+i) = -gama_QI[i];
    // # Quarentenado assintomático
        j = Mat::Qa*NEA + i;
        V(j, j) = p->mu_eq(i) + p->gama_RQA(i);
        V(j, Mat::A*NEA+i) = -gama_QA[i];
    // # Hospitalizado
        j = Mat::H*NEA + i;
        V(j, j) = p->mu_eq(i) + p->gama_HR(i) + p->mu_cov(i);
        V(j, Mat::I*NEA+i) = -p->gama_H(i);
        V(j, Mat::Qi*NEA+i) = -p->gama_HQI(i);
    }

    // std::cout << "writing to file" << std::endl;
    // std::ofstream file1("F.c.csv");
    // file1 << F.format(CSVFormat);
    // file1.close();
    // std::ofstream file2("V.c.csv");
    // file2 << V.format(CSVFormat);
    // file2.close();

    Eigen::Matrix<double, 6*NEA, 6*NEA> next_generation_matrix = F * V.inverse();
    double R0 = spectral_radius(next_generation_matrix);
    
    return R0;
}

Output run_model(int model, Eigen::Matrix<double, NA, NEA>  y0, StaticParameters *p, DynamicParameters *dp, bool verbose = false) {    
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
        if (verbose){
            std::cout << "xI day " << day << "values : " << dp->vectorPhaseParameters.at(i)->xI << std::endl;
            std::cout << "xA day " << day << "values : " << dp->vectorPhaseParameters.at(i)->xA << std::endl;
        }
    
        read_array_per_day(params.xI, dp->vectorPhaseParameters.at(i)->xI, day);
        read_array_per_day(params.xA, dp->vectorPhaseParameters.at(i)->xA, day);
        read_matrix_per_day(params.beta, dp->vectorPhaseParameters.at(i)->beta, day);
    }

    if (verbose){
        view_array(params.Lambda, "Lambda");
        view_array(params.mu_eq, "mu_eq");
        view_array(params.alpha, "alpha");
        view_array(params.ksi, "ksi");
        view_array(params.rho, "rho");
        view_array(params.phi, "phi");
        view_array(params.eta, "eta");
        view_array(params.a, "a");
        view_array(params.mu_cov, "mu_cov");
        view_array(params.theta, "theta");
        view_array(params.gama_A, "gama_A");
        view_array(params.gama_H, "gama_H");
        view_array(params.gama_HR, "gama_HR");
        view_array(params.gama_RI, "gama_RI");
        view_array(params.gama_RA, "gama_RA");
        view_array(params.gama_RQI, "gama_RQI");
        view_array(params.gama_RQA, "gama_RQA");
        view_array(params.gama_HQI, "gama_HQI");
        view_array(params.Tc, "Tc");
        view_array(params.Tlc, "Tlc");
    }

    
    for (int k = 0; k < NEA; k++)
	{
		for (int day = 0; day < MAX_DAYS; day++)
		{
            params.gama[day][k] = p->gama(k);
			params.gama_QI[day][k] = params.xI[day][k] * (params.gama_H[k] + params.gama_RI[k]) / (1 - params.xI[day][k]);
			params.gama_QA[day][k] = params.xA[day][k] * params.gama_RA[k] / (1 - params.xA[day][k]);
		}
	}
    for (int k = 0; k < NEA; k++)
	{
		for (int day = 0; day < MAX_DAYS; day++)
		{
            
            // std::cout << k << " " << day << " " << params.xI[day][k] << " " << params.gama_QI[day][k] << std::endl;
            params.gama[day][k] = p->gama(k);
			params.gama_QI[day][k] = params.xI[day][k] * (params.gama_H[k] + params.gama_RI[k]) / (1 - params.xI[day][k]);
			params.gama_QA[day][k] = params.xA[day][k] * params.gama_RA[k] / (1 - params.xA[day][k]);
		}
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
        .def_readwrite("xI", &PhaseParameters::xI)
        .def_readwrite("beta", &PhaseParameters::beta)
        .def_readwrite("xA", &PhaseParameters::xA);
    
    py::class_<DynamicParameters>(m, "DynamicParameters")
        .def(py::init<Eigen::Matrix<double, NEA, 1> ,
            Eigen::Matrix<double, NEA, 1> ,
            Eigen::Matrix<double, NEA, NEA> >())
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
        .def_readwrite("gama", &StaticParameters::gama)
        .def_readwrite("gama_A", &StaticParameters::gama_A)
        .def_readwrite("gama_H", &StaticParameters::gama_H)
        .def_readwrite("gama_HR", &StaticParameters::gama_HR)
        .def_readwrite("gama_RI", &StaticParameters::gama_RI)
        .def_readwrite("gama_RA", &StaticParameters::gama_RA)
        .def_readwrite("gama_RQI", &StaticParameters::gama_RQI)
        .def_readwrite("gama_RQA", &StaticParameters::gama_RQA)
        .def_readwrite("gama_HQI", &StaticParameters::gama_HQI)
        .def_readwrite("Tc", &StaticParameters::Tc)
        .def_readwrite("Tlc", &StaticParameters::Tlc)
        .def_readwrite("rel_pop", &StaticParameters::rel_pop);;

    m.def("model", &run_model, "A function which adds two numbers");
    m.def("calculateR0", &calculateR0, "A function which adds two numbers");

}