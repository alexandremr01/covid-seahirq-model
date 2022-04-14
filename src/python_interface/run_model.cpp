#include "run_model.h"

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


    for (unsigned int i =0; i < dp->days.size(); i++){
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
