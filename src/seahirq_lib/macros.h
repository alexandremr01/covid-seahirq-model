#pragma once


#ifndef MAX_DAYS_REAL
#define MAX_DAYS_REAL 200
#endif

// Numeros de faixas etarias
#ifndef NEA
#define NEA 16
#endif

// Numero de variaveis
#ifndef NA
#define NA 13
#endif

#ifndef MAX_DAYS
#define MAX_DAYS 400
#endif

struct ScenarioParameters;
    
typedef void (*DerivFunc)(double, double[][NA], double[][NA], ScenarioParameters *);

struct Model {
    enum Enum {
        SIR,
        SEIR,
        SEAIR,
        SEAHIRQ
    };
};
