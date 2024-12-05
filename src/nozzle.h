#ifndef NOZZLE_H
#define NOZZLE_H

#include <math.h>

typedef enum ExpansionFlavor
{
    SUBSONIC,
    SONIC_THROAT,
    NORM_SHK_EXIT,
    NORM_SHK_IN,
    OBL_SHK,
    IDEAL,
    EXP_FAN,
} ExpansionFlavor;

typedef struct ExpSubsonic
{
    int dummy;
} ExpSubsonic;

typedef struct ExpSonicThroat
{
    int dummy;
} ExpSonicThroat;

typedef struct ExpNormShkExit
{
    double shock_loc;
} ExpNormShkExit;

typedef struct ExpNormShkIn
{
    double shock_loc;
} ExpNormShkIn;

typedef struct ExpOblShk
{
    double shk_angle;
} ExpOblShk;

typedef struct ExpIdeal
{
    int dummy;
} ExpIdeal;

typedef struct ExpFan
{
    double turn_angle;
} ExpFan;

typedef union ExpansionAny
{
    ExpSubsonic subsonic;
    ExpSonicThroat sonic_throat;
    ExpNormShkExit norm_shk_exit;
    ExpNormShkIn norm_shk_in;
    ExpOblShk obl_shk;
    ExpIdeal ideal;
    ExpFan expfan;
} ExpansionAny;

typedef struct Expansion
{
    ExpansionFlavor flavor;
    ExpansionAny props;
} Expansion;

double M2_norm_shock(double M1, double gamma);

double Me_norm_shock_in(double gm, double p0, double p_b, double Ae_At);

double rho_norm_shock(double M1, double gamma);

double T_norm_shock(double M1, double gamma);

double p_norm_shock(double M1, double gamma);

double p0_norm_shock(double M1, double gamma);

double calc_C_F(double p0, double p_e, double p_a, double Ae_At, double gm);

double calc_C_F_corr(double p0, double p_e, double p_a, double Ae_At, double gm,
                     double lm_div);

double calc_cstar(double T0, double Rspec, double gm);

#endif // NOZZLE_H
