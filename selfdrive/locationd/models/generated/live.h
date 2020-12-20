/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7224527735921008559);
void inv_err_fun(double *nom_x, double *true_x, double *out_6250972390253074679);
void H_mod_fun(double *state, double *out_466223495476770526);
void f_fun(double *state, double dt, double *out_1767157372309052851);
void F_fun(double *state, double dt, double *out_5023570274078175190);
void h_3(double *state, double *unused, double *out_7140167188094692363);
void H_3(double *state, double *unused, double *out_236230850389828639);
void h_4(double *state, double *unused, double *out_1671790582997955301);
void H_4(double *state, double *unused, double *out_660516074330181767);
void h_9(double *state, double *unused, double *out_6753794233430078957);
void H_9(double *state, double *unused, double *out_1311201906635528488);
void h_10(double *state, double *unused, double *out_3697551042535373529);
void H_10(double *state, double *unused, double *out_7585822146060278075);
void h_12(double *state, double *unused, double *out_4775480389959497143);
void H_12(double *state, double *unused, double *out_2662331185030725080);
void h_31(double *state, double *unused, double *out_3977751566529858516);
void H_31(double *state, double *unused, double *out_8065806821902453616);
void h_32(double *state, double *unused, double *out_8152387222118460780);
void H_32(double *state, double *unused, double *out_3048835086191338114);
void h_13(double *state, double *unused, double *out_4049783251278300914);
void H_13(double *state, double *unused, double *out_1020065197045277602);
void h_14(double *state, double *unused, double *out_6753794233430078957);
void H_14(double *state, double *unused, double *out_1311201906635528488);
void h_19(double *state, double *unused, double *out_6093875439514080789);
void H_19(double *state, double *unused, double *out_7947782488440415268);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);