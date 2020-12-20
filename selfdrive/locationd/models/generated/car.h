/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_9032984207453160067);
void inv_err_fun(double *nom_x, double *true_x, double *out_4725152637874738267);
void H_mod_fun(double *state, double *out_2450607807844962768);
void f_fun(double *state, double dt, double *out_4768208189780400912);
void F_fun(double *state, double dt, double *out_7115390955048770172);
void h_25(double *state, double *unused, double *out_6402560618806271123);
void H_25(double *state, double *unused, double *out_9160002418099297048);
void h_24(double *state, double *unused, double *out_4412050316000987318);
void H_24(double *state, double *unused, double *out_988171895309078388);
void h_30(double *state, double *unused, double *out_4692970059448074351);
void H_30(double *state, double *unused, double *out_5038717930877977310);
void h_26(double *state, double *unused, double *out_5550176116699538588);
void H_26(double *state, double *unused, double *out_8551243072007920864);
void h_27(double *state, double *unused, double *out_8988385170581227712);
void H_27(double *state, double *unused, double *out_4708392490421332356);
void h_29(double *state, double *unused, double *out_424647924512771714);
void H_29(double *state, double *unused, double *out_8381368370826071752);
void h_28(double *state, double *unused, double *out_4205285541325701366);
void H_28(double *state, double *unused, double *out_5121870455338762450);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
