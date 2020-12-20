
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_9032984207453160067) {
   out_9032984207453160067[0] = delta_x[0] + nom_x[0];
   out_9032984207453160067[1] = delta_x[1] + nom_x[1];
   out_9032984207453160067[2] = delta_x[2] + nom_x[2];
   out_9032984207453160067[3] = delta_x[3] + nom_x[3];
   out_9032984207453160067[4] = delta_x[4] + nom_x[4];
   out_9032984207453160067[5] = delta_x[5] + nom_x[5];
   out_9032984207453160067[6] = delta_x[6] + nom_x[6];
   out_9032984207453160067[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4725152637874738267) {
   out_4725152637874738267[0] = -nom_x[0] + true_x[0];
   out_4725152637874738267[1] = -nom_x[1] + true_x[1];
   out_4725152637874738267[2] = -nom_x[2] + true_x[2];
   out_4725152637874738267[3] = -nom_x[3] + true_x[3];
   out_4725152637874738267[4] = -nom_x[4] + true_x[4];
   out_4725152637874738267[5] = -nom_x[5] + true_x[5];
   out_4725152637874738267[6] = -nom_x[6] + true_x[6];
   out_4725152637874738267[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2450607807844962768) {
   out_2450607807844962768[0] = 1.0;
   out_2450607807844962768[1] = 0.0;
   out_2450607807844962768[2] = 0.0;
   out_2450607807844962768[3] = 0.0;
   out_2450607807844962768[4] = 0.0;
   out_2450607807844962768[5] = 0.0;
   out_2450607807844962768[6] = 0.0;
   out_2450607807844962768[7] = 0.0;
   out_2450607807844962768[8] = 0.0;
   out_2450607807844962768[9] = 1.0;
   out_2450607807844962768[10] = 0.0;
   out_2450607807844962768[11] = 0.0;
   out_2450607807844962768[12] = 0.0;
   out_2450607807844962768[13] = 0.0;
   out_2450607807844962768[14] = 0.0;
   out_2450607807844962768[15] = 0.0;
   out_2450607807844962768[16] = 0.0;
   out_2450607807844962768[17] = 0.0;
   out_2450607807844962768[18] = 1.0;
   out_2450607807844962768[19] = 0.0;
   out_2450607807844962768[20] = 0.0;
   out_2450607807844962768[21] = 0.0;
   out_2450607807844962768[22] = 0.0;
   out_2450607807844962768[23] = 0.0;
   out_2450607807844962768[24] = 0.0;
   out_2450607807844962768[25] = 0.0;
   out_2450607807844962768[26] = 0.0;
   out_2450607807844962768[27] = 1.0;
   out_2450607807844962768[28] = 0.0;
   out_2450607807844962768[29] = 0.0;
   out_2450607807844962768[30] = 0.0;
   out_2450607807844962768[31] = 0.0;
   out_2450607807844962768[32] = 0.0;
   out_2450607807844962768[33] = 0.0;
   out_2450607807844962768[34] = 0.0;
   out_2450607807844962768[35] = 0.0;
   out_2450607807844962768[36] = 1.0;
   out_2450607807844962768[37] = 0.0;
   out_2450607807844962768[38] = 0.0;
   out_2450607807844962768[39] = 0.0;
   out_2450607807844962768[40] = 0.0;
   out_2450607807844962768[41] = 0.0;
   out_2450607807844962768[42] = 0.0;
   out_2450607807844962768[43] = 0.0;
   out_2450607807844962768[44] = 0.0;
   out_2450607807844962768[45] = 1.0;
   out_2450607807844962768[46] = 0.0;
   out_2450607807844962768[47] = 0.0;
   out_2450607807844962768[48] = 0.0;
   out_2450607807844962768[49] = 0.0;
   out_2450607807844962768[50] = 0.0;
   out_2450607807844962768[51] = 0.0;
   out_2450607807844962768[52] = 0.0;
   out_2450607807844962768[53] = 0.0;
   out_2450607807844962768[54] = 1.0;
   out_2450607807844962768[55] = 0.0;
   out_2450607807844962768[56] = 0.0;
   out_2450607807844962768[57] = 0.0;
   out_2450607807844962768[58] = 0.0;
   out_2450607807844962768[59] = 0.0;
   out_2450607807844962768[60] = 0.0;
   out_2450607807844962768[61] = 0.0;
   out_2450607807844962768[62] = 0.0;
   out_2450607807844962768[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_4768208189780400912) {
   out_4768208189780400912[0] = state[0];
   out_4768208189780400912[1] = state[1];
   out_4768208189780400912[2] = state[2];
   out_4768208189780400912[3] = state[3];
   out_4768208189780400912[4] = state[4];
   out_4768208189780400912[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4768208189780400912[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4768208189780400912[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7115390955048770172) {
   out_7115390955048770172[0] = 1;
   out_7115390955048770172[1] = 0;
   out_7115390955048770172[2] = 0;
   out_7115390955048770172[3] = 0;
   out_7115390955048770172[4] = 0;
   out_7115390955048770172[5] = 0;
   out_7115390955048770172[6] = 0;
   out_7115390955048770172[7] = 0;
   out_7115390955048770172[8] = 0;
   out_7115390955048770172[9] = 1;
   out_7115390955048770172[10] = 0;
   out_7115390955048770172[11] = 0;
   out_7115390955048770172[12] = 0;
   out_7115390955048770172[13] = 0;
   out_7115390955048770172[14] = 0;
   out_7115390955048770172[15] = 0;
   out_7115390955048770172[16] = 0;
   out_7115390955048770172[17] = 0;
   out_7115390955048770172[18] = 1;
   out_7115390955048770172[19] = 0;
   out_7115390955048770172[20] = 0;
   out_7115390955048770172[21] = 0;
   out_7115390955048770172[22] = 0;
   out_7115390955048770172[23] = 0;
   out_7115390955048770172[24] = 0;
   out_7115390955048770172[25] = 0;
   out_7115390955048770172[26] = 0;
   out_7115390955048770172[27] = 1;
   out_7115390955048770172[28] = 0;
   out_7115390955048770172[29] = 0;
   out_7115390955048770172[30] = 0;
   out_7115390955048770172[31] = 0;
   out_7115390955048770172[32] = 0;
   out_7115390955048770172[33] = 0;
   out_7115390955048770172[34] = 0;
   out_7115390955048770172[35] = 0;
   out_7115390955048770172[36] = 1;
   out_7115390955048770172[37] = 0;
   out_7115390955048770172[38] = 0;
   out_7115390955048770172[39] = 0;
   out_7115390955048770172[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7115390955048770172[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7115390955048770172[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7115390955048770172[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7115390955048770172[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7115390955048770172[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7115390955048770172[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7115390955048770172[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7115390955048770172[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7115390955048770172[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7115390955048770172[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7115390955048770172[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7115390955048770172[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7115390955048770172[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7115390955048770172[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7115390955048770172[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7115390955048770172[56] = 0;
   out_7115390955048770172[57] = 0;
   out_7115390955048770172[58] = 0;
   out_7115390955048770172[59] = 0;
   out_7115390955048770172[60] = 0;
   out_7115390955048770172[61] = 0;
   out_7115390955048770172[62] = 0;
   out_7115390955048770172[63] = 1;
}
void h_25(double *state, double *unused, double *out_6402560618806271123) {
   out_6402560618806271123[0] = state[6];
}
void H_25(double *state, double *unused, double *out_9160002418099297048) {
   out_9160002418099297048[0] = 0;
   out_9160002418099297048[1] = 0;
   out_9160002418099297048[2] = 0;
   out_9160002418099297048[3] = 0;
   out_9160002418099297048[4] = 0;
   out_9160002418099297048[5] = 0;
   out_9160002418099297048[6] = 1;
   out_9160002418099297048[7] = 0;
}
void h_24(double *state, double *unused, double *out_4412050316000987318) {
   out_4412050316000987318[0] = state[4];
   out_4412050316000987318[1] = state[5];
}
void H_24(double *state, double *unused, double *out_988171895309078388) {
   out_988171895309078388[0] = 0;
   out_988171895309078388[1] = 0;
   out_988171895309078388[2] = 0;
   out_988171895309078388[3] = 0;
   out_988171895309078388[4] = 1;
   out_988171895309078388[5] = 0;
   out_988171895309078388[6] = 0;
   out_988171895309078388[7] = 0;
   out_988171895309078388[8] = 0;
   out_988171895309078388[9] = 0;
   out_988171895309078388[10] = 0;
   out_988171895309078388[11] = 0;
   out_988171895309078388[12] = 0;
   out_988171895309078388[13] = 1;
   out_988171895309078388[14] = 0;
   out_988171895309078388[15] = 0;
}
void h_30(double *state, double *unused, double *out_4692970059448074351) {
   out_4692970059448074351[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5038717930877977310) {
   out_5038717930877977310[0] = 0;
   out_5038717930877977310[1] = 0;
   out_5038717930877977310[2] = 0;
   out_5038717930877977310[3] = 0;
   out_5038717930877977310[4] = 1;
   out_5038717930877977310[5] = 0;
   out_5038717930877977310[6] = 0;
   out_5038717930877977310[7] = 0;
}
void h_26(double *state, double *unused, double *out_5550176116699538588) {
   out_5550176116699538588[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8551243072007920864) {
   out_8551243072007920864[0] = 0;
   out_8551243072007920864[1] = 0;
   out_8551243072007920864[2] = 0;
   out_8551243072007920864[3] = 0;
   out_8551243072007920864[4] = 0;
   out_8551243072007920864[5] = 0;
   out_8551243072007920864[6] = 0;
   out_8551243072007920864[7] = 1;
}
void h_27(double *state, double *unused, double *out_8988385170581227712) {
   out_8988385170581227712[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4708392490421332356) {
   out_4708392490421332356[0] = 0;
   out_4708392490421332356[1] = 0;
   out_4708392490421332356[2] = 0;
   out_4708392490421332356[3] = 1;
   out_4708392490421332356[4] = 0;
   out_4708392490421332356[5] = 0;
   out_4708392490421332356[6] = 0;
   out_4708392490421332356[7] = 0;
}
void h_29(double *state, double *unused, double *out_424647924512771714) {
   out_424647924512771714[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8381368370826071752) {
   out_8381368370826071752[0] = 0;
   out_8381368370826071752[1] = 1;
   out_8381368370826071752[2] = 0;
   out_8381368370826071752[3] = 0;
   out_8381368370826071752[4] = 0;
   out_8381368370826071752[5] = 0;
   out_8381368370826071752[6] = 0;
   out_8381368370826071752[7] = 0;
}
void h_28(double *state, double *unused, double *out_4205285541325701366) {
   out_4205285541325701366[0] = state[5];
   out_4205285541325701366[1] = state[6];
}
void H_28(double *state, double *unused, double *out_5121870455338762450) {
   out_5121870455338762450[0] = 0;
   out_5121870455338762450[1] = 0;
   out_5121870455338762450[2] = 0;
   out_5121870455338762450[3] = 0;
   out_5121870455338762450[4] = 0;
   out_5121870455338762450[5] = 1;
   out_5121870455338762450[6] = 0;
   out_5121870455338762450[7] = 0;
   out_5121870455338762450[8] = 0;
   out_5121870455338762450[9] = 0;
   out_5121870455338762450[10] = 0;
   out_5121870455338762450[11] = 0;
   out_5121870455338762450[12] = 0;
   out_5121870455338762450[13] = 0;
   out_5121870455338762450[14] = 1;
   out_5121870455338762450[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
