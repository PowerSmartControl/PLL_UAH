
#if defined(_WIN32) || defined(_WIN64)
  #define _USE_MATH_DEFINES
   #pragma warning (disable : 4068)
#endif

#include <math.h>
#include <float.h>
#include <stdint.h>
#include "cos_table.h"

#define CHANNELS_PER_ADC  8
#define N_ADCS            3
#define N_CHANNELS        (CHANNELS_PER_ADC*N_ADCS)

// Application constants:
#define TS ((const float) 9.77e-6)
#define FDSC ((const float) 50)
#define FPI ((const float) 20)

#define ADQ_GAIN_SENS_1    (0.2969) ///< ADQGAINXX * Sens Gain
#define ADQ_GAIN_SENS_9    (0.2965) ///< ADQGAINXX * Sens Gain
#define ADQ_GAIN_SENS_17   (0.2964) ///< ADQGAINXX * Sens Gain

#define ADQ_OFFSET_SENS_1  (-1.022) ///< Sensor offset = (VREF * Sens Gain)
#define ADQ_OFFSET_SENS_9  (-0.6827) ///< Sensor offset = (VREF * Sens Gain)
#define ADQ_OFFSET_SENS_17 (1.7784) ///< Sensor offset = (VREF * Sens Gain)

//-----------------------------------------------

#define ABC_ALPHA_BETA_11   ((const float) sqrt(2.0/3))
#define ABC_ALPHA_BETA_12   ((const float) (-0.5*sqrt(2.0/3)))

#define ABC_ALPHA_BETA_22   ((const float) (M_SQRT2/2))

#define abc2alphaBeta(a,b,c,alphaBeta)  \
  alphaBeta[0] = ABC_ALPHA_BETA_11*a + ABC_ALPHA_BETA_12*(b + c); \
  alphaBeta[1] = ABC_ALPHA_BETA_22*(b - c);

#define alphaBeta2dqPos(alphaBeta, theta, dq)  \
  dq[0] = cos2(theta)*alphaBeta[0] + sin2(theta)*alphaBeta[1];  \
  dq[1] = cos2(theta)*alphaBeta[1] - sin2(theta)*alphaBeta[0];

#define alphaBeta2dqNeg(alphaBeta, theta, dq)  \
  dq[0] = cos2(theta)*alphaBeta[0] - sin2(theta)*alphaBeta[1];  \
  dq[1] = sin2(theta)*alphaBeta[0] + cos2(theta)*alphaBeta[1];


#define COS_TABLE_SIZE (512)

#define M_PI_2_INT      (1*COS_TABLE_SIZE)
#define M_PI_INT        (2*COS_TABLE_SIZE)
#define M_3_PI_2_INT    (3*COS_TABLE_SIZE)
#define M_2_PI_INT      (4*COS_TABLE_SIZE)

#define M_PI_2_FP       ((1*COS_TABLE_SIZE)<<8)
#define M_PI_FP         ((2*COS_TABLE_SIZE)<<8)
#define M_3_PI_2_FP     ((3*COS_TABLE_SIZE)<<8)
#define M_2_PI_FP       ((4*COS_TABLE_SIZE)<<8)

#define radf2i(x) ((int) ((x) * (const float) ((M_2_PI_INT)/(2*M_PI))))

#define fmod2pi(angle)                   \
  if (angle >= ((const float) (2*M_PI))) \
    angle -= ((const float) (2*M_PI));   \
  if (angle < 0)                         \
    angle += ((const float) (2*M_PI));

#define imod2pi(angle)                                \
  if ((int) angle >= M_2_PI_INT)                      \
    angle = (uint16_t) ((int) angle - (M_2_PI_INT));  \
  if ((int) angle < 0)                                \
    angle = (uint16_t) ((int) angle + (M_2_PI_INT));

#define pow2(x) (x*x)
    
inline float cos2(uint16_t x) {
  float y;
  uint16_t op1, op2, x2;
  float sign;
  uint16_t i;

  x2 = x;

  if(0 <= x2 && x2 < M_PI_2_INT) {
    op1 = x2;
    op2 = 0;
    sign = 1.0;
  }
  else if(M_PI_2_INT <= x2 && x2 < M_PI_INT) {
    op1 = M_PI_INT-1;
    op2 = x2;
    sign = -1.0;
  }
  else if(M_PI_INT <= x2 && x2 < M_3_PI_2_INT) {
    op1 = x2;
    op2 = M_PI_INT;
    sign = -1.0;
  }
  else {
    op1 = M_2_PI_INT-1;
    op2 = x2;
    sign = 1.0;
  }

  i = op1 - op2;
  y = sign * cos_table[i];

  return y;
}

inline float sin2(uint16_t x) {
  float y;
  uint16_t op1, op2, x2;
  float sign;
  uint16_t i, j;

  x2 = x;

  if(0 <= x2 && x2 < M_PI_2_INT) {
    op1 = x2;
    op2 = 0;
    sign = 1.0;
  }
  else if(M_PI_2_INT <= x2 && x2 < M_PI_INT) {
    op1 = M_PI_INT-1;
    op2 = x2;
    sign = 1.0;
  }
  else if(M_PI_INT <= x2 && x2 < M_3_PI_2_INT) {
    op1 = x2;
    op2 = M_PI_INT;
    sign = -1.0;
  }
  else {
    op1 = M_2_PI_INT-1;
    op2 = x2;
    sign = -1.0;
  }

  i = op1 - op2;
  j = (COS_TABLE_SIZE-1) - i;
  y = sign * cos_table[j];

  return y;
}

    
typedef enum{
  EGA = 0,
  EGB = 1,
  EGC = 2,
  EDAT_SIZE = 3
} edat;
    
//-----------------------------------------------


#if 1

inline const float calc_kp(const float epsilon, const float wn, const float ts) {
  const float const1 = 2*(1-exp(-epsilon*wn*ts)*cos(wn*ts*sqrt(1-pow(epsilon, 2))));
  const float kp = const1/ts;

  return kp;
}

inline const float calc_ki(const float epsilon, const float wn, const float ts) {
  const float kp = calc_kp(epsilon, wn, ts);
  const float const1 = 2*(1-exp(-epsilon*wn*ts)*cos(wn*ts*sqrt(1-pow(epsilon, 2))));
  const float alfa = (1-exp(-2*epsilon*wn*ts))/const1;
  const float ki = (1-alfa)*kp/ts;
  
  return ki;
}

// Member varibles & const

const float k_adc = 1000.0;
const float ts = TS;
const float wn_dsc = (const float) (2*M_PI*FDSC);
const float kdsc = (const float) 1.4142;
const float wn_pi = (const float) (2*M_PI*FPI);
const float kp = calc_kp((float) 0.707, wn_pi, ts);
const float ki = calc_ki((float) 0.707, wn_pi, ts);

// top
float theta;
float theta_dly = 0.0;
float dq[2] = {0.0, 0.0};
float mod_edq = FLT_MIN;
bool sequence = 0;

// dsc (Variables SOGI_QSG)
float alfa = (float) 0.0;
float qalfa = (float) 0.0;
float beta = (float) 0.0;
float qbeta = (float) 0.0;

// pi
float integralPll = 0.0;  // Error modificado del PI del SPLL



#ifdef COMMON_SCALE_FACTOR

void scale(const short int adc[N_CHANNELS], float eabc[3])
{
#pragma HLS INLINE
  const unsigned char adc_n_bits = 12;
  const float scale_factor = (float) ((1.0/(1<<(adc_n_bits-1)))*k_adc);

  scale_loop : for (uint8_t i=0; i<3; i++)
  {
#pragma HLS PIPELINE II=1
    eabc[i] = adc[i*8+1] * scale_factor;
  }
}

#else

void scale(const short int adc[N_CHANNELS], float eabc[3])
{
#pragma HLS INLINE
  const float gain[3] = {ADQ_GAIN_SENS_1, ADQ_GAIN_SENS_9, ADQ_GAIN_SENS_17};
  const float offset[3] = {ADQ_OFFSET_SENS_1, ADQ_OFFSET_SENS_9, ADQ_OFFSET_SENS_17};
  
  scale_loop : for (uint8_t i=0; i<3; i++)
  {
#pragma HLS PIPELINE II=1
    eabc[i] = (adc[i*8+1] * gain[i])+ offset[i];    
  }
}

#endif


inline void sogi(float in1, const float ts, const float wn, float *v, float *qv)
{
  float error_qv, pre_v, pre_qv;

  error_qv = in1 - (*qv);
  pre_v = error_qv * (wn*ts);
  pre_qv = (*v) * (wn*ts);

  (*v) += pre_v;
  (*qv) += pre_qv;
}

void dsc(const float eabc[3], uint16_t itheta,
    float eAlphaBeta[2], float eAlphaBetaPos[2], float eAlphaBetaNeg[2],
    float edq[2], float edqPos[2], float edqNeg[2],
    float *mod_edq, bool *sequence)
{
#pragma HLS INLINE

  float modEdq_pow2, modEdqPos_pow2, modEdqNeg_pow2;
  float sogi1_in1, sogi2_in1;

  // TRANSFORMACION abc -> AlphaBeta
  abc2alphaBeta(eabc[EGA],eabc[EGB],eabc[EGC],eAlphaBeta);

  // SOGI_QSG para alfa
  sogi1_in1 = ((float) (0.5*kdsc)*eAlphaBeta[0])-(alfa*kdsc);
  sogi(sogi1_in1, ts, wn_dsc, &alfa, &qalfa);

  // SOGI_QSG para beta
  sogi2_in1 = ((float) (0.5*kdsc)*eAlphaBeta[1])-(beta*kdsc);
  sogi(sogi2_in1, ts, wn_dsc, &beta, &qbeta);

  eAlphaBetaPos[0] =  alfa - qbeta;
  eAlphaBetaPos[1] = qalfa +  beta;
  modEdqPos_pow2 = pow2(eAlphaBetaPos[0]) + pow2(eAlphaBetaPos[1]);

  eAlphaBetaNeg[0] = alfa + qbeta;
  eAlphaBetaNeg[1] = beta - qalfa;
  modEdqNeg_pow2 = pow2(eAlphaBetaNeg[0]) + pow2(eAlphaBetaNeg[1]);

  alphaBeta2dqPos(eAlphaBetaPos, itheta, edqPos);
  alphaBeta2dqNeg(eAlphaBetaNeg, itheta, edqNeg);

  if(modEdqPos_pow2 >= modEdqNeg_pow2)
  {
    modEdq_pow2 = modEdqPos_pow2;
    edq[0] = edqPos[0];
    edq[1] = edqPos[1];
    (*sequence) = false;
  }
  else
  {
    modEdq_pow2 = modEdqNeg_pow2;
    edq[0] = edqNeg[0];
    edq[1] = edqNeg[1];
    (*sequence) = true;
  }

  (*mod_edq) = sqrtf(modEdq_pow2);
}


void pi(const float edq[2], float mod_edq, bool sequence, float *theta_dly, uint16_t *itheta, float *omega)
{
#pragma HLS INLINE

  // Variables de acumulación. Por tratarse de variables de acumulación se definen static
  // Error modificado del PI del SPLL
  // Variables intermedias
  float errorPllpu;
  float sign;

  // Calculo del error
  if(sequence) {
    sign = kp;
  }
  else {
    sign = -kp;
  }

  if(mod_edq == 0) 
    mod_edq = FLT_MIN;
  errorPllpu = (sign * edq[0])/mod_edq;

  // PI
  (*omega) = errorPllpu + integralPll; // + FEEDFORWARD_PLL;
  // Error de integración añadiendo la ganancia antiwindup
  integralPll += ((ts*ki/kp)*(errorPllpu));

  // Calculo de la theta, integrando omega.
  theta = (*theta_dly) + (ts * (*omega));

  // Se acota el valor de theta entre 0 y 2pi.
  *itheta = radf2i(theta);
  imod2pi(*itheta);
  fmod2pi(theta);
  (*theta_dly) = (theta);
}



void unscale(const float alpha_beta[2], const float alpha_beta_p[2], const float alpha_beta_n[2],
  const float edq[2], const float edq_p[2], const float edq_n[2],
  float mod_edq, float theta, float omega, float sequence,
  float result[16])
{
#pragma HLS INLINE

  result[0] = alpha_beta[0];    //  ret->ealpha = alphabeta[0];
  result[1] = alpha_beta_p[0];  //  ret->ealpha_p = alphabeta_p[0];
  result[2] = alpha_beta_n[0];  //  ret->ealpha_n = alphabeta_n[0];

  result[3] = alpha_beta[1];    //  ret->ebeta = alphabeta[1];
  result[4] = alpha_beta_p[1];  //  ret->ebeta_p = alphabeta_p[1];
  result[5] = alpha_beta_n[1];  //  ret->ebeta_n = alphabeta_n[1];

  result[6] = edq[0];           //  ret->ed = edq[0];
  result[7] = edq_p[0];         //  ret->ed_p = edq_p[0];
  result[8] = edq_n[0];         //  ret->ed_n = edq_n[0];

  result[9] = edq[1];           //  ret->eq = edq[1];
  result[10] = edq_p[1];        //  ret->eq_p = edq_p[1];
  result[11] = edq_n[1];        //  ret->eq_n = edq_n[1];

  result[12] = theta;           //  ret->theta = theta;
  result[13] = omega;           //  ret->omega = omega;
  result[14] = mod_edq;         //  ret->modEdq = modEdq;
  result[15] = sequence;        //  ret->sequence = sequence;
}



void cpll(short int adc[N_CHANNELS], float result[16]) {
#pragma HLS INLINE

  float eabc[3];
  
  float alpha_beta[2];
  float alpha_beta_p[2];
  float alpha_beta_n[2];
  
  float dq_p[2];
  float dq_n[2];
  
  uint16_t itheta;
  float omega;

  scale(adc, eabc);

  pi(dq, mod_edq, sequence, &theta_dly, &itheta, &omega);

  dsc(eabc, itheta, 
    alpha_beta, alpha_beta_p, alpha_beta_n,
    dq, dq_p, dq_n,
    &mod_edq, &sequence);

  unscale(
    alpha_beta, alpha_beta_p, alpha_beta_n,
    dq, dq_p, dq_n,
    mod_edq, theta, omega, sequence,
    result);
}






#endif



void pll_dummy(short int adc[N_CHANNELS], float dout[16]) {
  int i;
  
  for(i = 0; i < 16; i++) {
//    dout[i] = (float) (adc[0] + adc[1] + i*adc[2]);
    dout[i] = (float) (adc[i]);
  }    
}


void pll_core(short int adc[N_CHANNELS], float dout[16]) {
// Para tener dos datos por acceso:
//#pragma HLS ARRAY_RESHAPE variable=adc block factor=2 dim=1


//#pragma HLS RESOURCE variable=adc core=SPRAMD
#pragma HLS RESOURCE variable=adc  core=RAM_1P_BRAM
#pragma HLS RESOURCE variable=dout core=RAM_1P_BRAM

// NOTE:
// The bram interface mode is functional identical to the ap_memory interface.
// The only difference is how the ports are implemented when the design is
// used in Vivado IP Integrator:
// * An ap_memory interface is displayed as multiple and separate ports.
// * A bram interface is displayed as a single grouped port which can be
// connected to a Xilinx block RAM using a single point-to-point connection.
#if 1
#pragma HLS INTERFACE ap_memory port=adc
#pragma HLS INTERFACE ap_memory port=dout
#else
#pragma HLS INTERFACE bram port=adc
#pragma HLS INTERFACE s_axilite port=dout bundle=BUS_A
#pragma HLS INTERFACE s_axilite port=return bundle=BUS_A
#endif

  
#ifdef USE_DUMMY
  pll_dummy(adc, dout);  
#else
  cpll(adc, dout);  
#endif

}
