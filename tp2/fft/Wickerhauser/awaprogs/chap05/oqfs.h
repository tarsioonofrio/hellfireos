/*
 * 
 * Declare arrays of orthogonal quadrature filter coefficients.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef OQFS_HDR_ALREADY_INCLUDED
#define OQFS_HDR_ALREADY_INCLUDED

#include "real.h"

#define SR2  (1.4142135623730950488) /* sqrt(2.0)  */
#define SR3  (1.7320508075688772935) /* sqrt(3.0)  */
#define SR10 (3.1622776601683793320) /* sqrt(10.0) */
#define SR15 (3.8729833462074168852) /* sqrt(15.0) */
#define A  (2.6613644236006609279) /* (1/4)sqrt(2)[1+sqrt(10)+sqrt(5+2sqrt(10))] */
#define B  (0.2818103350856762928) /* 3.0*0.25/A */


/****************************************************************
 * Orthogonal quadrature mirror filters
 *
 * Sources:
 *
 *     Gregory Beylkin
 *     Ronald R. Coifman
 *     Ingrid Daubechies
 *     P. P. Vaidyanathan
 */

/******************* Beylkin 18 ********************/

static real
  b18soqf[18] = {
    9.93057653743539270E-02,
    4.24215360812961410E-01,
    6.99825214056600590E-01, /* Maximum absolute value */
    4.49718251149468670E-01,
    -1.10927598348234300E-01,
    -2.64497231446384820E-01,
    2.69003088036903200E-02,
    1.55538731877093800E-01,
    -1.75207462665296490E-02,
    -8.85436306229248350E-02,
    1.96798660443221200E-02,
    4.29163872741922730E-02,
    -1.74604086960288290E-02,
    -1.43658079688526110E-02,
    1.00404118446319900E-02,
    1.48423478247234610E-03,
    -2.73603162625860610E-03,
    6.40485328521245350E-04 
    };
static int
  b18salpha = 0,
  b18somega = 17;

static real
  b18doqf[18] = {
    6.40485328521245350E-04,
    2.73603162625860610E-03,
    1.48423478247234610E-03,
    -1.00404118446319900E-02,
    -1.43658079688526110E-02,
    1.74604086960288290E-02,
    4.29163872741922730E-02,
    -1.96798660443221200E-02,
    -8.85436306229248350E-02,
    1.75207462665296490E-02,
    1.55538731877093800E-01,
    -2.69003088036903200E-02,
    -2.64497231446384820E-01,
    1.10927598348234300E-01,
    4.49718251149468670E-01,
    -6.99825214056600590E-01, /* Maximum absolute value */
    4.24215360812961410E-01,
    -9.93057653743539270E-02
    };
static int
  b18dalpha = 0,
  b18domega = 17;


/******************* Coifman 6 ********************/


static real
  c06soqf[6] = {
    -7.2732619512526 E-02,
    3.3789766245748 E-01,
    8.5257202021160 E-01,
    3.8486484686486 E-01,
    -7.2732619512526 E-02,
    -1.5655728135792 E-02
    };
static int
  c06salpha = 0,
  c06somega = 5;

static real
  c06doqf[6] = {
    1.5655728135792 E-02,
    -7.2732619512526 E-02,
    -3.8486484686486 E-01,
    8.5257202021160 E-01,
    -3.3789766245748 E-01,
    -7.2732619512526 E-02
  };
static int
  c06dalpha = 0,
  c06domega = 5;

/******************* Coifman 12 ********************/

static real
  c12soqf[12] = {
    1.63873364631797850E-02,
    -4.14649367819664850E-02,
    -6.73725547222998740E-02,
    3.86110066823092900E-01,
    8.12723635449606130E-01, /* Maximum absolute value */
    4.17005184423777600E-01,
    -7.64885990782645940E-02,
    -5.94344186464712400E-02,
    2.36801719468767500E-02,
    5.61143481936598850E-03,
    -1.82320887091009920E-03,
    -7.20549445368115120E-04
    };
static int
  c12salpha = 0,
  c12somega = 11;

static real
  c12doqf[12] = {
    -7.20549445368115120E-04,
    1.82320887091009920E-03,
    5.61143481936598850E-03,
    -2.36801719468767500E-02,
    -5.94344186464712400E-02,
    7.64885990782645940E-02,
    4.17005184423777600E-01,
    -8.12723635449606130E-01, /* Maximum absolute value */
    3.86110066823092900E-01,
    6.73725547222998740E-02,
    -4.14649367819664850E-02,
    -1.63873364631797850E-02
    };
static int
  c12dalpha = 0,
  c12domega = 11;

/******************* Coifman 18 ********************/

static real
  c18soqf[18] = {
    -3.79351286437787590E-03,
    7.78259642567078690E-03,
    2.34526961421191030E-02,
    -6.57719112814312280E-02,
    -6.11233900029556980E-02,
    4.05176902409616790E-01,
    7.93777222625620340E-01, /* Maximum absolute value */
    4.28483476377618690E-01,
    -7.17998216191705900E-02,
    -8.23019271063202830E-02,
    3.45550275733444640E-02,
    1.58805448636159010E-02,
    -9.00797613673228960E-03,
    -2.57451768812796920E-03,
    1.11751877082696180E-03,
    4.66216959820144030E-04,
    -7.09833025057049280E-05,
    -3.45997731974026950E-05
    };
static int
  c18salpha = 0,
  c18somega = 17;

static real
  c18doqf[18] = {
    -3.45997731974026950E-05,
    7.09833025057049280E-05,
    4.66216959820144030E-04,
    -1.11751877082696180E-03,
    -2.57451768812796920E-03,
    9.00797613673228960E-03,
    1.58805448636159010E-02,
    -3.45550275733444640E-02,
    -8.23019271063202830E-02,
    7.17998216191705900E-02,
    4.28483476377618690E-01,
    -7.93777222625620340E-01, /* Maximum absolute value */
    4.05176902409616790E-01,
    6.11233900029556980E-02,
    -6.57719112814312280E-02,
    -2.34526961421191030E-02,
    7.78259642567078690E-03,
    3.79351286437787590E-03
    };
static int
  c18dalpha = 0,
  c18domega = 17;

/******************* Coifman 24 ********************/

static real
  c24soqf[24] = {
    8.92313668220275710E-04,
    -1.62949201311084900E-03,
    -7.34616632765623490E-03,
    1.60689439640692360E-02,
    2.66823001556288040E-02,
    -8.12666996803130540E-02,
    -5.60773133164719500E-02,
    4.15308407030430150E-01,
    7.82238930920498790E-01, /* Maximum absolute value */
    4.34386056491468390E-01,
    -6.66274742630007520E-02,
    -9.62204420335636970E-02,
    3.93344271229132190E-02,
    2.50822618451469330E-02,
    -1.52117315272391490E-02,
    -5.65828668594603800E-03,
    3.75143615692490270E-03,
    1.26656192867951870E-03,
    -5.89020756811437840E-04,
    -2.59974552319421750E-04,
    6.23390338657646180E-05,
    3.12298760780433580E-05,
    -3.25968044485761290E-06,
    -1.78498455869993380E-06
    };
static int
  c24salpha = 0,
  c24somega = 23;

static real
  c24doqf[24] = {
    -1.78498455869993380E-06,
    3.25968044485761290E-06,
    3.12298760780433580E-05,
    -6.23390338657646180E-05,
    -2.59974552319421750E-04,
    5.89020756811437840E-04,
    1.26656192867951870E-03,
    -3.75143615692490270E-03,
    -5.65828668594603800E-03,
    1.52117315272391490E-02,
    2.50822618451469330E-02,
    -3.93344271229132190E-02,
    -9.62204420335636970E-02,
    6.66274742630007520E-02,
    4.34386056491468390E-01,
    -7.82238930920498790E-01, /* Maximum absolute value */
    4.15308407030430150E-01,
    5.60773133164719500E-02,
    -8.12666996803130540E-02,
    -2.66823001556288040E-02,
    1.60689439640692360E-02,
    7.34616632765623490E-03,
    -1.62949201311084900E-03,
    -8.92313668220275710E-04
    };
static int
  c24dalpha = 0,
  c24domega = 23;

/******************* Coifman 30 ********************/

static real
  c30soqf[30] = {
    -2.12080863336306810E-04,
    3.58589677255698600E-04,
    2.17823630484128470E-03,
    -4.15935878160399350E-03,
    -1.01311175380455940E-02,
    2.34081567615927950E-02,
    2.81680290621414970E-02,
    -9.19200105488064130E-02,
    -5.20431632162377390E-02,
    4.21566206728765440E-01,
    7.74289603740284550E-01, /* Maximum absolute value. */
    4.37991626228364130E-01,
    -6.20359639056089690E-02,
    -1.05574208705835340E-01,
    4.12892087407341690E-02,
    3.26835742832495350E-02,
    -1.97617790117239590E-02,
    -9.16423115304622680E-03,
    6.76418541866332000E-03,
    2.43337320922405380E-03,
    -1.66286376908581340E-03,
    -6.38131296151377520E-04,
    3.02259519791840680E-04,
    1.40541148901077230E-04,
    -4.13404844919568560E-05,
    -2.13150140622449170E-05,
    3.73459674967156050E-06,
    2.06380639023316330E-06,
    -1.67408293749300630E-07,
    -9.51579170468293560E-08
    };
static int
  c30salpha = 0,
  c30somega = 29;

static real
  c30doqf[30] = {
    -9.51579170468293560E-08,
    1.67408293749300630E-07,
    2.06380639023316330E-06,
    -3.73459674967156050E-06,
    -2.13150140622449170E-05,
    4.13404844919568560E-05,
    1.40541148901077230E-04,
    -3.02259519791840680E-04,
    -6.38131296151377520E-04,
    1.66286376908581340E-03,
    2.43337320922405380E-03,
    -6.76418541866332000E-03,
    -9.16423115304622680E-03,
    1.97617790117239590E-02,
    3.26835742832495350E-02,
    -4.12892087407341690E-02,
    -1.05574208705835340E-01,
    6.20359639056089690E-02,
    4.37991626228364130E-01,
    -7.74289603740284550E-01, /* Maximum absolute value. */
    4.21566206728765440E-01,
    5.20431632162377390E-02,
    -9.19200105488064130E-02,
    -2.81680290621414970E-02,
    2.34081567615927950E-02,
    1.01311175380455940E-02,
    -4.15935878160399350E-03,
    -2.17823630484128470E-03,
    3.58589677255698600E-04,
    2.12080863336306810E-04
    };
static int
  c30dalpha = 0,
  c30domega = 29;

/******************* Daubechies 2 ********************/
static real
  d02soqf[2] = {
    (0.5)*SR2,			/* 0.707106781186547, */
    (0.5)*SR2,			/* 0.707106781186547  */
  };
static int
  d02salpha = 0,
  d02somega = 1;

static real
  d02doqf[2] = {
    (0.5)*SR2,			/* 0.707106781186547, */
    (-0.5)*SR2,			/* -0.707106781186547  */
  };
static int
  d02dalpha = 0,
  d02domega = 1;

/******************* Daubechies 4 ********************/

static real
  d04soqf[4] = {
    (1.0+SR3)/(4.0*SR2),	/*  4.82962913144534160E-01, */
    (3.0+SR3)/(4.0*SR2),	/*  8.36516303737807940E-01, */
    (3.0-SR3)/(4.0*SR2),	/*  2.24143868042013390E-01, */
    (1.0-SR3)/(4.0*SR2)		/* -1.29409522551260370E-01 */
    };
static int
  d04salpha = 0,
  d04somega = 3;

static real
  d04doqf[4] = {
    (1.0-SR3)/(4.0*SR2),	/* -1.29409522551260370E-01, */
    (SR3-3.0)/(4.0*SR2),	/* -2.24143868042013390E-01, */
    (3.0+SR3)/(4.0*SR2),	/*  8.36516303737807940E-01, */
    -(1.0+SR3)/(4.0*SR2)	/* -4.82962913144534160E-01 */
    };
static int
  d04dalpha = 0,
  d04domega = 3;

/******************* Daubechies 6 ********************/

static real
  d06soqf[6] = {
    0.125*A,			/*   3.32670552950082630E-01, */
    0.125*(SR2+2.0*A-B),	/*   8.06891509311092550E-01, */
    0.125*(3.0*SR2-2.0*B),	/*   4.59877502118491540E-01, */
    0.125*(3.0*SR2-2.0*A),	/*  -1.35011020010254580E-01, */
    0.125*(SR2+2.0*B-A),	/*  -8.54412738820266580E-02, */
    0.125*B			/*   3.52262918857095330E-02  */
    };
static int
  d06salpha = 0,
  d06somega = 5;

static real
  d06doqf[6] = {
    0.125*B,			/*   3.52262918857095330E-02, */
    0.125*(A-SR2-2.0*B),	/*   8.54412738820266580E-02, */
    0.125*(3.0*SR2-2.0*A),	/*  -1.35011020010254580E-01, */
    0.125*(2.0*B-3.0*SR2),	/*  -4.59877502118491540E-01, */
    0.125*(SR2+2.0*A-B),	/*   8.06891509311092550E-01, */
    -0.125*A			/*  -3.32670552950082630E-01  */
    };
static int
  d06dalpha = 0,
  d06domega = 5;

/******************* Daubechies 8 ********************/

static real
  d08soqf[8] = {
    2.30377813309000010E-01,
    7.14846570553000050E-01, /* Maximum absolute value */
    6.30880767930000030E-01,
    -2.79837694169999990E-02,
    -1.87034811718999990E-01,
    3.08413818359999990E-02,
    3.28830116670000010E-02,
    -1.05974017850000000E-02
    };
static int
  d08salpha = 0,
  d08somega = 7;

static real
  d08doqf[8] = {
    -1.05974017850000000E-02,
    -3.28830116670000010E-02,
    3.08413818359999990E-02,
    1.87034811718999990E-01,
    -2.79837694169999990E-02,
    -6.30880767930000030E-01,
    7.14846570553000050E-01, /* Maximum absolute value */
    -2.30377813309000010E-01
    };
static int
  d08dalpha = 0,
  d08domega = 7;

/******************* Daubechies 10 ********************/

static real
  d10soqf[10] = {
    1.60102397974000000E-01,
    6.03829269797000020E-01,
    7.24308528437999980E-01, /* Maximum absolute value. */
    1.38428145901000000E-01,
    -2.42294887066000000E-01,
    -3.22448695850000020E-02,
    7.75714938400000050E-02,
    -6.24149021300000020E-03,
    -1.25807519990000000E-02,
    3.33572528500000010E-03
    };
static int
  d10salpha = 0,
  d10somega = 9;

static real
  d10doqf[10] = {
    3.33572528500000010E-03,
    1.25807519990000000E-02,
    -6.24149021300000020E-03,
    -7.75714938400000050E-02,
    -3.22448695850000020E-02,
    2.42294887066000000E-01,
    1.38428145901000000E-01,
    -7.24308528437999980E-01, /* Maximum absolute value. */
    6.03829269797000020E-01,
    -1.60102397974000000E-01
    };
static int
  d10dalpha = 0,
  d10domega = 9;

/******************* Daubechies 12 ********************/

static real
  d12soqf[12] = {
    1.11540743350000000E-01,
    4.94623890397999980E-01,
    7.51133908021000000E-01, /* Maximum absolute value. */
    3.15250351709000010E-01,
    -2.26264693965000010E-01,
    -1.29766867567000010E-01,
    9.75016055869999950E-02,
    2.75228655299999990E-02,
    -3.15820393180000010E-02,
    5.53842200999999980E-04,
    4.77725751100000020E-03,
    -1.07730108500000000E-03
    };
static int
  d12salpha = 0,
  d12somega = 11;

static real
  d12doqf[12] = {
    -1.07730108500000000E-03,
    -4.77725751100000020E-03,
    5.53842200999999980E-04,
    3.15820393180000010E-02,
    2.75228655299999990E-02,
    -9.75016055869999950E-02,
    -1.29766867567000010E-01,
    2.26264693965000010E-01,
    3.15250351709000010E-01,
    -7.51133908021000000E-01, /* Maximum absolute value. */
    4.94623890397999980E-01,
    -1.11540743350000000E-01
    };
static int
  d12dalpha = 0,
  d12domega = 11;

/******************* Daubechies 14 ********************/

static real
  d14soqf[14] = {
    7.78520540849999970E-02,
    3.96539319482000000E-01,
    7.29132090845999950E-01, /* Maximum absolute value. */
    4.69782287405000000E-01,
    -1.43906003928999990E-01,
    -2.24036184993999990E-01,
    7.13092192669999990E-02,
    8.06126091510000060E-02,
    -3.80299369350000010E-02,
    -1.65745416310000000E-02,
    1.25509985560000000E-02,
    4.29577973000000010E-04,
    -1.80164070400000000E-03,
    3.53713800000000020E-04
    };
static int
  d14salpha = 0,
  d14somega = 13;

static real
  d14doqf[14] = {
    3.53713800000000020E-04,
    1.80164070400000000E-03,
    4.29577973000000010E-04,
    -1.25509985560000000E-02,
    -1.65745416310000000E-02,
    3.80299369350000010E-02,
    8.06126091510000060E-02,
    -7.13092192669999990E-02,
    -2.24036184993999990E-01,
    1.43906003928999990E-01,
    4.69782287405000000E-01,
    -7.29132090845999950E-01, /* Maximum absolute value. */
    3.96539319482000000E-01,
    -7.78520540849999970E-02
    };
static int
  d14dalpha = 0,
  d14domega = 13;

/******************* Daubechies 16 ********************/

static real
  d16soqf[16] = {
    5.44158422430000010E-02,
    3.12871590914000020E-01,
    6.75630736296999990E-01, /* Maximum absolute value. */
    5.85354683654000010E-01,
    -1.58291052559999990E-02,
    -2.84015542961999990E-01,
    4.72484573999999990E-04,
    1.28747426619999990E-01,
    -1.73693010020000010E-02,
    -4.40882539310000000E-02,
    1.39810279170000000E-02,
    8.74609404700000050E-03,
    -4.87035299299999960E-03,
    -3.91740372999999990E-04,
    6.75449405999999950E-04,
    -1.17476784000000000E-04
    };
static int
  d16salpha = 0,
  d16somega = 15;

static real
  d16doqf[16] = {
    -1.17476784000000000E-04,
    -6.75449405999999950E-04,
    -3.91740372999999990E-04,
    4.87035299299999960E-03,
    8.74609404700000050E-03,
    -1.39810279170000000E-02,
    -4.40882539310000000E-02,
    1.73693010020000010E-02,
    1.28747426619999990E-01,
    -4.72484573999999990E-04,
    -2.84015542961999990E-01,
    1.58291052559999990E-02,
    5.85354683654000010E-01,
    -6.75630736296999990E-01, /* Maximum absolute value. */
    3.12871590914000020E-01,
    -5.44158422430000010E-02
    };
static int
  d16dalpha = 0,
  d16domega = 15;

/******************* Daubechies 18 ********************/

static real
  d18soqf[18] = {
    3.80779473639999980E-02,
    2.43834674613000010E-01,
    6.04823123690000020E-01,
    6.57288078050999980E-01, /* Maximum absolute value */
    1.33197385824999990E-01,
    -2.93273783279000000E-01,
    -9.68407832229999930E-02,
    1.48540749337999990E-01,
    3.07256814790000010E-02,
    -6.76328290610000020E-02,
    2.50947114999999980E-04,
    2.23616621239999990E-02,
    -4.72320475800000040E-03,
    -4.28150368199999970E-03,
    1.84764688300000000E-03,
    2.30385764000000010E-04,
    -2.51963189000000020E-04,
    3.93473200000000030E-05
    };
static int
  d18salpha = 0,
  d18somega = 17;

static real
  d18doqf[18] = {
    3.93473200000000030E-05,
    2.51963189000000020E-04,
    2.30385764000000010E-04,
    -1.84764688300000000E-03,
    -4.28150368199999970E-03,
    4.72320475800000040E-03,
    2.23616621239999990E-02,
    -2.50947114999999980E-04,
    -6.76328290610000020E-02,
    -3.07256814790000010E-02,
    1.48540749337999990E-01,
    9.68407832229999930E-02,
    -2.93273783279000000E-01,
    -1.33197385824999990E-01,
    6.57288078050999980E-01, /* Maximum absolute value */
    -6.04823123690000020E-01,
    2.43834674613000010E-01,
    -3.80779473639999980E-02
    };
static int
  d18dalpha = 0,
  d18domega = 17;

/******************* Daubechies 20 ********************/

static real
  d20soqf[20] = {
    2.66700579010000010E-02,
    1.88176800078000000E-01,
    5.27201188931999960E-01,
    6.88459039454000000E-01, /* Maximum absolute value */
    2.81172343661000020E-01,
    -2.49846424326999990E-01,
    -1.95946274376999990E-01,
    1.27369340336000000E-01,
    9.30573646040000060E-02,
    -7.13941471659999970E-02,
    -2.94575368219999990E-02,
    3.32126740589999970E-02,
    3.60655356700000010E-03,
    -1.07331754830000000E-02,
    1.39535174700000000E-03,
    1.99240529500000020E-03,
    -6.85856695000000030E-04,
    -1.16466855000000000E-04,
    9.35886700000000050E-05,
    -1.32642030000000010E-05
    };
static int
  d20salpha = 0,
  d20somega = 19;

static real
  d20doqf[20] = {
    -1.32642030000000010E-05,
    -9.35886700000000050E-05,
    -1.16466855000000000E-04,
    6.85856695000000030E-04,
    1.99240529500000020E-03,
    -1.39535174700000000E-03,
    -1.07331754830000000E-02,
    -3.60655356700000010E-03,
    3.32126740589999970E-02,
    2.94575368219999990E-02,
    -7.13941471659999970E-02,
    -9.30573646040000060E-02,
    1.27369340336000000E-01,
    1.95946274376999990E-01,
    -2.49846424326999990E-01,
    -2.81172343661000020E-01,
    6.88459039454000000E-01, /* Maximum absolute value */
    -5.27201188931999960E-01,
    1.88176800078000000E-01,
    -2.66700579010000010E-02
    };
static int
  d20dalpha = 0,
  d20domega = 19;

/******************* Vaidyanathan 24 ********************/

static real
  v24soqf[24] = {
    -6.29061181907475230E-05,
    3.43631904821029190E-04,
    -4.53956619637219290E-04,
    -9.44897136321949270E-04,
    2.84383454683556460E-03,
    7.08137504052444710E-04,
    -8.83910340861387800E-03,
    3.15384705589700400E-03,
    1.96872150100727140E-02,
    -1.48534480052300990E-02,
    -3.54703986072834530E-02,
    3.87426192934114400E-02,
    5.58925236913735480E-02,
    -7.77097509019694100E-02,
    -8.39288843661128300E-02,
    1.31971661416977720E-01,
    1.35084227129481260E-01,
    -1.94450471766478170E-01,
    -2.63494802488459910E-01,
    2.01612161775308660E-01,
    6.35601059872214940E-01, /* Maximum absolute value. */
    5.72797793210734320E-01,
    2.50184129504662180E-01,
    4.57993341109767180E-02
    };
static int
  v24salpha = 0,
  v24somega = 23;

static real
  v24doqf[24] = {
    4.57993341109767180E-02,
    -2.50184129504662180E-01,
    5.72797793210734320E-01,
    -6.35601059872214940E-01, /* Maximum absolute value. */
    2.01612161775308660E-01,
    2.63494802488459910E-01,
    -1.94450471766478170E-01,
    -1.35084227129481260E-01,
    1.31971661416977720E-01,
    8.39288843661128300E-02,
    -7.77097509019694100E-02,
    -5.58925236913735480E-02,
    3.87426192934114400E-02,
    3.54703986072834530E-02,
    -1.48534480052300990E-02,
    -1.96872150100727140E-02,
    3.15384705589700400E-03,
    8.83910340861387800E-03,
    7.08137504052444710E-04,
    -2.84383454683556460E-03,
    -9.44897136321949270E-04,
    4.53956619637219290E-04,
    3.43631904821029190E-04,
    6.29061181907475230E-05
    };
static int
  v24dalpha = 0,
  v24domega = 23;


#undef SR2 
#undef SR3 
#undef SR10
#undef SR15
#undef A 
#undef B 

#endif				/* OQFS_HDR_ALREADY_INCLUDED */
