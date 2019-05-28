#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define bPI 3.1415926535898
#define MU 3.986005e14
#define bOMEGAE84 7.2921151467e-5

typedef struct
{
    long double crs;
    long double delta_n;
    long double M0;
    long double cuc;
    long double e;
    long double cus;
    long double sqrt_a;
    long double TOE;
    long double cic;
    long double omega0;
    long double cis;
    long double i0;
    long double crc;
    long double omega_e;
    long double omega_point;
    long double i_point;
    long double WN;
}data;

void location(data *pdata)
{
    long double n,n0;
    long double M;
    long double t;
    long double E;
    long double f;
    long double u_dd;
    long double r_dd;
    long double delta_u, delta_r, delta_i;
    long double u, r, i;
    long double x, y,omega;
    long double X, Y, Z;
    n0 = sqrt(MU) / pow(pdata->sqrt_a, 3);
    n = n0 + pdata->delta_n;
    M = pdata->M0 + n * (t - pdata->TOE);
    E = M;
    for (int i = 0; i < 7; i++)
        E = M + pdata->e*sin(E);
    f = atan2(sqrt(1 - pdata->e*pdata->e)*sin(E), cos(E) - pdata->e);
    u_dd = pdata->omega_e + f;
    r_dd = pdata->sqrt_a*pdata->sqrt_a*(1 - pdata->e*cos(E));
    delta_u = pdata->cuc*cos(2 * u_dd) + pdata->cus*sin(2 * u_dd);
    delta_r = pdata->crc*cos(2 * u_dd) + pdata->crs*sin(2 * u_dd);
    delta_i = pdata->cic*cos(2 * u_dd) + pdata->cis*sin(2 * u_dd);
    u = u_dd + delta_u;
    r = r_dd + delta_r;
    i = pdata->i0 + pdata->i_point*(t - pdata->TOE) + delta_i;
    x = r * cos(u);
    y = r * sin(u);
    omega = pdata->omega0 + (pdata->omega_point - pdata->omega_e)*(t - pdata->TOE) - pdata->omega_e*pdata->TOE;
    X = x * cos(omega) - y * cos(i)*sin(omega);
    Y = x * sin(omega) - y * cos(i)*cos(omega);
    Z = y * sin(omega);
    printf("X=%f\nY=%f\nZ=%f\n",X,Y,Z);
}
void main()
{
    data *pdata;
    *pdata = { -1.721875000000E+01,4.502687555010E-09,1.413760604187E+00,-7.990747690201E-07,7.598234573379E-03,
        1.118145883083E-05,5.153709835052E+03,4.536000000000E+05,-1.303851604462E-08,-1.095067942661E-01,1.527369022369E-07,
        9.571235745530E-01, 1.640000000000E+02, -2.656176299285E+00, -8.037477650349E-09,-5.193073455211E-10 };
    location(pdata);
    system("pause");
}

/*long double roota = 5153.79589081;
    long double toe = 93600.0;
    long double m0 = 1.05827953357;
    long double e = 0.00223578442819;
    long double delta_n = 0.465376527657e-08;
    long double smallomega = 2.06374037770;
    long double cus = 0.177137553692e-05;
    long double cuc = 0.457651913166e-05;
    long double crs = 88.6875000000;
    long double crc = 344.968750000;
    long double cis = -0.856816768646e-07;
    long double cic = 0.651925802231e-07;
    long double idot = 0.342514267094e-09;
    long double i0 = 0.961685061380;
    long double bigomega0 = 1.64046615454;
    long double earthrate = bOMEGAE84;
    long double bigomegadot = -0.856928551657e-08;
    long double t = 86400.00;
    long double A;
    long double n0, n;
    long double tk;
    long double mk, ek, vk, tak, ik, omegak, phik, uk, rk;
    long double corr_u, corr_r, corr_i;
    long double xpk, ypk;
    long double xk, yk, zk;
    long double mkdot, ekdot, takdot, ukdot, ikdot, rkdot, omegakdot;
    long double xpkdot, ypkdot;
    long double xkdot, ykdot, zkdot;
    long iter;
    A = roota * roota; //roota is thesquareroot of A
    n0 = sqrt(MU/ (A*A*A)); //bGM84 is what theICD-200 calls Greek mu
    tk = t - toe; //t is the time of the pos. & vel. request.
    n = n0 + delta_n;
    mk = m0 + n * tk;
    mkdot = n;
    ek = mk;
    for (iter = 0; iter < 7; iter++) ek = mk + e * sin(ek); //Overkill for small e
    ekdot = mkdot / (1.0 - e * cos(ek));
    //In the line, below, tak is the trueanomaly (which is nu in the ICD-200).
    tak = atan2(sqrt(1.0 - e * e)*sin(ek), cos(ek) - e);
    takdot = sin(ek)*ekdot*(1.0 + e * cos(tak)) / (sin(tak)*(1.0 - e * cos(ek)));
    phik = tak + smallomega;
    corr_u = cus * sin(2.0*phik) + cuc * cos(2.0*phik);
    corr_r = crs * sin(2.0*phik) + crc * cos(2.0*phik);
    corr_i = cis * sin(2.0*phik) + cic * cos(2.0*phik);
    uk = phik + corr_u;
    rk = A * (1.0 - e * cos(ek)) + corr_r;
    ik = i0 + idot * tk + corr_i;
    ukdot = takdot + 2.0*(cus*cos(2.0*uk) - cuc * sin(2.0*uk))*takdot;
    rkdot = A * e*sin(ek)*n / (1.0 - e * cos(ek)) + 2.0*(crs*cos(2.0*uk) - crc * sin(2.0*uk))*takdot;
    ikdot = idot + (cis*cos(2.0*uk) - cic * sin(2.0*uk))*2.0*takdot;
    xpk = rk * cos(uk);
    ypk = rk * sin(uk);
    xpkdot = rkdot * cos(uk) - ypk * ukdot;
    ypkdot = rkdot * sin(uk) + xpk * ukdot;
    omegak = bigomega0 + (bigomegadot - earthrate)*tk - earthrate * toe;
    omegakdot = (bigomegadot - earthrate);
    xk = xpk * cos(omegak) - ypk * sin(omegak)*cos(ik);
    yk = xpk * sin(omegak) + ypk * cos(omegak)*cos(ik);
    zk = ypk * sin(ik);
    xkdot = (xpkdot - ypk * cos(ik)*omegakdot)*cos(omegak)
        - (xpk*omegakdot + ypkdot * cos(ik) - ypk * sin(ik)*ikdot)*sin(omegak);
    ykdot = (xpkdot - ypk * cos(ik)*omegakdot)*sin(omegak)
        + (xpk*omegakdot + ypkdot * cos(ik) - ypk * sin(ik)*ikdot)*cos(omegak);
    zkdot = ypkdot * sin(ik) + ypk * cos(ik)*ikdot;
    //Results follow.
    printf("BCpos: t, xk, yk, zk: %9.3Lf %21.11Lf %21.11Lf %21.11Lf\n", t, xk, yk, zk);
    //BCpos: t, xk, yk, zk: 86400.000 -12611434.19782218519 -13413103.97797041226 19062913.07357876760
        printf("BCvel: t, Vxk, Vyk, Vzk: %9.3Lf %16.10Lf %16.10Lf %16.10Lf\n", t, xkdot, ykdot, zkdot);
    //BCvel: t, Vxk, Vyk, Vzk: 86400.000 266.2803795674 -2424.7683468482 -1529.7620784616
    //Use the positions at 86400.000+0.005 and 86400.000-0.005 fornumerical computation check.
    //Perfect agreement is precluded becausewehave alimited precision machine.
        printf("xdotnumerical.01: %20.12lf\n", (-12611432.86642217014 - (-12611435.52922596496)) / 0.01);
    //xdotnumerical.01: 266.280379332602
    printf("ydotnumerical.01: %20.12lf\n", (-13413116.10180993562 - (-13413091.85412646334)) / 0.01);
    //ydotnumerical.01: -2424.768347293139
    printf("zdotnumerical.01: %20.12lf\n", (19062905.42476327563 - (19062920.72238405509)) / 0.01);
    //zdotnumerical.01: -1529.762077704072 */
