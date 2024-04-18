/*
molDyn.c: In function ‘dscal’:
molDyn.c:83:19: error: ‘j’ not specified in enclosing ‘parallel’
   83 |             sx[j] *= sa;
      |                   ^~
molDyn.c:80:9: error: enclosing ‘parallel’
   80 | #pragma omp parallel for default(none) shared(n, sa, sx, incx) private(i)
      |         ^~~
molDyn.c: In function ‘forces’:
molDyn.c:175:22: error: ‘epot’ not specified in enclosing ‘parallel’
  175 |                 epot += (rrd6 - rrd3);
      |                      ^~
molDyn.c:138:9: error: enclosing ‘parallel’
  138 | #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez)
      |         ^~~
molDyn.c:177:21: error: ‘vir’ not specified in enclosing ‘parallel’
  177 |                 vir -= rd * r148;
      |                     ^~
molDyn.c:138:9: error: enclosing ‘parallel’
  138 | #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez)
      |         ^~~
molDyn.c: In function ‘mxwell’:
molDyn.c:352:9: error: ‘n3’ not specified in enclosing ‘parallel’
  352 | #pragma omp parallel for default(none) shared(vh, npart, tscale) private(i, s, v1, v2, r)
      |         ^~~
molDyn.c:352:9: error: enclosing ‘parallel’
molDyn.c:367:9: error: ‘n3’ not specified in enclosing ‘parallel’
  367 | #pragma omp parallel for default(none) shared(vh, npart, h) private(i)
      |         ^~~
molDyn.c:367:9: error: enclosing ‘parallel’
molDyn.c:369:12: error: ‘sp’ not specified in enclosing ‘parallel’
  369 |         sp += vh[i];
      |            ^~
molDyn.c:367:9: error: enclosing ‘parallel’
  367 | #pragma omp parallel for default(none) shared(vh, npart, h) private(i)
      |         ^~~
molDyn.c:371:9: error: ‘n3’ not specified in enclosing ‘parallel’
  371 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:371:9: error: enclosing ‘parallel’
molDyn.c:375:14: error: ‘ekin’ not specified in enclosing ‘parallel’
  375 |         ekin += vh[i] * vh[i];
      |              ^~
molDyn.c:371:9: error: enclosing ‘parallel’
  371 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:379:9: error: ‘n3’ not specified in enclosing ‘parallel’
  379 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:379:9: error: enclosing ‘parallel’
molDyn.c:383:9: error: ‘n3’ not specified in enclosing ‘parallel’
  383 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:383:9: error: enclosing ‘parallel’
molDyn.c:387:14: error: ‘ekin’ not specified in enclosing ‘parallel’
  387 |         ekin += vh[i] * vh[i];
      |              ^~
molDyn.c:383:9: error: enclosing ‘parallel’
  383 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:391:9: error: ‘n3’ not specified in enclosing ‘parallel’
  391 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:391:9: error: enclosing ‘parallel’
molDyn.c:395:9: error: ‘n3’ not specified in enclosing ‘parallel’
  395 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
molDyn.c:395:9: error: enclosing ‘parallel’
molDyn.c:399:14: error: ‘ekin’ not specified in enclosing ‘parallel’
  399 |         ekin += vh[i] * vh[i];
      |              ^~
molDyn.c:395:9: error: enclosing ‘parallel’
  395 | #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
      |         ^~~
*/
// #include <stdio.h>
// #include <math.h>
// #include <time.h>
// #include <omp.h>

// void dfill(int, double, double[], int);
// void domove(int, double[], double[], double[], double);
// void dscal(int, double, double[], int);
// void fcc(double[], int, int, double);
// void forces(int, double[], double[], double, double);
// double mkekin(int, double[], double[], double, double);
// void mxwell(double[], int, double, double);
// void prnout(int, double, double, double, double, double, double, int, double);
// double velavg(int, double[], double, double);
// double secnds(void);

// /*
//  * Variable declarations
//  */

// double epot;
// double vir;
// double count;

// /*
//  * function dfill: initializes double precision array to scalar value
//  */
// void dfill(int n, double val, double a[], int ia)
// {
//     int i;
//     for (i = 0; i < (n - 1) * ia + 1; i += ia)
//         a[i] = val;
// }

// /*
//  * Move particles
//  */
// void domove(int n3, double x[], double vh[], double f[], double side)
// {
//     int i;

// #pragma omp parallel for default(none) shared(x, vh, f, side, n3) private(i)
//     for (i = 0; i < n3; i++)
//     {
//         x[i] += vh[i] + f[i];
//         /*
//          *  Periodic boundary conditions
//          */
//         if (x[i] < 0.0)
//             x[i] += side;
//         if (x[i] > side)
//             x[i] -= side;
//         /*
//          *  Partial velocity updates
//          */
//         vh[i] += f[i];
//         /*
//          *  Initialize forces for the next iteration
//          */
//         f[i] = 0.0;
//     }
// }

// /*
//  * Scales an array
//  */
// void dscal(int n, double sa, double sx[], int incx)
// {
//     int i, j;

//     if (incx == 1)
//     {
// #pragma omp parallel for default(none) shared(n, sa, sx) private(i)
//         for (i = 0; i < n; i++)
//             sx[i] *= sa;
//     }
//     else
//     {
//         j = 0;
// #pragma omp parallel for default(none) shared(n, sa, sx, incx) private(i, j)
//         for (i = 0; i < n; i++)
//         {
//             sx[j] *= sa;
//             j += incx;
//         }
//     }
// }

// /*
//  * Generate fcc lattice for atoms inside the box
//  */
// void fcc(double x[], int npart, int mm, double a)
// {
//     int ijk = 0;
//     int i, j, k, lg;

// #pragma omp parallel for default(none) shared(x, mm, a, ijk) private(i, j, k, lg)
//     for (lg = 0; lg < 2; lg++)
//         for (i = 0; i < mm; i++)
//             for (j = 0; j < mm; j++)
//                 for (k = 0; k < mm; k++)
//                 {
//                     x[ijk] = i * a + lg * a * 0.5;
//                     x[ijk + 1] = j * a + lg * a * 0.5;
//                     x[ijk + 2] = k * a;
//                     ijk += 3;
//                 }

// #pragma omp parallel for default(none) shared(x, mm, a, ijk) private(i, j, k, lg)
//     for (lg = 1; lg < 3; lg++)
//         for (i = 0; i < mm; i++)
//             for (j = 0; j < mm; j++)
//                 for (k = 0; k < mm; k++)
//                 {
//                     x[ijk] = i * a + (2 - lg) * a * 0.5;
//                     x[ijk + 1] = j * a + (lg - 1) * a * 0.5;
//                     x[ijk + 2] = k * a + a * 0.5;
//                     ijk += 3;
//                 }
// }

// /*
//  * Compute forces and accumulate the virial and the potential
//  */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     int i, j;
//     double sideh, rcoffs;
//     double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
//     double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
//     double forcex, forcey, forcez;

//     vir = 0.0;
//     epot = 0.0;
//     sideh = 0.5 * side;
//     rcoffs = rcoff * rcoff;

// #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh, epot, vir) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         xi = x[i];
//         yi = x[i + 1];
//         zi = x[i + 2];
//         fxi = 0.0;
//         fyi = 0.0;
//         fzi = 0.0;

//         for (j = i + 3; j < npart * 3; j += 3)
//         {
//             xx = xi - x[j];
//             yy = yi - x[j + 1];
//             zz = zi - x[j + 2];
//             if (xx < -sideh)
//                 xx += side;
//             if (xx > sideh)
//                 xx -= side;
//             if (yy < -sideh)
//                 yy += side;
//             if (yy > sideh)
//                 yy -= side;
//             if (zz < -sideh)
//                 zz += side;
//             if (zz > sideh)
//                 zz -= side;
//             rd = xx * xx + yy * yy + zz * zz;

//             if (rd <= rcoffs)
//             {
//                 rrd = 1.0 / rd;
//                 rrd2 = rrd * rrd;
//                 rrd3 = rrd2 * rrd;
//                 rrd4 = rrd2 * rrd2;
//                 rrd6 = rrd2 * rrd4;
//                 rrd7 = rrd6 * rrd;
//                 epot += (rrd6 - rrd3);
//                 r148 = rrd7 - 0.5 * rrd4;
//                 vir -= rd * r148;
//                 forcex = xx * r148;
//                 fxi += forcex;
//                 f[j] -= forcex;
//                 forcey = yy * r148;
//                 fyi += forcey;
//                 f[j + 1] -= forcey;
//                 forcez = zz * r148;
//                 fzi += forcez;
//                 f[j + 2] -= forcez;
//             }
//         }
//         f[i] += fxi;
//         f[i + 1] += fyi;
//         f[i + 2] += fzi;
//     }
// }

// /*
//  * Main program: Molecular Dynamics simulation.
//  */
// int main()
// {
//     int mm = 15;
//     int npart = 4 * mm * mm * mm;
//     int move;
//     double x[npart * 3], vh[npart * 3], f[npart * 3];
//     double ekin;
//     double vel;
//     double sc;
//     double start, time;

//     /*
//      * Parameter definitions
//      */

//     double den = 0.83134;
//     double side = pow((double)npart / den, 0.3333333);
//     double tref = 0.722;
//     double rcoff = (double)mm / 4.0;
//     double h = 0.064;
//     int irep = 10;
//     int istop = 20;
//     int iprint = 5;
//     int movemx = 20;

//     double a = side / (double)mm;
//     double hsq = h * h;
//     double hsq2 = hsq * 0.5;
//     double tscale = 16.0 / ((double)npart - 1.0);
//     double vaver = 1.13 * sqrt(tref / 24.0);

//     /*
//      * Initial output
//      */

//     printf(" Molecular Dynamics Simulation example program\n");
//     printf(" ---------------------------------------------\n");
//     printf(" number of particles is ............ %6d\n", npart);
//     printf(" side length of the box is ......... %13.6f\n", side);
//     printf(" cut off is ........................ %13.6f\n", rcoff);
//     printf(" reduced temperature is ............ %13.6f\n", tref);
//     printf(" basic timestep is ................. %13.6f\n", h);
//     printf(" temperature scale interval ........ %6d\n", irep);
//     printf(" stop scaling at move .............. %6d\n", istop);
//     printf(" print interval .................... %6d\n", iprint);
//     printf(" total no. of steps ................ %6d\n", movemx);

//     /*
//      * Generate fcc lattice for atoms inside box
//      */
//     fcc(x, npart, mm, a);
//     /*
//      * Initialise velocities and forces (which are zero in fcc positions)
//      */
//     mxwell(vh, 3 * npart, h, tref);
//     dfill(3 * npart, 0.0, f, 1);
//     /*
//      * Start of md
//      */
//     printf("\n    i       ke         pe            e         temp   "
//            "   pres      vel      rp\n  -----  ----------  ----------"
//            "  ----------  --------  --------  --------  ----\n");

//     start = secnds();

//     for (move = 1; move <= movemx; move++)
//     {

//         /*
//          * Move the particles and partially update velocities
//          */
//         domove(3 * npart, x, vh, f, side);

//         /*
//          * Compute forces in the new positions and accumulate the virial
//          * and potential energy.
//          */
//         forces(npart, x, f, side, rcoff);

//         /*
//          * Scale forces, complete update of velocities and compute k.e.
//          */
//         ekin = mkekin(npart, f, vh, hsq2, hsq);

//         /*
//          * Average the velocity and temperature scale if desired
//          */
//         vel = velavg(npart, vh, vaver, h);
//         if (move < istop && fmod(move, irep) == 0)
//         {
//             sc = sqrt(tref / (tscale * ekin));
//             dscal(3 * npart, sc, vh, 1);
//             ekin = tref / tscale;
//         }

//         /*
//          * Sum to get full potential energy and virial
//          */
//         if (fmod(move, iprint) == 0)
//             prnout(move, ekin, epot, tscale, vir, vel, count, npart, den);
//     }

//     time = secnds() - start;

//     printf("Elapsed time =  %f\n", (float)time);
// }

// time_t starttime = 0;

// double secnds()
// {

//     return omp_get_wtime();
// }

// #include <stdio.h>
// /*
//  * Scale forces, update velocities and compute K.E.
//  */
// double mkekin(int npart, double f[], double vh[], double hsq2, double hsq)
// {
//     int i;
//     double sum = 0.0, ekin;

// #pragma omp parallel for default(none) shared(npart, f, vh, hsq2) reduction(+:sum) private(i)
//     for (i = 0; i < 3 * npart; i++)
//     {
//         f[i] *= hsq2;
//         vh[i] += f[i];
//         sum += vh[i] * vh[i];
//     }
//     ekin = sum / hsq;

//     return (ekin);
// }

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>

// void srand48(long);
// double drand48(void);
// /*
//  * Sample Maxwell distribution at temperature tref
//  */
// void mxwell(double vh[], int n3, double h, double tref)
// {
//     int i;
//     int npart = n3 / 3;
//     double r, tscale, v1, v2, s, ekin = 0.0, sp = 0.0, sc;

//     srand48(4711);
//     tscale = 16.0 / ((double)npart - 1.0);

// #pragma omp parallel for default(none) shared(vh, npart, tscale) private(i, s, v1, v2, r)
//     for (i = 0; i < n3; i += 2)
//     {
//         s = 2.0;
//         while (s >= 1.0)
//         {
//             v1 = 2.0 * drand48() - 1.0;
//             v2 = 2.0 * drand48() - 1.0;
//             s = v1 * v1 + v2 * v2;
//         }
//         r = sqrt(-2.0 * log(s) / s);
//         vh[i] = v1 * r;
//         vh[i + 1] = v2 * r;
//     }

// #pragma omp parallel for default(none) shared(vh, npart, h) private(i)
//     for (i = 0; i < n3; i += 3)
//         sp += vh[i];
//     sp /= (double)npart;
// #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
//     for (i = 0; i < n3; i += 3)
//     {
//         vh[i] -= sp;
//         ekin += vh[i] * vh[i];
//     }

//     sp = 0.0;
// #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
//     for (i = 1; i < n3; i += 3)
//         sp += vh[i];
//     sp /= (double)npart;
// #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
//     for (i = 1; i < n3; i += 3)
//     {
//         vh[i] -= sp;
//         ekin += vh[i] * vh[i];
//     }

//     sp = 0.0;
// #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
//     for (i = 2; i < n3; i += 3)
//         sp += vh[i];
//     sp /= (double)npart;
// #pragma omp parallel for default(none) shared(vh, npart, sp) private(i)
//     for (i = 2; i < n3; i += 3)
//     {
//         vh[i] -= sp;
//         ekin += vh[i] * vh[i];
//     }

//     sc = h * sqrt(tref / (tscale * ekin));
// #pragma omp parallel for default(none) shared(vh, n3, sc) private(i)
//     for (i = 0; i < n3; i++)
//         vh[i] *= sc;
// }

// /*
//  * Print out interesting information at current timestep
//  */
// void prnout(int move, double ekin, double epot, double tscale, double vir,
//             double vel, double count, int npart, double den)
// {
//     double ek, etot, temp, pres, rp;

//     ek = 24.0 * ekin;
//     epot *= 4.0;
//     etot = ek + epot;
//     temp = tscale * ekin;
//     pres = den * 16.0 * (ekin - vir) / (double)npart;
//     vel /= (double)npart;
//     rp = (count / (double)npart) * 100.0;
//     printf(" %6d%12.4f%12.4f%12.4f%10.4f%10.4f%10.4f%6.1f\n",
//            move, ek, epot, etot, temp, pres, vel, rp);
// }

// /*
//  * Compute average velocity
//  */
// double velavg(int npart, double vh[], double vaver, double h)
// {
//     int i;
//     double vaverh = vaver * h;
//     double vel = 0.0;
//     double sq;
//     extern double count;

//     count = 0.0;
// #pragma omp parallel for default(none) shared(npart, vh, vaverh, h, count) private(i, sq) reduction(+:vel)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         sq = sqrt(vh[i] * vh[i] + vh[i + 1] * vh[i + 1] + vh[i + 2] * vh[i + 2]);
//         if (sq > vaverh)
//             count++;
//         vel += sq;
//     }
//     vel /= h;

//     return (vel);
// }

// #include <stdio.h>
// #include <math.h>
// #include <time.h>
// #include <omp.h>

// void dfill(int, double, double[], int);
// void domove(int, double[], double[], double[], double);
// void dscal(int, double, double[], int);
// void fcc(double[], int, int, double);
// void forces(int, double[], double[], double, double);
// double mkekin(int, double[], double[], double, double);
// void mxwell(double[], int, double, double);
// void prnout(int, double, double, double, double, double, double, int, double);
// double velavg(int, double[], double, double);
// double secnds(void);

// /*
//  * Variable declarations
//  */

// double epot;
// double vir;
// double count;

// /*
//  * function dfill: initializes double precision array to scalar value
//  */
// void dfill(int n, double val, double a[], int ia)
// {
//     int i;
//     for (i = 0; i < (n - 1) * ia + 1; i += ia)
//         a[i] = val;
// }

// /*
//  * Move particles
//  */
// void domove(int n3, double x[], double vh[], double f[], double side)
// {
//     int i;

// #pragma omp parallel for default(none) shared(x, vh, f, side, n3) private(i)
//     for (i = 0; i < n3; i++)
//     {
//         x[i] += vh[i] + f[i];
//         /*
//          *  Periodic boundary conditions
//          */
//         if (x[i] < 0.0)
//             x[i] += side;
//         if (x[i] > side)
//             x[i] -= side;
//         /*
//          *  Partial velocity updates
//          */
//         vh[i] += f[i];
//         /*
//          *  Initialize forces for the next iteration
//          */
//         f[i] = 0.0;
//     }
// }

// /*
//  * Scales an array
//  */
// void dscal(int n, double sa, double sx[], int incx)
// {
//     int i, j;

//     if (incx == 1)
//     {
// #pragma omp parallel for default(none) shared(n, sa, sx) private(i)
//         for (i = 0; i < n; i++)
//             sx[i] *= sa;
//     }
//     else
//     {
// #pragma omp parallel for default(none) shared(n, sa, sx, incx) private(i, j)
//         for (i = 0; i < n; i += incx)
//         {
//             j = i;
//             sx[j] *= sa;
//         }
//     }
// }

// /*
//  * Generate fcc lattice for atoms inside the box
//  */
// // void fcc(double x[], int npart, int mm, double a)
// // {
// //     int ijk = 0;
// //     int i, j, k, lg;

// // #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg)
// //     for (lg = 0; lg < 2; lg++)
// //         for (i = 0; i < mm; i++)
// //             for (j = 0; j < mm; j++)
// //                 for (k = 0; k < mm; k++)
// //                 {
// //                     x[ijk] = i * a + lg * a * 0.5;
// //                     x[ijk + 1] = j * a + lg * a * 0.5;
// //                     x[ijk + 2] = k * a;
// //                     ijk += 3;
// //                 }

// // #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg)
// //     for (lg = 1; lg < 3; lg++)
// //         for (i = 0; i < mm; i++)
// //             for (j = 0; j < mm; j++)
// //                 for (k = 0; k < mm; k++)
// //                 {
// //                     x[ijk] = i * a + (2 - lg) * a * 0.5;
// //                     x[ijk + 1] = j * a + (lg - 1) * a * 0.5;
// //                     x[ijk + 2] = k * a + a * 0.5;
// //                     ijk += 3;
// //                 }
// // }

// void fcc(double x[], int npart, int mm, double a)
// {
//     int i, j, k, lg, ijk;

// #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg, ijk)
//     for (i = 0; i < mm; i++)
//     {
//         for (j = 0; j < mm; j++)
//         {
//             for (k = 0; k < mm; k++)
//             {
//                 lg = i + j + k;
//                 ijk = 3 * (i * mm * mm + j * mm + k);
//                 if (lg % 2 == 0)
//                 {
//                     x[ijk] = i * a;
//                     x[ijk + 1] = j * a;
//                     x[ijk + 2] = k * a;
//                 }
//                 else
//                 {
//                     x[ijk] = i * a + lg * a * 0.5;
//                     x[ijk + 1] = j * a + (2 - lg) * a * 0.5;
//                     x[ijk + 2] = k * a + (2 - lg) * a * 0.5;
//                 }
//             }
//         }
//     }
// }

// /*
//  * Compute forces and accumulate the virial and the potential
//  */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     int i, j;
//     double sideh, rcoffs;
//     double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
//     double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
//     double forcex, forcey, forcez;

//     vir = 0.0;
//     epot = 0.0;
//     sideh = 0.5 * side;
//     rcoffs = rcoff * rcoff;

// #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez) reduction(+:epot, vir)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         xi = x[i];
//         yi = x[i + 1];
//         zi = x[i + 2];
//         fxi = 0.0;
//         fyi = 0.0;
//         fzi = 0.0;

//         for (j = i + 3; j < npart * 3; j += 3)
//         {
//             xx = xi - x[j];
//             yy = yi - x[j + 1];
//             zz = zi - x[j + 2];
//             if (xx < -sideh)
//                 xx += side;
//             if (xx > sideh)
//                 xx -= side;
//             if (yy < -sideh)
//                 yy += side;
//             if (yy > sideh)
//                 yy -= side;
//             if (zz < -sideh)
//                 zz += side;
//             if (zz > sideh)
//                 zz -= side;
//             rd = xx * xx + yy * yy + zz * zz;

//             if (rd <= rcoffs)
//             {
//                 rrd = 1.0 / rd;
//                 rrd2 = rrd * rrd;
//                 rrd3 = rrd2 * rrd;
//                 rrd4 = rrd2 * rrd2;
//                 rrd6 = rrd2 * rrd4;
//                 rrd7 = rrd6 * rrd;
//                 epot += (rrd6 - rrd3);
//                 r148 = rrd7 - 0.5 * rrd4;
//                 vir -= rd * r148;
//                 forcex = xx * r148;
//                 fxi += forcex;
//                 f[j] -= forcex;
//                 forcey = yy * r148;
//                 fyi += forcey;
//                 f[j + 1] -= forcey;
//                 forcez = zz * r148;
//                 fzi += forcez;
//                 f[j + 2] -= forcez;
//             }
//         }
//         f[i] += fxi;
//         f[i + 1] += fyi;
//         f[i + 2] += fzi;
//     }
// }

// /*
//  * Main program: Molecular Dynamics simulation.
//  */
// int main()
// {
//     int mm = 15;
//     int npart = 4 * mm * mm * mm;
//     int move;
//     double x[npart * 3], vh[npart * 3], f[npart * 3];
//     double ekin;
//     double vel;
//     double sc;
//     double start, time;

//     /*
//      * Parameter definitions
//      */

//     double den = 0.83134;
//     double side = pow((double)npart / den, 0.3333333);
//     double tref = 0.722;
//     double rcoff = (double)mm / 4.0;
//     double h = 0.064;
//     int irep = 10;
//     int istop = 20;
//     int iprint = 5;
//     int movemx = 20;

//     double a = side / (double)mm;
//     double hsq = h * h;
//     double hsq2 = hsq * 0.5;
//     double tscale = 16.0 / ((double)npart - 1.0);
//     double vaver = 1.13 * sqrt(tref / 24.0);

//     /*
//      * Initial output
//      */

//     printf(" Molecular Dynamics Simulation example program\n");
//     printf(" ---------------------------------------------\n");
//     printf(" number of particles is ............ %6d\n", npart);
//     printf(" side length of the box is ......... %13.6f\n", side);
//     printf(" cut off is ........................ %13.6f\n", rcoff);
//     printf(" reduced temperature is ............ %13.6f\n", tref);
//     printf(" basic timestep is ................. %13.6f\n", h);
//     printf(" temperature scale interval ........ %6d\n", irep);
//     printf(" stop scaling at move .............. %6d\n", istop);
//     printf(" print interval .................... %6d\n", iprint);
//     printf(" total no. of steps ................ %6d\n", movemx);

//     /*
//      * Generate fcc lattice for atoms inside box
//      */
//     fcc(x, npart, mm, a);
//     /*
//      * Initialise velocities and forces (which are zero in fcc positions)
//      */
//     mxwell(vh, 3 * npart, h, tref);
//     dfill(3 * npart, 0.0, f, 1);
//     /*
//      * Start of md
//      */
//     printf("\n    i       ke         pe            e         temp   "
//            "   pres      vel      rp\n  -----  ----------  ----------"
//            "  ----------  --------  --------  --------  ----\n");

//     start = secnds();

//     for (move = 1; move <= movemx; move++)
//     {

//         /*
//          * Move the particles and partially update velocities
//          */
//         domove(3 * npart, x, vh, f, side);

//         /*
//          * Compute forces in the new positions and accumulate the virial
//          * and potential energy.
//          */
//         forces(npart, x, f, side, rcoff);

//         /*
//          * Scale forces, complete update of velocities and compute k.e.
//          */
//         ekin = mkekin(npart, f, vh, hsq2, hsq);

//         /*
//          * Average the velocity and temperature scale if desired
//          */
//         vel = velavg(npart, vh, vaver, h);
//         if (move < istop && fmod(move, irep) == 0)
//         {
//             sc = sqrt(tref / (tscale * ekin));
//             dscal(3 * npart, sc, vh, 1);
//             ekin = tref / tscale;
//         }

//         /*
//          * Sum to get full potential energy and virial
//          */
//         if (fmod(move, iprint) == 0)
//             prnout(move, ekin, epot, tscale, vir, vel, count, npart, den);
//     }

//     time = secnds() - start;

//     printf("Elapsed time =  %f\n", (float)time);
// }

// double secnds()
// {

//     return omp_get_wtime();
// }

// #include <stdlib.h>
// #include <math.h>

// /*
//  * Scale forces, update velocities and compute K.E.
//  */
// // double mkekin(int npart, double f[], double vh[], double hsq2, double hsq)
// // {
// //     int i;
// //     double sum = 0.0, ekin;

// // #pragma omp parallel for default(none) shared(npart, f, vh) private(i) reduction(+ \
// //                                                                                  : sum)
// //     for (i = 0; i < npart * 3; i++)
// //     {
// //         vh[i] += hsq2 * f[i];
// //         sum += vh[i] * vh[i];
// //     }
// //     ekin = sum * hsq;
// //     return ekin;
// // }

// double mkekin(int npart, double f[], double vh[], double hsq2, double hsq)
// {
//     int i;
//     double sum = 0.0, ekin;

// #pragma omp parallel for default(none) shared(npart, f, vh, hsq2) private(i) reduction(+ \
//                                                                                        : sum)
//     for (i = 0; i < npart * 3; i++)
//     {
//         vh[i] += hsq2 * f[i];
//         sum += vh[i] * vh[i];
//     }
//     ekin = sum * hsq;
//     return ekin;
// }

// /*
//  * Generate a Maxwellian distribution of velocities
//  */
// void mxwell(double vh[], int npart, double h, double tref)
// {
//     int i;
//     double ts, sp, v1, v2, r;

//     sp = 0.0;

//     for (i = 0; i < npart; i += 2)
//     {
//         do
//         {
//             v1 = 2.0 * rand() / RAND_MAX - 1.0;
//             v2 = 2.0 * rand() / RAND_MAX - 1.0;
//             r = v1 * v1 + v2 * v2;
//         } while (r >= 1.0);

//         ts = sqrt(-2.0 * log(r) / r);
//         sp += v1;
//         vh[i] = v1 * ts * sqrt(tref);
//         if ((i + 1) < npart)
//         {
//             sp += v2;
//             vh[i + 1] = v2 * ts * sqrt(tref);
//         }
//     }

//     /*
//      * Check average speed
//      */
//     sp /= (double)npart;
//     count = 0.0;

//     for (i = 0; i < npart; i++)
//         vh[i] -= sp;
// }

// /*
//  * Print out the state of the simulation
//  */
// void prnout(int move, double ekin, double epot, double tscale, double vir, double vel, double count, int npart, double den)
// {
//     double pres, rp;

//     pres = den * (ekin - 0.5 * vir) / npart;
//     rp = count / (double)(npart * 3);
//     printf("%6d %10.5f %12.5f %12.5f %9.5f %9.5f %9.5f %5.3f\n",
//            move, ekin, epot, (ekin + epot) / npart, tscale * ekin, pres, vel, rp);
// }

// /*
//  * Calculate average velocity
//  */
// double velavg(int npart, double vh[], double vaver, double h)
// {
//     int i;
//     double vel = 0.0;

//     for (i = 0; i < npart; i++)
//         vel += vh[i] * vh[i];

//     vel = sqrt(vel * h);
//     count += fabs(vel - vaver);
//     return vel;
// }

// #include <stdio.h>
// #include <math.h>
// #include <time.h>
// #include <omp.h>

// void dfill(int, double, double[], int);
// void domove(int, double[], double[], double[], double);
// void dscal(int, double, double[], int);
// void fcc(double[], int, int, double);
// void forces(int, double[], double[], double, double);
// double mkekin(int, double[], double[], double, double);
// void mxwell(double[], int, double, double);
// void prnout(int, double, double, double, double, double, double, int, double);
// double velavg(int, double[], double, double);
// double secnds(void);

// /*
//  * Variable declarations
//  */

// double epot;
// double vir;
// double count;

// omp_lock_t lock_epot;
// omp_lock_t lock_vir;
// omp_lock_t lock_count;

// /*
//  * function dfill: initializes double precision array to scalar value
//  */
// void dfill(int n, double val, double a[], int ia)
// {
//     int i;
//     for (i = 0; i < (n - 1) * ia + 1; i += ia)
//         a[i] = val;
// }

// /*
//  * Move particles
//  */
// void domove(int n3, double x[], double vh[], double f[], double side)
// {
//     int i;

// #pragma omp parallel for default(none) shared(x, vh, f, side, n3) private(i)
//     for (i = 0; i < n3; i++)
//     {
//         x[i] += vh[i] + f[i];
//         /*
//          *  Periodic boundary conditions
//          */
//         if (x[i] < 0.0)
//             x[i] += side;
//         if (x[i] > side)
//             x[i] -= side;
//         /*
//          *  Partial velocity updates
//          */
//         vh[i] += f[i];
//         /*
//          *  Initialize forces for the next iteration
//          */
//         f[i] = 0.0;
//     }
// }

// /*
//  * Scales an array
//  */
// void dscal(int n, double sa, double sx[], int incx)
// {
//     int i, j;

//     if (incx == 1)
//     {
// #pragma omp parallel for default(none) shared(n, sa, sx) private(i)
//         for (i = 0; i < n; i++)
//             sx[i] *= sa;
//     }
//     else
//     {
//         j = 0;
// #pragma omp parallel for default(none) shared(n, sa, sx, incx) private(i)
//         for (i = 0; i < n; i++)
//         {
//             sx[j] *= sa;
//             j += incx;
//         }
//     }
// }

// /*
//  * Generate fcc lattice for atoms inside the box
//  */
// void fcc(double x[], int npart, int mm, double a)
// {
//     int ijk = 0;
//     int i, j, k, lg;

// #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg, ijk)
//     for (lg = 0; lg < 2; lg++)
//         for (i = 0; i < mm; i++)
//             for (j = 0; j < mm; j++)
//                 for (k = 0; k < mm; k++)
//                 {
//                     ijk = 3 * (i * mm * mm + j * mm + k) + lg * 3 * mm;
//                     x[ijk] = i * a + lg * a * 0.5;
//                     x[ijk + 1] = j * a + lg * a * 0.5;
//                     x[ijk + 2] = k * a;
//                 }

// #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg, ijk)
//     for (lg = 1; lg < 3; lg++)
//         for (i = 0; i < mm; i++)
//             for (j = 0; j < mm; j++)
//                 for (k = 0; k < mm; k++)
//                 {
//                     ijk = 3 * (i * mm * mm + j * mm + k) + lg * 3 * mm;
//                     x[ijk] = i * a + (2 - lg) * a * 0.5;
//                     x[ijk + 1] = j * a + (lg - 1) * a * 0.5;
//                     x[ijk + 2] = k * a + a * 0.5;
//                 }
// }

// /*
//  * Compute forces and accumulate the virial and the potential
//  */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     int i, j;
//     double sideh, rcoffs;
//     double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
//     double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
//     double forcex, forcey, forcez;

//     vir = 0.0;
//     epot = 0.0;
//     sideh = 0.5 * side;
//     rcoffs = rcoff * rcoff;

// #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez) reduction(+:epot, vir)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         xi = x[i];
//         yi = x[i + 1];
//         zi = x[i + 2];
//         fxi = 0.0;
//         fyi = 0.0;
//         fzi = 0.0;

//         for (j = i + 3; j < npart * 3; j += 3)
//         {
//             xx = xi - x[j];
//             yy = yi - x[j + 1];
//             zz = zi - x[j + 2];
//             if (xx < -sideh)
//                 xx += side;
//             if (xx > sideh)
//                 xx -= side;
//             if (yy < -sideh)
//                 yy += side;
//             if (yy > sideh)
//                 yy -= side;
//             if (zz < -sideh)
//                 zz += side;
//             if (zz > sideh)
//                 zz -= side;
//             rd = xx * xx + yy * yy + zz * zz;

//             if (rd <= rcoffs)
//             {
//                 rrd = 1.0 / rd;
//                 rrd2 = rrd * rrd;
//                 rrd3 = rrd2 * rrd;
//                 rrd4 = rrd2 * rrd2;
//                 rrd6 = rrd2 * rrd4;
//                 rrd7 = rrd6 * rrd;
//                 r148 = rrd7 - 0.5 * rrd6;
//                 epot += rrd6 - rrd3;
//                 forcex = xx * r148;
//                 forcey = yy * r148;
//                 forcez = zz * r148;
//                 fxi += forcex;
//                 fyi += forcey;
//                 fzi += forcez;
//                 f[j] -= forcex;
//                 f[j + 1] -= forcey;
//                 f[j + 2] -= forcez;
//                 vir += rd * r148;
//             }
//         }
//         f[i] += fxi;
//         f[i + 1] += fyi;
//         f[i + 2] += fzi;
//     }
// }

// /*
//  * function mkekin: computes kinetic energy
//  */
// double mkekin(int npart, double vh[], double hsq2, double sp)
// {
//     int i;
//     double ekin = 0.0;

// #pragma omp parallel for default(none) shared(npart, vh, hsq2, sp) private(i) reduction(+ \
//                                                                                   : ekin)
//     for (i = 0; i < npart * 3; i++)
//     {
//         vh[i] += hsq2 * sp;
//         ekin += vh[i] * vh[i];
//     }
//     return ekin;
// }

// /*
//  * function mxwell: compute maxwellian velocity distribution
//  */
// void mxwell(double vh[], int npart, double h, double tscale)
// {
//     int i;

// #pragma omp parallel for default(none) shared(vh, npart, h, tscale) private(i)
//     for (i = 0; i < npart * 3; i++)
//         vh[i] = h * vh[i] * sqrt(tscale);
// }

// /*
//  * function prnout: prints out results
//  */
// void prnout(int natoms, double time, double ekin, double sp, double t, double pres, double rp, int nblok, double r)
// {
//     printf("%8d %14.6lf %14.6lf %14.6lf %14.6lf %14.6lf %14.6lf %8d %14.6lf\n",
//            natoms, time, ekin, sp / (3.0 * natoms), t, pres, rp, nblok, r);
// }

// /*
//  * function velavg: computes average velocities
//  */
// double velavg(int npart, double vh[], double h, double vel)
// {
//     int i;
//     double vavg = 0.0;

// #pragma omp parallel for default(none) shared(npart, vh, h) private(i) reduction(+ \
//                                                                               : vavg)
//     for (i = 0; i < npart * 3; i++)
//         vavg += vh[i] * vh[i];

//     vavg = sqrt(vavg / (3 * npart));
//     vel += h * (2.0 * count - 3.0) / (2.0 * count);
//     return vavg;
// }

// /*
//  * function secnds: returns elapsed time
//  */
// double secnds(void)
// {
//     struct timespec tp;
//     double etime;

//     clock_gettime(CLOCK_MONOTONIC, &tp);
//     etime = (double)tp.tv_sec + (double)tp.tv_nsec * 1.0e-9;
//     return etime;
// }

// int main(int argc, char *argv[])
// {
//     double side, rcoff, rcoffs, hsq, hsq2, tscale, vaver, sp;
//     double ekin, t, pres, rp, time, stime, velp;
//     double *x, *v, *f;
//     int natoms, nsteps, npart, nblok, n, mm, i;
//     int niter = 100;

//     natoms = 108;
//     side = 16.8645;
//     nsteps = 10;
//     rcoff = side / 2.0;
//     hsq = 0.03189;
//     hsq2 = hsq * 0.5;
//     tscale = 0.9;
//     nblok = 10;
//     rcoffs = rcoff * rcoff;

//     npart = 4 * pow((natoms / 4), 1.0 / 3.0);
//     x = (double *)malloc(3 * natoms * sizeof(double));
//     v = (double *)malloc(3 * natoms * sizeof(double));
//     f = (double *)malloc(3 * natoms * sizeof(double));

//     omp_init_lock(&lock_epot);
//     omp_init_lock(&lock_vir);
//     omp_init_lock(&lock_count);

//     dfill(natoms * 3, 0.0, f, 1);
//     fcc(x, natoms, npart, side);
//     mxwell(v, natoms, hsq, tscale);

//     time = secnds();

//     for (n = 1; n <= nsteps; n++)
//     {
//         forces(natoms, x, f, side, rcoff);
//         ekin = mkekin(natoms, v, hsq2, 1.0);
//         t = tscale * ekin;
//         pres = natoms * t / side + vir / (3.0 * side);
//         sp = 0.0;
//         for (i = 0; i < natoms * 3; i++)
//             sp += v[i] * v[i];
//         sp /= natoms;
//         sp = (sp - 1.0) * tscale;
//         stime = secnds() - time;

//         prnout(natoms, stime, ekin, sp, t, pres, rp, nblok, r);

//         if (n % nblok == 0)
//         {
//             velp = velavg(natoms, v, hsq, 0.0);
//             vaver += velp;
//             printf("Average velocity is: %f\n", velp);
//             count++;
//         }
//         domove(natoms * 3, x, v, f, side);
//     }

//     vaver /= niter;
//     printf("Average velocity over all steps is: %f\n", vaver);

//     free(x);
//     free(v);
//     free(f);

//     return 0;
// }

// #include <stdio.h>
// #include <math.h>
// #include <time.h>
// #include <stdlib.h> // Added for malloc and free
// #include <omp.h>

// void dfill(int, double, double[], int);
// void domove(int, double[], double[], double[], double);
// void dscal(int, double, double[], int);
// void fcc(double[], int, int, double);
// void forces(int, double[], double[], double, double);
// double mkekin(int npart, double vh[], double hsq2, double sp);
// void mxwell(double[], int, double, double);
// void prnout(int, double, double, double, double, double, int);
// double velavg(int, double[], double);
// double secnds(void);

// /*
//  * Variable declarations
//  */

// double epot;
// double vir;
// double count;

// omp_lock_t lock_epot;
// omp_lock_t lock_vir;
// omp_lock_t lock_count;

// /*
//  * function dfill: initializes double precision array to scalar value
//  */
// void dfill(int n, double val, double a[], int ia)
// {
//     int i;
//     for (i = 0; i < (n - 1) * ia + 1; i += ia)
//         a[i] = val;
// }

// /*
//  * Move particles
//  */
// void domove(int n3, double x[], double vh[], double f[], double side)
// {
//     int i;

// #pragma omp parallel for default(none) shared(x, vh, f, side, n3) private(i)
//     for (i = 0; i < n3; i++)
//     {
//         x[i] += vh[i] + f[i];
//         /*
//          *  Periodic boundary conditions
//          */
//         if (x[i] < 0.0)
//             x[i] += side;
//         if (x[i] > side)
//             x[i] -= side;
//         /*
//          *  Partial velocity updates
//          */
//         vh[i] += f[i];
//         /*
//          *  Initialize forces for the next iteration
//          */
//         f[i] = 0.0;
//     }
// }

// /*
//  * Scales an array
//  */
// void dscal(int n, double sa, double sx[], int incx)
// {
//     int i, j;

//     if (incx == 1)
//     {
// #pragma omp parallel for default(none) shared(n, sa, sx) private(i)
//         for (i = 0; i < n; i++)
//             sx[i] *= sa;
//     }
//     else
//     {
//         j = 0;
// #pragma omp parallel for default(none) shared(n, sa, sx, incx, j) private(i)
//         for (i = 0; i < n; i++)
//         {
//             sx[j] *= sa;
//             j += incx;
//         }
//     }
// }

// /*
//  * Generate fcc lattice for atoms inside the box
//  */
// void fcc(double x[], int npart, int mm, double a)
// {
//     int ijk = 0;
//     int i, j, k, lg;

// #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg, ijk)
//     for (lg = 0; lg < 2; lg++)
//         for (i = 0; i < mm; i++)
//             for (j = 0; j < mm; j++)
//                 for (k = 0; k < mm; k++)
//                 {
//                     ijk = 3 * (i * mm * mm + j * mm + k) + lg * 3 * mm;
//                     x[ijk] = i * a + lg * a * 0.5;
//                     x[ijk + 1] = j * a + lg * a * 0.5;
//                     x[ijk + 2] = k * a;
//                 }

// #pragma omp parallel for default(none) shared(x, mm, a) private(i, j, k, lg, ijk)
//     for (lg = 1; lg < 3; lg++)
//         for (i = 0; i < mm; i++)
//             for (j = 0; j < mm; j++)
//                 for (k = 0; k < mm; k++)
//                 {
//                     ijk = 3 * (i * mm * mm + j * mm + k) + lg * 3 * mm;
//                     x[ijk] = i * a + (2 - lg) * a * 0.5;
//                     x[ijk + 1] = j * a + (lg - 1) * a * 0.5;
//                     x[ijk + 2] = k * a + a * 0.5;
//                 }
// }

// /*
//  * Compute forces and accumulate the virial and the potential
//  */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     int i, j;
//     double sideh, rcoffs;
//     double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
//     double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
//     double forcex, forcey, forcez;

//     vir = 0.0;
//     epot = 0.0;
//     sideh = 0.5 * side;
//     rcoffs = rcoff * rcoff;

// #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez) reduction(+:epot, vir)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         xi = x[i];
//         yi = x[i + 1];
//         zi = x[i + 2];
//         fxi = 0.0;
//         fyi = 0.0;
//         fzi = 0.0;

//         for (j = i + 3; j < npart * 3; j += 3)
//         {
//             xx = xi - x[j];
//             yy = yi - x[j + 1];
//             zz = zi - x[j + 2];
//             /*
//              *  minimum image convention
//              */
//             if (xx > sideh)
//                 xx -= side;
//             if (yy > sideh)
//                 yy -= side;
//             if (zz > sideh)
//                 zz -= side;
//             if (xx < -sideh)
//                 xx += side;
//             if (yy < -sideh)
//                 yy += side;
//             if (zz < -sideh)
//                 zz += side;
//             rd = xx * xx + yy * yy + zz * zz;
//             if (rd <= rcoffs)
//             {
//                 rrd = 1.0 / rd;
//                 rrd2 = rrd * rrd;
//                 rrd3 = rrd2 * rrd;
//                 rrd4 = rrd2 * rrd2;
//                 rrd6 = rrd2 * rrd4;
//                 rrd7 = rrd6 * rrd;
//                 r148 = rrd7 - 0.5 * rrd6;
//                 epot += rrd6 - rrd3;
//                 forcex = xx * r148;
//                 forcey = yy * r148;
//                 forcez = zz * r148;
//                 fxi += forcex;
//                 fyi += forcey;
//                 fzi += forcez;
//                 f[j] -= forcex;
//                 f[j + 1] -= forcey;
//                 f[j + 2] -= forcez;
//                 vir += rd * r148;
//             }
//         }
//         f[i] += fxi;
//         f[i + 1] += fyi;
//         f[i + 2] += fzi;
//     }
// }

// /*
//  * function mkekin: computes kinetic energy
//  */
// double mkekin(int npart, double vh[], double hsq2, double sp)
// {
//     int i;
//     double ekin = 0.0;

// #pragma omp parallel for default(none) shared(npart, vh, hsq2, sp) private(i) reduction(+ \
//                                                                                   : ekin)
//     for (i = 0; i < npart * 3; i++)
//     {
//         vh[i] += hsq2 * sp;
//         ekin += vh[i] * vh[i];
//     }
//     return ekin;
// }

// /*
//  * function mxwell: compute maxwellian velocity distribution
//  */
// void mxwell(double vh[], int npart, double h, double tscale)
// {
//     int i;

// #pragma omp parallel for default(none) shared(vh, npart, h, tscale) private(i)
//     for (i = 0; i < npart * 3; i++)
//         vh[i] = h * vh[i] * sqrt(tscale);
// }

// /*
//  * function prnout: prints out results
//  */
// void prnout(int natoms, double time, double ekin, double sp, double t, double pres, int nblok)
// {
//     printf("%8d %14.6lf %14.6lf %14.6lf %14.6lf %14.6lf %8d\n",
//            natoms, time, ekin, sp / (3.0 * natoms), t, pres, nblok);
// }

// /*
//  * function velavg: computes average velocities
//  */
// double velavg(int npart, double vh[], double h)
// {
//     int i;
//     double vavg = 0.0;

// #pragma omp parallel for default(none) shared(npart, vh, h) private(i) reduction(+ \
//                                                                               : vavg)
//     for (i = 0; i < npart * 3; i++)
//         vavg += vh[i] * vh[i];

//     vavg = sqrt(vavg / (3 * npart));
//     return vavg;
// }

// /*
//  * function secnds: returns elapsed time
//  */
// double secnds(void)
// {
//     struct timespec tp;
//     double etime;

//     clock_gettime(CLOCK_MONOTONIC, &tp);
//     etime = (double)tp.tv_sec + (double)tp.tv_nsec * 1.0e-9;
//     return etime;
// }

// int main(int argc, char *argv[])
// {
//     double side, rcoff, hsq, hsq2, tscale, vaver, sp;
//     double ekin, t, pres, time, stime, velp;
//     double *x, *v, *f;
//     int natoms, nsteps, npart, nblok, n;

//     natoms = 108;
//     side = 16.8645;
//     nsteps = 10;
//     rcoff = side / 2.0;
//     hsq = 0.03189;
//     hsq2 = hsq * 0.5;
//     tscale = 0.9;
//     nblok = 10;

//     npart = 4 * pow((natoms / 4), 1.0 / 3.0);
//     x = (double *)malloc(3 * natoms * sizeof(double));
//     v = (double *)malloc(3 * natoms * sizeof(double));
//     f = (double *)malloc(3 * natoms * sizeof(double));

//     omp_init_lock(&lock_epot);
//     omp_init_lock(&lock_vir);
//     omp_init_lock(&lock_count);

//     dfill(natoms * 3, 0.0, f, 1);
//     fcc(x, natoms, npart, side);
//     mxwell(v, natoms, hsq, tscale);

//     time = secnds();

//     for (n = 1; n <= nsteps; n++)
//     {
//         forces(natoms, x, f, side, rcoff);
//         ekin = mkekin(natoms, v, hsq2, 1.0);
//         t = tscale * ekin;
//         pres = natoms * t / side + vir / (3.0 * side);
//         sp = 0.0;
//         for (int i = 0; i < natoms * 3; i++)
//             sp += v[i] * v[i];
//         sp /= natoms;
//         sp = (sp - 1.0) * tscale;
//         stime = secnds() - time;

//         prnout(natoms, stime, ekin, sp, t, pres, nblok);

//         if (n % nblok == 0)
//         {
//             velp = velavg(natoms, v, hsq);
//             vaver += velp;
//             printf("Average velocity is: %f\n", velp);
//         }
//         domove(natoms * 3, x, v, f, side);
//     }

//     vaver /= nsteps;
//     printf("Average velocity over all steps is: %f\n", vaver);

//     free(x);
//     free(v);
//     free(f);

//     return 0;
// }

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

/*
 *  Function declarations
 */

void dfill(int, double, double[], int);

void domove(int, double[], double[], double[], double);

void dscal(int, double, double[], int);

void fcc(double[], int, int, double);

void forces(int, double[], double[], double, double);

double
mkekin(int, double[], double[], double, double);

void mxwell(double[], int, double, double);

void prnout(int, double, double, double, double, double, double, int, double);

double
velavg(int, double[], double, double);

double
secnds(void);

/*
 *  Variable declarations
 */

double epot;
double vir;
double count;

/*
 *  function dfill : intialises double precision array to scalar value
 */
  void
  dfill(int n, double val, double a[], int ia){
    int i;

    for (i=0; i<(n-1)*ia+1; i+=ia)
      a[i] = val;
  }

/*
 *  Move particles
 */
  void
  domove(int n3, double x[], double vh[], double f[], double side){
    int i;

    for (i=0; i<n3; i++) {
      x[i] += vh[i]+f[i];
  /*
   *  Periodic boundary conditions
   */
      if (x[i] < 0.0)  x[i] += side;
      if (x[i] > side) x[i] -= side;
  /*
   *  Partial velocity updates
   */
      vh[i] += f[i];
  /*
   *  Initialise forces for the next iteration
   */
      f[i] = 0.0;
    }
  }

/*
 *  Scales an array
 */
  void
  dscal(int n,double sa,double sx[], int incx){
    int i,j;

    if (incx == 1) {
      for (i=0; i<n; i++)
        sx[i] *= sa;
    } else {
      j = 0;
      for (i=0; i<n; i++) {
        sx[j] *= sa;
        j += incx;
      }
    }
  }

/*
 *   Generate fcc lattice for atoms inside the box
 */
  void
  fcc(double x[], int npart, int mm, double a){
    int ijk=0;
    int i,j,k,lg;

    for (lg=0; lg<2; lg++)
      for (i=0; i<mm; i++)
        for (j=0; j<mm; j++)
          for (k=0; k<mm; k++) {
            x[ijk]   = i*a+lg*a*0.5;
            x[ijk+1] = j*a+lg*a*0.5;
            x[ijk+2] = k*a;
            ijk += 3;
          }

    for (lg=1; lg<3; lg++)
      for (i=0; i<mm; i++)
        for (j=0; j<mm; j++)
          for (k=0; k<mm; k++) {
            x[ijk]   = i*a+(2-lg)*a*0.5;
            x[ijk+1] = j*a+(lg-1)*a*0.5;
            x[ijk+2] = k*a+a*0.5;
            ijk += 3;
          }

  }

/*
 *  Compute forces and accumulate the virial and the potential
 */
extern double epot, vir;

/*
 * Compute forces and accumulate the virial and the potential
 */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     int i, j;
//     double sideh, rcoffs;
//     double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
//     double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
//     double forcex, forcey, forcez;

//     vir = 0.0;
//     epot = 0.0;
//     sideh = 0.5 * side;
//     rcoffs = rcoff * rcoff;

// #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez) reduction(+:epot,vir)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         xi = x[i];
//         yi = x[i + 1];
//         zi = x[i + 2];
//         fxi = 0.0;
//         fyi = 0.0;
//         fzi = 0.0;

//         for (j = i + 3; j < npart * 3; j += 3)
//         {
//             xx = xi - x[j];
//             yy = yi - x[j + 1];
//             zz = zi - x[j + 2];
//             if (xx < -sideh)
//                 xx += side;
//             if (xx > sideh)
//                 xx -= side;
//             if (yy < -sideh)
//                 yy += side;
//             if (yy > sideh)
//                 yy -= side;
//             if (zz < -sideh)
//                 zz += side;
//             if (zz > sideh)
//                 zz -= side;
//             rd = xx * xx + yy * yy + zz * zz;

//             if (rd <= rcoffs)
//             {
//                 rrd = 1.0 / rd;
//                 rrd2 = rrd * rrd;
//                 rrd3 = rrd2 * rrd;
//                 rrd4 = rrd2 * rrd2;
//                 rrd6 = rrd2 * rrd4;
//                 rrd7 = rrd6 * rrd;
//                 epot += (rrd6 - rrd3);
//                 r148 = rrd7 - 0.5 * rrd4;
//                 vir -= rd * r148;
//                 forcex = xx * r148;
//                 fxi += forcex;
//                 f[j] -= forcex;
//                 forcey = yy * r148;
//                 fyi += forcey;
//                 f[j + 1] -= forcey;
//                 forcez = zz * r148;
//                 fzi += forcez;
//                 f[j + 2] -= forcez;
//             }
//         }
//         f[i] += fxi;
//         f[i + 1] += fyi;
//         f[i + 2] += fzi;
//     }
// }


// /*
//  * Compute forces and accumulate the virial and the potential
//  */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     int i, j;
//     double sideh, rcoffs;
//     double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
//     double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
//     double forcex, forcey, forcez;
//     double local_epot = 0.0;
//     double local_vir = 0.0;

//     sideh = 0.5 * side;
//     rcoffs = rcoff * rcoff;

// #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) private(i, j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, forcex, forcey, forcez) reduction(+:local_epot, local_vir)
//     for (i = 0; i < npart * 3; i += 3)
//     {
//         xi = x[i];
//         yi = x[i + 1];
//         zi = x[i + 2];
//         fxi = 0.0;
//         fyi = 0.0;
//         fzi = 0.0;

//         for (j = i + 3; j < npart * 3; j += 3)
//         {
//             xx = xi - x[j];
//             yy = yi - x[j + 1];
//             zz = zi - x[j + 2];
//             if (xx < -sideh)
//                 xx += side;
//             if (xx > sideh)
//                 xx -= side;
//             if (yy < -sideh)
//                 yy += side;
//             if (yy > sideh)
//                 yy -= side;
//             if (zz < -sideh)
//                 zz += side;
//             if (zz > sideh)
//                 zz -= side;
//             rd = xx * xx + yy * yy + zz * zz;

//             if (rd <= rcoffs)
//             {
//                 rrd = 1.0 / rd;
//                 rrd2 = rrd * rrd;
//                 rrd3 = rrd2 * rrd;
//                 rrd4 = rrd2 * rrd2;
//                 rrd6 = rrd2 * rrd4;
//                 rrd7 = rrd6 * rrd;
//                 local_epot += (rrd6 - rrd3);
//                 r148 = rrd7 - 0.5 * rrd4;
//                 local_vir -= rd * r148;
//                 forcex = xx * r148;
//                 fxi += forcex;
//                 forcey = yy * r148;
//                 fyi += forcey;
//                 forcez = zz * r148;
//                 fzi += forcez;
//                 #pragma omp atomic update
//                 f[j] -= forcex;
//                 #pragma omp atomic update
//                 f[j + 1] -= forcey;
//                 #pragma omp atomic update
//                 f[j + 2] -= forcez;
//             }
//         }
//         #pragma omp atomic update
//         f[i] += fxi;
//         #pragma omp atomic update
//         f[i + 1] += fyi;
//         #pragma omp atomic update
//         f[i + 2] += fzi;
//     }
//     #pragma omp atomic update
//     epot += local_epot;
//     #pragma omp atomic update
//     vir += local_vir;
// }

/*
 * Compute forces and accumulate the virial and the potential
 */
void forces(int npart, double x[], double f[], double side, double rcoff)
{
    double sideh = 0.5 * side;
    double rcoffs = rcoff * rcoff;
    double local_epot = 0.0;
    double local_vir = 0.0;

    #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh, epot, vir) reduction(+:local_epot, local_vir) schedule(static)
    for (int i = 0; i < npart * 3; i += 3)
    {
        double local_fxi = 0.0;
        double local_fyi = 0.0;
        double local_fzi = 0.0;

        for (int j = i + 3; j < npart * 3; j += 3)
        {
            double xx = x[i] - x[j];
            double yy = x[i + 1] - x[j + 1];
            double zz = x[i + 2] - x[j + 2];

            // Apply periodic boundary conditions
            xx -= side * round(xx / side);
            yy -= side * round(yy / side);
            zz -= side * round(zz / side);

            double rd = xx * xx + yy * yy + zz * zz;

            if (rd <= rcoffs)
            {
                double rrd = 1.0 / rd;
                double rrd2 = rrd * rrd;
                double rrd3 = rrd2 * rrd;
                double rrd4 = rrd2 * rrd2;
                double rrd6 = rrd2 * rrd4;
                double rrd7 = rrd6 * rrd;
                local_epot += (rrd6 - rrd3);
                double r148 = rrd7 - 0.5 * rrd4;
                local_vir -= rd * r148;

                double forcex = xx * r148;
                double forcey = yy * r148;
                double forcez = zz * r148;

                local_fxi += forcex;
                local_fyi += forcey;
                local_fzi += forcez;

                #pragma omp atomic update
                f[j] -= forcex;
                #pragma omp atomic update
                f[j + 1] -= forcey;
                #pragma omp atomic update
                f[j + 2] -= forcez;
            }
        }

        #pragma omp atomic update
        local_epot += local_epot;
        #pragma omp atomic update
        local_vir += local_vir;

        #pragma omp atomic update
        epot += local_epot;
        #pragma omp atomic update
        vir += local_vir;

        #pragma omp atomic update
        f[i] += local_fxi;
        #pragma omp atomic update
        f[i + 1] += local_fyi;
        #pragma omp atomic update
        f[i + 2] += local_fzi;
    }
}

// Unsuccessful
// /*
//  * Compute forces and accumulate the virial and the potential
//  */
// void forces(int npart, double x[], double f[], double side, double rcoff)
// {
//     double sideh = 0.5 * side;
//     double rcoffs = rcoff * rcoff;
//     double local_epot = 0.0;
//     double local_vir = 0.0;

//     #pragma omp parallel for default(none) shared(npart, x, f, side, rcoffs, sideh) reduction(+:local_epot, local_vir) schedule(static)
//     for (int i = 0; i < npart; i++)
//     {
//         double local_fxi = 0.0;
//         double local_fyi = 0.0;
//         double local_fzi = 0.0;

//         for (int j = 0; j < npart; j++)
//         {
//             if (i != j)
//             {
//                 double xx = x[3 * i] - x[3 * j];
//                 double yy = x[3 * i + 1] - x[3 * j + 1];
//                 double zz = x[3 * i + 2] - x[3 * j + 2];

//                 // Apply periodic boundary conditions
//                 xx -= side * round(xx / side);
//                 yy -= side * round(yy / side);
//                 zz -= side * round(zz / side);

//                 double rd = xx * xx + yy * yy + zz * zz;

//                 if (rd <= rcoffs)
//                 {
//                     double rrd = 1.0 / rd;
//                     double rrd2 = rrd * rrd;
//                     double rrd3 = rrd2 * rrd;
//                     double rrd4 = rrd2 * rrd2;
//                     double rrd6 = rrd2 * rrd4;
//                     double rrd7 = rrd6 * rrd;
//                     local_epot += (rrd6 - rrd3);
//                     double r148 = rrd7 - 0.5 * rrd4;
//                     local_vir -= rd * r148;

//                     double forcex = xx * r148;
//                     double forcey = yy * r148;
//                     double forcez = zz * r148;

//                     local_fxi += forcex;
//                     local_fyi += forcey;
//                     local_fzi += forcez;
//                 }
//             }
//         }

//         #pragma omp atomic update
//         epot += local_epot;
//         #pragma omp atomic update
//         vir += local_vir;

//         #pragma omp atomic update
//         f[3 * i] += local_fxi;
//         #pragma omp atomic update
//         f[3 * i + 1] += local_fyi;
//         #pragma omp atomic update
//         f[3 * i + 2] += local_fzi;
//     }
// }



/*
 *  Main program : Molecular Dynamics simulation.
 */
int main(int arc, char **argv)
{
  int mm = atoi(argv[1]); //15;
  int npart = 4 * mm * mm * mm;
  int move;
  double x[npart * 3], vh[npart * 3], f[npart * 3];
  double ekin;
  double vel;
  double sc;
  double start, time;

  /*
   *  Parameter definitions
   */

  double den = 0.83134;
  double side = pow((double)npart / den, 0.3333333);
  double tref = 0.722;
  double rcoff = (double)mm / 4.0;
  double h = 0.064;
  int irep = 10;
  int istop = 20;
  int iprint = 5;
  int movemx = 20;

  double a = side / (double)mm;
  double hsq = h * h;
  double hsq2 = hsq * 0.5;
  double tscale = 16.0 / ((double)npart - 1.0);
  double vaver = 1.13 * sqrt(tref / 24.0);

  /*
   *  Initial output
   */

  printf(" Molecular Dynamics Simulation example program\n");
  printf(" ---------------------------------------------\n");
  printf(" number of particles is ............ %6d\n", npart);
  printf(" side length of the box is ......... %13.6f\n", side);
  printf(" cut off is ........................ %13.6f\n", rcoff);
  printf(" reduced temperature is ............ %13.6f\n", tref);
  printf(" basic timestep is ................. %13.6f\n", h);
  printf(" temperature scale interval ........ %6d\n", irep);
  printf(" stop scaling at move .............. %6d\n", istop);
  printf(" print interval .................... %6d\n", iprint);
  printf(" total no. of steps ................ %6d\n", movemx);

  /*
   *  Generate fcc lattice for atoms inside box
   */
  fcc(x, npart, mm, a);
  /*
   *  Initialise velocities and forces (which are zero in fcc positions)
   */
  mxwell(vh, 3 * npart, h, tref);
  dfill(3 * npart, 0.0, f, 1);
  /*
   *  Start of md
   */
  printf("\n    i       ke         pe            e         temp   "
         "   pres      vel      rp\n  -----  ----------  ----------"
         "  ----------  --------  --------  --------  ----\n");

  start = secnds();

  for (move = 1; move <= movemx; move++)
  {

    /*
     *  Move the particles and partially update velocities
     */
    domove(3 * npart, x, vh, f, side);

    /*
     *  Compute forces in the new positions and accumulate the virial
     *  and potential energy.
     */
    forces(npart, x, f, side, rcoff);

    /*
     *  Scale forces, complete update of velocities and compute k.e.
     */
    ekin = mkekin(npart, f, vh, hsq2, hsq);

    /*
     *  Average the velocity and temperature scale if desired
     */
    vel = velavg(npart, vh, vaver, h);
    if (move < istop && fmod(move, irep) == 0)
    {
      sc = sqrt(tref / (tscale * ekin));
      dscal(3 * npart, sc, vh, 1);
      ekin = tref / tscale;
    }

    /*
     *  Sum to get full potential energy and virial
     */
    if (fmod(move, iprint) == 0)
      prnout(move, ekin, epot, tscale, vir, vel, count, npart, den);
  }

  time = secnds() - start;

  printf("Elapsed time =  %f\n", (float)time);
}

time_t starttime = 0;

double secnds()
{

  return omp_get_wtime();
}
#include <stdio.h>
/*
 *  Scale forces, update velocities and compute K.E.
 */
  double
  mkekin(int npart, double f[], double vh[], double hsq2, double hsq){
    int i;
    double sum=0.0, ekin;

    for (i=0; i<3*npart; i++) {
      f[i]*=hsq2;
      vh[i]+=f[i];
      sum+=vh[i]*vh[i];
    }
    ekin=sum/hsq;

    return(ekin);
  }
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

  void srand48(long);
  double drand48(void);
/*
 *  Sample Maxwell distribution at temperature tref
 */
  void
  mxwell(double vh[], int n3, double h, double tref){
    int i;
    int npart=n3/3;
    double r, tscale, v1, v2, s, ekin=0.0, sp=0.0, sc;
    
    srand48(4711);
    tscale=16.0/((double)npart-1.0);

    for (i=0; i<n3; i+=2) {
      s=2.0;
      while (s>=1.0) {
        v1=2.0*drand48()-1.0;
        v2=2.0*drand48()-1.0;
        s=v1*v1+v2*v2;
      }
      r=sqrt(-2.0*log(s)/s);
      vh[i]=v1*r;
      vh[i+1]=v2*r;
    }

    for (i=0; i<n3; i+=3) sp+=vh[i];
    sp/=(double)npart;
    for(i=0; i<n3; i+=3) {
      vh[i]-=sp;
      ekin+=vh[i]*vh[i];
    }

    sp=0.0;
    for (i=1; i<n3; i+=3) sp+=vh[i];
    sp/=(double)npart;
    for(i=1; i<n3; i+=3) {
      vh[i]-=sp;
      ekin+=vh[i]*vh[i];
    }

    sp=0.0;
    for (i=2; i<n3; i+=3) sp+=vh[i];
    sp/=(double)npart;
    for(i=2; i<n3; i+=3) {
      vh[i]-=sp;
      ekin+=vh[i]*vh[i];
    }

    sc=h*sqrt(tref/(tscale*ekin));
    for (i=0; i<n3; i++) vh[i]*=sc;
  }

/*
 *   Print out interesting information at current timestep
 */
  void
  prnout(int move, double ekin, double epot, double tscale, double vir,
         double vel, double count, int npart, double den){
    double ek, etot, temp, pres, rp;

    ek=24.0*ekin;
    epot*=4.0;
    etot=ek+epot;
    temp=tscale*ekin;
    pres=den*16.0*(ekin-vir)/(double)npart;
    vel/=(double)npart;
    rp=(count/(double)npart)*100.0;
    printf(" %6d%12.4f%12.4f%12.4f%10.4f%10.4f%10.4f%6.1f\n",
           move,ek,epot,etot,temp,pres,vel,rp);

  }

/*
 *  Compute average velocity
 */
  double
  velavg(int npart, double vh[], double vaver, double h){
    int i;
    double vaverh=vaver*h;
    double vel=0.0;
    double sq;
    extern double count;

    count=0.0;
    for (i=0; i<npart*3; i+=3){
      sq=sqrt(vh[i]*vh[i]+vh[i+1]*vh[i+1]+vh[i+2]*vh[i+2]);
      if (sq>vaverh) count++;
      vel+=sq;
    }
    vel/=h;

    return(vel);
  }
