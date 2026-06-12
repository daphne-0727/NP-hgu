/*----------------------------------------------------------------\
@ C-Tutorial by Young-Keun Kim - Handong Global University

Author           : SSS LAB
Created          : 05-03-2021
Modified         : 08-19-2022
Language/ver     : C++ in MSVS2022

Description      : myNP_tutorial.h
/----------------------------------------------------------------*/

#pragma once
#define		_MY_NP_H
#define		PI		3.14159265358979323846264338327950288419716939937510582

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// power function
extern double power(double _x, int N);

//  Taylor series approximation for sin(x) (input unit: [rad])
extern double sinTaylor(double _x);

// Taylor series approximation for sin(x) (input unit: [deg])
extern double sindTaylor(double _x);

// factorial function
extern double factorial(double N);

extern void printVec(double* vec, int row);



/*==========================================================================*/
/*                    NONLINEAR		   		*/
/*==========================================================================*/
/* Bisection Method
	_a      : initial value #1
	_b      : initial value #2
	_tol   : tolerance
*/
double bisection(double func(double x), double _a0, double _b0, double _tol);
/* Newton-Raphson Method without Passing a Function  */
double newtonRaphson(double func(double x), double dfunc(double x), double _x0, double _tol);


/*==========================================================================*/
/*                  Interpolation	   		*/
/*==========================================================================*/
double linearInterp(double* x_arr, double* y_arr, int n, double xq);

/*==========================================================================*/
/*                   Differentiation	   		*/
/*==========================================================================*/
/* First order Differentiation from discrete dataset points*/
void gradient1D(double dydx[], double _x[], double _y[], int m);

/* First order Differentiation from a function*/
void gradientFunc(double dydx[], double func(const double x), double _x[], int m);

/* Second order Differentiation from discrete dataset points*/
void acceleration(double dy2dx2[], double x[], double y[], int m);


/*==========================================================================*/
/*                 Integration	   		*/
/*==========================================================================*/
//mixof simpson13 and simpson38
double simpson(double x[], double y[], int m);

//receiving function
double integral(double func(const double x), double a, double b, int n);

/*==========================================================================*/
/*                 ODE	   		*/
/*==========================================================================*/
void odeEU(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init);
void odeRK2(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init);
void odeRK3(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y0);
void sys2RK2(double y1[], double y2[], void odeFuncSys(double dYdt[], const double t, const double Y[]), const double t0, const double tf, const double h, const double y1_init, const double y2_init);

/*------------------------------------------------------------------------------------------*/
/*			Linear Regression							 		*/
/*------------------------------------------------------------------------------------------*/
// Calculates coefficients of least squares regression - Line
void linearRegression(double z_opt[], double xdata[], double ydata[], int dataN);
// Calculates coefficients of least squares regression - Nth order polynomial
void polyFit(double vecZ[], double vecX[], double vecY[], int dataL, int orderN);
// expRegression - Line
void expRegression(double z_opt[], double xdata[], double ydata[], int dataN);