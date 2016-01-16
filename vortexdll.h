

#ifdef  NS2D_H
#else
#define NS2D_H _declspec(dllimport)
#endif

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <malloc.h>

using std::vector;
using namespace std;
extern "C" {
#include "../laspack/laspack.h"
}


#define _AFX_SECURE_NO_WARNINGS


double*  Make1DArray(int m);
double** Make2DArray(int m, int n);
void free2DArray(double **a, int m);
 class ns2d
{
public:
	ns2d();
	ns2d(int tnx, int tny);
	~ns2d();

	void fdm();
	void fdm_im();
	void fdm2();
	void fvm1_ex();
	void fvm2_ex();
	void fvm1_im();
	void fvm2_im();
	void solve();
	void solve(int);
	void output();
	void output(char *);
	void output(double);
	void setDt(double);
	void setRe(double);
	void setSorwei(double);
	void setDim(int, int);
	void help();
	void init();
	void setT(double);
	void setT0(double);
	void setOutT(double);
	void isCon(bool);
	void setMethod(int);
	void loadData(char*);
	void setOutName(char*);

private:
	int nx, ny, T0;
	double Re, sorwei, Reinv;
	double ** u, ** v, ** vor, ** phi, **p;
	double dx, dy, dxinv, dsinv, dyinv, dx2inv, dy2inv;
	double dxinvt, dyinvt, dx2invt, dy2invt;
	double dt, dxdy, dxt, dyt, Rey, Rex, dxdyinv, t00;
	void initParameters();
	void boundCondition();
	void streamEq();
	void streamEq2();
	void constructStfMatrix();
	void constructVorMatrix1();
	void constructVorMatrix2();
	void constructVorMatrix3();
	void presure();
	void velocity();
	bool isStrM;
	int method;
	bool continuous;
	double totalT, outT;
	char* outputName;
	bool outName;
	QMatrix stfMatrix, vorMatrix;
	Vector  stfBp, stfX, vorX, vorBp;
};
 /*
void main2();
void ns2d_fdm(ns2d * n);
void ns2d_fdm_im(ns2d * n);
void ns2d_fdm2(ns2d * n);
void ns2d_fvm1_ex(ns2d * n);
void ns2d_fvm2_ex(ns2d * n);
void ns2d_fvm1_im(ns2d * n);
void ns2d_fvm2_im(ns2d * n);
void ns2d_solve0(ns2d * n);
void ns2d_solve1(ns2d * n, int);
void ns2d_output0(ns2d * n);
void ns2d_output1(ns2d * n, char *);
void ns2d_output2(ns2d * n, double);
void ns2d_setDt(ns2d * n, double);
void ns2d_setRe(ns2d * n, double);
void ns2d_setSorwei(ns2d * n, double);
void ns2d_setDim(ns2d * n, int, int);
void ns2d_help(ns2d * n);
void ns2d_init(ns2d * n);
void ns2d_setT(ns2d * n, double);
void ns2d_setT0(ns2d * n, double);
void ns2d_setOutT(ns2d * n, double);
void ns2d_isCon(ns2d * n, bool);
void ns2d_setMethod(ns2d * n, int);
void ns2d_loadData(ns2d * n, char*);
void ns2d_setOutName(ns2d * n, char*);*/

extern "C" _declspec(dllexport) void main2();
extern "C" _declspec(dllexport) ns2d* ns2d_new();
extern "C" _declspec(dllexport) void ns2d_fdm(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_fdm_im(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_fdm2(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_fvm1_ex(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_fvm2_ex(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_fvm1_im(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_fvm2_im(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_solve0(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_solve1(ns2d * n, int);
extern "C" _declspec(dllexport) void ns2d_output0(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_output1(ns2d * n, char *);
extern "C" _declspec(dllexport) void ns2d_output2(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_setDt(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_setRe(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_setSorwei(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_setDim(ns2d * n, int, int);
extern "C" _declspec(dllexport) void ns2d_help(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_init(ns2d * n);
extern "C" _declspec(dllexport) void ns2d_setT(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_setT0(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_setOutT(ns2d * n, double);
extern "C" _declspec(dllexport) void ns2d_isCon(ns2d * n, bool);
extern "C" _declspec(dllexport) void ns2d_setMethod(ns2d * n, int);
extern "C" _declspec(dllexport) void ns2d_loadData(ns2d * n, char*);
extern "C" _declspec(dllexport) void ns2d_setOutName(ns2d * n, char*);
