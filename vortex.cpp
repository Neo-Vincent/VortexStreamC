#pragma warning(disable:4996)
#pragma warning(disable : 4996)
#define _AFX_SECURE_NO_WARNINGS
#define NS2D_H     _declspec(dllexport);
#include "vortexdll.h"

double*  Make1DArray(int m)
{
	double *a;
	a = (double*)malloc(sizeof(double)*m);
	return a;
}

double** Make2DArray(int m, int n)
{
	int i;
	double **a;
	a = (double**)malloc(sizeof(double*)*m);
	for (i = 0; i<m; i++)
		a[i] = (double*)malloc(sizeof(double)*n);
	return a;
}

void free2DArray(double **a, int m)
{
	int i;
	for (i = 0; i<m; i++)
		free(a[i]);
	free(a);
}

ns2d::ns2d() {
	nx = 100;
	ny = 100;
	initParameters();
}
ns2d::ns2d(int tnx, int tny) {
	nx = tnx;
	ny = tny;
	initParameters();
}

void ns2d::initParameters() {
	isStrM = false;
	Re = 400.;
	sorwei = 1.8;
	dt = 1.e-3;
	//dt = 1.0;
	method = 0;
	continuous = False;
	totalT = 100.;
	outT = 0.5;
	T0 = 0;
	t00 = 0.;
	outName = False;
	init();
}
void ns2d::init() {
	u = Make2DArray(nx + 1, ny + 1);
	v = Make2DArray(nx + 1, ny + 1);
	vor = Make2DArray(nx + 1, ny + 1);
	p = Make2DArray(nx + 1, ny + 1);
	phi = Make2DArray(nx + 1, ny + 1);
	dx = 1. / nx;
	dy = 1. / ny;
	dx2inv = 1. / (dx*dx);
	dy2inv = 1. / (dy*dy);
	dsinv = 0.5 / (dx2inv + dy2inv);
	dxinv = 1. / dx;
	dyinv = 1. / dy;
	Reinv = 1. / Re;
	dx2invt = dt*dx2inv;
	dy2invt = dt*dy2inv;
	dxinvt = dt*dxinv;
	dyinvt = dx*dyinv;
	dxdy = dx*dy;
	dxt = dx*dt;
	dyt = dy*dt;
	Rey = Reinv*dy*dt / dx;
	Rex = Reinv*dx*dt / dy;
	dxdyinv = dxinv*dyinv;
	T0 = t00 / dt;
	int i, j;
	// init flows
	for (i = 0; i <= nx; i++)
		for (j = 0; j <= ny; j++)
		{
			u[i][j] = 0.;
			v[i][j] = 0.;
			vor[i][j] = 0.;
			p[i][j] = 1.;
			phi[i][j] = 0.;
		}

	for (i = 1; i < nx; i++) {
		u[i][ny] = 1.;
	}
}

ns2d::~ns2d() {
	//Q_Destr(&stfMatrix);
	//V_Destr(&stfBp);
	//V_Destr(&stfX);
	//Q_Destr(&vorMatrix);
	//V_Destr(&vorBp);
	//V_Destr(&vorX);
	free2DArray(u, nx + 1);
	free2DArray(v, nx + 1);
	free2DArray(vor, nx + 1);
	free2DArray(p, nx + 1);
	free2DArray(phi, nx + 1);
}


void ns2d::boundCondition() {
	int i, j;
	for (i = 1; i<nx; i++) {
		vor[i][0] = dyinv*(u[i][1] - u[i][0]);
		vor[i][ny] = dyinv*(u[i][ny] - u[i][ny - 1]);
	}
	for (j = 1; j<ny; j++) {
		vor[0][j] = -dxinv*(v[1][j] - v[0][j]);
		vor[nx][j] = -dxinv*(v[nx][j] - v[nx - 1][j]);
	}
}
void ns2d::streamEq() {
	int i, j, iter;
	double Res, Res1, tmp, ttt;
	// stream function, solving using SOR iteration
	for (iter = 0; iter<500; iter++)
	{
		Res = 0.;
		for (i = 1; i < nx; i++) {
			for (j = 1; j < ny; j++)
			{
				tmp = phi[i][j];
				ttt = dsinv*(dx2inv*(phi[i + 1][j] + phi[i - 1][j]) +
					dy2inv*(phi[i][j + 1] + phi[i][j - 1]) - vor[i][j]); // Gaussian-Seidel iteration;
				phi[i][j] = (1. - sorwei)*phi[i][j] + sorwei*ttt;
				Res += fabs(phi[i][j] - tmp);
			}
		}
		if (iter == 0)
			Res1 = Res;
		else if (Res / Res1<1.e-2)
			break;
	}
}
void ns2d::streamEq2() {
	if (!isStrM) {
		constructStfMatrix();
	}
	int i, j;
	SetRTCAccuracy(1.e-3);
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {
			if (i == 0 || j == 0 || i == nx || j == ny) {
				// right side value
				V_SetCmp(&stfBp, i*(nx + 1) + j + 1, phi[i][j]);
			}
			else {
				V_SetCmp(&stfBp, i*(nx + 1) + j + 1, vor[i][j]);
			}
		}
	}
	CGIter(&stfMatrix, &stfX, &stfBp, 400, SSORPrecond, 1.0);
	//GMRESIter(&stfMatrix, &stfX, &stfBp, 400, SSORPrecond, 1.0);
	for (i = 0; i < nx; i++) {
		for (j = 0; j <= ny; j++) {

			phi[i][j] = stfX.Cmp[i*(nx + 1) + j + 1];
		}
	}
	Q_Destr(&stfMatrix);
	V_Destr(&stfBp);
	V_Destr(&stfX);
}
void ns2d::constructStfMatrix() {
	isStrM = True;
	int size, i, j;
	size = (nx + 1)*(ny + 1);
	double a;
	a = -1. / dsinv;
	// stf solution vector
	V_Constr(&stfX, "stf", size, Normal, True);
	V_SetAllCmp(&stfX, 1.);
	Q_Constr(&stfMatrix, "stf eq matrix", size, True, Rowws, Normal, True);
	V_Constr(&stfBp, "rightH", size, Normal, True);
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {
			if (i == 0 || j == 0 || i == nx || j == ny) {
				// right side value
				V_SetCmp(&stfBp, i*(nx + 1) + j + 1, phi[i][j]);
				// i row, 1 non-zero values
				Q_SetLen(&stfMatrix, i*(nx + 1) + j + 1, 1);
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 0, i*(nx + 1) + j + 1, 1.);
			}
			else {
				/*
				V_SetCmp(&stfBp, i*(nx + 1) + j + 1, vor[i][j]);
				Q_SetLen(&stfMatrix, i*(nx + 1) + j + 1, 5); // i*(nx + 1) + j + 1 row, 5 non-zero values
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 0, (i - 1)*(nx + 1) + j + 1, dx2inv);
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 1,  i     *(nx + 1) + j    , dy2inv);
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 2, (i    )*(nx + 1) + j + 1, a     );
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 3, (i    )*(nx + 1) + j + 2, dy2inv);
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 4, (i + 1)*(nx + 1) + j + 1, dx2inv);
				*/
				V_SetCmp(&stfBp, i*(nx + 1) + j + 1, vor[i][j]);
				Q_SetLen(&stfMatrix, i*(nx + 1) + j + 1, 3); // i*(nx + 1) + j + 1 row, 5 non-zero values
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 0, (i)*(nx + 1) + j + 1, a);
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 1, (i)*(nx + 1) + j + 2, dy2inv);
				Q_SetEntry(&stfMatrix, i*(nx + 1) + j + 1, 2, (i + 1)*(nx + 1) + j + 1, dx2inv);
			}
		}
	}
}
void ns2d::velocity() {
	int i, j;
	for (i = 1; i < nx; i++) {
		for (j = 1; j < ny; j++)
		{
			u[i][j] = 0.5* dyinv * (phi[i][j + 1] - phi[i][j - 1]);
			v[i][j] = -0.5* dxinv * (phi[i + 1][j] - phi[i - 1][j]);
		}
	}
}
void ns2d::presure() {
	int i, j, iter;
	double Res, Res1, tmp, ttt;
	// stream function, solving using SOR iteration
	for (iter = 0; iter<500; iter++)
	{
		Res = 0.;
		for (i = 1; i < nx; i++) {
			for (j = 1; j < ny; j++)
			{
				tmp = p[i][j];
				ttt = dsinv*(dx2inv*(p[i + 1][j] + p[i - 1][j]) +
					dy2inv*(p[i][j + 1] + p[i][j - 1]) -
					((u[i + 1][j] - u[i - 1][j])*(v[i][j + 1] - v[i][j - 1])*dxdyinv -
						(u[i][j + 1] - u[i][j - 1])*(v[i + 1][j] - v[i - 1][j])*dxdyinv)); // Gaussian-Seidel iteration;
				p[i][j] = (1. - sorwei)*p[i][j] + sorwei*ttt;
				Res += fabs(p[i][j] - tmp);
			}
		}
		if (iter == 0)
			Res1 = Res;
		else if (Res / Res1<1.e-2)
			break;
	}
}
void ns2d::solve()
{
	init();
	printf("initialization OK!\n");
	int n;
	int N = totalT / dt + 0.5;
	int outN = outT / dt + 05;
	int onepercent = 0.01*N;
	if (onepercent == 0) {
		onepercent = 0.1*N;
	}
	if (onepercent == 0) {
		onepercent = N;
	}
	printf("Re=%f,T=%f,dt=%f,method=%d\n", Re, totalT, dt, method);
	for (n = T0; n<N + 1; n++)
	{

		switch (method)
		{
		case 0:
			fdm();
			break;
		case 1:
			fdm2();
			break;
		case 2:
			fvm1_ex();
			break;
		case 3:
			fvm2_ex();
			break;
		case 4:

			fvm1_im();
			break;
		case 5:
			fvm2_im();
			break;
		default:
			fvm2_im();
			break;
		}

		if (n%onepercent == 0) {
			printf("%2.1f%% completed, t=%f \n", float(n) / N*100., dt*n);
		}
		if (n % outN == 0) {
			if (continuous) {
				output(dt*n);
			}
			else
				if (outName) {
					output(outputName);
				}
				else
					output();
		}
	}
	if (outName) {
		output(outputName);
	}
	else
		output();
}
void ns2d::solve(int Method)
{

	printf("initialization OK!\n");
	int n;
	for (n = 0; n<500000; n++)
	{

		switch (Method)
		{
		case 0:
			fdm();
			break;
		case 1:
			fdm2();
			break;
		case 2:
			fvm1_ex();
			break;
		case 3:
			fvm2_ex();
			break;
		case 4:

			fvm1_im();
			break;
		case 5:

			fvm2_im();
			break;
		default:
			fdm();
			break;
		}
		//fvm1_im();
		//fvm1_ex();
		//fvm2_ex();
		printf("n=%d \n", n);
		printf("vor=%f", vor[70][70]);
		if (n % 50 == 0)
			printf("n=%d \n", n);
		if (n % 500 == 0) {
			output();
			// return 0;
		}
	}
}



void ns2d::fdm() {
	//copy from the Li's code
	int i, j;
	double fxl, fxr, fyl, fyr, flux, fluy;
	boundCondition();
	// vortex evolution, solve using first order Euler method
	for (i = 1; i < nx; i++) {
		for (j = 1; j < ny; j++)
		{
			fxr = (vor[i + 1][j] - vor[i][j]) * dxinv;
			fxl = (vor[i][j] - vor[i - 1][j]) * dxinv;
			fyr = (vor[i][j + 1] - vor[i][j]) * dyinv;
			fyl = (vor[i][j] - vor[i][j - 1]) * dyinv;
			flux = -0.5*((u[i][j] - fabs(u[i][j]))*fxr + (u[i][j] + fabs(u[i][j]))*fxl);
			fluy = -0.5*((v[i][j] - fabs(v[i][j]))*fyr + (v[i][j] + fabs(v[i][j]))*fyl);

			vor[i][j] = vor[i][j] + dt*(flux + fluy +
				1. / Re*(dx2inv*(vor[i + 1][j] - 2.*vor[i][j] + vor[i - 1][j]) +
					dy2inv*(vor[i][j + 1] - 2.*vor[i][j] + vor[i][j - 1])));
		}
	}
	streamEq();
	velocity();
}
void ns2d::fdm_im() {
	//copy from the Li's code
	int i, j;
	double fxl, fxr, fyl, fyr, flux, fluy;
	boundCondition();
	// vortex evolution, solve using first order Euler method
	for (i = 1; i < nx; i++) {
		for (j = 1; j < ny; j++)
		{
			fxr = (vor[i + 1][j] - vor[i][j]) * dxinv;
			fxl = (vor[i][j] - vor[i - 1][j]) * dxinv;
			fyr = (vor[i][j + 1] - vor[i][j]) * dyinv;
			fyl = (vor[i][j] - vor[i][j - 1]) * dyinv;
			flux = -0.5*((u[i][j] - fabs(u[i][j]))*fxr + (u[i][j] + fabs(u[i][j]))*fxl);
			fluy = -0.5*((v[i][j] - fabs(v[i][j]))*fyr + (v[i][j] + fabs(v[i][j]))*fyl);

			vor[i][j] = vor[i][j] + dt*(flux + fluy +
				1. / Re*(dx2inv*(vor[i + 1][j] - 2.*vor[i][j] + vor[i - 1][j]) +
					dy2inv*(vor[i][j + 1] - 2.*vor[i][j] + vor[i][j - 1])));
		}
	}
	streamEq();
	velocity();
}
void ns2d::fdm2() {
	//copy from the Li's code
	int i, j;
	double fxl, fxr, fyl, fyr, flux, fluy;
	boundCondition();
	// vortex evolution, solve using first order Euler method
	for (i = 1; i < nx; i++) {
		for (j = 1; j < ny; j++)
		{
			fxr = (vor[i + 1][j] - vor[i][j]) * dxinv;
			fxl = (vor[i][j] - vor[i - 1][j]) * dxinv;
			fyr = (vor[i][j + 1] - vor[i][j]) * dyinv;
			fyl = (vor[i][j] - vor[i][j - 1]) * dyinv;
			flux = -0.5*((u[i][j] - fabs(u[i][j]))*fxr + (u[i][j] + fabs(u[i][j]))*fxl);
			fluy = -0.5*((v[i][j] - fabs(v[i][j]))*fyr + (v[i][j] + fabs(v[i][j]))*fyl);

			vor[i][j] = vor[i][j] + dt*(flux + fluy +
				1. / Re*(dx2inv*(vor[i + 1][j] - 2.*vor[i][j] + vor[i - 1][j]) +
					dy2inv*(vor[i][j + 1] - 2.*vor[i][j] + vor[i][j - 1])));
		}
	}
	streamEq2();
	velocity();
}
void ns2d::fvm1_ex() {
	int i, j;
	boundCondition();
	double ue, uw, aue, auw, uomega;
	double vn, vs, avn, avs, vomega;
	double right;

	for (i = 1; i < nx; i++) {
		for (j = 1; j < ny; j++)
		{
			ue = (u[i + 1][j] + u[i][j])*0.5;
			uw = (u[i - 1][j] + u[i][j])*0.5;
			aue = fabs(ue);
			auw = fabs(uw);
			uomega = -dxinv*0.5*((vor[i][j] * (ue - uw + aue + auw) + vor[i + 1][j] * (ue - aue) - vor[i - 1][j] * (uw + auw)));
			vn = (v[i][j + 1] + v[i][j])*0.5;
			vs = (v[i][j - 1] + v[i][j])*0.5;
			avn = fabs(vn);
			avs = fabs(vs);
			vomega = -dyinv*0.5*((vor[i][j] * (vn - vs + avn + avs) + vor[i][j + 1] * (vn - avs) - vor[i][j - 1] * (vs + avs)));
			right = 1. / Re*((vor[i + 1][j] + vor[i - 1][j] - 2 * vor[i][j])*dx2inv + (vor[i][j + 1] + vor[i][j - 1] - 2 * vor[i][j])*dy2inv);
			vor[i][j] = dt*(uomega + vomega + right) + vor[i][j];
		}
	}
	streamEq();
	velocity();
}
void ns2d::fvm2_ex() {
	int i, j;
	boundCondition();
	double ue, uw, aue, auw, uomega;
	double vn, vs, avn, avs, vomega;
	double right;

	for (i = 1; i < nx; i++) {
		for (j = 1; j < ny; j++)
		{
			ue = (u[i + 1][j] + u[i][j])*0.5;
			uw = (u[i - 1][j] + u[i][j])*0.5;
			aue = fabs(ue);
			auw = fabs(uw);
			vn = (v[i][j + 1] + v[i][j])*0.5;
			vs = (v[i][j - 1] + v[i][j])*0.5;
			avn = fabs(vn);
			avs = fabs(vs);
			if (i == 1 || j == 1 || i == nx - 1 || j == ny - 1) {
				uomega = -dxinv*0.5*((vor[i][j] * (ue - uw + aue + auw) + vor[i + 1][j] * (ue - aue) - vor[i - 1][j] * (uw + auw)));
				vomega = -dyinv*0.5*((vor[i][j] * (vn - vs + avn + avs) + vor[i][j + 1] * (vn - avs) - vor[i][j - 1] * (vs + avs)));
				right = 1. / Re*((vor[i + 1][j] + vor[i - 1][j] - 2 * vor[i][j])*dx2inv + (vor[i][j + 1] + vor[i][j - 1] - 2 * vor[i][j])*dy2inv);
				vor[i][j] = dt*(uomega + vomega + right) + vor[i][j];
			}
			else {
				uomega = -dxinv*(0.75*vor[i][j] * (ue - uw + aue + auw) + 0.25* vor[i + 1][j] * (3 * ue - 3 * aue + uw - auw)
					- 0.25*vor[i - 1][j] * (3 * uw + 3 * auw + ue + aue)
					+ 0.25*vor[i - 2][j] * (uw + auw) - 0.25*vor[i + 2][j] * (ue - aue));
				vomega = -dyinv*(0.75*vor[i][j] * (vn - vs + avn + avs) + 0.25* vor[i][j + 1] * (3 * vn - 3 * avn + vs - avs)
					- 0.25*vor[i][j - 1] * (3 * vs + 3 * avs + vn + avn)
					+ 0.25*vor[i][j - 2] * (vs + avs) - 0.25*vor[i][j + 2] * (vn - avn));
				right = 1. / Re*((vor[i + 1][j] + vor[i - 1][j] - 2 * vor[i][j])*dx2inv + (vor[i][j + 1] + vor[i][j - 1] - 2 * vor[i][j])*dy2inv);
				vor[i][j] = dt*(uomega + vomega + right) + vor[i][j];
			}

		}
	}
	streamEq();
	velocity();
}
void ns2d::constructVorMatrix1() {
	int size, i, j;
	size = (nx + 1)*(ny + 1);
	double ue, uw, vn, vs, aue, auw, avn, avs;
	double a, b, c, d, e;
	// stf solution vector
	V_Constr(&vorX, "vor", size, Normal, True);
	V_SetAllCmp(&vorX, 1.);
	Q_Constr(&vorMatrix, "vor eq matrix", size, False, Rowws, Normal, True);
	V_Constr(&vorBp, "rightH", size, Normal, True);
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {


			if (i == 0 || j == 0 || i == nx || j == ny) {
				// right side value
				V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j]);
				// i row, 1 non-zero values
				Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 1);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, i*(nx + 1) + j + 1, 1.);
			}
			else {
				ue = (u[i + 1][j] + u[i][j])*0.5;
				uw = (u[i - 1][j] + u[i][j])*0.5;
				aue = fabs(ue);
				auw = fabs(uw);
				vn = (v[i][j + 1] + v[i][j])*0.5;
				vs = (v[i][j - 1] + v[i][j])*0.5;
				avn = fabs(vn);
				avs = fabs(vs);
				a = -0.5*(uw + auw)*dyt - Rey;  //Ww
				b = -0.5*(vs + avs)*dxt - Rex;//Ws
				c = dxdy + 0.5*((ue + aue - uw + auw)*dyt + (vn + avn - vs + avs)*dxt) + 2.*(Rey + Rex); //wp
				d = 0.5*(vn - avn)*dxt - Rex;//wn
				e = 0.5*(ue - aue)*dyt - Rey;//we
				V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j] * dxdy);
				Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 5); // i*(nx + 1) + j + 1 row, 5 non-zero values
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, (i - 1)*(nx + 1) + j + 1, a);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 1, i     *(nx + 1) + j, b);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 2, (i)*(nx + 1) + j + 1, c);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 3, (i)*(nx + 1) + j + 2, d);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 4, (i + 1)*(nx + 1) + j + 1, e);
			}
		}
	}
}
void ns2d::fvm1_im() {

	int i, j;
	boundCondition();
	constructVorMatrix1();
	SetRTCAccuracy(1.e-3);
	//CGIter(&vorMatrix, &vorX, &vorBp, 400, SSORPrecond, 1.0);
	GMRESIter(&vorMatrix, &vorX, &vorBp, 400, SSORPrecond, 1.0);
	for (i = 0; i < nx; i++) {
		for (j = 0; j <= ny; j++) {
			vor[i][j] = vorX.Cmp[i*(nx + 1) + j + 1];
		}
	}
	Q_Destr(&vorMatrix);
	V_Destr(&vorBp);
	V_Destr(&vorX);
	streamEq();
	velocity();
}

void ns2d::constructVorMatrix2() {
	int size, i, j;
	size = (nx + 1)*(ny + 1);
	double ue, uw, vn, vs, aue, auw, avn, avs;
	double a, b, c, d, e, f, g, h, k;
	// stf solution vector
	V_Constr(&vorX, "vor", size, Normal, True);
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {
			V_SetCmp(&vorX, i*(nx + 1) + j + 1, vor[i][j] + 0.000001);
		}
	}
	//V_SetAllCmp(&vorX, 0.1);
	Q_Constr(&vorMatrix, "vor_matrix", size, False, Rowws, Normal, True);
	V_Constr(&vorBp, "rightH", size, Normal, True);
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {
			if (i == 0 || j == 0 || i == nx || j == ny) {
				// right side value
				V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j]);
				Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 1);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, i*(nx + 1) + j + 1, 1.);
			}
			else {
				ue = (u[i + 1][j] + u[i][j])*0.5;
				uw = (u[i - 1][j] + u[i][j])*0.5;
				aue = fabs(ue);
				auw = fabs(uw);
				vn = (v[i][j + 1] + v[i][j])*0.5;
				vs = (v[i][j - 1] + v[i][j])*0.5;
				avn = fabs(vn);
				avs = fabs(vs);
				if (i == 1 || j == 1 || i == nx - 1 || j == ny - 1) {
					a = -0.5*(uw + auw)*dyt - Rey;  //Ww
					b = -0.5*(vs + avs)*dxt - Rex;//Ws
					c = dxdy + 0.5*((ue + aue - uw + auw)*dyt + (vn + avn - vs + avs)*dxt) + 2.*(Rey + Rex); //wp
					d = 0.5*(vn - avn)*dxt - Rex;//wn
					e = 0.5*(ue - aue)*dyt - Rey;//we
					V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j] * dxdy);
					Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 5); // i*(nx + 1) + j + 1 row, 5 non-zero values
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, (i - 1)*(nx + 1) + j + 1, a);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 1, i      *(nx + 1) + j, b);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 2, i      *(nx + 1) + j + 1, c);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 3, i      *(nx + 1) + j + 2, d);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 4, (i + 1)*(nx + 1) + j + 1, e);
					//printf("11^_^!!\n");
				}
				else {

					a = 0.25*(uw + auw)*dyt;  //Www
					b = -Rey - 0.25*(ue + aue + 3. * uw + 3.*auw); //Ww
					c = 0.25*(vs + avs)*dxt;  //Wss
					d = -Rex - 0.25*(vn + avn + 3. * vs + 3.*avs); //Ws
					e = dxdy + 0.75*((ue - uw + aue + auw)*dyt + (vn - vs + avn + avs)*dxt) + 2.*(Rey + Rex); //Wp
					f = -Rex + 0.25*(3.*vn - 3.*avn + vs - avs)*dyt;//Wn
					g = 0.25*(-vn + avn)*dxt;//Wnn
					h = -Rey + 0.25*(3.*ue - 3.*aue + uw - auw)*dyt;//We
					k = 0.25*(-ue + aue)*dyt;//Wee

					V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j] * dxdy);
					Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 9); // i*(nx + 1) + j + 1 row, 5 non-zero values
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, (i - 2)*(nx + 1) + j + 1, a);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 1, (i - 1)*(nx + 1) + j + 1, b);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 2, i      *(nx + 1) + j - 1, c);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 3, i      *(nx + 1) + j, d);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 4, i      *(nx + 1) + j + 1, e);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 5, i      *(nx + 1) + j + 2, f);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 6, i      *(nx + 1) + j + 3, g);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 7, (i + 1)*(nx + 1) + j + 1, h);
					Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 8, (i + 2)*(nx + 1) + j + 1, k);

				}
			}
		}
	}
}
void ns2d::constructVorMatrix3() {
	fvm1_im();
	int size, i, j;
	size = (nx + 1)*(ny + 1);
	double ue, uw, vn, vs, aue, auw, avn, avs;
	double a, b, c, d, e, f, g, h, k;
	// stf solution vector
	V_Constr(&vorX, "vor", size, Normal, True);
	V_SetAllCmp(&vorX, 1.);
	Q_Constr(&vorMatrix, "vor eq matrix", size, False, Rowws, Normal, True);
	V_Constr(&vorBp, "rightH", size, Normal, True);
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++) {


			if (i == 0 || j == 0 || i == nx || j == ny || i == 1 || j == 1 || i == nx - 1 || j == ny - 1) {
				//if(True){
				// right side value
				V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j]);
				Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 1);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, i*(nx + 1) + j + 1, 1.);
			}
			else {
				ue = (u[i + 1][j] + u[i][j])*0.5;
				uw = (u[i - 1][j] + u[i][j])*0.5;
				aue = fabs(ue);
				auw = fabs(uw);
				vn = (v[i][j + 1] + v[i][j])*0.5;
				vs = (v[i][j - 1] + v[i][j])*0.5;
				avn = fabs(vn);
				avs = fabs(vs);

				a = 0.25*(uw + auw)*dyt;  //Www
				b = -Rey - 0.25*(ue + aue + 3. * uw + 3.*auw); //Ww
				c = 0.25*(vs + avs)*dxt;  //Wss
				d = -Rex - 0.25*(vn + avn + 3. * vs + 3.*avs); //Ws
				e = dxdy + 0.75*((ue - uw + aue + auw)*dyt + (vn - vs + avn + avs)*dxt) + 2.*(Rey + Rex); //Wp
				f = -Rex + 0.25*(3.*vn - 3.*avn + vs - avs)*dyt;//Wn
				g = 0.25*(-vn + avn)*dxt;//Wnn
				h = -Rey + 0.25*(3.*ue - 3.*aue + uw - auw)*dyt;//We
				k = 0.25*(-ue + aue)*dyt;//Wee

				V_SetCmp(&vorBp, i*(nx + 1) + j + 1, vor[i][j] * dxdy);
				Q_SetLen(&vorMatrix, i*(nx + 1) + j + 1, 9); // i*(nx + 1) + j + 1 row, 5 non-zero values
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 0, (i - 2)*(nx + 1) + j + 1, a);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 1, (i - 1)*(nx + 1) + j + 1, b);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 2, (i)*(nx + 1) + j - 1, c);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 3, (i)*(nx + 1) + j, d);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 4, (i)*(nx + 1) + j + 1, e);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 5, (i)*(nx + 1) + j + 2, f);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 6, (i)*(nx + 1) + j + 3, g);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 7, (i + 1)*(nx + 1) + j + 1, h);
				Q_SetEntry(&vorMatrix, i*(nx + 1) + j + 1, 8, (i + 2)*(nx + 1) + j + 1, k);

			}
		}
	}
}
void ns2d::fvm2_im() {
	int i, j;
	boundCondition();
	constructVorMatrix2();
	SetRTCAccuracy(1.e-7);
	//CGIter(&vorMatrix, &vorX, &vorBp, 400, SSORPrecond, 1.0);
	//GMRESIter(&vorMatrix, &vorX, &vorBp, 8000, JacobiPrecond,1.0);
	//QMRIter(&vorMatrix, &vorX, &vorBp, 1000, JacobiPrecond, 0.000000001);
	GMRESIter(&vorMatrix, &vorX, &vorBp, 1000, JacobiPrecond, 1.e-10);
	int Iter = GetLastNoIter();
	double Res = GetLastAccuracy();
	cout << "iter=" << Iter << endl;
	cout << "res=" << Res << endl;


	for (i = 0; i < nx; i++) {
		for (j = 0; j <= ny; j++) {
			vor[i][j] = vorX.Cmp[i*(nx + 1) + j + 1];
		}
	}
	streamEq();
	velocity();
	Q_Destr(&vorMatrix);
	V_Destr(&vorBp);
	V_Destr(&vorX);

}


void ns2d::output()
{
	int i, j;
	FILE *fp;
	fp = fopen("out_ns.dat", "w+");
	fprintf(fp, " variables=\"x\",\"y\",\"u\",\"v\",\"vorticity\",\"stream line\",\"P\" \n");
	fprintf(fp, " zone i=%d   j=%d f=point \n", nx + 1, ny + 1);
	presure();
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++)
		{
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f \n ", dx*i, dy*j, u[i][j], v[i][j], vor[i][j], phi[i][j], p[i][j]);
		}
	}
	fclose(fp);
}
void ns2d::output(double outDt)
{
	int i, j;
	FILE *fp;
	char filename[40];
	sprintf(filename, "out_ns_%6.2f.dat", outDt);
	fp = fopen(filename, "w+");
	fprintf(fp, " variables=\"x\",\"y\",\"u\",\"v\",\"vorticity\",\"stream line\",\"P\" \n");
	fprintf(fp, " zone i=%d   j=%d f=point \n", nx + 1, ny + 1);
	presure();
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++)
		{
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f \n ", dx*i, dy*j, u[i][j], v[i][j], vor[i][j], phi[i][j], p[i][j]);
		}
	}
	fclose(fp);
}
void ns2d::output(char * filename)
{
	int i, j;
	FILE *fp;
	fp = fopen(filename, "w+");
	fprintf(fp, " variables=\"x\",\"y\",\"u\",\"v\",\"vorticity\",\"stream line\",\"P\" \n");
	fprintf(fp, " zone i=%d   j=%d f=point \n", nx + 1, ny + 1);
	presure();
	for (i = 0; i <= nx; i++) {
		for (j = 0; j <= ny; j++)
		{
			fprintf(fp, "%f,%f,%f,%f,%f,%f,%f \n ", dx*i, dy*j, u[i][j], v[i][j], vor[i][j], phi[i][j], p[i][j]);
		}
	}
	fclose(fp);
}



void ns2d::loadData(char* filename) {
	FILE *file;
	char buffer[512];
	char temp[10];
	char* t = 0;
	file = fopen(filename, "r");
	int i, j;
	double xij, yij;
	double uij, vij, vorij, phiij, pij;
	int tnx, tny;
	char aa[10], bb[10], cc[10];
	fgets(buffer, sizeof(buffer), file);
	fgets(buffer, sizeof(buffer), file);
	i = -1;
	j = 0;
	for (int k = 0; k <= sizeof(buffer); k++) {
		if (buffer[k] == '=') {
			i = 0;
			k++;
		}
		if (i >= 0) {
			temp[i] = buffer[k];
			i++;
		}
		if (buffer[k] == ' ') {
			if (i > 0) {
				if (j == 0) {
					tnx = atoi(temp);
					j++;
				}

				if (j == 1) {
					tny = atoi(temp);
					break;
				}
			}
			i = -1;
		}
	}
	nx = tnx;
	ny = tny;
	init();
	printf("loading file\n");
	while (!feof(file))
	{
		fgets(buffer, sizeof(buffer), file);
		t = buffer;
		sscanf(t, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &xij, &yij, &uij, &vij, &vorij, &phiij, &pij);
		i = xij / dx + 0.5;
		j = yij / dy + 0.5;
		u[i][j] = uij;
		v[i][j] = vij;
		phi[i][j] = phiij;
		vor[i][j] = vorij;
		p[i][j] = pij;
	}
	fclose(file);
}

void ns2d::setDim(int tnx, int tny) {
	nx = tnx;
	ny = tny;
}
void ns2d::setDt(double Dt) {
	dt = Dt;
}
void ns2d::setRe(double RE) {
	Re = RE;
}
void ns2d::setSorwei(double sw) {
	if (sw >= 1.0 && sw <= 2.0) {
		sorwei = sw;
	}
	else {
		printf("ERROR!! Overrelaxation Weight is out range! Please check!\nUse SORWEI=1.8 as the weight.\n");
	}
}
void ns2d::setT(double t) {
	totalT = t;
}
void ns2d::setT0(double t) {
	t00 = t;
}
void ns2d::setOutT(double t) {
	outT = t;
}
void ns2d::setMethod(int m) {
	method = m;
}
void ns2d::isCon(bool isc) {
	continuous = isc;
}
void ns2d::setOutName(char * outn) {
	outputName = outn;
	outName = True;
}
void ns2d::help() {
	printf("Usage:\n\
-h print this message\n\n\
-Re [number] set the Reynolds number, ex: -Re 800 means set Re=800\n\n\
-dt [number] set the time step, default value is 0.001s\n\n\
-Dim [number] [number] set the space mesh,\n ex -Dim 100 100 means meshing the space into a 100 X 100 matrix\n\n\
-t [number] set the time of simulation.\n\n\
-outT [number] set periode of print out put\n\n\
-a [yes/no] no option is print the last result, \n yes is print all result with a time periode T\n\n\
-load [filename] load history result file\n as initial data to compute the evaluation\n\
warming: for using this function, \n all paramerters should be well set, such dt,Re\n\n\
-t0 [number] set the initial time as t0 when we load a result file\n\n\
-o [filename] set the output filename\n\n\
-m [method] set the method of simulation, \
the availablee options are:\n\
fdm  finite diffrence method\n\
fvm1ex explicite finite volume method with 1er order of upwind scheme\n\
fvm2ex explicite finite volume method with 2nd order of upwind scheme\n\
fvm1im implicite finite volume method with 1er order of upwind scheme\n\
fvm2im explicite finite volume method with 2nd order of upwind scheme\n\
		");
}

bool argError(int t, int i) {
	if (i >= t) {
		printf("Argument ERROR!!");
		return False;
	}
	return True;
}
/*
int main(int argc, char *argv[])
{
ns2d vortex;
int filename = 0;
if (argc != 1) {
for (int i = 1; i < argc; i++) {
if (strcmp(argv[i],"-h")==0) {
vortex.help();
return 0;
}
if (strcmp(argv[i], "-Re")==0 && argError(argc,i)) {
double RE=atof(argv[i+1]);
vortex.setRe(RE);
}
if (strcmp(argv[i], "-dt")==0 && argError(argc, i)) {
double RE = atof(argv[i + 1]);
vortex.setDt(RE);
}
if (strcmp(argv[i], "-Dim")==0 && argError(argc, i+1)) {
int nx=atoi(argv[i+1]),ny=atoi(argv[i+2]);
vortex.setDim(nx,ny);
}
if (strcmp(argv[i], "-t")==0 && argError(argc, i)) {
double RE = atof(argv[i + 1]);
vortex.setT(RE);
}
if (strcmp(argv[i], "-t0") == 0 && argError(argc, i)) {
double RE = atof(argv[i + 1]);
vortex.setT0(RE);
}
if (strcmp(argv[i], "-outT")==0 && argError(argc, i)) {
double RE = atof(argv[i + 1]);
vortex.setOutT(RE);
}
if (strcmp(argv[i], "-a")==0 && argError(argc, i)) {
bool isCon = False;
if (strcmp(argv[i + 1], "no")==0) isCon = False;
if (strcmp(argv[i + 1], "yes" )==0) isCon = True;
vortex.isCon(isCon);
}
if (strcmp(argv[i], "-m")==0 && argError(argc, i)) {
int method = 0;
if (strcmp(argv[i + 1],"fdm"   )==0) method=0;
if (strcmp(argv[i + 1],"fvm1ex")==0) method = 2;
if (strcmp(argv[i + 1],"fvm2ex")==0) method = 3;
if (strcmp(argv[i + 1],"fvm1im")==0) method = 4;
if (strcmp(argv[i + 1],"fvm2im")==0) method = 5;
vortex.setMethod(method);
}
if (strcmp(argv[i], "-load") == 0 && argError(argc, i)) {
filename = i + 1;
//vortex.loadData(argv[i + 1]);
}
if (strcmp(argv[i], "-o") == 0 && argError(argc, i)) {
vortex.setOutName(argv[i + 1]);
//vortex.loadData(argv[i + 1]);
}
}
}
printf("\n");
vortex.init();
if (filename != 0) {
vortex.loadData(argv[filename]);
}
vortex.solve();
//0:fdm 1:fdm2 2:fvm_ex1 3:fvm_ex2 4:fvm_im1 5:fvm_im2
return 0;
}
*/
void main2() {
	ns2d vortex;
	vortex.solve(5);
}/**/

	ns2d * ns2d_new() {
		return new ns2d();
	}

	void ns2d_fdm(ns2d * n) {
		n->fdm();
	}
	void ns2d_fdm_im(ns2d * n) {
		n->fdm_im();
	}
	void ns2d_fdm2(ns2d * n) {
		n->fdm2();
	}
	void ns2d_fvm1_ex(ns2d * n) {
		n->fvm1_ex();
	}
	void ns2d_fvm2_ex(ns2d * n) {
		n->fvm2_ex();
	}
	void ns2d_fvm1_im(ns2d * n) {
		n->fvm1_im();
	}
	void ns2d_fvm2_im(ns2d * n) {
		n->fvm2_im();
	}
	void ns2d_solve0(ns2d * n) {
		n->solve();
	}
	void ns2d_solve1(ns2d * n, int a) {
		n->solve(a);
	}
	void ns2d_output0(ns2d * n) {
		n->output();
	}
	void ns2d_output1(ns2d * n,char * a) {
		n->output(a);
	}
	void ns2d_output2(ns2d * n, double a) {
		n->output(a);
	}
	void ns2d_setDt(ns2d * n, double  a) {
		n->setDt(a);
	}
	void ns2d_setRe(ns2d * n, double  a) {
		n->setRe(a);
	}
	void ns2d_setSorwei(ns2d * n, double a) {
		n->setSorwei(a);
	}
	void ns2d_setDim(ns2d * n, int a, int b) {
		n->setDim(a, b);
	}
	void ns2d_help(ns2d * n) {
		n->help();
	}
	void ns2d_init(ns2d * n) {
		n->init();
	}
	void ns2d_setT(ns2d * n, double a) {
		n->setT(a);
	}
	void ns2d_setT0(ns2d * n, double a) {
		n->setT0(a);
	}
	void ns2d_setOutT(ns2d * n, double a) {
		n->setOutT(a);
	}
	void ns2d_isCon(ns2d * n, bool a) {
		n->isCon(a);
	}
	void ns2d_setMethod(ns2d * n, int a) {
		n->setMethod(a);
	}
	void ns2d_loadData(ns2d * n, char* a) {
		n->loadData(a);
	}
	void ns2d_setOutName(ns2d * n, char* a) {
		n->setOutName(a);
	}