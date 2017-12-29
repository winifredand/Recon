#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include <Windows.h>
#include "reconst.h"

using namespace std;

double Likelihood(double x, double y, double z, double e, double *X, double *Y, double *Z, double *q)
{
	int i;
	double R = 1.0;
	double Rd[100], rd[100];
	double cos_theta[100];
	double f_cos[100];
	double ee[100];
	double u[100];
	double L = 0.0;
	double W;
	W = Effi_QE * Light_yield * S_pmt / 4 / PI;  //Other parameter's product;

	cout << "Optical model's parameters" << endl; 
	printf("i = i\tRd\t        rd\t        cos_theta\t f_cos\t          ee\t        u\t          L\t\n");
	for (i = 0; i < 100; i++)
	{
		Rd[i] = (x - X[i]) * (x - X[i]) + (y - Y[i]) * (y - Y[i]) + (z - Z[i]) * (z - Z[i]);
		rd[i] = sqrt(Rd[i]);
		if (rd[i] <= 0) break;
		cos_theta[i] = ((x - X[i]) * X[i] + (y - Y[i]) * Y[i] + (z - Z[i]) * Z[i]) / 2 / rd[i] / R;
		f_cos[i] = 0.999946 + cos_theta[i] * (0.101046 + cos_theta[i] * (-1.040140 + cos_theta[i] * (1.013810 + cos_theta[i] * -0.410953)));
		ee[i] = exp( - rd[i] / Lambda_a);
		u[i] = W * e * ee[i] * f_cos[i] / Rd[i];
		L += 2 * ((u[i] - q[i]) + log(q[i] / u[i]) * q[i]);

		printf("i = %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\n", i, Rd[i], rd[i], cos_theta[i], f_cos[i], ee[i], u[i], L);
	}
	system("pause");
	return L;
}

void GenerateInitData(double *xx, double *yy, double *zz, double *QQ)
{
	int i;
	double x, y, z, E;

	srand((unsigned)time(0));

	cout << "The PMT 's position :" << endl;

	for (i = 0; i < 100; i++)
	{
		x = rand() / (double)(RAND_MAX);
		y = rand() / (double)(RAND_MAX);
		z = rand() / (double)(RAND_MAX);
		E = rand() % 50;

		if ((x * x + y * y + z * z) == 1.0 && (x != y) || (z != y) || (x != z) && ((E > 20 && E < 50) || (E > 120 && E < 150)))
		{
			xx[i] = x;
			yy[i] = y;
			zz[i] = z;
			QQ[i] = E;
			printf("%lf	%lf	%lf	%lf\n", xx[i], yy[i], zz[i], QQ[i]);
		}
	}
}

void ReadData(double *x, double *y, double *z, double *Q)   //Read PMT's position and detected charge;
{
	FILE *fp;
	int i;

	if (ReadDataFromFile == 1)
	{
		GenerateInitData(x, y, z , Q);  //Find the features of 100 PMTs that are uniform distributed in ball's surface;
	}
	else
	{
		fp = fopen("data.inp", "r");
		if (fp == NULL) Error("Can not open ctl.inp file. Exit DEM.\n");
		for (i = 0; i < 100; i++)
		{
			fscanf(fp, "%lf	%lf	%lf	%lf", &x[i], &y[i], &z[i], &Q[i]);
		}
		fclose(fp);
	}

}

void Gradient(double *grad, double x, double y, double z, double e, double *x_PMT, double *y_PMT, double *z_PMT, double *Q)
{
	double delt_v = 1E-2;
	double x_up, x_down;
	double y_up, y_down;
	double z_up, z_down;
	double e_up, e_down;

	x_up = x + delt_v;
	y_up = y + delt_v;
	z_up = z + delt_v;
	e_up = e + delt_v;

	x_down = x - delt_v;
	y_down = y - delt_v;
	z_down = z - delt_v;
	e_down = e - delt_v;
	
	double L_x_up, L_x_down;
	double L_y_up, L_y_down;
	double L_z_up, L_z_down;
	double L_e_up, L_e_down;

	L_x_up = Likelihood(x_up, y, z, e, x_PMT, x_PMT,x_PMT, Q);
	L_x_down = Likelihood(x_down, y, z, e, x_PMT, x_PMT, x_PMT, Q);
	grad[0] = (L_x_up - L_x_down) / 2 / delt_v;

	L_y_up = Likelihood(x, y_up, z, e, x_PMT, x_PMT, x_PMT, Q);
	L_y_down = Likelihood(x, y_down, z, e, x_PMT, x_PMT, x_PMT, Q);
	grad[1] = (L_y_up - L_y_down) / 2 / delt_v;

	L_z_up = Likelihood(x, y, z_up, e, x_PMT, x_PMT, x_PMT, Q);
	L_z_down = Likelihood(x, y, z_down, e, x_PMT, x_PMT, x_PMT, Q);
	grad[2] = (L_z_up - L_z_down) / 2 / delt_v;

	L_e_up = Likelihood(x, y, z, e_up, x_PMT, x_PMT, x_PMT, Q);
	L_e_down = Likelihood(x, y, z, e_down, x_PMT, x_PMT, x_PMT, Q);
	grad[3] = (L_e_up - L_e_down) / 2 / delt_v;
}

void UpdateHesse(double(*H_old)[4], double (*H_new)[4], double *s, double *y)
{
	int i, j;
	double SS[4][4], sy = 0.0;
	for (i = 0; i < 4; i++)
	{
		sy += s[i] * y[i];
		for (j = 0; j < 4; j++){
			SS[i][j] = s[i] * s[j];
		}
	}
	
	double H42[4] = { 0. }, H43[4][4], H44[4][4];
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++){
			H42[i] += H_new[i][j] * y[j];
		}
	}
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++){
			H43[i][j] = H42[i] *y[j];
		}
	}

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++){
			H44[i][j] = H43[i][j] * H_new[i][j];
		}
	}

	double H31[4] = { 0. }, H32 = 0.0;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++){
			H31[i] += y[j] * H_new[j][i];
		}
	}
	for (i = 0; i < 4; i++)
	{
		H32 += H31[i] * y[i];
	}

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++){
			H_old[i][j] = H_new[i][j] + SS[i][j] / sy - H43[i][j] / H32;
		}
	}
}

/*
void CheckData(void)
{
	int i;
	FILE *fp;

	if ((fp = fopen("PMTdata.dat", "w")) == NULL)
		Error("\nCan not open file for debug_particle.\n");

	fprintf(fp, "(x,y,z) E)\n");
	for (i = 0; i < 100; i++)
	{
		fprintf(fp, "(%lf,%lf,%lf)\t", X_PMT[i], X_PMT[i], X_PMT[i]);
		fprintf(fp, "%lf\t", Q_expected[i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}
*/
int main()
{
	int i, j, k = 0;
	double x = 0.4, y = 0.5, z = 0.6, E = 2.506;  //Recontest data;
	double X_PMT[NUM];
	double Y_PMT[NUM];
	double Z_PMT[NUM];  //The position of PMT;
	double Q_expected[NUM];  //The expected charge; 

	ReadData(X_PMT, Y_PMT, Z_PMT, Q_expected);  //Read PMT's position and detected charge;
	//CheckData();  //Check the random data; 

	double x_new, x_old;
	double y_new, y_old;
	double z_new, z_old;
	double E_new, E_old;
	double Grad_new[4], Grad_old[4];
	double He_new[4][4], He_old[4][4];
	double Y[4], S[4];
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			if (i == j) He_old[i][j] = 1;
			He_old[i][j] = 0;
		}
	}
	
	x_old = x;
	y_old = y;
	z_old = z;
	E_old = E;

	Gradient(Grad_old, x_old, y_old, z_old, E_old, X_PMT, Y_PMT, Z_PMT, Q_expected);
	cout << Grad_old[0] << endl;
	system("pause");
	
	do
	{
		x_new = x_old;
		y_new = y_old;
		z_new = z_old;
		E_new = E_old;

		for (i = 0; i < 4; i++)
		{
			Grad_new[i] = Grad_old[i];
			for (j = 0; j < 4; j++)
			{
				He_new[i][j] = He_old[i][j];
			}
		}

		double HG[4];
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 4; j++)
			{
				HG[i] += He_new[i][j] * Grad_new[j];
			}
		}
		x_old = x_new - t * HG[0];
		y_old = y_new - t * HG[1];
		z_old = z_new - t * HG[2];
		E_old = E_new - t * HG[3];

		Gradient(Grad_old, x_old, y_old, z_old, E_old, X_PMT, Y_PMT, Z_PMT, Q_expected);

		for (i = 0; i < 4; i++)
		{
			Y[i] = Grad_old[i] - Grad_new[i];
		}

		S[0] = x_old - x_new;
		S[1] = y_old - y_new;
		S[2] = z_old - z_new;
		S[3] = E_old - E_new;

		UpdateHesse(He_old, He_new, S, Y);		

		k ++;
	} while (k <= 1);
	//while (abs(S[0]) < eps && abs(S[1]) < eps && abs(S[2]) < eps && abs(S[3]) < eps);

	cout << "The result: " << endl;
	cout << "(x,y,z) E :" << "\t" << x_old << "\t" << y_old
		<< "\t" << z_old << "\t" << E_old << endl;
	cout << "The K: " << k << endl;

	system("pause");
	return 0;
}
