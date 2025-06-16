#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <ctime>
#include <iomanip>

int main() {
	int N_p = 50; // количество частиц в одной €чейке
	int i_max = 2001; // количество узлов
	int N = 2 * N_p * (i_max - 1); // глобальное количество частиц
	double* E = new double[i_max + 1];
	double* x = new double[N];
	double* v = new double[N];
	double* rho = new double[i_max];
	double* phi = new double[i_max]; // потенциал (необходимо найти)
	double eps = 1.e-10;
	double Max = 0.; 
	double c;
	int counter = 0;
	double u_0 = 0.8;
	double* a = new double[i_max + 1]; // прогоночный коэффициент 
	double* b = new double[i_max + 1]; // прогоночный коэффициент
	double m_i = 1. / N_p;
	double x_max = 100.;
	double h = x_max / (i_max - 1.);
	int N_1 = N_p * (i_max - 1);//количество частиц до добавлени€
	double EE = 0.;
	double tau = 0.0125; // добавить u_0
	double T_end = 3000.;
	double par_1 = u_0 * tau;
	double par_2 = x_max / N_1;
	int n_end = T_end / tau;
	double MAX_1;
	double MAX_2;

	for (int i = 0; i < N; i++) {
		x[i] = 0.;
		v[i] = 0.;
	}

	for (int i = 0; i < N_1; i++) {
		x[i] = par_2 * (i + 0.5);
		v[i] = u_0;
	}

	double first = x[0];
	//	std::cout << x[0] << " " << x[N_1 - 1] << std::endl; 
	for (int i = 0; i < i_max; i++) { //возьмем разность логарифма и (x_max-x)
		phi[i] = 0.;
		//std::cout << phi[i] << std::endl;
	}
	//std::ofstream file_max("D:\\Sp\\Phi data\\phi_max.txt");
	for (int n = 0; n < T_end / tau; n++) {

		for (int i = 0; i < i_max; i++) {
			rho[i] = 0.; 
		}

		//	std::cout << "n=" << n<< "   N_1=" << N_1<< std::endl;
		double ss = 0.;
		for (int i = 0; i < N_1; i++) {
			double d = x[i] / h;
			int t = int(d);
			ss = double(d - t);
			rho[t] = rho[t] + m_i * (1. - ss);
			rho[t + 1] = rho[t + 1] + m_i * ss;

		}
		rho[0] = 2. * rho[0];
		rho[i_max - 1] = 2. * rho[i_max - 1];
		//  std::cout << n-1 << " rho[0]= " << rho[0] << std::endl;
		a[1] = 0.; //исправить гран услови€ дл потенциала
		phi[0] = 0.;
		b[1] = phi[0];
		//phi[i_max] = Right;
		counter = 0;
		while (Max > eps || counter == 0) {
			for (int i = 1; i < i_max; i++) {
				double s = exp(phi[i]);
				double s3 = h * h * s;
				a[i + 1] = 1. / (2. + s3 - a[i]);
				b[i + 1] = (b[i] - s3 * (1. - phi[i]) + h * h * rho[i]) * a[i + 1];

			}
			phi[i_max - 1] = (a[i_max] * b[i_max - 1] + b[i_max]) / (1 - a[i_max] * a[i_max - 1]);
			Max = 0.;
			double s1 = 0.;
			for (int i = i_max - 2; i > 0; i--) {
				double s2 = a[i + 1] * phi[i + 1] + b[i + 1];
				double s3 = abs(s2 - phi[i]);
				if (s3 > Max) Max = s3;
				phi[i] = s2;


				//for (int i = 1; i < i_max; i++) {
				//	phi[i] = phi_1[i];
				//}
				std::cout<<counter++<<std::endl;
			}

		}
		MAX_1 = 0.;
		MAX_2 = 0.;
		for (int i = 0; i < i_max; i++) {
			if (phi[i] > MAX_1) {
				MAX_1 = phi[i];
				MAX_2 = i * h;
			}
		}

		//file_max << tau * n << "  " << MAX_2 << "  " << MAX_1 << std::endl;

		for (int i = 1; i < i_max; i++) { //нужны гран услови€ дл€ E
			E[i] = (-phi[i] + phi[i - 1]) / h;
		}
		E[i_max] = -E[i_max - 1];
		E[0] = 0.;
		for (int i = 0; i < N_1; i++) {
			double d = 0.5 + x[i] / h;
			int t = int(d);
			double ss = double(d - t);
			EE = E[t] * (1. - ss) + E[t + 1] * ss;
			//EE=0.;
			v[i] = v[i] + tau * EE;
			x[i] = x[i] + tau * v[i];
			if (x[i] < 0.) {
				x[i] = -x[i];
				v[i] = -v[i];
			}
			if (x[i] > x_max) {
				x[i] = 2. * x_max - x[i];
				v[i] = -v[i];
			}

		} first = first + par_1;
		//	std::cout << first << std::endl;
		while (first > par_2) {
			//std::cout << N_1 << " " << last << " fefefefefe" << std::endl;
			N_1 = N_1 + 1;
			x[N_1 - 1] = first - par_2;
			//	std::cout << "n=" << n<< ", N_1="<< N_1<< ", x_last=" << x[N_1-1] << std::endl;
			v[N_1 - 1] = u_0;
			first = first - par_2;

		}


		if (n % 80 == 0) { // нормировать
			std::ofstream file_1;
			std::ofstream file_2;// сам файл
			char filename[10]; // временный буфер
			std::string path_1;
			std::string path_2;
			std::string path_t = ".txt";
			_itoa_s(n, filename, 10);
			path_1 = "D:\\Sp\\DiplomB08\\rho_E_phi"; // часть имени файла (оно будет посто€нным)
			path_1 += filename + path_t; // собираем путь и им€ дл€ нового файла
			path_2 = "D:\\Sp\\DiplomB08\\velocity"; // часть имени файла (оно будет посто€нным)
			path_2 += filename + path_t; // собираем путь и им€ дл€ нового файла
			file_1.open(path_1.c_str());
			for (int i = 0; i < i_max; i++) {
				file_1 << std::fixed << std::setprecision(20) << i * h << "  " << rho[i] << "  " << E[i] << "  " << phi[i] << std::endl;
			}
			file_1.close();
			file_2.open(path_2.c_str());
			for (int i = 0; i < N_1; i++) {
				file_2 << x[i] << "  " << v[i] << std::endl;
			}
			file_2.close();
		}
	}

			//file_1.open("D:\\Sp\\Phi data\Phi_Data_u_0_0.8.txt");
			//file_1 << MAX;
		//file_max.close();
	}



