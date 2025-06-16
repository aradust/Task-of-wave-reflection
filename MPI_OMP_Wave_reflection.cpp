#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstdio>
#include <iomanip>
#include <locale.h>
#include <vector>
#include <sstream>
#include <utility> // Для std::pair
#include <algorithm>
#include <random>
#include <ctime>
#include <chrono>


int N_p = 500; // количество частиц в одной ячейке
int i_max = 251; // количество узлов
int N_1 = N_p * (i_max - 1);//изначальное количество частиц до добавления
int N = 2 * N_p * (i_max - 1); // глобальное количество частиц
std::vector<double> E(i_max, 0); // электрическое поле
std::vector<double> x_e(N, 0);//координата электронов
std::vector<double> x_i(N, 0);//координата ионов
std::vector<double> v_e(N, 0);//скорость электронов
std::vector<double> v_i(N, 0);//скорость ионов
std::vector<double> rho_i(i_max + 1, 0.);//плотность ионов
std::vector<double> rho_e(i_max + 1, 0.);//плотность электронов
std::vector<double> phi(i_max + 1, 0.);
double x_max = 100.;//длина области расчета
//исходное расстояние между частицами
double m_e = 1.;//масса электронов
double m_i = 1. / N_p;
double Ion_m = m_e * 1600;//масса ионов
double h = x_max / (i_max - 1.);//шаг по пространству
double tau = 0.025; // шаг по времени
double Time = 30;//время расчета
int tau_quantity = Time / tau; //количество шагов расчёта
double u_0 = 0.3 / 40.; //скорость потока
double par_1 = u_0 * tau;
double par_2 = x_max / N_1;
double x0 = 0.;
int P = 0.;
double X_0 = 0.;
double X_1 = 0.;
double static first_ion = 0.;
double static first_el = 0.;
//std::ofstream file_first("first.txt");





/*inline void Density(std::vector<double>& Particle_rho_, std::vector<double>& Particle_x_) {
	double ss = 0.;
	for (int i = 0; i < N_1; i++) {
		double d = 0.5 + Particle_x_.at(i) / h;
		int t = int(d);
		ss = double(d - t);
		Particle_rho_.at(t) = Particle_rho_.at(t) + (1. - ss);
		Particle_rho_.at(t + 1) = Particle_rho_.at(t + 1) + ss;

	}
	Particle_rho_.at(1) = Particle_rho_.at(1) + Particle_rho_.at(0);
	Particle_rho_.at(i_max - 1) = Particle_rho_.at(i_max - 1) + Particle_rho_.at(i_max);

}*/


inline void piksr2(int n, std::vector<double>& arr, std::vector<double>& brr)
{
	int i, j;
	float a, b;
	for (j = 2; j <= n; ++j) {
		a = arr[j];
		b = brr[j];
		i = j - 1;
		while (i > 0 && arr[i] > a) {
			arr[i + 1] = arr[i];
			brr[i + 1] = brr[i];
			--i;
		}
		arr[i + 1] = a;
		brr[i + 1] = b;
	}
}

inline void sort2(int n, std::vector<double>& arr, std::vector<double>& brr)
{
	int  i, ir, j, jstack, k, l, istack[50];
	double a, b, temp;
	jstack = 0;
	l = 1;
	ir = n;
label_1: if (ir - l < 7)
{
	for (j = l + 1; j <= ir; ++j)
	{
		a = arr[j];
		b = brr[j];
		for (i = j - 1; i >= 1; --i)
		{
			if (arr[i] <= a)  goto label_2;
			arr[i + 1] = arr[i];
			brr[i + 1] = brr[i];
		}
		i = l - 1;
	label_2: arr[i + 1] = a;
		brr[i + 1] = b;
	}
	if (jstack == 0) return;
	ir = istack[jstack];
	l = istack[jstack - 1];
	jstack = jstack - 2;
}
else
{
	k = (l + ir) / 2;
	temp = arr[k];
	arr[k] = arr[l + 1];
	arr[l + 1] = temp;
	temp = brr[k];
	brr[k] = brr[l + 1];
	brr[l + 1] = temp;
	if (arr[l] > arr[ir])
	{
		temp = arr[l];
		arr[l] = arr[ir];
		arr[ir] = temp;
		temp = brr[l];
		brr[l] = brr[ir];
		brr[ir] = temp;
	}
	if (arr[l + 1] > arr[ir])
	{
		temp = arr[l + 1];
		arr[l + 1] = arr[ir];
		arr[ir] = temp;
		temp = brr[l + 1];
		brr[l + 1] = brr[ir];
		brr[ir] = temp;
	}
	if (arr[l] > arr[l + 1])
	{
		temp = arr[l];
		arr[l] = arr[l + 1];
		arr[l + 1] = temp;
		temp = brr[l];
		brr[l] = brr[l + 1];
		brr[l + 1] = temp;

	}
	i = l + 1;
	j = ir;
	a = arr[l + 1];
	b = brr[l + 1];
label_3: i = i + 1;
	if (arr[i] < a) goto label_3;
label_4:  j = j - 1;
	if (arr[j] > a) goto label_4;
	if (j < i) goto label_5;
	temp = arr[i];
	arr[i] = arr[j];
	arr[j] = temp;
	temp = brr[i];
	brr[i] = brr[j];
	brr[j] = temp;
	goto label_3;
label_5: arr[l + 1] = arr[j];
	arr[j] = a;
	brr[l + 1] = brr[j];
	brr[j] = b;

	jstack = jstack + 2;

	if (jstack > 50) printf("NSTACK too small in sort");
	if (ir - i + 1 >= j - l)
	{
		istack[jstack] = ir;
		istack[jstack - 1] = i;
		ir = j - 1;
	}
	else
	{
		istack[jstack] = j - 1;
		istack[jstack - 1] = l;
		l = i;
	}
}
goto label_1;
}






inline void Quick_Sort(std::vector<double>& x_, std::vector<double>& v_, int N_1_) {






	// Создаем вектор пар (координата, скорость)
	std::vector<std::pair<double, double>> combined;

	for (int i = 0; i < N_1_; ++i) {
		combined.emplace_back(x_[i], v_[i]);
	}

	// Сортируем вектор пар по первой координате
	std::sort(combined.begin(), combined.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
		return a.first < b.first; // Сравниваем по координате
		});


	// Разделяем отсортированные пары обратно на два массива
	for (int i = 0; i < N_1_; ++i) {
		x_[i] = combined[i].first;
		v_[i] = combined[i].second;
	}

}


inline double Maxwell() {
	double Electron_velocity_min_ = -3.;
	double Electron_velocity_max_ = 3.;
	std::random_device rd;
	std::mt19937 gen(rd());
	while (true) {
		std::uniform_real_distribution<> dist_1(Electron_velocity_min_, Electron_velocity_max_);
		std::uniform_real_distribution<> dist_2(0, 1);
		double p_1 = dist_1(gen);
		double p_2 = dist_2(gen);
		double val = -(p_1 * p_1) / 2.;
		double velocity;
		if (p_2 < exp(val)) {
			return p_1;

		}
	}
}



inline void Many_process_Density(std::vector<double>& x_, std::vector<double>& rho_1_, int myrank_, const std::string& str, int N_1_) {
	// Локальный массив для глобальной плотности (i_max + 1 узлов)
	std::vector<double> rho_(i_max + 1, 0.0);

	// Получаем количество потоков OpenMP
	int num_threads = omp_get_max_threads();

	// Создаем локальные буферы плотности для каждого потока
	std::vector<std::vector<double>> local_rho(num_threads, std::vector<double>(i_max + 1, 0.0));

	// Параллельное вычисление локальных плотностей
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		auto& local = local_rho[tid];

#pragma omp for
		for (int i = 0; i < N_1_; ++i) {
			double d = 0.5 + x_[i] / h;
			int t = static_cast<int>(d);
			double ss = d - t;

			if (t >= 0 && t + 1 <= i_max) {
				local[t] += (1.0 - ss);
				local[t + 1] += ss;
			}
		}
	}

	// Суммируем результаты из всех потоков в общий буфер rho_
#pragma omp parallel for
	for (int t = 0; t <= i_max; ++t) {
		for (int tid = 0; tid < num_threads; ++tid) {
			rho_[t] += local_rho[tid][t];
		}
	}

	// ⛔ ВАЖНО: MPI не потокобезопасен по умолчанию, вызываем Allreduce ВНЕ OpenMP
	MPI_Allreduce(rho_.data(), rho_1_.data(), i_max + 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Граничные условия
	rho_1_[1] += rho_1_[0];
	rho_1_[i_max - 1] += rho_1_[i_max];
	rho_1_[0] = rho_1_[1];
	rho_1_[i_max] = rho_1_[i_max - 1];
}
 



/*
for (int i = 0; i < N_1_; i++) {
		double d = 0.5 + x_[i] / h;
		int t = int(d);
		ss = double(d - t);
		rho_[t] = rho_[t] + (1. - ss);
		rho_[t + 1] = rho_[t + 1] + ss;

	}
*/


inline void Move(std::vector<double>& x_, std::vector<double>& v_, int NPE_, int myrank_, double M, double V_t, double X_0_, double X_1_, int n, double& first_, std::vector<double>& E_) {

	MPI_Status status;
	int i = 0;
	int j = 0;
	int k = j;
	int r_dif = 0;
	int l_dif = 0;
	double buf = 0.;
	int from_left = 0;
	int from_right = 0;
	std::vector<double> r_buf;
	std::vector<double> l_buf;
	std::vector<double> From_left, From_right;



	double EE = 0.;
	for (int i = 0; i < N_1; ++i) {
		double d = x_.at(i) / h;
		int t = int(d);
		double ss = double(d - t);
		EE = E_.at(t) * (1. - ss) + E_.at(t + 1) * ss;

		v_.at(i) = v_.at(i) + tau / M * EE;
		x_.at(i) = x_.at(i) + tau * v_.at(i);
		if (x_.at(i) < 0.) {
			x_.at(i) = -x_.at(i);
			v_.at(i) = -v_.at(i);
		}
		if (x_.at(i) > x_max) {
			x_.at(i) = 2. * x_max - x_.at(i);
			v_.at(i) = -v_.at(i);
		}


	}

	if (myrank_ == 0) {
		//	file_first<<n<<" "<<"initial first "<<first_<<" "<<std::endl;
		first_ = first_ + par_1;
		while (first_ >= par_2) {

			N_1 = N_1 + 1;
			x_.at(N_1 - 1) = first_ - par_2;
			//			file_first<<"first= "<<first_<<" "<<n<<"  x_0="<<x_[0]<<" x_N= "<<x_[N_1-1]<<std::endl;
			//			file_first<<"   "<<std::endl;
			v_.at(N_1 - 1) = u_0 + V_t * Maxwell();

			first_ = first_ - par_2;

		}
	}

	sort2(N_1 - 1, x_, v_);


	if (myrank_ != NPE_ - 1) {
		for (j = N_1 - 1; j >= 0; --j) {
			if (x_[j] < X_1_) {
				k = j + 1;
				break;
			}
		}



		r_dif = N_1 - k;
		if (r_dif != 0) {
			r_buf.resize(2 * r_dif);
			for (i = 0; i < r_dif; ++i) {
				r_buf.at(i) = x_[k + i];
				r_buf.at(i + r_dif) = v_[k + i];
			}
		}

	}


	if (myrank_ != 0) {
		for (j = 0; j < N_1; ++j) {
			if (x_[j] > X_0_) {
				k = j - 1;
				break;
			}
		}


		l_dif = k + 1;
		l_buf.resize(2 * l_dif);
		for (i = 0; i < l_dif; ++i) {
			l_buf[i] = x_[i];
			l_buf[i + l_dif] = v_[i];
		}


	}


	if (myrank_ != NPE_ - 1) {



		MPI_Send(&r_dif, 1, MPI_INT, myrank_ + 1, 0, MPI_COMM_WORLD);
		if (r_dif != 0) {
			MPI_Send(r_buf.data(), 2 * r_dif, MPI_DOUBLE, myrank_ + 1, 0, MPI_COMM_WORLD);

		}
	}


	if (myrank_ != 0) {
		MPI_Recv(&from_left, 1, MPI_INT, myrank_ - 1, 0, MPI_COMM_WORLD, &status);

		if (from_left != 0) {
			From_left.resize(2 * from_left);
			MPI_Recv(From_left.data(), 2 * from_left, MPI_DOUBLE, myrank_ - 1, 0, MPI_COMM_WORLD, &status);

		}
	}

	if (myrank_ != 0) {
		MPI_Send(&l_dif, 1, MPI_INT, myrank_ - 1, 0, MPI_COMM_WORLD);
		if (l_dif != 0) {
			MPI_Send(l_buf.data(), 2 * l_dif, MPI_DOUBLE, myrank_ - 1, 0, MPI_COMM_WORLD);

		}
	}


	if (myrank_ != NPE_ - 1) {
		MPI_Recv(&from_right, 1, MPI_INT, myrank_ + 1, 0, MPI_COMM_WORLD, &status);
		if (from_right != 0) {
			From_right.resize(2 * from_right);
			MPI_Recv(From_right.data(), 2 * from_right, MPI_DOUBLE, myrank_ + 1, 0, MPI_COMM_WORLD, &status);
		}
	}
	for (i = 0; i < N_1 - l_dif - r_dif; ++i) {
		x_[i] = x_[l_dif + i];
		v_[i] = v_[l_dif + i];
	}

	for (i = 0; i < from_left; ++i) {
		x_[i + N_1 - l_dif - r_dif] = From_left[i];
		v_[i + N_1 - l_dif - r_dif] = From_left[i + from_left];
	}

	for (i = 0; i < from_right; ++i) {
		k = i + N_1 - l_dif - r_dif + from_left;
		x_[k] = From_right[i];
		v_[k] = From_right[i + from_right];
	}
	N_1 = N_1 - l_dif - r_dif + from_left + from_right;


}



int main(int argc, char** argv) {
	auto start = std::chrono::high_resolution_clock::now();
	MPI_Init(&argc, &argv);
	int myrank;
	int NPE;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &NPE);
	MPI_Status status;
	X_0 = x_max * myrank / NPE;
	X_1 = x_max / NPE * (myrank + 1);
	std::ofstream file_1;
	std::ofstream file_2;// сам файл
	P = (i_max - 1) / NPE;
	x0 = myrank * P * h;
	N_1 = N_1 / NPE;
	int N_1e = N_1;
	int N_1i = N_1;





	//распределяем частицы
	for (int i = 0; i < N_1; i++) { //? границы 
		x_e.at(i) = X_0 + par_2 * (i + 0.5);
		x_i.at(i) = X_0 + par_2 * (i + 0.5);
		v_e.at(i) = u_0 + Maxwell();
		v_i.at(i) = u_0;
	}


	first_el = x_e[0];
	first_ion = x_i[0];

	//организуем цикл по времени

	for (int n = 0; n < tau_quantity + 1; ++n) {

		//std::cout<<"n "<<n<<std::endl;
				//находим плотности ионов и электронов
		rho_e.assign(rho_e.size(), 0.);
		rho_i.assign(rho_i.size(), 0.);
		Many_process_Density(x_e, rho_e, myrank, "e", N_1e);
		Many_process_Density(x_i, rho_i, myrank, "i", N_1i);

		for (int i = 0; i < i_max + 1; i++) {
			rho_e.at(i) = rho_e.at(i) / N_p;
			rho_i.at(i) = rho_i.at(i) / N_p;
		}


		E.at(i_max - 1) = 0.;
		for (int i = i_max - 2; i >= 0; --i) {
			E.at(i) = E.at(i + 1) - h * (rho_i.at(i + 1) - rho_e.at(i + 1));

		}

		for (int i = 0; i < i_max; ++i) phi.at(i + 1) = phi.at(i) - h * E.at(i);

		N_1 = N_1e;

		Move(x_e, v_e, NPE, myrank, -m_e, 1., X_0, X_1, n, first_el, E);
		N_1e = N_1;
		N_1 = N_1i;
		Move(x_i, v_i, NPE, myrank, Ion_m, 0., X_0, X_1, n, first_ion, E);
		N_1i = N_1;

	/*	if (n % 5000
			== 0) {
			std::ostringstream filename_2;
			std::ostringstream filename_El;
			filename_El << "D://Sp//OMP//Field_And_Density" << myrank << "_" << n << ".txt";
			filename_2 << "D://Sp//OMP//Coordinate_And_Velocity" << myrank << "_" << n << ".txt";
			std::ofstream file_21(filename_2.str());

			std::ofstream file_El(filename_El.str());

			for (int i = 0; i < i_max; ++i) {
				file_El << std::fixed << h * i << " " << E.at(i) << " " << rho_i.at(i) << " " << rho_e.at(i) << " " << phi.at(i) << std::endl;
			}
			file_El.close();

			for (int i = 0; i < N_1; ++i) {
				file_21 << std::fixed << x_e.at(i) << " " << v_e.at(i) << " " << x_i.at(i) << " " << v_i.at(i) << std::endl;
			}
			file_21.close();
		}


		*/
	}
	MPI_Finalize();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;  // Разница
	std::cout << "time c++: " << duration.count() << " sec" << std::endl;
	return 0;
}
