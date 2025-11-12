#include"../include/opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo³añ funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };

		solution x1(x0);
		solution x2(x0 + d);

		x1.fit_fun(ff, ud1, ud2);
		x2.fit_fun(ff, ud1, ud2);

		if (solution::f_calls >= Nmax)
		{
			p[0] = m2d(x1.x);
			p[1] = m2d(x2.x);
			if (p[0] > p[1]) std::swap(p[0], p[1]);
			return p;
		}

		if (m2d(x1.y) <= m2d(x2.y))
		{
			d = -d;
			std::swap(x1, x2);
		}

		while (true)
		{
			solution x3(m2d(x2.x) + d);
			x3.fit_fun(ff, ud1, ud2);

			//Sprawdzenie, czy "przeskoczyliśmy" minimum
			if (m2d(x3.y) > m2d(x2.y))
			{
				// Znaleziono przedział. Minimum jest między x1 a x3.
				p[0] = m2d(x1.x);
				p[1] = m2d(x3.x);
				// Porządkowanie wyniku, aby p[0] < p[1]
				if (p[0] > p[1]) std::swap(p[0], p[1]);
				break;
			}

			// Sprawdzenie warunku stopu
			if (solution::f_calls >= Nmax)
			{
				// Zwracamy najlepszy znaleziony dotąd przedział
				p[0] = m2d(x1.x);
				p[1] = m2d(x2.x);
				if (p[0] > p[1]) std::swap(p[0], p[1]);
				break;
			}

			// Przejście do następnej iteracji
			x1 = x2;
			x2 = x3;
			d *= alpha; // Zwiększenie kroku
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution A(a), B(b), C, D;
		vector<double> fib_seq = {0.0, 1.0};
		while (fib_seq.back() < (B.x(0) - A.x(0)) / epsilon) {
			fib_seq.push_back(fib_seq[fib_seq.size() - 1] + fib_seq[fib_seq.size() - 2]);
		}
		int k = fib_seq.size()-1;
		//cout << "Liczba iteracji: " << k <<" "<<fib_seq.back()<<" "<<(B.x(0) - A.x(0)) / epsilon<< endl;
		C.x(0) = B.x(0) - (fib_seq[k - 1] / fib_seq[k]) * (B.x(0) - A.x(0));
		D.x(0) = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun(ff);
		D.fit_fun(ff);

		for (int i = 0; i < k - 3; i++) {
			if (C.y(0) < D.y(0)) {
				B.x = D.x;
			}
			else {
				A.x = C.x;
			}
			C.x(0) = B.x(0) - (fib_seq[k - i - 2] / fib_seq[k - i - 1]) * (B.x(0) - A.x(0));
			D.x(0) = A.x(0) + B.x(0) - C.x(0);
			C.fit_fun(ff);
			D.fit_fun(ff);
		}

		Xopt.x = C.x(0);
		Xopt.y = C.y(0);
		Xopt.flag = 1;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution A(a), B(b), C, D, D1;
		C.x = (a + b) / 2;
		A.fit_fun(ff);
		B.fit_fun(ff);
		C.fit_fun(ff);
		double l, m;
		int i = 0;
		while (true) {
			l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
			m = (A.y(0) * (B.x(0) - C.x(0))) + (B.y(0) * (C.x(0) - A.x(0))) + (C.y(0) * (A.x(0) - B.x(0)));
			if (m <= 0) {
				Xopt.flag = 0;
				break;
			}
			D1.x = D.x;
			D.x = l / (2 * m);
			D.fit_fun(ff);
			if (A.x(0) < D.x(0) && D.x(0) < C.x(0))
			{
				if (D.y(0) < C.y(0))
				{
					B.x = C.x;
					C.x = D.x;
					B.fit_fun(ff);
					C.fit_fun(ff);
				}
				else
					A.x = D.x;
				A.fit_fun(ff);
			}
			else if (C.x(0) < D.x(0) && D.x(0) < B.x(0))
			{
				if (D.y(0) < C.y(0))
				{
					A.x = C.x;
					C.x = D.x;
					A.fit_fun(ff);
					C.fit_fun(ff);
				}
				else
					B.x = D.x;
				B.fit_fun(ff);
			}
			else {
				std::cerr << "Punkt D jest poza zakresem." << std::endl;
				Xopt.flag = 0;
				break;
			}
			i++;
			if (i > Nmax) {
				std::cerr << "Osiagnieto limit iteracji petli." << std::endl;
				Xopt.flag = 0;
				break;
			}
			if (B.x(0) - A.x(0) < epsilon || abs(D.x(0) - D1.x(0)) <= gamma)
			{
				break;
			}
		}
		D.fit_fun(ff);
		Xopt.x = D.x(0);
		Xopt.y = D.y(0);
		Xopt.flag = 1;
		return Xopt;

	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}


solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
		Xopt.clear_calls();
		solution XB, XB_temp, X(x0);
		do {
			XB.x = X.x;
			XB.fit_fun(ff);
			X = HJ_trial(ff, XB, s);
			if (X.y < XB.y) {
				do {
					XB_temp = XB;
					XB = X;
					X.x = 2 * XB.x - XB_temp.x;
					X.fit_fun(ff);
					X = HJ_trial(ff, X, s);
					if (X.f_calls > Nmax) {
						Xopt.flag = 0;
						return Xopt;
					}
				} while (X.y < XB.y);
				X = XB;
			}
			else {
				s = alpha * s;
			}
			if (X.f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
		} while (s > epsilon);

		Xopt.x = XB.x;
		Xopt.y = XB.y;
		Xopt.flag = 1;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod 
		solution X;
		int n = get_len(XB.x);
		matrix E(n, n);
		for (int j = 0; j < n; j++) {
			E(j, j) = 1.0;
		}
		for (int j = 0; j < n; j++) {
			X.x = XB.x + (s * E[j]);
			X.fit_fun(ff);
			if (X.y < XB.y) {
				XB = X;
			}
			else {
				X.x = XB.x - (s * E[j]);
				X.fit_fun(ff);
				if (X.y < XB.y) {
					XB = X;
				}
			}
		}
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji
		solution::clear_calls();
		//int i = 0;
		int n = get_len(x0);
		matrix d(n, n);
		for (int j = 0; j < n; j++) {
			d(j, j) = 1.0;
		}
		matrix lamda(n, 1), p(n, 1), s(s0);
		solution XB, X_temp;
		XB.x = x0;
		XB.fit_fun(ff);

		while (true) {

			//zacznijmy szukać rozwiązania zgodnie z zadaną bazą 
			for (int j = 0; j < n; j++) {
				X_temp.x = XB.x + (s(j) * d[j]);
				X_temp.fit_fun(ff);
				if (X_temp.y(0) < XB.y(0)) {
					XB = X_temp;
					lamda(j) += s(j);
					s(j) *= alpha;
				}
				else {
					p(j) = p(j) + 1;
					s(j) *= -beta;
				}
			}

			//czy któryś z kroków pogorszył sytuacje?
			bool change = true;
			for (int j = 0; j < n; j++) {
				if (p(j) == 0 || lamda(j) == 0)
				{
					change = false;
					break;
				}

			}

			//zmiana bazy jeżeli jest to konieczne
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = lamda(i);

				Q = d * Q;
				v = Q[0] / norm(Q[0]);
				d.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + (trans(Q[i]) * d[j]) * d[j];
					v = Q[i] - temp;
					d.set_col(v, i);
				}
				s = s0;
				lamda = matrix(n, 1);
				p = matrix(n, 1);
			}

			//czy wynik jest wystarczająco dokładny?
			double max_s = abs(s(0));
			for (int i = 1; i < n; ++i) {
				if (max_s < abs(s(i))) {
					max_s = abs(s(i));
				}
			}
			if (max_s < epsilon || solution::f_calls > Nmax) {
				return XB;
			}
		}
		

	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}



solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
