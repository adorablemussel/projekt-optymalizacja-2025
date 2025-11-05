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
		solution X, XB, Xold;

		XB.x = x0;
		XB.y = ff(x0, ud1, ud2);
		solution::f_calls++;

		int fcalls = 1;

		do
		{
			X = HJ_trial(ff, XB, s, ud1, ud2);
			fcalls += (int)X.f_calls;

			if (X.y(0, 0) < XB.y(0, 0))
			{
				do
				{
					Xold = XB;
					XB = X;

					//ruch wzorcowy
					X.x = 2 * XB.x - Xold.x;
					X.y = ff(X.x, ud1, ud2);
					solution::f_calls++;
					fcalls++;

					//sprawdzamy wokół nowego punktu
					X = HJ_trial(ff, X, s, ud1, ud2);
					fcalls++;

					if (fcalls > Nmax)
					{
						XB.flag = -1; // jak przekroczono limit
						return XB;
					}
				} while (X.y(0, 0) < XB.y(0, 0));

				X = XB; // cofnięcie do najlepszego punktu
			}
			else
			{
				s *= alpha; // zmniejszenie kroku
			}

			if (fcalls > Nmax)
			{
				XB.flag = -1;
				return XB;
			}

		} while (s > epsilon);

		XB.flag = 1; // sukces
		Xopt = XB;
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
		//Tu wpisz kod funkcji
		int n = 2; 				//liczba współrzędnych
		matrix X_try = XB.x;
		matrix Y_try;
		
		for (int j = 0; j < n; j++){	//j to aktualna współrzędna czyli 0 = poziom, 1 = pion
			// dla +s
			X_try = XB.x;
			X_try(0, j) += s;
			Y_try = ff(X_try, ud1, ud2);
			XB.f_calls++;

			if(Y_try < XB.y(0,0)){
				XB.x = X_try;
				XB.y = Y_try;
			} 
			else {

				// dla -s
				X_try = XB.x; // reset
				X_try(0, j) -= s;
				Y_try = ff(X_try, ud1, ud2);
				XB.f_calls++;
				
				if(Y_try < XB.y(0,0)){
					XB.x = X_try;
					XB.y = Y_try;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
