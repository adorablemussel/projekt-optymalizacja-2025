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
		solution XB, tmp, X(x0);
		do {
			XB.x = X.x;
			XB.fit_fun(ff);
			X = HJ_probuj(ff, XB, s);
			if (X.y < XB.y) {
				do {
					tmp = XB;
					XB = X;
					X.x = 2 * XB.x - tmp.x;
					X.fit_fun(ff);
					X = HJ_probuj(ff, X, s);
					if (solution::f_calls > Nmax) {
						Xopt.flag = 0;
						return Xopt;
					}
				} while (X.y < XB.y);
				X = XB;
			}
			else {
				s = alpha * s;
			}
			if (solution::f_calls > Nmax) {
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

solution HJ_probuj(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod 
		solution X;
		int n = get_len(XB.x);
		matrix E(n, n);
		for (int i = 0; i < n; i++) {
			E(i, i) = 1.0;
		}
		for (int i = 0; i < n; i++) {
			X.x = XB.x + (s * E[i]);
			X.fit_fun(ff);
			if (X.y < XB.y) {
				XB = X;
			}
			else {
				X.x = XB.x - (s * E[i]);
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
		int n = get_len(x0);
		matrix d(n, n);
		for (int i = 0; i < n; i++) {
			d(i, i) = 1.0;
		}
		matrix lambda(n, 1), p(n, 1), s(s0);
		solution XB, X_temp;
		XB.x = x0;
		XB.fit_fun(ff);

		while (true) {
			for (int i = 0; i < n; i++) {
				X_temp.x = XB.x + (s(i) * d[i]);
				X_temp.fit_fun(ff);
				if (X_temp.y(0) < XB.y(0)) {
					XB = X_temp;
					lambda(i) += s(i);
					s(i) *= alpha;
				}
				else {
					p(i) = p(i) + 1;
					s(i) *= -beta;
				}
			}
			bool change = true;
			for (int i = 0; i < n; i++) {
				if (p(i) == 0 || lambda(i) == 0)
				{
					change = false;
					break;
				}

			}
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = lambda(i);

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
				lambda = matrix(n, 1);
				p = matrix(n, 1);
			}
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
		solution XB;
		XB.x = x0;
		XB.fit_fun(ff, ud1, ud2);

		solution XT;
		XT = XB;

		double s = 2.5; //początkowy rozmiar simpleksu
		double alpha = 1.0; //odbicie
		double beta = 0.5; //zwężenie
		double gamma = 2.0; //ekspansja
		double delta = 0.5; //redukcja

		do
		{
			XT.x = XB.x;
			XT = sym_NM(ff, XB.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
			c *= dc;

			if (solution::f_calls > Nmax)
			{
				XT.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

			if (norm(XT.x - XB.x) < epsilon)
				break;

			XB = XT;
		} while (true);

		return XT;
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
        // Pobranie wymiaru problemu
        int dims = get_len(x0);
        
        // Inicjalizacja simpleksu (wektora rozwiązań)
        std::vector<solution> S(dims + 1);
        
        // Punkt startowy
        S[0].x = x0;
        S[0].fit_fun(ff, ud1, ud2);

        // Generowanie pozostałych wierzchołków simpleksu
        matrix basis_vector(dims, 1);
        for (int i = 1; i <= dims; ++i)
        {
            basis_vector = 0.0;
            
            matrix offset(dims, 1);
            offset(i - 1, 0) = 1.0; // Ustawienie 1 na przekątnej/pozycji

			// Skalowanie wektora bazowego do długości boku simpleksu
            S[i].x = S[0].x +  offset * s;
            S[i].fit_fun(ff, ud1, ud2);
        }

        // Zmienne indeksów
        int idx_best = 0;
        int idx_worst = 0;

        while (true)
        {
            // Znalezienie najlepszego i najgorszego punktu
            idx_best = 0;
            idx_worst = 0;
            
            for (size_t i = 1; i < S.size(); ++i)
            {
                if (S[i].y < S[idx_best].y) idx_best = i;
                if (S[i].y > S[idx_worst].y) idx_worst = i;
            }

            // Sprawdzenie warunku stopu (rozmiar simpleksu)
            double max_dist = 0.0;
            for (const auto& point : S)
            {
                double current_dist = norm(S[idx_best].x - point.x);
                if (current_dist > max_dist) max_dist = current_dist;
            }

            if (max_dist < epsilon) break;
            if (solution::f_calls > Nmax)
            {
                S[idx_best].flag = 0;
                throw std::string("Maximum amount of f_calls reached!");
            }

            // Wyznaczenie środka ciężkości (Centroid)
            matrix centroid = S[0].x;
            for(size_t i = 1; i < S.size(); ++i) centroid = centroid + S[i].x;
            
            centroid = centroid - S[idx_worst].x;
            centroid = centroid / (double)(S.size() - 1); // Dzielimy przez n (ilość punktów bez najgorszego)
            centroid = matrix(dims, 1); // Zerowanie
            for(size_t i=0; i<S.size(); ++i) {
                if(i != idx_worst) centroid = centroid + S[i].x;
            }
            centroid = centroid / (double)S.size();


            // Odbicie
            solution P_refl;
            P_refl.x = centroid + alpha * (centroid - S[idx_worst].x);
            P_refl.fit_fun(ff, ud1, ud2);

            if (P_refl.y < S[idx_best].y)
            {
                // Ekspansja
                solution P_exp;
                P_exp.x = centroid + gamma * (P_refl.x - centroid);
                P_exp.fit_fun(ff, ud1, ud2);

                if (P_exp.y < P_refl.y)
                    S[idx_worst] = P_exp;
                else
                    S[idx_worst] = P_refl;
            }
            else
            {
                if (S[idx_best].y <= P_refl.y && P_refl.y < S[idx_worst].y)
                {
                    // Akceptacja odbicia, brak ekspansji
                    S[idx_worst] = P_refl;
                }
                else
                {
                    // Zawężenie
                    solution P_contr;
                    P_contr.x = centroid + beta * (S[idx_worst].x - centroid);
                    P_contr.fit_fun(ff, ud1, ud2);

                    if (P_contr.y >= S[idx_worst].y)
                    {
                        // Redukcja globalna
                        for (size_t i = 0; i < S.size(); ++i)
                        {
                            if (i == idx_best) continue;
                            S[i].x = delta * (S[i].x + S[idx_best].x);
                            S[i].fit_fun(ff, ud1, ud2);
                        }
                    }
                    else
                    {
                        S[idx_worst] = P_contr;
                    }
                }
            }
        }

        return S[idx_best];
    }
    catch (const std::string& msg)
    {
        throw "solution sym_NM(...):\n" + msg;
    }
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
		Xopt.x = x0;
        // pomocnicze zmienne
        matrix x_prev;
        matrix d;
        matrix g;

        while (true)
        {
            // zapamiętaj poprzedni punkt
            x_prev = Xopt.x;
            // sprawdzenie limitu wywołań funkcji
            if (Xopt.f_calls >= Nmax)
			{
                //throw string("Przekroczono maksymalna liczbe wywolan funkcji");
				Xopt.flag = -1;   // brak zbieżności
    			return Xopt;
			}


            // gradient w aktualnym punkcie
            g = Xopt.grad(gf, ud1, ud2);

            // kierunek najszybszego spadku
            d = (-1.0) * g;

            // aktualizacja punktu
            Xopt.x = Xopt.x + h0 * d;

			Xopt.fit_fun(ff, ud1, ud2);

            // sprawdzenie warunku stopu
            if (norm(Xopt.x - x_prev) < epsilon)
                break;
        }

        Xopt.flag = 0;
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
		 // punkt startowy
        Xopt.x = x0;

        matrix x_prev;
        matrix g_prev;
        matrix g_curr;
        matrix d;

        // gradient początkowy
        g_curr = Xopt.grad(gf, ud1, ud2);

        // kierunek początkowy
        d = (-1.0) * g_curr;

        while (true)
        {
            x_prev = Xopt.x;
            g_prev = g_curr;
            // kontrola liczby wywołań funkcji
            if (Xopt.f_calls >= Nmax)
			{
				//throw string("Przekroczono maksymalna liczbe wywolan funkcji");
				Xopt.flag = -1;   // brak zbieżności
    			return Xopt;
			}

            // krok
            Xopt.x = Xopt.x + h0 * d;
			Xopt.fit_fun(ff, ud1, ud2);

            // nowy gradient
            g_curr = Xopt.grad(gf, ud1, ud2);

            // współczynnik beta (Fletcher–Reeves)
            double beta =
                pow(norm(g_curr), 2.0) /
                pow(norm(g_prev), 2.0);

            // nowy kierunek
            d = (-1.0) * g_curr + beta * d;

            // warunek stopu
            if (norm(Xopt.x - x_prev) < epsilon)
                break;
        }

        Xopt.flag = 0;

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
		// punkt startowy
        Xopt.x = x0;

        matrix x_prev;
        matrix g;
        matrix H;
        matrix d;

        while (true)
        {
            // zapamiętaj poprzedni punkt
            x_prev = Xopt.x;

            // kontrola liczby wywołań funkcji
            if (Xopt.f_calls >= Nmax)
			{
				//throw string("Przekroczono maksymalna liczbe wywolan funkcji");
				Xopt.flag = -1;   // brak zbieżności
    			return Xopt;
			}

            // gradient
            g = Xopt.grad(gf, ud1, ud2);

            // hesjan
            H = Xopt.hess(Hf, ud1, ud2);

            // kierunek Newtona: d = - H^{-1} * g
            d = (-1.0) * (inv(H) * g);

            // aktualizacja punktu
            Xopt.x = Xopt.x + h0 * d;
			Xopt.fit_fun(ff, ud1, ud2);


            // warunek stopu
            if (norm(Xopt.x - x_prev) < epsilon)
                break;
        }

        Xopt.flag = 0;
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
		const double alpha = (sqrt(5.0) - 1.0) / 2.0;

        double a_i = a;
        double b_i = b;

        double c = b_i - alpha * (b_i - a_i);
        double d = a_i + alpha * (b_i - a_i);

        int i = 0;

        Xopt.x = matrix(1, 1, c);
        matrix fc = Xopt.fit_fun(ff, ud1, ud2);

        Xopt.x = matrix(1, 1, d);
        matrix fd = Xopt.fit_fun(ff, ud1, ud2);

		while ((b_i - a_i) > epsilon)
        {
            if (Xopt.f_calls >= Nmax)
                throw string("Przekroczono maksymalna liczbe wywolan funkcji");

            if (fc(0) < fd(0))
            {
                b_i = d;
                d = c;
                fd = fc;

                c = b_i - alpha * (b_i - a_i);
                Xopt.x = matrix(1, 1, c);
                fc = Xopt.fit_fun(ff, ud1, ud2);
            }
            else
            {
                a_i = c;
                c = d;
                fc = fd;

                d = a_i + alpha * (b_i - a_i);
                Xopt.x = matrix(1, 1, d);
                fd = Xopt.fit_fun(ff, ud1, ud2);
            }

            i++;
        }

        double h_opt = 0.5 * (a_i + b_i);

        Xopt.x = matrix(1, 1, h_opt);
        Xopt.flag = 0;

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

		int n = get_len(x0);
		matrix D = ident_mat(n), A(n, 2);
		solution X, P, h;
		X.x = x0;
		double* ab;

		while (true)
		{
			P = X;

			for (int i = 0; i < n; ++i)
			{
				A.set_col(P.x, 0);
				A.set_col(D[i], 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
				P.x = P.x + h.x * D[i];
			}

			if (norm(P.x - X.x) < epsilon)
			{
				Xopt = X;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				break;
			}

			if (solution::f_calls > Nmax)
			{
				Xopt = X;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}

			for (int i = 0; i < n - 1; ++i)
				D.set_col(D[i + 1], i);

			D.set_col(P.x - X.x, n - 1);
			A.set_col(P.x, 0);
			A.set_col(D[n - 1], 1);
			ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
			h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
			X.x = P.x + h.x * D[n - 1];
		}

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
