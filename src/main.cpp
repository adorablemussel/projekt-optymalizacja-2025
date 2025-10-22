/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include "opt_alg.h"


std::string double_to_string_comma(double x, int precision = 6) {
    std::ostringstream oss;
    oss << x;
    std::string s = oss.str();
    for (auto& c : s) {
        if (c == '.') c = ',';
    }
    return s;
}

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		// lab0();
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok³adnoœæ
	int Nmax = 10000;										// maksymalna liczba wywo³añ funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dok³adne rozwi¹zanie optymalne
	solution opt;											// rozwi¹zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dok³adnoœæ
	epsilon = 1e-2;											// maksymalna liczba wywo³añ funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad³a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	ofstream Sout("data/results/symulacja_lab0.csv");					// definiujemy strumieñ do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
}

void lab1()
{
	std::ofstream Sout("data/results/symulacja_lab1.csv");
	//zadanie teoretyczne

	double* res = new double[2] { 0.f, 0.f };
	matrix lb(-100), ub(100);
	double x0, d = 5, alpha = 1.5;
	double epsilon = 1e-2;
	double gamma = 1e-6;
	int Nmax = 10000;
	solution wynik,wynik2;
	for(int j=0;j<3;j++){
		alpha = (rand() % 200+100)/100.0;
		Sout << "\n\n\nWspolczynnik ekspansji: " << alpha << '\n'
		<< "Lp.;x(0);a;b;fcalls_Eksp;x*_Fib;y*_Fib;fcalls_Fib;Minimum lokalne/globalne;x*_Lagr;y*_Lagr;fcalls_Lagr;Minimum lokalne/globalne\n\n";
		for (int i = 0; i < 100; i++)
		{
			solution::clear_calls();
			x0 = rand() % 200-99;
			res = expansion(ff1T, x0, d, alpha, Nmax,lb,ub);
			int callsExp = solution::f_calls;
			solution::clear_calls();
			wynik = fib(ff1T, res[0], res[1], epsilon);
			int callsFib = wynik.f_calls;
			solution::clear_calls();
			wynik2 = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
			int callsLag = wynik2.f_calls;
			string czyGlobalne = (wynik.x(0) > 53.0) ? "globalne" : "lokalne";
			cout << i + 1 << ';'
			<< x0 << ';'
			<< double_to_string_comma(res[0]) << ';'
			<< double_to_string_comma(res[1]) << ';'
			<< callsExp<< ';'
			<< wynik.x << ' '
			<< wynik.y<< ' '
			<< callsFib<< ';'
			<< czyGlobalne << ';'
			<< wynik2.x << ' '
			<< wynik2.y << ' '
			<< callsLag<< ';'
			<< czyGlobalne << '\n';

			Sout << i + 1 << ';'
			<< x0 << ';'
			<< double_to_string_comma(res[0]) << ';'
			<< double_to_string_comma(res[1]) << ';'
			<< callsExp<< ';'
			<< wynik.x << ' '
			<< wynik.y<< ' '
			<< callsFib<< ';'
			<< czyGlobalne << ';'
			<< wynik2.x << ' '
			<< wynik2.y << ' '
			<< callsLag<< ';'
			<< czyGlobalne << '\n';
		}
	}
	Sout.close();

	//Funkcja testowa
	
	double* p = new double[2];
	p = expansion(ff1T, 53.0, 1.0, 2.0, Nmax, lb, ub);
	solution opt;
	opt = fib(ff1T, p[0],p[1], epsilon, lb, ub);
	cout << opt << endl << endl;
	matrix Y0 = matrix(3,1);
	Y0(0)=5.0;	//objętość zbiornika a
	Y0(1)=1.0;	//objętość zbiornika b
	Y0(2)=20.0;	//temperatura zbiornika b
	matrix ud1(1, 1), ud2;
	ud1(0) = 0.002;
	matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, ud1, ud2);
	ofstream Sout2("data/results/symulacja_lab1.csv");		// definiujemy strumieñ do pliku .csv
	Sout2 << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout2.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
	solution xd=0.005;

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
