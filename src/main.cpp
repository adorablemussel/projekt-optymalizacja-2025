/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include "../include/opt_alg.h"
#include <ctime>


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
		srand((unsigned)time(nullptr));

		//lab0();
		//lab1();
		lab2();
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
	//funkcja testowa
	std::ofstream Sout1("data/results/symulacja_lab1_testowa.csv");

	double* res = new double[2] { 0.f, 0.f };
	matrix lb(-100), ub(100);
	double x0;
	double d = 5.0;
	double alpha;
	double epsilon = 1e-2;
	double gamma = 1e-4;
	int Nmax = 10000;
	solution wynikFib, wynikLag;

	for(int j=0;j<3;j++){
		alpha = (rand() %100+100)/100.0;
		Sout1 << "\n\n\nWspolczynnik ekspansji: " << alpha << '\n'
		<< "Lp.;x(0);a;b;fcalls_Eksp;x*_Fib;y*_Fib;fcalls_Fib;Minimum lokalne/globalne;x*_Lagr;y*_Lagr;fcalls_Lagr;Minimum lokalne/globalne\n\n";
		
		for (int i = 0; i < 100; i++)
		{
			solution::clear_calls();
			x0 = rand() % 201-100;
			res = expansion(ff1T, x0, d, alpha, Nmax, lb, ub);
			int callsExp = solution::f_calls;
			
			solution::clear_calls();
			wynikFib = fib(ff1T, res[0], res[1], epsilon);
			int callsFib = wynikFib.f_calls;
			
			solution::clear_calls();
			wynikLag = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
			int callsLag = wynikLag.f_calls;
			string czyGlobalne = (wynikFib.x(0) > 53.0) ? "globalne" : "lokalne";

			// cout << i + 1 << ';'
			// << x0 << ';'
			// << double_to_string_comma(res[0]) << ';'
			// << double_to_string_comma(res[1]) << ';'
			// << callsExp<< ';'
			// << wynikFib.x << ' '
			// << wynikFib.y << ' '
			// << callsFib<< ';'
			// << czyGlobalne << ';'
			// << wynikLag.x << ' '
			// << wynikLag.y << ' '
			// << callsLag<< ';'
			// << czyGlobalne << '\n';

			Sout1 << i + 1 << ';'
			<< x0 << ';'
			<< double_to_string_comma(res[0]) << ';'
			<< double_to_string_comma(res[1]) << ';'
			<< callsExp<< ';'
			<< wynikFib.x << ' '
			<< wynikFib.y << ' '
			<< callsFib<< ';'
			<< czyGlobalne << ';'
			<< wynikLag.x << ' '
			<< wynikLag.y << ' '
			<< callsLag<< ';'
			<< czyGlobalne << '\n';
		}
	}
	Sout1.close();

	
	//problem rzeczywisty
	epsilon = 1e-5;
	gamma = 1e-6;
	double a_start = 0.0001;// przedział zrobiony w m^2 dlatego epsilon i gamma takie duże żeby precyzja nie uproszczała za bardzo wyniku
    double b_start = 0.01;
	Nmax = 500;
	solution::clear_calls();
	solution opt = fib(ff1R, a_start, b_start, epsilon);// optymalizacja fibonacim
	cout << opt << endl << endl;
	solution::clear_calls();
	solution opt1 = lag(ff1R, a_start, b_start, epsilon, gamma, Nmax);// optymalizacja lag
	cout << opt1 << endl << endl;
	matrix Y0 = matrix(3,1);
	Y0(0)=5.0;	//objętość zbiornika a
	Y0(1)=1.0;	//objętość zbiornika b
	Y0(2)=20.0;	 //temperatura zbiornika b
	matrix ud1(1, 1);
	ud1(0) = m2d(opt.x);
	cout<<m2d(opt.x)<<endl;
	cout<<m2d(opt.y) + 50.0 << endl;
	matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, ud1);
	ofstream Sout2("data/results/symulacja_lab1_rzeczywisty.csv");		// definiujemy strumieñ do pliku .csv
	Sout2 << "t(s); V(A); V(B); temp(B);\n\n";
	Sout2 << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout2.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();

	ud1(0) = m2d(opt1.x);
	cout<<m2d(opt1.x)<<endl;
	cout<<m2d(opt1.y) + 50.0 << endl;
	matrix* Y1 = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, ud1);
	ofstream Sout3("data/results/symulacja_lab1_rzeczywisty1.csv");		// definiujemy strumieñ do pliku .csv
	Sout3 << "t(s); V(A); V(B); temp(B);\n\n";
	Sout3 << hcat(Y1[0], Y1[1]);								// zapisyjemy wyniki w pliku
	Sout3.close();											// zamykamy strumieñ
	Y1[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y1[1].~matrix();
	
}

void lab2()
{
	// testy
	matrix x = matrix(1, 2, 0.0);
	double x_1 = -1.0;				//pierwsza współrzędna (oś x)
	double x_2 = 0.5;				//druga współrzędna (oś y)
	x(0, 0) = x_1;
	x(0, 1) = x_2;

	cout << "x = " << x << endl;
	cout << "f(x) = " << ff2T(x, matrix(), matrix()) << endl;

	double testS = 0.2;				//testowy krok
	matrix ud1, ud2; //puste

	solution xb;					//x początkowy (z dwoma współrzędnymi i wartoscia)
	xb.x = x;
	xb.y = ff2T(x, ud1, ud2);

	solution testHJ;
	testHJ = HJ_trial(ff2T, xb, testS, ud1, ud2);

	cout << "Wynik HJ_trial:" << endl;
	cout << "x = " << testHJ.x << endl;
	cout << "y = " << testHJ.y << endl;
	cout << "f_calls = " << solution::f_calls << endl;
	// koniec testów


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
