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
#include <clocale>

class MySeparator : public std::numpunct<char> {
protected:
    char do_decimal_point() const override {
        return ','; // TU ZMIENIASZ: wpisz '.' dla kropki lub ',' dla przecinka
    }
};

std::locale my_loc(std::locale(), new MySeparator());

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
    	// Narzucamy ją na std::cout
		std::locale::global(my_loc);
    	std::cout.imbue(my_loc);
		//lab0();
		//lab1();
		//lab2();
		lab3();
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

	// testowa funkcja celu
	srand(time(NULL));
	std::ofstream Sout("symulacja_lab3_testowa.csv");

	matrix X;
	double step = 0.0001, alpha = 0.8, beta = 0.1, epsilon = 0.0001;
	double a, b;
	int Nmax = 1000;

	// zadanie teoretyczne
	Sout << "Lp." << ";" <<"x1(0)" << ";" << "x2(0)" << ";" << "x1* (hooke)" << ";" << "x2* (hooke)" << ";" << "y (hooke)" << ";" << "f_calls" << ";" << "x1* (rosen)" << ";" << "x2* (rosen)" << ";" << "y (rosen)" << ";" << "f_calls" << ";";
	for (int i = 1; i <= 100; i++)
	{
		a = ((rand() % 200) / 100.0) - 1;
		b = ((rand() % 200) / 100.0) - 1;
		alpha = 0.8;
		X = matrix(2, new double[2] {a, b});
		solution hooke = HJ(ff2T, X, step, alpha, epsilon, Nmax);
		Sout << i << ";" << a << ";" << b << ";" <<hooke.x(0) << ";" << hooke.x(1) << ";" << hooke.y << ";" << solution::f_calls << ";";
		//cout << hooke;

		
		alpha = 1.8;
		matrix Step = matrix(2, new double[2] { step, step});
		solution rosen = Rosen(ff2T, X, Step, alpha, beta, epsilon, Nmax);
		Sout << rosen.x(0) << ";" << rosen.x(1) << ";"<< rosen.y << ";" << solution::f_calls << "\n";
		//cout << rosen;
	}

	//problem rzeczywisty

	matrix x(2, 1, 5);
	cout << ff2R(x);
	X = matrix(2, new double[2] {5, 5});
	solution wynikHJ = HJ(ff2R, X, step, alpha, epsilon, Nmax);
	cout << wynikHJ;

	alpha = 1.8;
	matrix Step = matrix(2, new double[2] { step, step});
	solution wynikR = Rosen(ff2R, X, Step, alpha, beta, epsilon, Nmax);
	cout << wynikR;

	//symulacja 
	double t0 = 0.0;
	double tend = 100.0;
	double dt = 0.1;

	matrix Y0(2, 1); 
	Y0(0) = 0.0;     
	Y0(1) = 0.0;     

	matrix ud3(2, 1);
	ud3(0) = 3.14;   
	ud3(1) = 0.0;    

	matrix ud4(2, 1);
	ud4(0) = wynikHJ.x(0);     
	//ud4(0) = wynikR.x(0);     
	ud4(1) = wynikHJ.x(1);     
	//ud4(1) = wynikR.x(1);     

	matrix* result = solve_ode(df2, t0, dt, tend, Y0, ud3, ud4);

	int n = get_len(result[0]);
	cout << "Czas\tKat\tPredkosc katowa" << endl;
	for (int i = 0; i < n; i++) {
		cout << result[0](i) << "\t"      
			<< result[1](i, 0) << "\t"  
			<< result[1](i, 1) << endl; 
	}

	ofstream file("symulacja_lab2_rzeczywisty.csv");
	file << "Czas;Kat;Predkosc katowa\n";
	for (int i = 0; i < n; ++i) {
		file << result[0](i) << ";" << result[1](i, 0) << ";" << result[1](i, 1) << "\n";
	}
	file.close();

}

void lab3()
{
	srand(time(NULL));
	// funkcja testowa
	double epsilon = 1E-3;
	int Nmax = 10000;
	double c_in = 100;
	double dc_in = 0.2;
	double c_out = 1.0;
	double dc_out = 1.5;

	ofstream Sout("symulacja_lab3_testowa.csv");
	Sout << "x0_1;x0_2;x1_out;x2_out;norm_out;y_out;f_calls_out;x1_in;x2_in;norm_in;y_in;f_calls_in\n";
	
	solution testowa;
	matrix a = matrix(4.0);
	matrix test_x0{};

	for (int i = 0; i < 3; ++i) {
		if (i == 0)
			a = matrix(4.0);
		else if (i == 1)
			a = matrix(4.4934);
		else
			a = matrix(5.0);

		for (int j = 0; j < 100; ++j) {
			double x0_1 = 1.5 + (double)rand() / (double)RAND_MAX * (5.5-1.5);
            double x0_2 = 1.5 + (double)rand() / (double)RAND_MAX * (5.5-1.5);

			test_x0 = matrix(2, new double[2] {x0_1, x0_2});
			Sout << j+1 <<" ;"
			<<test_x0(0) << " ;"
			<< test_x0(1) << " ;";

			// Zewnętrzne rozwiązanie
			testowa = pen(ff3T_out, test_x0, c_out, dc_out, epsilon, Nmax, a);
			//cout << testowa;
			Sout << testowa.x(0)<< " ;"
			<< testowa.x(1) << " ;"
			<< sqrt(pow(testowa.x(0), 2) + pow(testowa.x(1), 2))<< " ;"
			<< testowa.y << " "
			<< testowa.f_calls << " ;";
			solution::clear_calls();

			// Wewnętrzne rozwiązanie
			testowa = pen(ff3T_in, test_x0, c_in, dc_in, epsilon, Nmax, a);
			//cout << testowa;
			Sout << testowa.x(0) << " ;" 
			<< testowa.x(1) << " ;"
			<< sqrt(pow(testowa.x(0), 2) + pow(testowa.x(1), 2)) << " ;"
			<< testowa.y << " "
			<< testowa.f_calls << " \n";
			solution::clear_calls();
		}
	}
	
	Sout.close();

	std::ofstream Sout2("symulacja_lab3_rzeczywista.csv");
	// Nagłówki w pliku CSV
	Sout2 << "time;x_position;y_position\n";

	//Dane zadania
	matrix ud1 = matrix(5, new double[5] {
		0.47, //C
		1.2, //rho
		0.12, //r
		0.6, //m
		9.81 //g
});

	//Początkowe wartości szukania minimum
	matrix x0 = matrix(2, new double[2] {5.0, -5.0});
	
	solution opt = pen(ff3R, x0, c_out, dc_out, epsilon, Nmax, ud1);
	//v0x(0) = 5.0, w(0) = -5.0
	//v0x*, w*, xend*, fcalls - wynik w cout
	cout << opt << "\n";

	matrix Y0(4, new double[4] {0.0, opt.x(0), 100, 0});
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, opt.x(1));

	for (int i = 0; i < get_len(Y[0]); ++i) {
		Sout2 << Y[0](i) << " ;"  // t
			<< Y[1](i, 0) << " ;"  // x
			<< Y[1](i, 2) << " \n";  // y
	}
	Sout2.close();
	delete[] Y;
	

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
