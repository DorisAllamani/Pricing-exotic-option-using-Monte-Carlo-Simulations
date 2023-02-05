#define _GLIBCXX_USE_CXX11_ABI 0/1
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>

using namespace std; //abrivation for using namespace standart

 
class monteCarlo 
{
    public:  // public variables ( private variables, protected variables)
    int  nScenarios, optionType;
    double S_t, r, sigma, strike1, strike2, strike3,  t, T;

    //constructer // the same name as the class// needed for the creation of the object
    public:
    monteCarlo(int nScenarios, double S_t, double T, double t, double r, double sigma, double strike1, double strike2, double strike3, double optionType)
    {
        this-> nScenarios = nScenarios; //Number of simulation
        this->S_t = S_t; // Stock price
        this->T = T; // maturity
        this->t= t; // interval 
        this->sigma = sigma; // volatility
        this->r = r; // interest rate
        this->strike1 = strike1; 
        this->strike2 = strike2;
        this->strike3 = strike3;
        this->optionType = optionType;
    }
   //method returning a pointr arrow
    double  *sample(int nScenarios, double S_t, double T, int t, double r, double sigma)
    {
        double *randomS_T = new double[nScenarios]; // Vector pointer for the result of our different prices

        std::default_random_engine generator(time(0));  // to use the normal_distrib below
       
        for (int i = 0; i < nScenarios; i++) 
        {
            std::normal_distribution<double> varGen1(0, 1); // normal distribution
                randomS_T[i] = S_t *exp((r -  0.5 * sigma * sigma )* (T - t))*exp( sigma * sqrt(T - t) * varGen1(generator));    
            //Option price at different guassian variable generate     
        }    
        return randomS_T;
    }
   
    double *run(int nScenarios, double S_T, double T, int t, double r, double sigma, double strike1, double strike2, double strike3, int optionType)
    {
        double sum = 0, vol = 0;
        double *a = sample(nScenarios, S_t, T, t, r, sigma); //Generate vector price
        double *b = new double[nScenarios]; 
        double *results = new double[2];
        
        for(int i = 1; i < nScenarios; i++)
        {
            if( optionType == 1) //bearSpread
            { 
                b[i] = exp(-r* (T - t)) * (std::max(a[i] - strike1, 0.0) + std::min(strike3 - a[i], 0.0));
                //Option price for a BearSpread
            }
            else if (optionType == 2)
            { //bullSpread
                b[i] = exp(-r* (T - t)) * (std::max(a[i] - strike3, 0.0) + std::min(strike1 - a[i], 0.0));
                // Option price for a BullSpread
            }
            else
            { //butterfly 
                b[i] =  exp(-r* (T - t)) * (std::max(a[i] - strike3, 0.0) + 2 * std::min(strike2 - a[i], 0.0) + std::max(a[i] - strike1, 0.0));
                //Option price for a BUTTERFLY strategy
            }            
        }
        for (int i =0; i < nScenarios; i++)
        {
            sum = sum + b[i];
        }
        results[0] = sum / nScenarios; //Price estimator
        
        for (int i =0; i< nScenarios; i++)
        {
            vol = vol + pow(b[i] - (results[0]), 2); // sigma estimator
        }
        results[1] = sqrt(vol / (nScenarios -1));
        return results;
       }
    
       
    double show(int nScenarios, double S_t, double T, int t, double r, double sigma, double strike1, double strike2, double strike3, int optionType)
    {
        double *a = sample(nScenarios, S_t, T, t, r, sigma);
        double *b = run(nScenarios, S_t, T, t, r, sigma, strike1, strike2, strike3, optionType);
        
        for (int i = 0; i < nScenarios; i++)
        {
            std::cout << "The generated stock price is: " << a[i] << std::endl;
        }
        std::cout << "\n**********************************************"<< std::endl;
        std::cout << "\n The price of the option is : " << b[0] << std::endl;
        std::cout << " The standard deviation is: " << b[1] << std::endl;
        std::cout << " The number of simulations: " << nScenarios << std::endl;
        std::cout << "\n**********************************************"<< std::endl;

        return 0;    
    }   
};

// END OF CLASS

// An approximation to the cumulative distribution function for the standard normal distribution
// Joshi approximation 


double norm_pdf(const double x) {
  return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

double norm_cdf(const double x) {
  double k = 1.0/(1.0 + 0.2316419*x);
  double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

  if (x >= 0.0) {
    return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
  } else {
    return 1.0 - norm_cdf(-x);
  }
}

double ft_greek (std::string s, double S, double K, double r, double sigma, double T)
{
    double d1;
    double d2;
    //double result;

    //By the Black and Schole model

    d1 = (log(S/K) + (r + sigma*sigma/2)*T) / (sigma * sqrt(T)); 
    d2 = d1 - sigma * sqrt(T);

    if (s == "delta")
        return(norm_cdf(d1));
    if (s == "gamma")
        return(norm_pdf(d1) / (S * sigma * sqrt(T)));
    if (s == "theta")
        return( -(S * norm_pdf(d1) * sigma) / (2 * sqrt(T)) - r * K * exp(-r * T) * norm_cdf(d2));
    if (s == "vega")
        return(S *norm_pdf(d1) * sqrt(T));
    if (s== "rho")
        return(K * T * exp(-r * T) * norm_cdf(d2));
    else 
        return(0);
}   



void show_greek(double finale_1,double finale_2,double finale_3,double optionType) // Display the greek according to what is demand
{
    double finale;

    if (optionType <1) 
{
    finale = finale_1 - finale_3;
    std::cout << "The result of computing the greek is: " << finale<< std::endl;
} 
else if (optionType <2) 
{
    finale = finale_3 - finale_1;
    std::cout << "The result of computing the greek is: " << finale<< std::endl;
} 
else 
{
  finale = finale_3 + finale_1 - 2* finale_2;
  std::cout << "The result of computing the greek is: " << finale << std::endl;
}
}

int main()
{
    // Importation of the data
    string name_csv;
    name_csv = "AZN.csv"; //AstraZeneca PLC Stoc price - Nasdaq
    std::cout << "  Hello user, this programs implements the price of 3 different options for the Astra Zeneca. \n";
    std::cout << "\n*********************************************************************************"<< std::endl;
    
    //Creating 2d array from csv file
    vector<vector<string >>  data; // 
    vector<string> row;
    string line, word;
 
    fstream file (name_csv, ios::in); //File Stream open the file 
    
    if(file.is_open())
    {
        while(getline(file, line))
        {
            row.clear();
            stringstream str(line);
            while(getline(str, word, ','))
                row.push_back(word);
                data.push_back(row);
        }
    }
    else
    {
        std::cout<<" Problem opening the file\n";
    }
    // computation with the data of our volatility, the mean and the interest rate
    float sum = 0, mean, stdev;
    
    for(int i = 0; i < 251; i++) 
    {
        sum = sum + std::stof(data[i][0]);
    }
    mean = sum / 251;

    std::cout << "Mean = " << mean;
    std::cout<<"\n";

    for(int m = 0; m < 251; m++) 
    {
        sum = sum + pow((std::stof(data[m][0])-mean), 2); // stoff convert string to float
    }

    stdev = sqrt((sum / (251-1))); //standard deviation given the data
    std::cout << " Standard deviation = " << stdev;
    std::cout<<"\n";

    float bond = 0.03456, inflation = 0.0327, r; 
    //10Y government bonds
    //long-term inflation rate
    r = (1 + bond) / (1 + inflation) - 1; //risk-free interest rate

    std::cout << " Risk-free interest rate = " << r;
    std::cout<<"\n";

    int nScenarios, optionType;
    double S_t = 67; //current price
    double T, t;
    double   sigma = stdev/S_t, strike1 = 72, strike2 = 67, strike3 = 62; // Explained in the PDF
   
    std::cout << "\n Please type 1 for a bear spread call:  ";
    std::cout << "\n Please type 2 for a bull spread call:  ";
    std::cout << "\n Please type 3 for a butterfly call:    ";
    std::cin >> optionType;
    while (optionType != 1 && optionType != 2 && optionType != 3)
    {
        std::cout << "Please type 1 or 2 or 3 to choose the option type:    ";
        std::cin >> optionType;
    }
    std::cout << "Please type the maturity T of the option in years:  ";
    std::cin >> T;
    while (T <= 0 || T > 2.0)
    {
        std::cout << "Please type a positive number and smaller than 2 :    ";
        std::cin >> T;
    }
    std::cout << "Please type the time when you want to price the option:  ";
    std::cin >> t;
    while (t < 0 || t >= T)
    {
        std::cout << "Please type a positive number and smaller than T :    ";
        std::cin >> t;
    }
      
    std::cout << "Please type the number of simulations for the Monte Carlo method of pricing an option:  ";
    std::cin >> nScenarios;

    while (nScenarios < 1000)
    {
        std::cout << "Please type an integer larger than 1000 :    ";
        std::cin >> nScenarios;
    }   
    monteCarlo mC(nScenarios, S_t, T, t, r, sigma, strike1, strike2, strike3, optionType); // Monte Carlo class
    //mC.sample(nScenarios, S_t, T, t, r, sigma);
    //mC.run(nScenarios, S_t, T, t, r, sigma, strike1, strike2, strike3,  optionType); // Simulation 
    mC.show(nScenarios, S_t, T, t, r, sigma, strike1, strike2, strike3, optionType); // Results 
//GREEKS PART
double* S_T;
S_T = mC.sample(nScenarios, S_t, T, t, r, sigma); //Price for each scenario
string greek;
double finale = 0, finale_1 = 0, finale_2 = 0, finale_3 =0 ;


while (greek != "no")
{
    std::cout << "Which Greek do you want to know ? (delta, gamma, vega, theta, rho, no):   ";
    std::cin >> greek;

    for (int i =0; i<nScenarios; i++)
    {
    finale_1 = finale_1 + ft_greek(greek,S_T[i],strike1,r,sigma,T);
    finale_2 = finale_2+ ft_greek(greek,S_T[i],strike2,r,sigma,T);
    finale_3 =finale_3 + ft_greek(greek,S_T[i],strike3,r,sigma,T);
    }
    show_greek(finale_1/nScenarios, finale_2/nScenarios,finale_3/nScenarios,optionType);
    finale_1 = 0;
    finale_2 =0;
    finale_3 = 0;
}

    std::cout << "\n \n \n \n \n\n\n\n\n\n\n\n\n\n********************************************************************************* \nThank you -- EKATARINA & DORIS & INES & KARINE" ;

return 0;
}

