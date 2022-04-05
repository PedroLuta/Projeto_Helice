#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace::std;



vector<double> linspace(double start, double end, int num)
{

  std::vector<double> linspaced;

  if (num == 0) return linspaced; 
  if (num == 1) 
    {
      linspaced.push_back(end);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

void print_vector(std::vector<double> vec)
{
  std::cout << "size: " << vec.size() << std::endl;
  for (double d : vec)
    std::cout << d << " ";
  std::cout << std::endl;
}







bool check_valid(double h, double s, double m, double aa, double bb, double cc, double dd, double ee)
{
    double f_prime_h = 4*aa*(pow(h, 3)) + 3*bb*(pow(h, 2)) + 2*cc*h + dd;
    double f_prime_s = 4*aa*(pow(s, 3)) + 3*bb*(pow(s, 2)) + 2*cc*s + dd;
    double f_prime_1 = 4*aa*(pow(1, 3)) + 3*bb*(pow(1, 2)) + 2*cc*1 + dd;
    double f_2prime_m = 12*aa*(pow(m, 2)) + 6*bb*m + 2*cc;
    double f_2prime_s = 12*aa*(pow(s, 2)) + 6*bb*s + 2*cc;
    double f_2prime_1 = 12*aa*(pow(1, 2)) + 6*bb*1 + 2*cc;
    double f_h = aa*(pow(h, 4)) + bb*(pow(h, 3)) + cc*(pow(h, 2)) + dd*h + ee;

    if ((f_prime_h < 0) || (f_prime_s > 0) || (f_prime_1 > 0) || (f_2prime_m > 0) || (f_2prime_s > 0) || (f_2prime_1 > 0) || (f_h < 0)) return false;
    else return true;
}



bool check_integral(double aa, double bb, double cc, double dd, double ee, double h, double wanted, double tol, double &integral_placeholder)
{
    double integral = aa*(1 - (pow(h,5)))/5 + bb*(1 - (pow(h,4)))/4 + cc*(1 - (pow(h,3)))/3 + dd*(1 - (pow(h,2)))/2 + ee*(1 - h);
    //cout << integral << endl;
    integral_placeholder = integral;
    if (abs(wanted - integral) < tol) return true;
    else return false;
}


void solve_linear_theory1(double h, double m, double s, double Croot, double Cmax, double Cat, double Ctip, double &aa, double &bb, double &cc, double &dd, double &ee )
{
    double AM[][5] = {{pow(h,4), pow(h,3), pow(h,2), h, 1.0}, {pow(m,4), pow(m,3), pow(m,2) , m, 1.0}, {pow(s,4), pow(s,3), pow(s,2), s, 1.0}, {1.0, 1.0, 1.0, 1.0, 1.0}, {4*(pow(m,3)), 3*(pow(m,2)), 2*m, 1.0, 0.0}};
    double BM[5] = {Croot, Cmax, Cat, Ctip, 0};
    double fdScaler, crScaler;

    for(int fd = 0; fd < 5; fd++)
    {
        fdScaler = AM[fd][fd];
        if (fdScaler == 0) continue;
        for(int j = 0; j < 5; j++) AM[fd][j] = AM[fd][j]/fdScaler;
        BM[fd] = BM[fd]/fdScaler;

        for(int i = 0; i < 5; i++)
        {
            if(i == fd) continue;
            crScaler = AM[i][fd];
            for(int j = 0; j < 5; j++) AM[i][j] = AM[i][j] - crScaler*AM[fd][j];
            BM[i] = BM[i] - crScaler*BM[fd];
        }
    }
    bool check = check_valid(h, s, m, BM[0], BM[1], BM[2], BM[3], BM[4]);
    if(check){
        aa = BM[0];
        bb = BM[1];
        cc = BM[2];
        dd = BM[3];
        ee = BM[4];
    } 
    else aa = 1337;
}


double function(double aa, double bb, double cc, double dd, double ee, double xx)
{
    return aa*(pow(xx,4)) + bb*(pow(xx,3)) + cc*(pow(xx,2)) + dd*xx + ee;
}









int main(){
    ofstream outFile;

    double a, b, c, d, e;
    double solutions;
    bool integral_check, check; //, skip
    int ii, counter_skip;
    double Area_Wanted = 0.4;
    double Tolerance = 0.0005;
    double x_disc[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    int counter = 0;
    double array_current[] = {0, 0, 0, 0, 0};
    double array_remember[] = {0, 0, 0, 0, 0};
    double integral_placeholder;

    //vector<double> h = linspace(0.05, 0.25, 41); //44
    vector<double> m;
    vector<double> s;
    vector<double> Cmax = linspace(0.3, 0.6, 31); //39
    vector<double> Croot;
    vector<double> Cat;
    vector<double> Ctip;

    vector<double>::iterator it1, it2, it3, it4, it5, it6, it7;

    outFile.open("output_cpp.txt");
    for(int i = 0; i < 6; i++) outFile << x_disc[i] << " ";
    outFile << endl << endl;

    double h = 0.15;

    //for(it1 = h.begin(); it1 != h.end(); ++it1)
    //{
        //cout << *it1 << endl;
        m = linspace(h, 1, 86); //37
        for(it2 = m.begin(); it2 != m.end(); ++it2)
        {
            cout << "CHANGING M " << *it2 << endl;
            s = linspace(*it2, 1, max(1, min(49, (int)floor((1 - *it2)/0.005)))); //27
            for(it3 = s.begin(); it3 != s.end(); ++it3)
            {
                for(it4 = Cmax.begin(); it4 != Cmax.end(); ++it4)
                {
                    Croot = linspace(0.01, *it4, max(1, min(43, (int)floor((*it4 - 0.01)/0.005)))); //17
                    for(it5 = Croot.begin(); it5 != Croot.end(); ++it5)
                    {
                        Cat = linspace(*it5, *it4, max(1, min(41, (int)floor((*it4 - *it5)/0.005)))); //15
                        for(it6 = Cat.begin(); it6 != Cat.end(); ++it6)
                        {
                            Ctip = linspace(0, *it6, max(1, min(39, (int)floor((*it6)/0.005)))); //13
                            for(it7 = Ctip.begin(); it7 != Ctip.end(); ++it7)
                            {
                                solve_linear_theory1(h, *it2, *it3, *it4, *it5, *it6, *it7, a, b, c, d, e);
                                if(a == 1337) continue;
                                integral_check = false;
                                integral_check = check_integral(a, b, c, d, e, h, Area_Wanted, Tolerance, integral_placeholder); 
                                if (integral_check)
                                {
                                    for(ii = 0; ii < 6; ii++) array_current[ii] = function(a, b, c, d, e, x_disc[ii]);
                                    //skip = false;
                                    counter_skip = 0;
                                    for(ii = 0; ii < 6; ii++)
                                    {
                                        if(abs(array_current[ii] - array_remember[ii]) < 0.005) counter_skip++;
                                    }
                                    if(counter_skip == 6) continue;
                                    outFile << a << " " << b << " " << c << " " << d << " " << e << "     " << integral_placeholder << endl << endl; //" h:" << h <<
                                    //for(ii = 0; ii < 6; ii++) outFile << function(a, b, c, d, e, x_disc[ii]) << " ";
                                    for(ii = 0; ii < 6; ii++) array_remember[ii] = array_current[ii];
                                    //outFile << endl << endl;
                                    counter++;
                                    cout << "distributions found: " << counter << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    //}

outFile.close();


    return 0;
}