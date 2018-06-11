#include "MSc_Math_Finance_2018.hpp"
#include <iostream>
#include <string>
using namespace std;

namespace MSc_Math_Finance_2018
{
  int testProject(void)
  {
    string home_location=getenv("HOME");
    int exit_status = 0;
    cout << " My new project:: MSc_Math_Finance_2018 " << endl;
    cout << " Home location " << home_location << endl;
    return exit_status;
  }
}

