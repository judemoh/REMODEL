#include <iostream>
using namespace std;

int main() {
  int sales = 95000;
  double state_tax = (0.04*sales);
  double county_tax = (0.02*sales);
  cout << "State tax: " << state_tax << endl
      << "County tax: " << county_tax << endl
      << "Sales" << sales << endl;

  return 0;
}