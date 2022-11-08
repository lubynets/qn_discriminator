#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <sstream>
#include <string>

template<typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

std::string StringBinNumber(int number) {
  if (number < 10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}

#endif//HELPER_H
