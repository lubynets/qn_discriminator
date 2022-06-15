#ifndef HELPER_H
#define HELPER_H

#include <sstream>
#include <iostream>
#include <string>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

#endif  //HELPER_H