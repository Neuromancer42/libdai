%module LibDAISWIGInterface
%include "std_vector.i"
%include "std_string.i"
%include "std_except.i"
namespace std {
        %template(DoubleVector) vector<double>;
}

%{
#include "libdai-swig-wrapper.h"
#include <iostream>
#include <map>
#include <string>
#include <vector>
%}

%include "libdai-swig-wrapper.h"
