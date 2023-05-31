/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "timeseries.hpp"
#include "subann.hpp"

namespace py = pybind11;


PYBIND11_MODULE(subann, m) {
    py::class_<TimeSeries>(m, "TimeSeries")
        .def(py::init<std::vector<std::vector<double>>>())
        .def(py::init<std::vector<std::vector<double>>, std::vector<int>>())
        .def("__getitem__", &TimeSeries::operator[])
        .def("values", &TimeSeries::values)
        .def("get_dim", &TimeSeries::get_dim)
        .def("size", &TimeSeries::size)
        .def("__iter__", [](const TimeSeries& ts) { return py::make_iterator(ts.begin(), ts.end()); },
            py::keep_alive<0, 1>())
        .def("projectAxis", &TimeSeries::projectAxis);
    ;
    py::class_<SubAnnDS>(m, "SubAnnDS")
        .def(py::init<TimeSeries&, double, double>())
        .def(py::init<TimeSeries&, double, double, double, bool>())
        .def("range_query", &SubAnnDS::range_query)
        
        .def("stats", &SubAnnDS::stats)
    ;


}
