BOOST_ROOT := ../../Boost
BOOST_INC := ${BOOST_ROOT}/include

eval: SplineEvaluation.cpp
	c++ -mavx -O3 -std=c++11 -I${BOOST_ROOT} SplineEvaluation.cpp -o eval
