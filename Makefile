BOOST_ROOT := ../../goost
BOOST_INC := ${BOOST_ROOT}/include

eval: SplineEvaluation.cpp
	c++ -mavx -O3 -I${BOOST_ROOT} SplineEvaluation.cpp -o eval

