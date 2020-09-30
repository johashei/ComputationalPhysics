#include "catch.hpp"
#include "EigenValueSolver.hpp"

// Tests for Jacobi_rotation:

TEST_CASE("Testing max max_offdiag"){
  int k,l;
  double maxnondiag;
  // Create a simple symmetric matrix:
  dmat A = {{-3,  7,  13  },{7,   113,-23 },{13,  -23,  5 }};
  // Initialize solver:
  Jacobi_rotation max_offdiag_tester;
  max_offdiag_tester.init(A,0,0); // tolerance and maxiter are not used here
  max_offdiag_tester.max_offdiag();
  // Get results:
  k = max_offdiag_tester.m_k;
  l = max_offdiag_tester.m_l;
  maxnondiag = max_offdiag_tester.m_maxnondiag;
  REQUIRE(k==1);
  REQUIRE(l==2);
  REQUIRE(maxnondiag==23);
}


TEST_CASE("Testing eigenvalues"){
  dmat A = {{ 32, -16,  0   },
            {-16, 32,   -16 },
            { 0,  -16,  32  }};
  Jacobi_rotation eigval_tester;
  eigval_tester.init(A,1e-8,10);
  eigval_tester.solve();
  vec eigvals = eigval_tester.Get_eigvals(0,2);
  REQUIRE(eigvals(2)==Approx(54.627417));
  REQUIRE(eigvals(0)==Approx(9.372583));
  REQUIRE(eigvals(1)==Approx(32.000000));
}

TEST_CASE("Testing eigenvectors"){
  dmat R;
  dmat A = {{ 32, -16,  0   },
            {-16, 32,   -16 },
            { 0,  -16,  32  }};
  Jacobi_rotation eigvec_tester;
  eigvec_tester.init(A,1e-8,10);
  eigvec_tester.solve();
  R = eigvec_tester.m_R;
  REQUIRE(R(0,0)==Approx(0.5));
  REQUIRE(R(1,0)==Approx(-0.707106781));
  REQUIRE(R(2,0)==Approx(0.5));
  REQUIRE(R(0,1)==Approx(0.5));
  REQUIRE(R(1,1)==Approx(0.707106781));
  REQUIRE(R(2,1)==Approx(0.5));
  REQUIRE(R(0,2)==Approx(-0.707106781));
  REQUIRE(R(1,2)==Approx(0.0));
  REQUIRE(R(2,2)==Approx(0.707106781));
}
