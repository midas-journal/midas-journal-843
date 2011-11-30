cmake_minimum_required(VERSION 2.8)

project( PointInCircle )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

# create the executable that will be used by all tests
add_executable( IsPointInCircle 
  IsPointInCircle.cxx
  predicates.c
  )
target_link_libraries( IsPointInCircle
  ${ITK_LIBRARIES}
  )

enable_testing()

# add tests to illustrate the different cases

# in test one we take one point of the square and move it AWAY from the center
# of the circle. The point is outside, the test fail as expected.
add_test( TestIsPointInCircle_1
  IsPointInCircle -0.2 -0.2 0
  )
set_tests_properties(
    TestIsPointInCircle_1 
    PROPERTIES WILL_FAIL "TRUE"
    )

# in test two we take one point of the square and move it TOWARD the center
# of the circle. The point is then inside, the test passes.
add_test( TestIsPointInCircle_2
  IsPointInCircle 0.2 0.2 0
  )

# in test 3 we move the point only sightly, the point is now outside, but
# rounding errors make the test pass. This is the example of the limitation 
# of float arithmetic
add_test( TestIsPointInCircle_3
  IsPointInCircle 0 -0.00000000000000002 0
  )

# in test 4, we move the point as in test 3 but we use exact geometrical
# predicate. This time, the test as the correct behavior.
add_test( TestIsPointInCircle_4
  IsPointInCircle 0 -0.00000000000000002 1
  )
set_tests_properties(
    TestIsPointInCircle_4 
    PROPERTIES WILL_FAIL "TRUE"
    )
