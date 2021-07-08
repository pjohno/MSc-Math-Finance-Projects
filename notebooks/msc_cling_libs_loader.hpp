
// add include and lib paths for this repository
#pragma cling add_include_path("../include")
#pragma cling add_library_path("../lib")
// add include and lib paths for MATH60082
#pragma cling add_include_path("../MATH60082/include")
#pragma cling add_library_path("../MATH60082/lib")

// load libraries
#pragma cling load("boost_system")
#pragma cling load("math60082")
#pragma cling load("mscmathfinanceprojects")
#pragma cling load("gsl")
#pragma cling load("cblas")

