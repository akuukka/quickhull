add_executable(QuickHullTests)
target_sourceS(QuickHullTests
    PUBLIC
        main.cpp
        QuickHullTests.cpp
        QuickHullTests.hpp
)
target_link_libraries(QuickHullTests
    PRIVATE
        QuickHull
)
add_test(NAME QuickHullTests
         COMMAND QuickHullTests)
